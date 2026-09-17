#include <XSContainer.h>
#include <XSModel/Data/Detector/MultiBinningResponse.h>
#include <XSModel/Data/Detector/Response.h>
#include <XSModel/Data/SpectralData.h>
#include <XSModel/GlobalContainer/DataContainer.h>
#include <XSUser/Global/XSGlobal.h>
#include <XSUser/UserInterface/xstcl.h>
#include <XSstreams.h>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct SpecConfig
{
  size_t specNum = 0;
  int region = 0;
  std::string obs;
  std::string mixGroup;
  std::string pha;
  std::string rmf;
  std::string arf;
};

struct RunConfig
{
  bool valid = true;
  bool packed = false;
  size_t nSources = 0;
  std::string modelExpr;
  std::vector<std::string> modelParams;
  std::string refObs;
  std::map<size_t, std::string> sourceModels;
  std::map<size_t, std::vector<std::string> > sourceParams;
  std::vector<SpecConfig> spectra;
  std::vector<std::string> arfs;
  std::map<std::pair<size_t, size_t>, std::string> arfMap;
};

std::string trim(const std::string& text)
{
  std::string::size_type begin = 0;
  while (begin < text.size() && std::isspace(static_cast<unsigned char>(text[begin]))) {
    ++begin;
  }
  std::string::size_type end = text.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
    --end;
  }
  return text.substr(begin, end - begin);
}

std::vector<std::string> splitWords(const std::string& text)
{
  std::istringstream input(text);
  std::vector<std::string> words;
  std::string word;
  while (input >> word) words.push_back(word);
  return words;
}

bool parsePositiveIndex(const std::string& text, size_t& value)
{
  std::string number = text;
  const std::string::size_type colon = number.find(':');
  if (colon != std::string::npos) number = number.substr(0, colon);
  if (number.size() > 1 && number[0] == 's') number = number.substr(1);
  if (number.empty()) return false;
  for (char c : number) {
    if (!std::isdigit(static_cast<unsigned char>(c))) return false;
  }
  value = static_cast<size_t>(std::stoul(number));
  return value > 0;
}

bool startsWith(const std::string& text, const std::string& prefix)
{
  return text.size() >= prefix.size() &&
         text.compare(0, prefix.size(), prefix) == 0;
}

bool numericToken(const std::string& text)
{
  if (text.empty()) return false;
  char* end = 0;
  std::strtod(text.c_str(), &end);
  return end && *end == '\0';
}

bool numericRow(const std::vector<std::string>& words)
{
  if (words.empty()) return false;
  for (const auto& word : words) {
    if (!numericToken(word)) return false;
  }
  return true;
}

std::string normalizeModelExpression(const std::string& text)
{
  std::string expr = trim(text);
  if (startsWith(expr, "xcrsarf*(") && expr.size() > 9 && expr.back() == ')') {
    expr = trim(expr.substr(9, expr.size() - 10));
  }
  return expr;
}

std::string compactModelExpression(const std::string& text)
{
  std::string compact;
  for (char c : text) {
    if (!std::isspace(static_cast<unsigned char>(c))) {
      compact += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
  }
  std::string::size_type pos = 0;
  while ((pos = compact.find("constant", pos)) != std::string::npos) {
    compact.replace(pos, 8, "const");
    pos += 5;
  }
  pos = 0;
  while ((pos = compact.find("powerlaw", pos)) != std::string::npos) {
    compact.replace(pos, 8, "po");
    pos += 2;
  }
  return compact;
}

bool isApecConst(const std::string& expr)
{
  return compactModelExpression(expr) == "apec*const";
}

bool isPowerlawConst(const std::string& expr)
{
  return compactModelExpression(expr) == "po*const";
}

bool isApecPowerlawConst(const std::string& expr)
{
  const std::string compact = compactModelExpression(expr);
  return compact == "(apec+po)*const" || compact == "(po+apec)*const";
}

std::string dirname(const std::string& path)
{
  const std::string::size_type pos = path.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return path.substr(0, pos);
}

std::string resolvePath(const std::string& baseDir, const std::string& path)
{
  if (path.empty() || path[0] == '/') return path;
  return baseDir + "/" + path;
}

std::string stripFitsSpecifier(const std::string& path)
{
  const std::string::size_type bracket = path.find('[');
  if (bracket == std::string::npos) return path;
  return path.substr(0, bracket);
}

bool readableFile(const std::string& path)
{
  std::ifstream input(stripFitsSpecifier(path).c_str());
  return input.good();
}

std::string tclQuote(const std::string& text)
{
  std::string out("{");
  for (char c : text) {
    if (c == '{' || c == '}' || c == '\\') out += '\\';
    out += c;
  }
  out += '}';
  return out;
}

RunConfig parseRunConfig(const std::string& path)
{
  RunConfig cfg;
  const std::string baseDir = dirname(path);
  size_t activeModel = static_cast<size_t>(-1);
  std::ifstream input(path.c_str());
  if (!input) {
    tcout << "xcrsarfrun: cannot open " << path << std::endl;
    return cfg;
  }

  std::string line;
  while (std::getline(input, line)) {
    const std::string::size_type comment = line.find('#');
    if (comment != std::string::npos) line = line.substr(0, comment);
    line = trim(line);
    if (line.empty()) continue;

    const std::vector<std::string> words = splitWords(line);
    if (words.empty()) continue;
    if (numericRow(words)) {
      if (activeModel == static_cast<size_t>(-1)) {
        tcout << "xcrsarfrun: parameter row has no preceding model line: "
              << line << std::endl;
        cfg.valid = false;
        return cfg;
      }
      if (activeModel == 0) {
        cfg.modelParams.push_back(line);
      } else {
        cfg.sourceParams[activeModel].push_back(line);
      }
    } else if (words[0] == "mode" && words.size() >= 2) {
      activeModel = static_cast<size_t>(-1);
      cfg.packed = words[1] == "packed";
    } else if (words[0] == "packed" && words.size() >= 2) {
      activeModel = static_cast<size_t>(-1);
      cfg.packed = words[1] == "yes" || words[1] == "true" || words[1] == "1";
    } else if (words[0] == "sources" && words.size() >= 2) {
      activeModel = static_cast<size_t>(-1);
      cfg.nSources = static_cast<size_t>(std::stoi(words[1]));
    } else if (words[0] == "model" && words.size() >= 2) {
      size_t sourceNum = 0;
      if (words.size() >= 3 && parsePositiveIndex(words[1], sourceNum)) {
        cfg.sourceModels[sourceNum] =
            normalizeModelExpression(line.substr(line.find(words[2])));
        cfg.sourceParams[sourceNum].clear();
        activeModel = sourceNum;
      } else {
        cfg.modelExpr = normalizeModelExpression(line.substr(line.find(words[1])));
        cfg.modelParams.clear();
        activeModel = 0;
      }
    } else if (words[0] == "source_model" && words.size() >= 3) {
      size_t sourceNum = 0;
      if (!parsePositiveIndex(words[1], sourceNum)) {
        tcout << "xcrsarfrun: invalid source_model index: " << words[1]
              << std::endl;
        cfg.valid = false;
        return cfg;
      }
      cfg.sourceModels[sourceNum] =
          normalizeModelExpression(line.substr(line.find(words[2])));
      cfg.sourceParams[sourceNum].clear();
      activeModel = sourceNum;
    } else if (words[0] == "reference_obs" && words.size() >= 2) {
      activeModel = static_cast<size_t>(-1);
      cfg.refObs = words[1];
    } else if (words[0] == "spectrum" && words.size() >= 7) {
      activeModel = static_cast<size_t>(-1);
      SpecConfig spec;
      spec.specNum = static_cast<size_t>(std::stoi(words[1]));
      spec.region = std::stoi(words[2]);
      spec.obs = words[3];
      spec.pha = resolvePath(baseDir, words[4]);
      spec.rmf = resolvePath(baseDir, words[5]);
      spec.arf = resolvePath(baseDir, words[6]);
      spec.mixGroup = words.size() >= 8 ? words[7] : spec.obs;
      cfg.spectra.push_back(spec);
    } else if (words[0] == "arf" && words.size() >= 4) {
      activeModel = static_cast<size_t>(-1);
      const size_t sourceRegion = static_cast<size_t>(std::stoi(words[1]));
      const size_t targetSpec = static_cast<size_t>(std::stoi(words[2]));
      const std::string arfPath = resolvePath(baseDir, words[3]);
      cfg.arfs.push_back(arfPath);
      cfg.arfMap[std::make_pair(sourceRegion, targetSpec)] = arfPath;
    }
  }

  if (!cfg.nSources) {
    for (const auto& spec : cfg.spectra) {
      if (spec.region > static_cast<int>(cfg.nSources)) {
        cfg.nSources = static_cast<size_t>(spec.region);
      }
    }
  }
  if (cfg.refObs.empty() && !cfg.spectra.empty()) cfg.refObs = cfg.spectra[0].obs;
  return cfg;
}

int evalLine(Tcl_Interp* interp, const std::string& command)
{
  const int status = Tcl_Eval(interp, command.c_str());
  if (status != TCL_OK) {
    tcout << "xcrsarfrun: failed command: " << command << std::endl;
  }
  return status;
}

int checkFiles(const RunConfig& cfg)
{
  for (const auto& spec : cfg.spectra) {
    if (!readableFile(spec.pha)) {
      tcout << "xcrsarfrun: cannot read PHA for spectrum " << spec.specNum
            << ": " << spec.pha << std::endl;
      return TCL_ERROR;
    }
    if (!readableFile(spec.rmf)) {
      tcout << "xcrsarfrun: cannot read RMF for spectrum " << spec.specNum
            << ": " << spec.rmf << std::endl;
      return TCL_ERROR;
    }
    if (!readableFile(spec.arf)) {
      tcout << "xcrsarfrun: cannot read diagonal ARF for spectrum "
            << spec.specNum << ": " << spec.arf << std::endl;
      return TCL_ERROR;
    }
  }
  for (const auto& arf : cfg.arfs) {
    if (!readableFile(arf)) {
      tcout << "xcrsarfrun: cannot read cross ARF: " << arf << std::endl;
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

std::string modelForSource(const RunConfig& cfg, size_t sourceNum)
{
  const std::map<size_t, std::string>::const_iterator found =
      cfg.sourceModels.find(sourceNum);
  if (found != cfg.sourceModels.end()) return found->second;
  if (cfg.modelExpr.empty() && sourceNum != 1) {
    const std::map<size_t, std::string>::const_iterator first =
        cfg.sourceModels.find(1);
    if (first != cfg.sourceModels.end()) return first->second;
  }
  return cfg.modelExpr;
}

std::vector<std::string> paramsForSource(const RunConfig& cfg, size_t sourceNum)
{
  const std::map<size_t, std::vector<std::string> >::const_iterator found =
      cfg.sourceParams.find(sourceNum);
  if (found != cfg.sourceParams.end()) return found->second;
  if (cfg.modelParams.empty() && sourceNum != 1) {
    const std::map<size_t, std::vector<std::string> >::const_iterator first =
        cfg.sourceParams.find(1);
    if (first != cfg.sourceParams.end()) return first->second;
  }
  return cfg.modelParams;
}

bool sourceHasSpectrum(const RunConfig& cfg, size_t sourceNum)
{
  for (const auto& spec : cfg.spectra) {
    if (spec.region == static_cast<int>(sourceNum)) return true;
  }
  return false;
}

bool usesPackedApecPowerlawSuperset(const RunConfig& cfg)
{
  bool hasPowerlaw = false;
  bool hasSuperset = false;
  for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
    if (!sourceHasSpectrum(cfg, sourceNum)) continue;
    const std::string expr = modelForSource(cfg, sourceNum);
    if (isPowerlawConst(expr)) hasPowerlaw = true;
    else if (isApecPowerlawConst(expr)) hasSuperset = true;
  }
  return hasPowerlaw || hasSuperset || isApecPowerlawConst(cfg.modelExpr);
}

std::string packedModelExpression(const RunConfig& cfg)
{
  if (usesPackedApecPowerlawSuperset(cfg)) return "(apec+po)*const";
  return modelForSource(cfg, 1);
}

std::vector<std::string> inactiveApecRows()
{
  std::vector<std::string> rows;
  rows.push_back("1.0 -1 1.0 1.0 25.0 25.0");
  rows.push_back("0.4 -1 0.1 0.1 1.0 1.0");
  rows.push_back("0.0546 -1 0.01 0.01 0.30 0.30");
  rows.push_back("0.0 -1 0.0 0.0 1e-1 1e-1");
  return rows;
}

std::vector<std::string> inactivePowerlawRows()
{
  std::vector<std::string> rows;
  rows.push_back("0.0 -1 0.0 0.0 3.0 3.0");
  rows.push_back("0.0 -1 0.0 0.0 1e-1 1e-1");
  return rows;
}

std::vector<std::string> packedParamsForSource(const RunConfig& cfg,
                                               size_t sourceNum)
{
  const std::string expr = modelForSource(cfg, sourceNum);
  const std::vector<std::string> params = paramsForSource(cfg, sourceNum);
  if (!usesPackedApecPowerlawSuperset(cfg)) return params;
  if (isApecPowerlawConst(expr)) return params;

  std::vector<std::string> packed;
  if (isApecConst(expr)) {
    if (params.size() < 5) return params;
    packed.insert(packed.end(), params.begin(), params.begin() + 4);
    const std::vector<std::string> inactive = inactivePowerlawRows();
    packed.insert(packed.end(), inactive.begin(), inactive.end());
    packed.push_back(params[4]);
    return packed;
  }
  if (isPowerlawConst(expr)) {
    if (params.size() < 3) return params;
    const std::vector<std::string> inactive = inactiveApecRows();
    packed.insert(packed.end(), inactive.begin(), inactive.end());
    packed.push_back(params[0]);
    packed.push_back(params[1]);
    packed.push_back(params[2]);
    return packed;
  }
  return params;
}

int checkModels(const RunConfig& cfg)
{
  if (cfg.packed) {
    if (!sourceHasSpectrum(cfg, 1)) {
      tcout << "xcrsarfrun: packed mode source 1 must have a spectrum row "
               "because source slot 1 is reserved for the xpack model."
            << std::endl;
      return TCL_ERROR;
    }
    const bool useSuperset = usesPackedApecPowerlawSuperset(cfg);
    const std::vector<std::string> templateParams = paramsForSource(cfg, 1);
    for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
      const std::string expr = modelForSource(cfg, sourceNum);
      if (expr.empty()) {
        tcout << "xcrsarfrun: missing model for source " << sourceNum
              << ". Add 'model 1 EXPR' as the packed default, a bare "
                 "'model EXPR' line, or an explicit 'model " << sourceNum
              << " EXPR'." << std::endl;
        return TCL_ERROR;
      }
      const std::vector<std::string> params = paramsForSource(cfg, sourceNum);
      if (params.empty()) {
        tcout << "xcrsarfrun: missing parameter rows for source " << sourceNum
              << "." << std::endl;
        return TCL_ERROR;
      }
      if (!sourceHasSpectrum(cfg, sourceNum)) continue;
      if (useSuperset) {
        size_t expected = 0;
        if (isApecConst(expr)) expected = 5;
        else if (isPowerlawConst(expr)) expected = 3;
        else if (isApecPowerlawConst(expr)) expected = 7;
        if (!expected) {
          tcout << "xcrsarfrun: packed mode can auto-combine only "
                   "apec*const, po*const, and (apec+po)*const. Source "
                << sourceNum << " uses " << expr << "." << std::endl;
          return TCL_ERROR;
        }
        if (params.size() != expected) {
          tcout << "xcrsarfrun: source " << sourceNum << " uses " << expr
                << " and has " << params.size()
                << " parameter rows; expected " << expected << "."
                << std::endl;
          return TCL_ERROR;
        }
      } else if (params.size() != templateParams.size()) {
        tcout << "xcrsarfrun: packed mode source " << sourceNum
              << " has " << params.size() << " parameter rows, but source 1 "
              << "has " << templateParams.size() << "." << std::endl;
        return TCL_ERROR;
      }
    }
    for (const auto& params : cfg.sourceParams) {
      if (params.first > cfg.nSources) {
        tcout << "xcrsarfrun: parameter block for source " << params.first
              << " is outside sources " << cfg.nSources << "." << std::endl;
        return TCL_ERROR;
      }
    }
    for (const auto& model : cfg.sourceModels) {
      if (model.first > cfg.nSources) {
        tcout << "xcrsarfrun: source model " << model.first
              << " is outside sources " << cfg.nSources << "." << std::endl;
        return TCL_ERROR;
      }
    }
    return TCL_OK;
  }

  for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
    if (modelForSource(cfg, sourceNum).empty()) {
      tcout << "xcrsarfrun: missing model for source " << sourceNum
            << ". Add a default 'model EXPR' line or 'model " << sourceNum
            << " EXPR'." << std::endl;
      return TCL_ERROR;
    }
  }
  for (const auto& model : cfg.sourceModels) {
    if (model.first > cfg.nSources) {
      tcout << "xcrsarfrun: source model " << model.first
            << " is outside sources " << cfg.nSources << "." << std::endl;
      return TCL_ERROR;
    }
  }
  for (const auto& params : cfg.sourceParams) {
    if (params.first > cfg.nSources) {
      tcout << "xcrsarfrun: parameter block for source " << params.first
            << " is outside sources " << cfg.nSources << "." << std::endl;
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

std::string buildModelCommand(const RunConfig& cfg, size_t sourceNum)
{
  std::ostringstream cmd;
  cmd << "model " << sourceNum << ":s" << sourceNum
      << " xcrsarf*(" << modelForSource(cfg, sourceNum) << ")";

  const std::vector<std::string> params = paramsForSource(cfg, sourceNum);
  for (const auto& row : params) {
    cmd << " & " << row;
  }
  cmd << " & /*";
  return cmd.str();
}

std::vector<size_t> pointOnlySources(const RunConfig& cfg)
{
  std::vector<size_t> sources;
  for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
    if (!sourceHasSpectrum(cfg, sourceNum)) sources.push_back(sourceNum);
  }
  return sources;
}

std::string buildPointSourceModelCommand(const RunConfig& cfg, size_t sourceNum,
                                         size_t sourceSlot)
{
  const std::vector<std::string> params = paramsForSource(cfg, sourceNum);
  const size_t perGroup = params.size();
  const size_t packedPerGroup =
      packedParamsForSource(cfg, static_cast<size_t>(cfg.spectra[0].region)).size();
  std::map<std::string, size_t> obsAnchor;
  for (const auto& spec : cfg.spectra) {
    obsAnchor.insert(std::make_pair(spec.obs, spec.specNum));
  }
  const size_t firstSpec = cfg.spectra.front().specNum;

  std::ostringstream cmd;
  cmd << "model " << sourceSlot << ":s" << sourceNum << " "
      << modelForSource(cfg, sourceNum);
  for (const auto& spec : cfg.spectra) {
    const size_t anchorSpec = obsAnchor[spec.obs];
    const size_t packedConstantPar = anchorSpec * packedPerGroup;
    for (size_t row = 0; row < params.size(); ++row) {
      const bool isConstant = row + 1 == perGroup;
      if (isConstant) {
        cmd << " & = xpack:p" << packedConstantPar;
      } else if (spec.specNum == firstSpec) {
        cmd << " & " << params[row];
      } else {
        cmd << " & = s" << sourceNum << ":p" << (row + 1);
      }
    }
  }
  cmd << " & /*";
  return cmd.str();
}

int evalCommand(Tcl_Interp* interp, const char* command,
                const std::string& index, const std::string& filename);

std::string buildPackedModelCommand(const RunConfig& cfg)
{
  const size_t perGroup =
      packedParamsForSource(cfg, static_cast<size_t>(cfg.spectra[0].region)).size();
  std::map<int, size_t> refSpecForRegion;
  std::map<std::string, size_t> obsAnchor;
  for (const auto& spec : cfg.spectra) {
    if (spec.obs == cfg.refObs) {
      refSpecForRegion.insert(std::make_pair(spec.region, spec.specNum));
    }
    obsAnchor.insert(std::make_pair(spec.obs, spec.specNum));
  }

  std::ostringstream cmd;
  cmd << "model 1:xpack xcrsarf*(" << packedModelExpression(cfg) << ")";
  for (const auto& spec : cfg.spectra) {
    const std::vector<std::string> params =
        packedParamsForSource(cfg, static_cast<size_t>(spec.region));
    const std::map<int, size_t>::const_iterator ref =
        refSpecForRegion.find(spec.region);
    const size_t refFirstPar =
        ref == refSpecForRegion.end() ? 0 : (ref->second - 1) * perGroup + 1;
    const size_t anchorSpec = obsAnchor[spec.obs];
    const size_t anchorPar = anchorSpec * perGroup;

    for (size_t row = 0; row < params.size(); ++row) {
      if (row + 1 == perGroup) {
        if (spec.specNum == anchorSpec) {
          cmd << " & " << params[row];
        } else {
          cmd << " & = xpack:p" << anchorPar;
        }
      } else if (ref != refSpecForRegion.end() && spec.specNum != ref->second) {
        cmd << " & = xpack:p" << (refFirstPar + row);
      } else {
        cmd << " & " << params[row];
      }
    }
  }
  cmd << " & /*";
  return cmd.str();
}

int setupPackedPointSourceResponses(Tcl_Interp* interp, const RunConfig& cfg,
                                    const std::vector<size_t>& pointSources)
{
  for (const size_t sourceNum : pointSources) {
    bool hasAnyArf = false;
    for (const auto& spec : cfg.spectra) {
      const std::map<std::pair<size_t, size_t>, std::string>::const_iterator found =
          cfg.arfMap.find(std::make_pair(sourceNum, spec.specNum));
      if (found == cfg.arfMap.end()) continue;

      hasAnyArf = true;
      std::ostringstream index;
      index << sourceNum << ":" << spec.specNum;
      if (evalCommand(interp, "response", index.str(), spec.rmf) != TCL_OK) {
        return TCL_ERROR;
      }
      if (evalCommand(interp, "arf", index.str(), found->second) != TCL_OK) {
        return TCL_ERROR;
      }
    }
    if (!hasAnyArf) {
      tcout << "xcrsarfrun: point/source-only source " << sourceNum
            << " has no 'arf " << sourceNum
            << " TARGET ARF' rows." << std::endl;
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}

int evalCommand(Tcl_Interp* interp, const char* command,
                const std::string& index, const std::string& filename)
{
  Tcl_Obj* objv[3];
  objv[0] = Tcl_NewStringObj(command, -1);
  objv[1] = Tcl_NewStringObj(index.c_str(), -1);
  objv[2] = Tcl_NewStringObj(filename.c_str(), -1);
  for (int i = 0; i < 3; ++i) Tcl_IncrRefCount(objv[i]);

  const int status = (std::string(command) == "response")
      ? XSGlobal::xsResponse(0, interp, 3, objv)
      : XSGlobal::xsArf(0, interp, 3, objv);

  for (int i = 0; i < 3; ++i) Tcl_DecrRefCount(objv[i]);
  return status;
}

int setupSlots(Tcl_Interp* interp, size_t nSources)
{
  const size_t nSpectra = XSContainer::datasets->numberOfSpectra();
  if (!nSpectra) {
    tcout << "xcrsarfsetup: no spectra are loaded." << std::endl;
    return TCL_ERROR;
  }
  if (nSources < 2) return TCL_OK;

  for (size_t specNum = 1; specNum <= nSpectra; ++specNum) {
    const SpectralData* spec = XSContainer::datasets->lookup(specNum);
    const Response* response = spec ? spec->detectorResponse(0) : 0;
    if (!response || !response->toMultiBinningResponse()) {
      tcout << "xcrsarfsetup: spectrum " << specNum
            << " has no source-1 RMF response to copy." << std::endl;
      return TCL_ERROR;
    }
    const MultiBinningResponse* mb = response->toMultiBinningResponse();
    const std::string rmf = mb->rmfName();
    const std::string arf = mb->arfName();
    if (rmf.empty() || arf.empty()) {
      tcout << "xcrsarfsetup: spectrum " << specNum
            << " needs source-1 RMF and ARF loaded first." << std::endl;
      return TCL_ERROR;
    }

    for (size_t sourceNum = 2; sourceNum <= nSources; ++sourceNum) {
      std::ostringstream index;
      index << sourceNum << ":" << specNum;
      if (evalCommand(interp, "response", index.str(), rmf) != TCL_OK) return TCL_ERROR;
      if (evalCommand(interp, "arf", index.str(), arf) != TCL_OK) return TCL_ERROR;
    }
  }
  return TCL_OK;
}

} // namespace

extern "C" int XCrossArfSetupCmd(ClientData, Tcl_Interp* interp, int objc,
                                 Tcl_Obj* CONST objv[])
{
  const size_t nSpectra = XSContainer::datasets->numberOfSpectra();
  if (!nSpectra) {
    tcout << "xcrsarfsetup: no spectra are loaded." << std::endl;
    return TCL_ERROR;
  }

  size_t nSources = nSpectra;
  if (objc > 2) {
    tcout << "Usage: xcrsarfsetup ?number_of_sources?" << std::endl;
    return TCL_ERROR;
  }
  if (objc == 2) {
    int parsed = 0;
    if (Tcl_GetIntFromObj(interp, objv[1], &parsed) != TCL_OK || parsed < 1) {
      tcout << "xcrsarfsetup: number_of_sources must be a positive integer."
            << std::endl;
      return TCL_ERROR;
    }
    nSources = static_cast<size_t>(parsed);
  }

  if (nSources < 2) {
    tcout << "xcrsarfsetup: only one source requested; nothing to add."
          << std::endl;
    return TCL_OK;
  }

  if (setupSlots(interp, nSources) != TCL_OK) return TCL_ERROR;

  tcout << "xcrsarfsetup: activated " << nSources << " source slots for "
        << nSpectra << " spectra using each spectrum's diagonal RMF/ARF."
        << std::endl;
  return TCL_OK;
}

extern "C" int XCrossArfRunCmd(ClientData, Tcl_Interp* interp, int objc,
                               Tcl_Obj* CONST objv[])
{
  if (objc < 2 || objc > 3) {
    tcout << "Usage: xcrsarfrun config.txt ?model_expression?" << std::endl;
    return TCL_ERROR;
  }

  const std::string configPath = Tcl_GetString(objv[1]);
  RunConfig cfg = parseRunConfig(configPath);
  if (!cfg.valid || cfg.spectra.empty() || !cfg.nSources) return TCL_ERROR;
  if (objc == 3) {
    cfg.modelExpr = Tcl_GetString(objv[2]);
    cfg.modelParams.clear();
    cfg.sourceModels.clear();
    cfg.sourceParams.clear();
  }
  if (checkModels(cfg) != TCL_OK) return TCL_ERROR;
  if (checkFiles(cfg) != TCL_OK) return TCL_ERROR;

  setenv("XCRSARF_CONFIG", configPath.c_str(), 1);

  std::ostringstream data;
  data << "data";
  for (const auto& spec : cfg.spectra) {
    data << " " << spec.specNum << ":" << spec.specNum
         << " " << tclQuote(spec.pha);
  }
  if (evalLine(interp, data.str()) != TCL_OK) return TCL_ERROR;

  for (const auto& spec : cfg.spectra) {
    std::ostringstream resp;
    resp << "response 1:" << spec.specNum << " " << tclQuote(spec.rmf);
    if (evalLine(interp, resp.str()) != TCL_OK) return TCL_ERROR;

    std::ostringstream arf;
    arf << "arf 1:" << spec.specNum << " " << tclQuote(spec.arf);
    if (evalLine(interp, arf.str()) != TCL_OK) return TCL_ERROR;
  }

  if (cfg.packed) {
    const std::vector<size_t> pointSources = pointOnlySources(cfg);
    if (setupPackedPointSourceResponses(interp, cfg, pointSources) != TCL_OK) {
      return TCL_ERROR;
    }

    if (evalLine(interp, buildPackedModelCommand(cfg)) != TCL_OK) {
      return TCL_ERROR;
    }

    for (const size_t sourceNum : pointSources) {
      if (evalLine(interp, buildPointSourceModelCommand(cfg, sourceNum,
                                                       sourceNum)) != TCL_OK) {
        return TCL_ERROR;
      }
    }

    if (evalLine(interp, "tclout modpar xpack") != TCL_OK) return TCL_ERROR;
    const char* val = Tcl_GetVar(interp, "xspec_tclout", TCL_GLOBAL_ONLY);
    if (!val) return TCL_ERROR;
    const size_t totalPars = static_cast<size_t>(std::stoul(val));
    const size_t perGroup = totalPars / cfg.spectra.size();
    if (!perGroup) return TCL_ERROR;

    const size_t expectedPerGroup =
        packedParamsForSource(cfg, static_cast<size_t>(cfg.spectra[0].region)).size();
    if (perGroup != expectedPerGroup) {
      tcout << "xcrsarfrun: packed model has " << perGroup
            << " parameters per data group, but the config supplies "
            << expectedPerGroup << " packed parameter rows." << std::endl;
      return TCL_ERROR;
    }
    tcout << "xcrsarfrun: packed mode loaded " << cfg.spectra.size()
          << " spectra, created one cross-ARF packed model, and added "
          << pointSources.size() << " point/source-only model(s)." << std::endl;
    return TCL_OK;
  }

  if (setupSlots(interp, cfg.nSources) != TCL_OK) return TCL_ERROR;

  for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
    if (evalLine(interp, buildModelCommand(cfg, sourceNum)) != TCL_OK) {
      return TCL_ERROR;
    }
  }

  std::map<std::string, size_t> obsAnchor;
  for (const auto& spec : cfg.spectra) {
    obsAnchor.insert(std::make_pair(spec.obs, spec.specNum));
  }
  const size_t refAnchorSpec = obsAnchor[cfg.refObs];
  for (size_t sourceNum = 1; sourceNum <= cfg.nSources; ++sourceNum) {
    std::ostringstream modparCmd;
    modparCmd << "tclout modpar s" << sourceNum;
    if (evalLine(interp, modparCmd.str()) != TCL_OK) return TCL_ERROR;
    const char* val = Tcl_GetVar(interp, "xspec_tclout", TCL_GLOBAL_ONLY);
    if (!val) continue;
    const size_t totalPars = static_cast<size_t>(std::stoul(val));
    const size_t perGroup = totalPars / cfg.spectra.size();
    if (!perGroup) continue;
    const std::vector<std::string> params = paramsForSource(cfg, sourceNum);
    if (params.size() > perGroup) {
      tcout << "xcrsarfrun: source " << sourceNum << " has " << params.size()
            << " parameter rows, but the model has only " << perGroup
            << " parameters per data group." << std::endl;
      return TCL_ERROR;
    }

    for (const auto& spec : cfg.spectra) {
      const size_t parNum = spec.specNum * perGroup;
      std::ostringstream untie;
      untie << "untie s" << sourceNum << ":" << parNum;
      if (evalLine(interp, untie.str()) != TCL_OK) return TCL_ERROR;

      if (spec.obs == cfg.refObs) {
        const size_t refPar = refAnchorSpec * perGroup;
        if (sourceNum == 1 && spec.specNum == refAnchorSpec) {
          std::ostringstream set;
          set << "newpar s1:" << refPar << " 1";
          if (evalLine(interp, set.str()) != TCL_OK) return TCL_ERROR;
          set.str("");
          set << "freeze s1:" << refPar;
          if (evalLine(interp, set.str()) != TCL_OK) return TCL_ERROR;
        } else {
          std::ostringstream tie;
          tie << "newpar s" << sourceNum << ":" << parNum
              << " = s1:" << refPar;
          if (evalLine(interp, tie.str()) != TCL_OK) return TCL_ERROR;
        }
      } else {
        const size_t anchorSpec = obsAnchor[spec.obs];
        const size_t anchorPar = anchorSpec * perGroup;
        if (sourceNum == 1 && spec.specNum == anchorSpec) {
          std::ostringstream set;
          set << "newpar s1:" << anchorPar << " 1";
          if (evalLine(interp, set.str()) != TCL_OK) return TCL_ERROR;
        } else {
          std::ostringstream tie;
          tie << "newpar s" << sourceNum << ":" << parNum
              << " = s1:" << anchorPar;
          if (evalLine(interp, tie.str()) != TCL_OK) return TCL_ERROR;
        }
      }
    }
  }

  tcout << "xcrsarfrun: loaded " << cfg.spectra.size() << " spectra, created "
        << cfg.nSources << " source models, and tied last-parameter constants by observation."
        << std::endl;
  return TCL_OK;
}
