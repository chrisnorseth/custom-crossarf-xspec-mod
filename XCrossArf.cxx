#include "XCrossArf.h"

#include <XSstreams.h>
#include <fitsio.h>

// from heasp
#include <SPutils.h>
#include <arf.h>

#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>
#include <vector>

namespace {

std::string stripFitsSpecifier(const std::string& path)
{
  const std::string::size_type bracket = path.find('[');
  if (bracket == std::string::npos) return path;
  return path.substr(0, bracket);
}

bool startsWith(const std::string& text, const std::string& prefix)
{
  return text.size() >= prefix.size() &&
         text.compare(0, prefix.size(), prefix) == 0;
}

std::string trim(const std::string& text)
{
  std::string::size_type begin = 0;
  while (begin < text.size() &&
         std::isspace(static_cast<unsigned char>(text[begin]))) {
    ++begin;
  }
  std::string::size_type end = text.size();
  while (end > begin &&
         std::isspace(static_cast<unsigned char>(text[end - 1]))) {
    --end;
  }
  return text.substr(begin, end - begin);
}

std::string cleanFitsString(const std::string& text)
{
  std::string value = trim(text);
  if (value.size() >= 2 && value.front() == '\'' && value.back() == '\'') {
    value = trim(value.substr(1, value.size() - 2));
  }
  return value;
}

std::vector<std::string> splitWords(const std::string& text)
{
  std::istringstream input(text);
  std::vector<std::string> words;
  std::string word;
  while (input >> word) words.push_back(word);
  return words;
}

std::string dirname(const std::string& path)
{
  const std::string::size_type pos = path.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return path.substr(0, pos);
}

bool isAbsolutePath(const std::string& path)
{
  return !path.empty() && path[0] == '/';
}

std::string resolvePath(const std::string& baseDir, const std::string& path)
{
  if (path.empty() || isAbsolutePath(path)) return path;
  return baseDir + "/" + path;
}

struct ConfigSpec
{
  int region = -1;
  std::string obs;
  std::string mixGroup;
  std::string pha;
  std::string rmf;
  std::string arf;
};

struct ConfigData
{
  bool loaded = false;
  bool packed = false;
  size_t nSources = 0;
  std::string modelExpr;
  std::map<size_t, ConfigSpec> spectra;
  std::map<std::pair<int, int>, std::string> arfs;
};

ConfigData readConfig()
{
  ConfigData cfg;
  const char* pathValue = std::getenv("XCRSARF_CONFIG");
  if (!pathValue || !*pathValue) return cfg;

  const std::string path(pathValue);
  std::ifstream input(path.c_str());
  if (!input) return cfg;

  cfg.loaded = true;
  const std::string baseDir = dirname(path);
  std::string line;
  while (std::getline(input, line)) {
    const std::string::size_type comment = line.find('#');
    if (comment != std::string::npos) line = line.substr(0, comment);
    line = trim(line);
    if (line.empty()) continue;

    const std::vector<std::string> words = splitWords(line);
    if (words.empty()) continue;
    if (words[0] == "mode" && words.size() >= 2) {
      cfg.packed = words[1] == "packed";
    } else if (words[0] == "packed" && words.size() >= 2) {
      cfg.packed = words[1] == "yes" || words[1] == "true" || words[1] == "1";
    } else if (words[0] == "sources" && words.size() >= 2) {
      cfg.nSources = static_cast<size_t>(std::stoi(words[1]));
    } else if (words[0] == "model" && words.size() >= 2) {
      std::string expr = line.substr(line.find(words[1]));
      cfg.modelExpr = trim(expr);
    } else if (words[0] == "spectrum" && words.size() >= 7) {
      ConfigSpec spec;
      const size_t specNum = static_cast<size_t>(std::stoi(words[1]));
      spec.region = std::stoi(words[2]);
      spec.obs = words[3];
      spec.pha = resolvePath(baseDir, words[4]);
      spec.rmf = resolvePath(baseDir, words[5]);
      spec.arf = resolvePath(baseDir, words[6]);
      spec.mixGroup = words.size() >= 8 ? words[7] : spec.obs;
      cfg.spectra[specNum] = spec;
    } else if (words[0] == "arf" && words.size() >= 4) {
      const int sourceRegion = std::stoi(words[1]);
      const int targetSpec = std::stoi(words[2]);
      cfg.arfs[std::make_pair(targetSpec, sourceRegion)] =
          resolvePath(baseDir, words[3]);
    }
  }
  return cfg;
}

} // namespace

XCrossArf::XCrossArf(const string& name)
  : MixUtility(name)
{
}

XCrossArf::~XCrossArf()
{
}

void XCrossArf::initialize(const std::vector<Real>& params,
			   const IntegerVector& spectrumNums,
			   const size_t sourceNum,
			   const std::string& modelName)
{
  if (!spectrumNums.size())
  {
     string msg("XCrossArf Model Error: No data assigned, mixing will not be performed.\n");
     throw YellowAlert(msg);
  }

  IntegerVector allSpectrumNums;
  for (size_t iSpec = 1; iSpec <= numberOfSpectra(); ++iSpec) {
    allSpectrumNums.push_back(static_cast<int>(iSpec));
  }

  const ConfigData config = readConfig();

  std::map<int, int> regionMap;
  for (auto& sNum : allSpectrumNums) {
    const int region = regionNumber(sNum);
    if (region > 0) regionMap[region] = sNum;
  }

  size_t nSources = config.nSources;
  if (!nSources) {
    nSources = regionMap.empty() ? 0 : static_cast<size_t>(regionMap.rbegin()->first);
  }
  if (sourceNum < 1 || sourceNum > nSources) {
    std::stringstream msg;
    msg << "XCrossArf: source number " << sourceNum
        << " does not correspond to a configured source/region.\n";
    throw YellowAlert(msg.str());
  }

  std::vector<IntegerVector> mixSets(1);
  mixSets[0] = allSpectrumNums;
  setSpecNumsMixSets(mixSets);

  const size_t nTargets = allSpectrumNums.size();
  std::vector<std::vector<RealArray> > mixFact(nTargets);
  RealArray mixEnergies;

  if (config.loaded && config.packed) {
    std::map<std::pair<std::string, int>, size_t> specForMixGroupRegion;
    for (const auto& specEntry : config.spectra) {
      specForMixGroupRegion[std::make_pair(specEntry.second.mixGroup,
                                           specEntry.second.region)] =
          specEntry.first;
    }

    for (size_t itarg = 0; itarg < nTargets; itarg++) {
      mixFact[itarg].resize(nTargets);
      const size_t targetSpec = allSpectrumNums[itarg];
      std::map<size_t, ConfigSpec>::const_iterator targetInfo =
          config.spectra.find(targetSpec);
      if (targetInfo == config.spectra.end()) continue;

      const string diagName = targetInfo->second.arf;
      if (diagName.empty()) {
        std::stringstream msg;
        msg << "XCrossArf: no diagonal ARF for spectrum " << targetSpec << "\n";
        throw YellowAlert(msg.str());
      }

      heasp::arf diagARF;
      if (diagARF.read(diagName)) {
        std::stringstream msg;
        msg << "XCrossArf: failed reading diagonal ARF " << diagName << "\n";
        throw YellowAlert(msg.str());
      }
      const std::vector<Real>& diagEff = diagARF.getEffArea();
      if (mixEnergies.size() == 0) {
        const std::vector<Real>& lowE = diagARF.getLowEnergy();
        const std::vector<Real>& highE = diagARF.getHighEnergy();
        mixEnergies.resize(lowE.size() + 1);
        mixEnergies[0] = lowE[0];
        for (size_t ie = 0; ie < highE.size(); ie++) mixEnergies[ie + 1] = highE[ie];
      }

      for (size_t sourceRegion = 1; sourceRegion <= nSources; ++sourceRegion) {
        std::map<std::pair<std::string, int>, size_t>::const_iterator inputSpec =
            specForMixGroupRegion.find(std::make_pair(targetInfo->second.mixGroup,
                                                      static_cast<int>(sourceRegion)));
        if (inputSpec == specForMixGroupRegion.end()) continue;

        string sourceName;
        std::map<std::pair<int, int>, std::string>::const_iterator found =
            config.arfs.find(std::make_pair(static_cast<int>(targetSpec),
                                            static_cast<int>(sourceRegion)));
        if (found != config.arfs.end()) sourceName = found->second;
        if (sourceName.empty() &&
            static_cast<int>(sourceRegion) == targetInfo->second.region) {
          sourceName = diagName;
        }
        if (sourceName.empty()) continue;

        heasp::arf sourceARF;
        if (sourceARF.read(sourceName)) {
          std::stringstream msg;
          msg << "XCrossArf: failed reading source ARF " << sourceName << "\n";
          throw YellowAlert(msg.str());
        }
        const std::vector<Real>& srcEff = sourceARF.getEffArea();
        if (diagEff.size() != srcEff.size()) {
          std::stringstream msg;
          msg << "XCrossArf: ARF bin count mismatch for " << sourceName
              << " and " << diagName << "\n";
          throw YellowAlert(msg.str());
        }

        size_t isrc = 0;
        for (; isrc < nTargets; ++isrc) {
          if (static_cast<size_t>(allSpectrumNums[isrc]) == inputSpec->second) break;
        }
        if (isrc == nTargets) continue;

        RealArray factor(srcEff.size());
        for (size_t ie = 0; ie < srcEff.size(); ie++) {
          factor[ie] = diagEff[ie] > 0.0 ? srcEff[ie] / diagEff[ie] : 0.0;
        }
        mixFact[itarg][isrc] = factor;
      }
    }

    tcout << "XCrossArf: packed mode using " << nSources
          << " source regions for " << nTargets << " target spectra."
          << std::endl;

    setMixingFactors(mixFact);
    setMixingEnergies(mixEnergies);
    verifyData();
    return;
  }

  for (size_t itarg = 0; itarg < nTargets; itarg++) {
    mixFact[itarg].resize(nTargets);
    const size_t targetSpec = allSpectrumNums[itarg];
    const int targetRegion = regionNumber(targetSpec);

    string diagName;
    if (config.loaded) {
      std::map<size_t, ConfigSpec>::const_iterator found =
          config.spectra.find(targetSpec);
      if (found != config.spectra.end()) diagName = found->second.arf;
    }
    if (diagName.empty()) diagName = arfName(targetSpec, 1);
    if (diagName.empty()) {
      std::stringstream msg;
      msg << "XCrossArf: no diagonal ARF for spectrum " << targetSpec << "\n";
      throw YellowAlert(msg.str());
    }

    string sourceName;
    if (config.loaded) {
      std::map<std::pair<int, int>, std::string>::const_iterator found =
          config.arfs.find(std::make_pair(static_cast<int>(targetSpec),
                                          static_cast<int>(sourceNum)));
      if (found != config.arfs.end()) sourceName = found->second;
      if (sourceName.empty() && static_cast<int>(sourceNum) == targetRegion) {
        sourceName = diagName;
      }
    } else {
      if (static_cast<int>(sourceNum) == targetRegion) {
        sourceName = diagName;
      } else {
        const std::map<string, string> xflt = xfltStrings(targetSpec);
        std::stringstream key;
        key << "CrossArfFrom" << sourceNum << "To" << targetRegion;
        std::map<string, string>::const_iterator found = xflt.find(key.str());
        if (found != xflt.end()) sourceName = found->second;
      }
    }
    if (sourceName.empty()) continue;

    heasp::arf diagARF;
    heasp::arf sourceARF;
    if (diagARF.read(diagName) || sourceARF.read(sourceName)) {
      std::stringstream msg;
      msg << "XCrossArf: failed reading ARF pair " << sourceName
          << " and " << diagName << "\n";
      throw YellowAlert(msg.str());
    }
    const std::vector<Real>& diagEff = diagARF.getEffArea();
    const std::vector<Real>& srcEff = sourceARF.getEffArea();
    if (diagEff.size() != srcEff.size()) {
      std::stringstream msg;
      msg << "XCrossArf: ARF bin count mismatch for " << sourceName
          << " and " << diagName << "\n";
      throw YellowAlert(msg.str());
    }

    if (mixEnergies.size() == 0) {
      const std::vector<Real>& lowE = sourceARF.getLowEnergy();
      const std::vector<Real>& highE = sourceARF.getHighEnergy();
      mixEnergies.resize(lowE.size() + 1);
      mixEnergies[0] = lowE[0];
      for (size_t ie = 0; ie < highE.size(); ie++) mixEnergies[ie + 1] = highE[ie];
    }

    RealArray factor(srcEff.size());
    for (size_t ie = 0; ie < srcEff.size(); ie++) {
      factor[ie] = diagEff[ie] > 0.0 ? srcEff[ie] / diagEff[ie] : 0.0;
    }
    mixFact[itarg][itarg] = factor;
  }

  tcout << "XCrossArf: source " << sourceNum
        << " using source-to-target ARFs for " << nTargets
        << " target spectra." << std::endl;

  setMixingFactors(mixFact);
  setMixingEnergies(mixEnergies);
  verifyData();
}

void XCrossArf::initializeForFit(const std::vector<Real>& params,
                                 bool paramsAreFrozen)
{
}

void XCrossArf::perform(const EnergyPointer& energy,
                        const std::vector<Real>& params,
                        GroupFluxContainer& flux)
{
  doMix(energy, flux);
}

void XCrossArf::verifyData()
{
}

std::map<std::pair<int, int>, string> XCrossArf::setArfFiles(
    const IntegerVector& spectrumNums, const size_t sourceNum)
{
  std::map<std::pair<int, int>, string> arffiles;

  for (auto& sNum : spectrumNums) {
    const int sourceRegion = regionNumber(sNum);
    if (sourceRegion <= 0) continue;

    string filename = arfName(sNum, sourceNum);
    if (filename.empty() && sourceNum != 1) filename = arfName(sNum, 1);
    arffiles[std::make_pair(sourceRegion, sourceRegion)] = filename;

    const std::map<string, string> xflt = xfltStrings(sNum);
    for (auto& xf : xflt) {
      const string& key = xf.first;
      if (startsWith(key, "CrossArfFrom")) {
	size_t tpt = key.find("To", 12);
	if (tpt == string::npos) {
	  string msg = "XFLT keyword from spectrum " + std::to_string(sNum)
	      + " has wrong format: cannot read from " + key + "\n";
	  throw YellowAlert(msg);
	}
	int sregnum, tregnum;
	try {
	  sregnum = stoi(key.substr(12, tpt - 12));
	} catch(...) {
	  string msg = "XFLT keyword from spectrum " + std::to_string(sNum)
	    + " has wrong format: cannot read integer from "
	    + key.substr(12, tpt - 12) + "\n";
	  throw YellowAlert(msg);
	}
	try {
	  tregnum = stoi(key.substr(tpt + 2, string::npos));
	} catch(...) {
	  string msg = "XFLT keyword from spectrum " + std::to_string(sNum)
	    + " has wrong format: cannot read integer from "
	    + key.substr(tpt + 2, string::npos) + "\n";
	  throw YellowAlert(msg);
	}
	arffiles[std::make_pair(tregnum, sregnum)] = xf.second;
      }
    }
  }

  return arffiles;
}

std::map<string, string> XCrossArf::xfltStrings(const size_t spectrumNumber) const
{
  std::map<string, string> values = xfltStringsFromHeader(spectrumNumber);
  const std::map<string, string> loadedValues = allXfltValuesStr(spectrumNumber);
  if (!values.empty()) {
    static std::set<size_t> reportedSpectra;
    if (values.size() > loadedValues.size() &&
        reportedSpectra.insert(spectrumNumber).second) {
      tcout << "XCrossArf: spectrum " << spectrumNumber << " using "
            << values.size() << " XFLT keywords read directly from PHA header"
            << " (" << loadedValues.size() << " loaded by XSPEC)." << std::endl;
    }
    return values;
  }
  return loadedValues;
}

std::map<string, string> XCrossArf::xfltStringsFromHeader(
    const size_t spectrumNumber) const
{
  std::map<string, string> values;
  const std::string phaPath = stripFitsSpecifier(fullPathName(spectrumNumber));

  fitsfile* fptr = 0;
  int status = 0;
  if (fits_open_file(&fptr, phaPath.c_str(), READONLY, &status)) {
    return values;
  }

  int hduType = 0;
  int moveStatus = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, const_cast<char*>("SPECTRUM"), 0,
                  &moveStatus);
  if (moveStatus) {
    moveStatus = 0;
    fits_movabs_hdu(fptr, 2, &hduType, &moveStatus);
  }
  if (moveStatus) {
    fits_close_file(fptr, &status);
    return values;
  }

  int nkeys = 0;
  if (fits_get_hdrspace(fptr, &nkeys, 0, &status)) {
    fits_close_file(fptr, &status);
    return values;
  }

  for (int ikey = 1; ikey <= nkeys; ++ikey) {
    char keyname[FLEN_KEYWORD] = {0};
    char keyval[FLEN_VALUE] = {0};
    char comment[FLEN_COMMENT] = {0};
    int keyStatus = 0;
    if (fits_read_keyn(fptr, ikey, keyname, keyval, comment, &keyStatus))
      continue;

    const std::string fitsKey(keyname);
    if (!startsWith(fitsKey, "XFLT")) continue;

    const std::string keyValue = cleanFitsString(keyval);
    const std::string::size_type colon = keyValue.find(':');
    if (colon == std::string::npos) continue;
    values[trim(keyValue.substr(0, colon))] = trim(keyValue.substr(colon + 1));
  }

  fits_close_file(fptr, &status);
  return values;
}

int XCrossArf::regionNumber(const size_t spectrumNumber) const
{
  const ConfigData config = readConfig();
  if (config.loaded) {
    std::map<size_t, ConfigSpec>::const_iterator found =
        config.spectra.find(spectrumNumber);
    if (found != config.spectra.end()) return found->second.region;
  }

  const std::map<string, string> xflt = xfltStrings(spectrumNumber);
  std::map<string, string>::const_iterator match = xflt.find("CrossArfRegion");
  if (match == xflt.end()) return -1;
  try {
    return std::stoi(match->second);
  } catch(...) {
    return -1;
  }
}
