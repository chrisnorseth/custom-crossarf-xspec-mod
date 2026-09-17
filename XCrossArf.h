#ifndef XCROSSARF_H
#define XCROSSARF_H 1

#include <map>
#include <string>
#include <utility>

#include <xsTypes.h>
#include <XSModel/Model/Component/MixUtility.h>

class XCrossArf : public MixUtility
{
 public:
  XCrossArf(const string& name);
  virtual ~XCrossArf();

  void initialize(const std::vector<Real>& params, const IntegerVector& specNums,
                  const size_t sourceNum, const std::string& modelName);
  void perform(const EnergyPointer& energy, const std::vector<Real>& params,
               GroupFluxContainer& flux);
  void initializeForFit(const std::vector<Real>& params, bool paramsAreFrozen);

 protected:
  void verifyData();

 private:
  XCrossArf(const XCrossArf& right);
  XCrossArf& operator=(const XCrossArf& right);

  std::map<std::pair<int, int>, string> setArfFiles(
      const IntegerVector& spectrumNums, const size_t sourceNum);
  std::map<string, string> xfltStrings(const size_t spectrumNumber) const;
  std::map<string, string> xfltStringsFromHeader(
      const size_t spectrumNumber) const;
  int regionNumber(const size_t spectrumNumber) const;
};

#endif
