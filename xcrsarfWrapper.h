//C++
#ifndef XCRSARFWRAPPER_H
#define XCRSARFWRAPPER_H

#include <xsTypes.h>
#include <XSFunctions/Utilities/XSCall.h>

class MixUtility;

class xcrsarfWrapper
{
 public:
  static void modFunction(const EnergyPointer& energyArray,
			  const std::vector<Real>& parameterValues,
			  GroupFluxContainer& flux,
			  MixUtility* mixGenerator,
			  const std::string& modelName);
  static MixUtility* createUtility();
};

template <>
void XSCall<xcrsarfWrapper>::operator()(
    const EnergyPointer& energyArray,
    const std::vector<Real>& parameterValues,
    GroupFluxContainer& flux,
    MixUtility* mixGenerator,
    const string& modelName) const;

template <>
MixUtility* XSCall<xcrsarfWrapper>::getUtilityObject() const;

#endif
