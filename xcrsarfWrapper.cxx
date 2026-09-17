#include <xsTypes.h>

#include "XCrossArf.h"
#include "xcrsarfWrapper.h"

void xcrsarfWrapper::modFunction(const EnergyPointer& energyArray,
			   const std::vector<Real>& parameterValues,
			   GroupFluxContainer& flux,
			   MixUtility* mixUtility,
			   const std::string& modelName)
{
  mixUtility->perform(energyArray, parameterValues, flux);
}

MixUtility* xcrsarfWrapper::createUtility()
{
   return new XCrossArf("xcrsarf");
}

template <>
void XSCall<xcrsarfWrapper>::operator()(
    const EnergyPointer& energyArray,
    const std::vector<Real>& parameterValues,
    GroupFluxContainer& flux,
    MixUtility* mixGenerator,
    const string& modelName) const
{
   xcrsarfWrapper::modFunction(energyArray, parameterValues, flux,
                               mixGenerator, modelName);
}

template <>
MixUtility* XSCall<xcrsarfWrapper>::getUtilityObject() const
{
   return xcrsarfWrapper::createUtility();
}
