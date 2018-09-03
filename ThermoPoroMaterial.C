// $Id$
//==============================================================================
//!
//! \file ThermoPoroMaterial.C
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for thermo-poro-elastic material models.
//!
//==============================================================================


#include "ThermoPoroMaterial.h"

#include "Functions.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"


  template<class T>
static bool propertyParse(PoroMaterial::FuncConstPair<T>& data,
                          const TiXmlElement* elem,
                          const std::string& attr,
                          const std::string& tag)
{
  std::string constant;
  if (utl::getAttribute(elem,attr.c_str(),constant)) {
    std::stringstream str;
    str << constant;
    str >> data.constant;
    return true;
  }

  const TiXmlElement* child = elem->FirstChildElement(tag);
  if (child) {
    IFEM::cout <<" ";
    std::string type;
    utl::getAttribute(child,"type",type,true);
    const TiXmlNode* aval;
    if ((aval = child->FirstChild()))
      data.function = data.parse(aval->Value(),type);

    return data.function != nullptr;
  }

 return false;
}


void ThermoPoroMaterial::parse(const TiXmlElement* elem)
{
  propertyParse(fheatcapacity, elem, "cpf", "fluidheatcapacity");
  propertyParse(sheatcapacity, elem, "cps", "solidheatcapacity");
  propertyParse(fconductivity, elem, "kappaf", "fluidconductivity");
  propertyParse(sconductivity, elem, "kappas", "solidconductivity");
  propertyParse(sexpansion, elem, "alphas", "solidexpansion");
  this->PoroMaterial::parse(elem);
}


void ThermoPoroMaterial::printLog() const
{
  this->PoroMaterial::printLog();
  IFEM::cout << "\tHeat capacities: "
             << "\n\t\tHeat capacity of Fluid, cpf= " << fheatcapacity.constant
             << "\n\t\tHeat capacity of Solid, cps= " << sheatcapacity.constant << std::endl;
  IFEM::cout << "\tThermal Conductivities: "
             << "\n\t\tThermal Conductivity of Fluid, kappaf= " << fconductivity.constant
             << "\n\t\tThermal Conductivity of Solid, kappas= " << sconductivity.constant
             << std::endl;
}


double ThermoPoroMaterial::getHeatCapacity(double T) const
{
  Vec3 X;
  return getPorosity(X)*getFluidDensity(X)*getFluidHeatCapacity(T) +
         (1.0-getPorosity(X))*getSolidDensity(X)*getSolidHeatCapacity(T);
}


double ThermoPoroMaterial::getFluidHeatCapacity (double T) const
{
  return fheatcapacity.evaluate(T);
}


double ThermoPoroMaterial::getSolidHeatCapacity (double T) const
{
  return sheatcapacity.evaluate(T);
}


double ThermoPoroMaterial::getFluidThermalConductivity(double T) const
{
  return fconductivity.evaluate(T);
}


double ThermoPoroMaterial::getSolidThermalConductivity(double T) const
{
  return sconductivity.evaluate(T);
}


double ThermoPoroMaterial::getThermalConductivity(double T) const
{
  Vec3 X;
  return pow(getFluidThermalConductivity(T),getPorosity(X))*
         pow(getSolidThermalConductivity(T),1.0-getPorosity(X));
}


double ThermoPoroMaterial::getSolidThermalExpansion(double T) const
{
  return sexpansion.evaluate(T);
}
