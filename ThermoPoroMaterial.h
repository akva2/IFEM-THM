// $Id$
//==============================================================================
//!
//! \file ThermoPoroMaterial.h
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for thermo-poro-elastic material models.
//!
//==============================================================================

#ifndef _THERMO_PORO_MATERIAL_H
#define _THERMO_PORO_MATERIAL_H

#include "Function.h"
#include "MatVec.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "PoroMaterial.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class ThermoPoroMaterial : public PoroMaterial
{
public:
  //! \brief Empty constructor.
  ThermoPoroMaterial() = default;

  //! \brief Empty destructor.
  virtual ~ThermoPoroMaterial() = default;

  //! \brief Parses material parementers from an XML element.
  void parse(const TiXmlElement*) override;

  //! \brief Prints out material parameters to the log stream.
  void printLog() const override;

  double getHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getFluidHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getSolidHeatCapacity(double T) const;
  //! \brief Evaluates the effective thermal conductivity at the current point
  double getThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the fluid at the current point
  double getFluidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the solid at the current point
  double getSolidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal expansion of the solid at the current point
  double getSolidThermalExpansion(double T) const;

protected:
  FuncConstPair<ScalarFunc> fheatcapacity; //!< Specific heat capacity for fluid
  FuncConstPair<ScalarFunc> sheatcapacity; //!< Specific heat capacity for solid
  FuncConstPair<ScalarFunc> fconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sexpansion;    //!< Thermal expansion
};

#endif
