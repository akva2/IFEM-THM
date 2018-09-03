// $Id$
//==============================================================================
//!
//! \file ThermoPoroElasticity.h
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for non-isothermal PoroElasticity problems
//!
//==============================================================================

#ifndef _THERMO_POROELASTICITY_H_
#define _THERMO_POROELASTICITY_H_

#include "Vec3.h"
#include "ElmMats.h"
#include "IntegrandBase.h"
#include "ThermoPoroMaterial.h"
#include "PoroElasticity.h"


/*!
 * \brief Class representing the integrand of a non-isothermal and fully
 	  coupled PoroElasticity problem
*/

class ThermoPoroElasticity : public PoroElasticity
{
  //! \brief Superclass for the ThermoPoroElasticity element matrices.
  //! \brief Default constructor.
  //! \param[in] ndof_displ Number of dofs in displacement
  //! \param[in] ndof_press Number of dofs in pressure
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  //! \param[in] dynamic Option for dynamic analysis
  //! (0: static analysis, 1: allocate M, 2: allocate M and D)
  //! \param[in] nsd Number of spatial dimensions
  class Mats : public PoroElasticity::Mats
  {
  public:
    //! \brief Default constructor.
    //! \param[in] ndof_displ Number of dofs in displacement
    //! \param[in] ndof_press Number of dofs in pressure
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    //! \param[in] dynamic Option for dynamic analysis
    //! (0: static analysis, 1: allocate M, 2: allocate M and D)
    //! \param[in] nsd Number of spatial dimensions
    Mats(size_t ndof_displ, size_t ndof_press, bool neumann, char dynamic, int nbasis, int nsd);

    //! \brief Empty destructor.
    virtual ~Mats() = default;

    //! \brief Returns the element level Newton matrix.
    const Matrix& getNewtonMatrix() const override;
    //! \brief Returns the element level right-hand-side vector.
    const Vector& getRHSVector() const override;

    size_t TT_K = 3; //!< Thermal stiffness matrix
    size_t uT_K = 5; //!< Displacement-temperature coupling
    size_t pT_K = 7; //!< Pressure-temperature coupling
    size_t TU_K = 8; //!< Temperature-displacement coupling
    size_t Tp_K = 9; //!< Temperature-pressure coupling
    size_t TT_M = 13; //!< Thermal mass matrix
  };

  //! \brief Class representing the element matrices for mixed formulation.
  class MixedElmMats : public Mats
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    MixedElmMats(size_t ndof_d, size_t ndof_p, bool neumann, char dyn, int nsd)
      : Mats(ndof_d, ndof_p, neumann, dyn, 2, nsd) {}
    //! \brief Empty destructor.
    virtual ~MixedElmMats() {}
  };

  //! \brief Class representing the element matrices for standard formulation.
  class StdElmMats : public Mats
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    StdElmMats(size_t ndof_d, size_t ndof_p, bool neumann, char dyn, int nsd)
      : Mats(ndof_d, ndof_p, neumann, dyn, 1, nsd) {}
    //! \brief Empty destructor.
    virtual ~StdElmMats() {}
  };

public:
  /*!
   * \brief Class for integrating Robin boundary conditions
  */
  class WeakDirichlet : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    WeakDirichlet(unsigned short int n);

    //! \brief Empty destructor.
    virtual ~WeakDirichlet() = default;

    //! \brief Defines the flux function
    void setFlux(RealFunc* f) { flux = f; }
    //! \brief Defines the temperature of the environment.
    void setEnvtTemperature(double envT) { Te = envT; }
    //! \brief Defines the thermal conductivity of the environment.
    void setEnvtConductivity(double envCond) { lambdae = envCond; }

    //! \brief Returns that this integrand has no interior contributions.
    bool hasInteriorTerms() const override { return false; }

    //! \brief Returns a local integral container for the given element.
    //! \param[in] nen1 Number of nodes on element for basis 1
    //! \param[in] nen2 Number of nodes on element for basis 2
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                    size_t, bool neumann) const override;

    //! \brief Initializes current element for numerical boundary integration (mixed)
    //! \param[in] MNPC1 Nodal point correspondence for basis 1
    //! \param[in] MNPC2 Nodal point correspondence for basis 2
    //! \param elmInt The local integral object for current element
    bool initElementBou(const std::vector<int>& MNPC,
                        const std::vector<size_t>& basis_sizes,
                        const std::vector<size_t>& elem_sizes,
                        LocalIntegral& elmInt) override;

    //! \brief Evaluates the integrands at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] time Parameters for nonlineat and time-dependent simulations
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at integration point
    bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                   const TimeDomain& time, const Vec3& X,
                   const Vec3& normal) const override;

  protected:
    unsigned short int nsd;       //!< Number of space dimensions
    RealFunc* flux;               //!< Flux function
    double Te;                    //!< Temperature of environment
    double lambdae;               //!< Thermal conductivity of environment
  };

	//! \brief The default constructor initializes all pointers to zero.
	//! \param[in] n Number of spatial dimensions
	//! \param[in] order The order of the time integration
	ThermoPoroElasticity(unsigned short int n = 2, int order = 1, bool stab = false);

	//! \brief Empty destructor. 
  virtual ~ThermoPoroElasticity() = default;

  //! \brief Defines the scaling values between U and P and U and T
  void setScalingValues(double scl_u_p, double scl_u_T)
    { scl1 = scl_u_p; scl2 = scl_u_T; }

  int getIntegrandType() const override
  { return SUPG ? ELEMENT_CORNERS | SECOND_DERIVATIVES : STANDARD; }

  //! \brief Returns a local integral container for the given element
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                  size_t, bool neumann) const override;

  //! \brief Initializes current element for numerical integration
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param elmInt The local integral object for current element
  bool initElement(const std::vector<int>& MNPC,
                   const std::vector<size_t>& basis_sizes,
                   const std::vector<size_t>& elem_sizes,
                   LocalIntegral& elmInt) override;

  //! \brief Initializes current element for numerical boundary integration (mixed)
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param elmInt The local integral object for current element
  bool initElementBou(const std::vector<int>& MNPC,
                      const std::vector<size_t>& basis_sizes,
                      const std::vector<size_t>& elem_sizes,
                      LocalIntegral& elmInt) override;

  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override;

  //! \brief Evaluates the integrand at a boundary point
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X,
                 const Vec3& normal) const override;

  //! \brief Finalizes the element quantities after the numerical integration
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t) override;

  //! \brief Returns the number of primary/secondary solution field components
  //! \param[in] fld Which field set to consider (1=primary,2=secondary)
  size_t getNoFields(int fld = 1) const override;

  //! \brief Returns the name of a primary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField1Name(size_t i, const char* prefix = 0) const override;

protected:
	unsigned short int nsd;    //!< Number of space dimensions (2)
  unsigned short int eS;     //!< Index to element load vector
  double gacc;               //!< Gravitational acceleration
  double scl1;               //!< Scaling value for displacement and pressure
  double scl2;               //!< Scaling value for displacement and temperature
  bool SUPG;                 //!< \e true to use SUPG stabilization for temperature
};

#endif  // _THERMO_POROELASTICITY_H_
