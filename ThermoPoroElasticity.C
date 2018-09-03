// $Id$
//==============================================================================
//!
//! \file ThermoPoroElasticity.C
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for non-isothermal PoroElasticity problems
//!
//==============================================================================

#include "ThermoPoroElasticity.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "EqualOrderOperators.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "StabilizationUtils.h"


//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  Uc = 0,         // Previous displacement
  Pc = 1,         // Previous pore pressure
  Tc = 2,         // Previous temperature
  Uo = 3,         // Current displacement
  Po = 4,         // Current pore pressure
  To = 5,         // Current temperature
  NSOL = 6
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fu = 1,
  Fp = 2,
  FT = 3,
  Fprev = 5,
  Fnow = 6,
  Fc = 6
};


ThermoPoroElasticity::Mats::Mats(size_t ndof_displ, size_t ndof_press,
                                 bool neumann, char dynamic,
                                 int nbasis, int nsd) :
  PoroElasticity::Mats(3, nbasis)
{
  up_Q = 4;
  pu_Q = 6;
  uu_M = 10;
  up_D = 11;
  pp_P = 12;
  NMAT = 13;

  rhsOnly = neumann;
  withLHS = !neumann;

  this->resize(NMAT, 7);

  this->redim(1, ndof_displ, nsd, 1);
  this->redim(2, ndof_press, 1, nbasis);
  this->redim(3, ndof_press, 1, nbasis);
  this->redimOffDiag(4,0);
  this->redimOffDiag(5,0);
  this->redimOffDiag(6,0);
  this->redimOffDiag(7,0);
  this->redimOffDiag(8,0);
  this->redimOffDiag(9,0);
  this->redimNewtonMat();

  if (withLHS) {
    if (dynamic > 0)
      A[uu_M].resize(A[uu_K].rows(), A[uu_K].rows());
    if (dynamic == 2)
      A[up_D].resize(A[up_Q].rows(), A[up_Q].cols());
    A[pp_P].resize(A[pp_S].rows(), A[pp_S].cols());
  }

  h = lambda = 0.0;
}


const Matrix& ThermoPoroElasticity::Mats::getNewtonMatrix() const
{
  const_cast<Matrix&>(A[pu_Q]).addBlock(A[up_Q], 1.0, 1, 1, true);
  return this->BlockElmMats::getNewtonMatrix();
}


const Vector& ThermoPoroElasticity::Mats::getRHSVector() const
{
  const Vector& result = this->BlockElmMats::getRHSVector();
  Vector& F = const_cast<Vector&>(b[0]);

  return result;
}


ThermoPoroElasticity::ThermoPoroElasticity(unsigned short int n, int order, bool stab) :
  PoroElasticity(n, true), nsd(n), gacc(9.81), SUPG(stab)
{
  primsol.resize(1+order);
  tracFld = nullptr;
  fluxFld = nullptr;
  eS = 1;
}


LocalIntegral* ThermoPoroElasticity::getLocalIntegral(const std::vector<size_t>& nen,
                                                      size_t, bool neumann) const
{
  MixedElmMats* result = new MixedElmMats(nen[0], nen[1], neumann, 0, nsd);
  result->b[Fc].resize(result->b[0].size());
  return result;
}


bool ThermoPoroElasticity::initElement(const std::vector<int>& MNPC,
                                       const std::vector<size_t>& elem_sizes,
                                       const std::vector<size_t>& basis_sizes,
                                       LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],nsd*basis_sizes[0],basis_sizes[0])|
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol[1],elmInt.vec[Uo])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** ThermoPoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool ThermoPoroElasticity::initElementBou(const std::vector<int>& MNPC,
                                          const std::vector<size_t>& elem_sizes,
                                          const std::vector<size_t>& basis_sizes,
                                          LocalIntegral& elmInt)
{
  return this->IntegrandBase::initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool ThermoPoroElasticity::evalIntMx(LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const TimeDomain& time, const Vec3& X) const
{
//  if (!this->PoroElasticity::evalIntMx(elmInt, fe, time, X))
//    return false;

  Mats& elMat = static_cast<Mats&>(elmInt);
  ThermoPoroMaterial* mat = dynamic_cast<ThermoPoroMaterial*>(material);

  if (!mat)
  {
    std::cerr << __FUNCTION__ << ":No material data." << std::endl;
    return false;
  }

  size_t i,j,k;

  Matrix Bmat, Cmat, CB;

  if(!this->formBmatrix(Bmat,fe.grad(1)))
    return false;

  SymmTensor sigma(nsd);
  double U;
  if(!mat->evaluate(Cmat,sigma,U,fe,X,sigma,sigma,0))
    return false;

  Vec3 hydCond = mat->getPermeability(X);
  double rhof = mat->getFluidDensity(X);
  double Ko = mat->getBulkMedium(X);
  double Ks = mat->getBulkSolid(X);
  double Kf = mat->getBulkFluid(X);
  double poro = mat->getPorosity(X);
  // Biot's coefficient
  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kf);
  // m vector for 2D
  Vector m(3);
  m(1) = m(2) = 1.0; m(3) = 0.0;

  // Integration of Cuu
  CB.multiply(Cmat,Bmat,false,false);
  CB *= -1.0 * fe.detJxW;
  elMat.A[elMat.uu_K].multiply(Bmat,CB,true,false,true);

  // Integration of Cup
  Matrix Cuptmp;
  Cuptmp.outer_product(m, fe.basis(2), false, scl1*alpha*fe.detJxW);

  elMat.A[elMat.up_Q].multiply(Bmat,Cuptmp,true,false,true);

  double T = elMat.vec[Tc].dot(fe.basis(2));
  double Tn = elMat.vec[To].dot(fe.basis(2));

  double rhoc = mat->getHeatCapacity(T);
  double lambda = mat->getThermalConductivity(T);

  // Integration of CuT
  Matrix BC, CuTtmp;
  BC.multiply(Bmat,Cmat,true,false);
  double alpha_s = mat->getSolidThermalExpansion(T);
  CuTtmp.outer_product(m, fe.basis(2), false, scl2*alpha_s/3.0*fe.detJxW);
  elMat.A[elMat.uT_K].multiply(BC,CuTtmp,false,false,true);

  // Integration of Cpp
  Matrix Cpp;
  Cpp.resize(fe.basis(2).size(),fe.basis(2).size());
  EqualOrderOperators::Weak::Mass(Cpp, fe, scl1*scl1*Minv, 2);

  // Integration of Kpp
  Matrix Kpp;
  Kpp.resize(fe.basis(2).size(),fe.basis(2).size());
  EqualOrderOperators::Weak::LaplacianCoeff(Kpp, fe, hydCond, scl1*scl1/(rhof*gacc), 2);

  elMat.A[elMat.pp_S] += Cpp;
  elMat.A[elMat.pp_S].add(Kpp,time.dt);

  // Integration of CTT;
  Matrix CTT;
  CTT.resize(fe.basis(2).size(),fe.basis(2).size());
  EqualOrderOperators::Weak::Mass(CTT, fe, scl2*scl2*rhoc, 2);

  // Integration of KTT
  double rhofcf = rhof * mat->getFluidHeatCapacity(T);
  Vector gradP;
  fe.grad(2).multiply(elMat.vec[Pc],gradP,true);
  Vec3 vel;
  for (k = 1; k <= nsd; k++)
    vel[k-1] = -1.0*gradP(k)*hydCond[k-1]/(rhof*gacc);

  Matrix KTT;
  KTT.resize(fe.basis(2).size(),fe.basis(2).size());

  EqualOrderOperators::Weak::Advection(KTT, fe, vel, scl2*scl2*rhofcf, 2);
  EqualOrderOperators::Weak::Laplacian(KTT, fe, scl2*scl2*lambda, false, 2);

  elMat.A[elMat.TT_K] += CTT;
  elMat.A[elMat.TT_K].add(KTT,time.dt);

  if (SUPG)
  {
    /*
    size_t ru = elMat.A[elMat.uu_K].rows();
    size_t rp = elMat.A[elMat.pp_S].rows();

    Matrix KTTs, CTTs;
    KTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    CTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    // Evaluate the element size
    double h = StabilizationUtils::getElementSize(fe.XC,nsd);
    // Calculate the Peclet number
    double alpha_e = h*fabs(nsd*vel[0])/(2*lambda);
    // Evaluate the SUPG stabilization parameter
    double taue;
    if (alpha_e >= 1.0)
      taue = h/(2*fabs(nsd*vel[0])) * (1/tanh(alpha_e) - 1/alpha_e);
    else
      taue = h * h / 12 / lambda;
    for (i = 1; i <= fe.basis(2).size(); i++) {
      for (j = 1; j <= fe.basis(2).size(); j++) {
        double a = 0.0, b = 0.0, c = 0.0;
        for (k = 1; k <= nsd; k++) {
          a += fe.grad(2)(i,k)*vel[k-1]*vel[k-1]*fe.grad(2)(j,k);
          b += fe.grad(2)(i,k)*vel[k-1]*fe.hess(2)(i,k,k);
          c += fe.grad(2)(i,k)*vel[k-1]*fe.basis(2)(j);
        }
        KTTs(i,j) += taue * rhoc * (rhoc * a + lambda * b) * fe.detJxW;
        CTTs(i,j) += taue * rhoc * c * fe.detJxW;
      }
    }
    elMat.A[elMat.TT_K] += CTTs;
    elMat.A[elMat.TT_K].add(KTTs,time.dt);
    for (i = 1; i <= rp; i++) {
      for (j = 1; j <= rp; j++) {
        size_t li = ru+2*i;
        size_t lj = ru+2*j;
        elMat.A[elMat.Kprev](li,lj) += CTTs(i,j);
      }
    }
    */
  }

  return true;
}


bool ThermoPoroElasticity::evalBouMx(LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const TimeDomain& time,
                                     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << " *** ThermoPoroElasticity::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr << " *** ThermoPoroElasticity::evalBouMx: No load vector." << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 tr2 = this->getTraction(Xt,normal);
  Xt.t -= time.dt;
  Vec3 tr1 = this->getTraction(Xt,normal);
  Vec3 dtr;
  dtr = tr2 - tr1;

  // Integration of vector Fu
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Fu](nsd*(i-1)+j) += -1.0*dtr[j-1]*fe.basis(1)(i)*fe.detJxW;

  // Integration of vector Fp
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[Fp](i) += 0.0*scl1; // TO DO

  // Integration of vector FT
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[FT](i) += 0.0*scl2; // TO DO

  return true;
}


bool ThermoPoroElasticity::finalizeElement(LocalIntegral& elmInt,
                                           const TimeDomain&, size_t)
{
  Mats& elMat = static_cast<Mats&>(elmInt);

  // Residual
  elMat.A[elMat.uu_K].multiply(elMat.vec[Uc], elMat.b[Fu], -1.0, 1.0);
  elMat.A[elMat.up_Q].multiply(elMat.vec[Pc], elMat.b[Fu], -1.0, 1.0);
  elMat.A[elMat.uT_K].multiply(elMat.vec[Tc], elMat.b[Fu], -1.0, 1.0);

  elMat.A[elMat.up_Q].multiply(elMat.vec[Uc], elMat.b[Fp], -1.0, 1.0, true); // pu
  elMat.A[elMat.pp_S].multiply(elMat.vec[Pc], elMat.b[Fp], -1.0, 1.0);

  elMat.A[elMat.TT_K].multiply(elMat.vec[Tc], elMat.b[FT], -1.0, 1.0);

  // Time stepping
  elMat.A[elMat.uu_K].multiply(elMat.vec[Uo], elMat.b[Fu], 1.0, 1.0);
  elMat.A[elMat.uT_K].multiply(elMat.vec[To], elMat.b[Fu], 1.0, 1.0);
  elMat.A[elMat.up_Q].multiply(elMat.vec[Po], elMat.b[Fu], 1.0, 1.0);
  elMat.A[elMat.up_Q].multiply(elMat.vec[Uo], elMat.b[Fp], 1.0, 1.0, true); // pu
  elMat.A[elMat.TT_K].multiply(elMat.vec[To], elMat.b[FT], 1.0, 1.0);

  return true;
}


size_t ThermoPoroElasticity::getNoFields(int fld) const
{
  return fld < 2 ? nsd+2 : this->PoroElasticity::getNoFields(fld);
}


std::string ThermoPoroElasticity::getField1Name(size_t i, const char* prefix) const
{
  static const char* s[5] = { "u_x", "u_y", "u_z", "p^w", "T" };
  if (i == 11) {
    if (nsd == 3)
      return "u_x&&u_y&&u_z";
    else
      return "u_x&&u_y";
  }
  if (i == 12)
    return "p^w&&T";

  if (nsd == 2 && i >= 2)
    ++i;

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


ThermoPoroElasticity::WeakDirichlet::WeakDirichlet(unsigned short int n) :
  nsd(n), flux(nullptr)
{
  primsol.resize(2);
}


bool ThermoPoroElasticity::WeakDirichlet::initElementBou(const std::vector<int>& MNPC,
                                                         const std::vector<size_t>& elem_sizes,
                                                         const std::vector<size_t>& basis_sizes,
                                                         LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],
                         nsd*basis_sizes[0],basis_sizes[0])|
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],
                         nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol[1],elmInt.vec[Uo])
           + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],
                         nsd*basis_sizes[0],basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],
                         nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** ThermoPoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


LocalIntegral* ThermoPoroElasticity::WeakDirichlet::getLocalIntegral(const std::vector<size_t>& nen,
                                                                     size_t,
                                                                     bool neumann) const
{
  return new MixedElmMats(nen[0], nen[1], neumann, 0, nsd);
}


bool ThermoPoroElasticity::WeakDirichlet::evalBouMx(LocalIntegral& elmInt,
                                                    const MxFiniteElement& fe,
                                                    const TimeDomain& time,
                                                    const Vec3& X,
                                                    const Vec3& normal) const
{
  return false;
  /*

  Mats& elMat = static_cast<Mats&>(elmInt);

  // Evaluate the Neumann heat flux on the boundary
  double qT = 0.0;
  if (flux)
    qT = (*flux)(X);

  size_t ru = elMat.A[elMat.uu_K].rows();

  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    for (size_t j = 1; j <= fe.basis(2).size(); j++)
    {
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[elMat.Know](li,lj)  += time.dt*fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
      elMat.A[elMat.TT_K](i,j) += time.dt*fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
    }
    elMat.b[FT](i) += time.dt*(qT*fe.basis(2)(i) + lambdae*Te*fe.basis(2)(i))*fe.detJxW;
  }

  Vector PoTo, PcTc, prevSol, currSol;
  PoTo.resize(2*elMat.vec[Po].size());
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Po].size(); i++) {
    PoTo(2*i-1) = elMat.vec[Po](i);
    PoTo(2*i  ) = elMat.vec[To](i);
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  elMat.b[Fc] = elMat.vec[Uc];
  elMat.b[Fc].insert(elMat.b[Fc].end(),PcTc.begin(),PcTc.end());

  return true;
  */
}
