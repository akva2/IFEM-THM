// $Id$
//==============================================================================
//!
//! \file main_THM.C
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Main program for the isogeometrics solver for thermo hydro mechanics
//!
//==============================================================================

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMThermoHydroMechanics.h"
#include "SIMargsBase.h"
#include "SIMSolver.h"
#include "ASMmxBase.h"
#include "Utilities.h"
#include "Profiler.h"
#include "VTF.h"

  template<class Dim>
int runSimulator(char* infile, bool supg)
{
  SIMThermoHydroMechanics<Dim> model(supg);
  SIMSolver<SIMThermoHydroMechanics<Dim>> solver(model);

  // Read input file
  if(!model.read(infile) || !solver.read(infile))
    return 1;

  // Configure finite element library
  if(!model.preprocess())
    return 2;

  // Setup integration
  model.setQuadratureRule(model.opt.nGauss[0],true);
  model.initSystem(model.opt.solver,1,1,false);
  model.init(solver.getTimePrm());
  model.setInitialConditions();
  model.setAssociatedRHS(0,0);
  model.setMode(SIM::DYNAMIC);

  if (model.opt.dumpHDF5(infile))
    solver.handleDataOutput(model.opt.hdf5, model.opt.saveInc,
                            model.opt.restartInc);

  if (!solver.solveProblem(infile,"Solving thermo-hydro-mechanics problem"))
    return 5;

  return 0;
}


/*!
  \brief Main program for the NURBS-based isogeometric THM solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -2D : Use two-parametric simulation driver
  \arg -supg : Use stream-line upwind stabilization. Suitable for quadratic or higher elements.
  \arg -mixed: Use a mixed (P_N/P_N-1) solver.
*/

int main(int argc, char ** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  std::vector<int> ignoredPatches;
  int i;
  char* infile = 0;
  bool supg = false;
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  SIMargsBase args("thermoporoelasticity");

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
  {
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass and args handled by args class
    else if (dummy.parseOldOptions(argc,argv,i))
      ; // Ignore the obsolete option
    else if (!strcasecmp(argv[i],"-supg"))
      supg = true;
    else if (!infile) {
      infile = argv[i];
      if (strcasestr(infile, ".xinp")) {
        if (!args.readXML(infile,false))
          return 1;
        i = 0;
      }
    } else
      std::cerr << "*** Unknown option ignored: " << argv[i] << std::endl;
  }

  if (!infile)
  {
    IFEM::cout << "Usage: " << argv[0]
      << " <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
      << " [-free] [-lag|-spec|-LR] [-1D|-2D] [-mixed] [-nGauss <n>]"
      << "\n       [-vtf <format> [-nviz <nviz>]"
      << " [-nu <nu> [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
      << "       [-eig <iop> [-nev <nev>] [-ncv <ncv>] [-shift <shf>]]\n"
      << "       [-ignore <p1> <p2> ...] [-fixDup]" << std::endl;
    return 0;
  }

  IFEM::cout << "\n >>> IFEM ThermoHydroMechanics Solver <<<"
    << "\n ======================================"
    << "\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout << " " << argv[i];
  IFEM::cout << "\n\n Input file: " << infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (args.dim == 3)
    return runSimulator<SIM3D>(infile,supg);
  else if (args.dim == 2)
    return runSimulator<SIM2D>(infile,supg);
  else
    return runSimulator<SIM1D>(infile,supg);

  return 1;
}
