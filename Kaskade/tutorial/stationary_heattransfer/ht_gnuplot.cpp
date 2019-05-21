/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "linalg/direct.hh"
#include "linalg/trivialpreconditioner.hh"
//#include "linalg/partialDirectPreconditioner.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/iluprecond.hh"      // PrecondType::ILUT, PrecondType::ILUK, PrecondType::ARMS
#include "linalg/iccprecond.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG, Euclid
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare, createUnitCube
#include "io/vtk.hh"
#include "io/gnuplot.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"

using namespace Kaskade;
#include "ht.hh"

int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer tutorial program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int  refinements = getParameter(pt, "refinements", 5),
       order       =  getParameter(pt, "order", 2),
       verbosity   = getParameter(pt, "verbosity", 1);
  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "output level (verbosity)       : " << verbosity << std::endl;

  int  direct, onlyLowerTriangle = false;
    
  DirectType directType;
//  IterateType iterateType = IterateType::CG;
  MatrixProperties property = MatrixProperties::SYMMETRIC;
  PrecondType precondType = PrecondType::NONE;
  std::string empty;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);
    
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  directType = static_cast<DirectType>(getParameter(pt, s, 0));

//  s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
//  iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  precondType = static_cast<PrecondType>(getParameter(pt, s, 0));

  property = MatrixProperties::SYMMETRIC;

  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || ( (precondType == PrecondType::ICC) && !direct ) )
  {
    onlyLowerTriangle = true;
    std::cout << 
      "Note: direct solver MUMPS/PARADISO or PrecondType::ICC preconditioner ===> onlyLowerTriangle is set to true!" 
      << std::endl;
  }

  boost::timer::cpu_timer gridTimer;
//   two-dimensional space: dim=2
  constexpr int dim=2;        
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  // using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = HeatFunctional<double,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;

  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " triangles, " << std::endl;
  std::cout << "      " << gridManager.grid().size(1) << " edges, " << std::endl;
  std::cout << "      " << gridManager.grid().size(2) << " points" << std::endl;
  std::cout << "computing time for generation of initial mesh: " << (double)(gridTimer.elapsed().user)/1e9 << "s\n";

  
  // construction of finite element space for the scalar solution T.
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
    
  Spaces spaces(&temperatureSpace);
    
  // construct variable list.
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
        
  std::string varNames[1] = { "u" };
    
  VariableSet variableSet(spaces,varNames);

  // construct variational functional
    
  double kappa = 1.0;
  double q = 1.0;
  Functional F(kappa,q);
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  std::cout << std::endl << "no of variables = " << nvars << std::endl;
  std::cout << "no of equations = " << neq   << std::endl;
  size_t dofs = variableSet.degreesOfFreedom(0,nvars);
  std::cout << "number of degrees of freedom = " << dofs   << std::endl;

  
  //construct Galerkin representation
  
  Assembler assembler(gridManager,spaces);
  VariableSet::VariableSet u(variableSet);
  VariableSet::VariableSet du(variableSet);

  size_t nnz = assembler.nnz(0,neq,0,nvars,onlyLowerTriangle);
  std::cout << "number of nonzero elements in the stiffness matrix: " << nnz << std::endl << std::endl;
  boost::timer::cpu_timer assembTimer;
  
  CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,neq>::init(spaces));
  solution = 0;
  
  assembler.assemble(linearization(F,u));
  CoefficientVectors rhs(assembler.rhs());
  AssembledGalerkinOperator<Assembler,0,neq,0,nvars> A(assembler, onlyLowerTriangle);
  MatrixAsTriplet<double> tri = A.get<MatrixAsTriplet<double> >();
  std::cout << "computing time for assemble: " << (double)(assembTimer.elapsed().user)/1e9 << "s\n";

//     for (k=0; k< nnz; k++)
//       {
//         printf("%3d %3d %e\n", tri.ridx[k], tri.cidx[k], tri.data[k]);
//       }

  if (direct)
  {
    boost::timer::cpu_timer directTimer;
    directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
    u.data = solution.data;
    std::cout << "computing time for direct solve: " << (double)(directTimer.elapsed().user)/1e9 << "s\n";
  }
  else
  {
    boost::timer::cpu_timer iteTimer;
    Dune::InverseOperatorResult res;
    const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
    int iteSteps = getParameter(pt, "solver.iteMax", 2000);
    double iteEps = getParameter(pt, "solver.iteEps", 1.0e-10);
    StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
    int lookAhead;
    switch (precondType)
    {
      case PrecondType::NONE:
      case PrecondType::HB:   lookAhead=50; break;
      default:                lookAhead=3; break;
    }
    lookAhead = getParameter(pt, "solver.lookAhead", lookAhead);
    termination.setLookAhead(lookAhead);

    switch (precondType)
    {
      case PrecondType::NONE:
      {
        std::cout << "selected preconditioner: NONE" << std::endl;
        TrivialPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > trivial;
        CG<LinearSpace,LinearSpace> cg(A,trivial,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ADDITIVESCHWARZ:
      {
        std::cout << "selected preconditioner: ADDITIVESCHWARZ" << std::endl;
        std::pair<size_t,size_t> idx = temperatureSpace.mapper().globalIndexRange(gridManager.grid().leafIndexSet().geomTypes(dim)[0]);
        AdditiveSchwarzPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > addschwarz(A,idx.first,idx.second,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,addschwarz,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUT:
      {
        std::cout << "selected preconditioner: ILUT" << std::endl;
//        std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int lfil = getParameter(pt, "solver.ILUT.lfil", 140);
        double dropTol = getParameter(pt, "solver.ILUT.dropTol", 0.01);
        ILUTPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > ilut(A,lfil,dropTol,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,ilut,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::BiCGSTABSolver<LinearSpace> cg(A,ilut,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUK:
      {
        std::cout << "selected preconditioner: ILUK" << std::endl;
//        std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int fill_lev = getParameter(pt, "solver.ILUK.fill_lev", 3);
        ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,fill_lev,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ARMS:
      {
        int lfil = getParameter(pt, "solver.ARMS.lfil", 140);
        int lev_reord = getParameter(pt, "solver.ARMS.lev_reord", 1);
        double dropTol = getParameter(pt, "solver.ARMS.dropTol", 0.01);
        double tolind = getParameter(pt, "solver.ARMS.tolind", 0.2);
        ARMSPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,lfil,dropTol,lev_reord,tolind,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ICC:
      {
        std::cout << "selected preconditioner: ICC" << std::endl;
        if (property != MatrixProperties::SYMMETRIC) 
        {
          std::cout << "PrecondType::ICC preconditioner of TAUCS lib has to be used with matrix.property==MatrixProperties::SYMMETRIC\n";
          std::cout << "i.e., call the executable with option --solver.property MatrixProperties::SYMMETRIC\n\n";
        }
        double dropTol = getParameter(pt, "solver.ICC.dropTol", 0.01);;
        ICCPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > icc(A,dropTol);
        CG<LinearSpace,LinearSpace> cg(A,icc,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ICC0:
      {
        std::cout << "selected preconditioner: ICC0" << std::endl;
        ICC_0Preconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > icc0(A);
        CG<LinearSpace,LinearSpace> cg(A,icc0,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::HB:
      {
        std::cout << "selected preconditioner: HB" << std::endl;
        HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type, AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type > hb(gridManager.grid());
        CG<LinearSpace,LinearSpace> cg(A,hb,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::BOOMERAMG:
      {
        int steps = getParameter(pt, "solver.BOOMERAMG.steps", iteSteps);
        int coarsentype = getParameter(pt, "solver.BOOMERAMG.coarsentype", 21);
        int interpoltype = getParameter(pt, "solver.BOOMERAMG.interpoltype", 0);
        int cycleType = getParameter(pt, "solver.BOOMERAMG.cycleType", 1);
        int relaxType = getParameter(pt, "solver.BOOMERAMG.relaxType", 3);
        int variant = getParameter(pt, "solver.BOOMERAMG.variant", 0);
        int overlap = getParameter(pt, "solver.BOOMERAMG.overlap", 1);
        double tol = getParameter(pt, "solver.BOOMERAMG.tol", iteEps);
        double strongThreshold = getParameter(pt, "solver.BOOMERAMG.strongThreshold", (dim==2)?0.25:0.6);
        BoomerAMG<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> >
                  boomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
                  strongThreshold,variant,overlap,1,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,boomerAMGPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::LoopSolver<LinearSpace> cg(A,boomerAMGPrecon,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::EUCLID:
      {
        std::cout << "selected preconditioner: EUCLID" << std::endl;
        int level      = getParameter(pt, "solver.EUCLID.level",1);
        double droptol = getParameter(pt, "solver.EUCLID.droptol",0.01);
        int printlevel = 0;
        if (verbosity>2) printlevel=verbosity-2;
        printlevel = getParameter(pt,"solver.EUCLID.printlevel",printlevel);
        int bj = getParameter(pt, "solver.EUCLID.bj",0);
        Euclid<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > EuclidPrecon(A,level,droptol,printlevel,bj,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,EuclidPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
     }
      break;
      case PrecondType::JACOBI:
      default:
      {
        std::cout << "selected preconditioner: JACOBI" << std::endl;
        JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobi(A,1.0);
        CG<LinearSpace,LinearSpace> cg(A,jacobi,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
    }
    solution *= -1.0;
    u.data = solution.data;
    
    std::cout << "iterative solve eps= " << iteEps << ": " 
              << (res.converged?"converged":"failed") << " after "
              << res.iterations << " steps, rate="
              << res.conv_rate << ", computing time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
  }
  
  // compute L2 norm of the solution
  boost::timer::cpu_timer outputTimer;
  L2Norm l2Norm;
  std::cout << "L2norm(solution) = " << l2Norm(boost::fusion::at_c<0>(u.data)) << std::endl;

  // output of solution in VTK format for visualization,
  // the data are written as ascii stream into file temperature.vtu,
  // possible is also binary
  writeVTKFile(u,"temperature",IoOptions().setOrder(order));
  std::cout << "graphical output finished, data in VTK format is written into file temperature.vtu \n";
  IoOptions gnuplotOptions{};
  //    gnuplotOptions.info = IoOptions::none; // or IoOptions::summary or IoOptions::detail
  writeGnuplotFile(u,"temperature",gnuplotOptions);
  std::cout << "graphical output finished, Gnuplot data are written into file temperature.data \n";
  
  // output of solution for Amira visualization,
  // the data are written in binary format into file temperature.am,
  // possible is also ascii
  // IoOptions options;
  // options.outputType = IoOptions::ascii;
  // LeafView leafGridView = gridManager.grid().leafGridView();
  // writeAMIRAFile(leafGridView,variableSet,u,"temperature",options);

  std::cout << "computing time for output: " << (double)(outputTimer.elapsed().user)/1e9 << "s\n";

  std::cout << "total computing time: " << (double)(totalTimer.elapsed().user)/1e9 << "s\n";
  std::cout << "End heat transfer tutorial program" << std::endl;
}
