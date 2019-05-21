/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                       */
/*  This file is part of the library KASKADE 7                 */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html         */
/*                                       */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin              */
/*                                       */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.  */
/*    see $KASKADE/academic.txt                        */
/*                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/embedded_errorest.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "linalg/direct.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/iluprecond.hh"      // PrecondType::ILUT, PrecondType::ILUK, PrecondType::ARMS
#include "linalg/iccprecond.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare, createUnitCube
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"

//#include "cubus.hh"
using namespace Kaskade;
#include "peaksource.hh"

#ifndef SPACEDIM
#define SPACEDIM 2
#endif

#if SPACEDIM==2
#define DEFAULT_REFINEMENTS 5
#else
#define DEFAULT_REFINEMENTS 2
#endif

int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer tutorial program using embedded error estimation" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int  refinements = getParameter(pt, "refinement", DEFAULT_REFINEMENTS),
       order       = getParameter(pt, "order", 2),
       verbosity   = getParameter(pt, "verbosity", 1);

  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order         : " << order << std::endl;
  std::cout << "output level (verbosity)     : " << verbosity << std::endl;

  DirectType directType;
//  IterateType iterateType = IterateType::CG;
  MatrixProperties property;
  PrecondType precondType = PrecondType::NONE;
  std::string empty;
  
  int direct, onlyLowerTriangle = false;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);

  // the user may select a value for solver.direct of 
  // the enumeration class {UMFPACK, PARDISO, MUMPS, SUPERLU, UMFPACK3264, UMFPACK64}
  // Remark: DirectType::PARDISO not yet available
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  directType = static_cast<DirectType>(getParameter(pt, s, 2));

  // the user may select a value for solver.iterate of 
  // the enumeration {CG, BICGSTAB, GMRES, PCG, APCG, SGS}
  // Remark: in this example only IterateType::CG is used.
//  s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
//  iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
   
  // the user may select a value for solver.preconditioner of 
  // the enumeration class {NONE, JACOBI, ILUT, ILUK, ARMS, ADDITIVESCHWARZ,
  //                  BOOMERAMG, EUCLID, SSOR, ICC0, ICC, ILUKS}
  // Remark: in this example only PrecondType::NONE,PrecondType::JACOBI, PrecondType::ICC, PrecondType::ICC0, PrecondType::BOOMERAMG are used.
  s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  precondType = static_cast<PrecondType>(getParameter(pt, s, 0));

  property = MatrixProperties::SYMMETRIC;
  std::cout << "discretization is symmetric" << std::endl;
  
  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || ( (precondType == PrecondType::ICC) && !direct ) )
  {
    onlyLowerTriangle = true;
    std::cout << 
      "Note: direct solver MUMPS/PARADISO or PrecondType::ICC preconditioner ===> onlyLowerTriangle is set to true!" 
      << std::endl;
  }

#if SPACEDIM==2
  //   two-dimensional space: dim=2
  constexpr int dim=2;        
  using Grid = Dune::UGGrid<dim>;
  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " triangles, " << std::endl;
#else
  //  three-dimensional space: dim=3
  constexpr int dim=3; 
  using Grid = Dune::UGGrid<dim>;
  GridManager<Grid> gridManager( createUnitCube<Grid>(0.5) );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " tetrahedra, " << std::endl;
  std::cout << "      " << gridManager.grid().size(1) << " triangles, " << std::endl;
#endif
  std::cout << "      " << gridManager.grid().size(dim-1) << " edges, " << std::endl;
  std::cout << "      " << gridManager.grid().size(dim) << " points" << std::endl;

  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  // alternative: using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,1,0> >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = PeaksourceFunctional<double,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;

  gridManager.setVerbosity(verbosity);
  gridManager.enforceConcurrentReads(true);
  
  // construction of finite element space for the scalar solution T
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),
               order);
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  std::string varNames[1] = { "T" };
  VariableSet variableSet(spaces,varNames);

  Functional F;
  
    //construct Galerkin representation

  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  std::cout << "no of variables = " << nvars << std::endl;
  std::cout << "no of equations = " << neq   << std::endl;
  
  Assembler assembler(gridManager,spaces);
  VariableSet::VariableSet xx(variableSet);
  
  size_t nnz  = assembler.nnz(0,neq,0,nvars,onlyLowerTriangle);
  size_t size = variableSet.degreesOfFreedom(0,nvars);
  if ( verbosity>0) std::cout << "init mesh: nnz = " << nnz << ", dof = " << size << std::endl;
  
  std::vector<std::pair<double,double> > tol(1);
  double atol = getParameter(pt, "solver.atol", 1.0e-5);
  double rtol = getParameter(pt, "solver.rtol", 1.0e-5);      
  tol[0] = std::make_pair(atol,rtol); 
  std::cout << std::endl << "Accuracy: atol = " << atol << ",  rtol = " << rtol << std::endl;


  bool accurate = true;
  int refSteps = -1;
  int iter=0;

  do {
    refSteps++;

    boost::timer::cpu_timer assembTimer;
    VariableSet::VariableSet x(variableSet);
    assembler.assemble(linearization(F,x));
    CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,neq>::init(spaces));
    solution = 0;
    CoefficientVectors rhs(assembler.rhs());
    AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, onlyLowerTriangle);
    MatrixAsTriplet<double> tri = A.get<MatrixAsTriplet<double> >();
    if ( verbosity>1) std::cout << "assemble: " << (double)assembTimer.elapsed().user/1e9 << "s\n";


    if (direct) {
      boost::timer::cpu_timer directTimer;
      directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
      x.data = solution.data;

      if ( verbosity>1) std::cout << "direct solve: " << (double)(directTimer.elapsed().user)/1e9 << "s\n";
    }
    else {
      //if ( verbosity>0) std::cout << "iterative solver: steps = " << iteSteps << ", eps = " << iteEps << std::endl;
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
          TrivialPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > trivial;
          CG<LinearSpace,LinearSpace> cg(A,trivial,defaultScalarProduct,termination,verbosity);
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
          ICCPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > icc(A,dropTol);
          CG<LinearSpace,LinearSpace> cg(A,icc,defaultScalarProduct,termination,verbosity);
          cg.apply(solution,rhs,res);
        }
        break;
        case PrecondType::ICC0:
        {
          std::cout << "selected preconditioner: ICC0" << std::endl;
          ICC_0Preconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > icc0(A);
          CG<LinearSpace,LinearSpace> cg(A,icc0,defaultScalarProduct,termination,verbosity);
          cg.apply(solution,rhs,res);
        }
        break;
        case PrecondType::HB:
        {
          std::cout << "selected preconditioner: HB" << std::endl;
          HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type, AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type > hb(gridManager.grid());
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
          BoomerAMG<AssembledGalerkinOperator<Assembler,0,1,0,1> >
          BoomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
          strongThreshold,variant,overlap,1,verbosity);
          CG<LinearSpace,LinearSpace> cg(A,BoomerAMGPrecon,defaultScalarProduct,termination,verbosity);
//          Dune::LoopSolver<LinearSpace> cg(A,BoomerAMGPrecon,iteEps,iteSteps,verbosity);
          cg.apply(solution,rhs,res);
        }
        break;
        case PrecondType::JACOBI:
        default:
        {
          JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > jacobi(A,1.0);
          CG<LinearSpace,LinearSpace> cg(A,jacobi,defaultScalarProduct,termination,verbosity);
          cg.apply(solution,rhs,res);
        }
        break;
      }
      solution *= -1.0;
      x.data = solution.data;
  
      if ( verbosity>0) std::cout << "iterative solve eps= " << iteEps << ": " 
          << (res.converged?"converged":"failed") << " after "
          << res.iterations << " steps, rate="
          << res.conv_rate << ", time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
    }
    
	// graphical output of solution
    std::ostringstream fn;
    fn << "graph/peak-grid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();

    // output of solution in VTK format for visualization,
    // the data are written as ascii stream into file temperature.vtu,
    // possible is also binary
    writeVTKFile(x,fn.str(),IoOptions().setOrder(order));
  
  // output of solution for Amira visualization,
  // the data are written in binary format into file temperature.am,
  // possible is also ascii
  //    IoOptions options;
  //    options.outputType = IoOptions::ascii;
  //    LeafView leafGridView = gridManager.grid().leafGridView();
  //    writeAMIRAFile(leafGridView,variableSet,x,fn.str(),options);



    VariableSet::VariableSet e = x;
    projectHierarchically(variableSet,e);
    e -= x;    
  
    accurate = embeddedErrorEstimator(variableSet,e,x,IdentityScaling(),tol,gridManager,verbosity);
    nnz = assembler.nnz(0,1,0,1,onlyLowerTriangle);;
    size_t size = variableSet.degreesOfFreedom(0,1);
    if ( verbosity>0) std::cout << "new mesh: nnz = " << nnz << ", dof = " << size << std::endl;
    
    //ridx.resize(nnz);
    //cidx.resize(nnz);
    //data.resize(nnz);
    //rhs.resize(size);
    //solution.resize(size);

    // VariableSet::VariableSet xx may be used beyond the do...while loop	
    xx.data = x.data;
    iter++; 
    if (iter>9) 
    {
      std::cout << "*** Maximum number of iterations exceeded ***" << std::endl;
      break;
    }
    
  }  while (!accurate); 
  
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End heat transfer (peak source) tutorial program" << std::endl;
}
