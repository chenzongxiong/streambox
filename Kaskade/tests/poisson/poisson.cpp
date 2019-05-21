/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

 /**
  * \defgroup tests Tests
  * \brief Classes and functions that test whether certain parts of Kaskade work as they should.
  */
 /**
 * @file
 * @ingroup tests
 * @brief  Test with simple poisson equation.
 * 
 * Testprogram for Kaskade: tests wether the error between calculated and exact solution behaves as expected. 
 * This test was built by changing the laplace example, such that the exact solution is known. 
 * Boundary conditions and f are adapted in poisson.hh. 
 * 
 * The test does not work with the following preconditioners: PrecondType::NONE, PrecondType::ARMS, PrecondType::ADDITIVESCHWARZ, PrecondType::HB. 
 * The test works with: PrecondType::JACOBI,  PrecondType::ILUT, PrecondType::ILUK, PrecondType::BOOMERAMG, PrecondType::EUCLID, PrecondType::ICC0, PrecondType::ICC and with direct solvers DirectType::MUMPS, DirectType::UMFPACK and DirectType::UMFPACK3264. 
 * Uses ContinuousLagrangeMapper. 
 */
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
#include "linalg/partialDirectPreconditioner.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/iluprecond.hh"      // PrecondType::ILUT, PrecondType::ILUK, PrecondType::ARMS
#include "linalg/iccprecond.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG, Euclid
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "utilities/enums.hh"
#include "utilities/kaskopt.hh"


using namespace Kaskade;
#include "poisson.hh"
#undef max


/**
 * @brief  This represents the known exact solution. 
 * It is \f$ u = \cos(3 \pi x)+\sin(4 \pi y) \f$.
 */
struct PoissonSolution
{
   using Scalar = double;
   static int const components = 1;
   using ValueType = Dune::FieldVector<Scalar,components>;

   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   ValueType value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
     Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
     double PI = 3.141592653589793238462;
     return std::cos(3*PI*x[0]) + std::sin(4*PI*x[1]);
   }
};


int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  boost::timer::cpu_timer totalTimer;
  std::cout << "Start test with poisson equation" << std::endl;
  
  int verbosityOpt = 0;
  bool dump = false; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);
  
  //true omits test with ansatz functions of order 5, so test will be faster
  bool quick = getParameter(pt,"quick",0);

  constexpr int dim=2;    
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = PoissonFunctional<double,VariableSet>;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using AssGalOp = AssembledGalerkinOperator<Assembler,0,neq,0,nvars>;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  
  std::cout << "dimension of space:             " << dim << std::endl;
  int refinements = 1;
  constexpr int maxRefSteps = 5;
  
  constexpr int orderLength = 4;
  //ansatz orders to be tested
  int order[orderLength] = {1,2,3,5};
  //error values for corresponding ansatz order and refinement
  double errors[orderLength][maxRefSteps] = {{0.226,0.075,0.021,5.3e-3,1.4e-3}, //order 1
               {0.15,0.007,4.6e-4,3.15e-5,2.24e-6}, // order 2
               {0.013,8.1e-4,4.96e-5,3.08e-6,1.92e-7}, //order 3
               {6.1e-4,7.34e-6,1.16e-7,1.82e-9,2.90e-11}}; //order 5
               
//   original values for L2-Norm of discretization error             
//   double errors[orderLength][maxRefSteps] = {{0.35,0.11,0.03,7.8e-3,2e-3}, //order 1
//               {0.22,0.01,6.9e-4,4.7e-5,3.3e-6}, // order 2
//               {0.02,1.2e-3,7.5e-5,4.7e-6,2.9e-7}, //order 3
//               {9.1e-4,1.1e-5,1.8e-7,2.7e-9,4.2e-11}}; //order 5
  // values for order of convergence for corresponding ansatz order
  double conOrder[orderLength] = {1.7,3.8,3.8,5.7};
  
  bool valid = true;
  std::stringstream message("Test succeeded", std::stringstream::out);
  
  int verbosity   = getParameter(pt, "verbosity", 1);
  // if true, then the test result will be written in a file
  bool result = getParameter(pt, "result",0);
  
  int  direct, onlyLowerTriangle = false;
    
  DirectType directType;
  MatrixProperties property = MatrixProperties::SYMMETRIC;
  PrecondType precondType = PrecondType::NONE;
  std::string empty;

  std::string s("names.type.");
  std::string solverTypeI = getParameter(pt, "solver.type", empty);
  std::string solverTypeII;
  s += solverTypeI;
  direct = getParameter(pt, s, 0);
    
  std::string directSolver = getParameter(pt, "solver.direct", empty);
  s = "names.direct." + directSolver;
  directType = static_cast<DirectType>(getParameter(pt, s, 0));

  std::string preconditioner = getParameter(pt, "solver.preconditioner", empty);
  s = "names.preconditioner." + preconditioner;
  precondType = static_cast<PrecondType>(getParameter(pt, s, 0));

  int blocks = getParameter(pt,"blocks",40);
  int nthreads = getParameter(pt,"threads",4);
  double rowBlockFactor = getParameter(pt,"rowBlockFactor",2.0);

  property = MatrixProperties::SYMMETRIC;

  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || (precondType == PrecondType::ICC) )
  {
    onlyLowerTriangle = true;
    if(verbosity > 0) {
      std::cout << 
  "Note: direct solver MUMPS/PARADISO or PrecondType::ICC preconditioner ===> onlyLowerTriangle is set to true!" 
  << std::endl;
    }
  }
  
  //iterate over the ansatz orders
  for(int i = 0;i < (orderLength-quick) && valid;i++) {
    if(verbosity > 0) {
      std::cout << "original mesh shall be refined: " << maxRefSteps + refinements << " times" << std::endl;
      std::cout << "discretization order          : " << order[i] << std::endl;
    }

    boost::timer::cpu_timer gridTimer;
    Dune::GridFactory<Grid> factory;

    // vertex coordinates v[0], v[1]
    Dune::FieldVector<double,dim> v;  
    v[0]=0; v[1]=0; factory.insertVertex(v);
    v[0]=1; v[1]=0; factory.insertVertex(v);
    v[0]=1; v[1]=1; factory.insertVertex(v);
    v[0]=0; v[1]=1; factory.insertVertex(v);
    // triangle defined by 3 vertex indices
    std::vector<unsigned int> vid(3);
    Dune::GeometryType gt(Dune::GeometryType::simplex,2);
    vid[0]=0; vid[1]=1; vid[2]=2; factory.insertElement(gt,vid);
    vid[0]=0; vid[1]=2; vid[2]=3; factory.insertElement(gt,vid);
    std::unique_ptr<Grid> grid( factory.createGrid() ) ;
    // the coarse grid will be refined once
    grid->globalRefine(refinements);
    if(verbosity > 0) {
      // some information on the refined mesh
      std::cout << std::endl << "Grid is read successfully" << std::endl;
      std::cout << "Grid: " << grid->size(0) << " triangles, " << std::endl;
      std::cout << "      " << grid->size(1) << " edges, " << std::endl;
      std::cout << "      " << grid->size(2) << " points" << std::endl;
      std::cout << "computing time for grid construction: " << (double)(gridTimer.elapsed().user)/1e9 << "s\n";
    }

    // a gridmanager is constructed 
    // as connector between geometric and algebraic information
    GridManager<Grid> gridManager(std::move(grid));

    
    // construction of finite element space for the scalar solution u.
    H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order[i]);
    
    Spaces spaces(&temperatureSpace);
    
    // construct variable list.
    // VariableDescription<int spaceId, int components, int Id>
    // spaceId: number of associated FEFunctionSpace
    // components: number of components in this variable
    // Id: number of this variable
    
    std::string varNames[1] = { "u" };
    
    VariableSet variableSet(spaces,varNames);

    // construct variational functional
    Functional F;
    if(verbosity > 0) {
      std::cout << std::endl << "no of variables = " << nvars << std::endl;
      std::cout << "no of equations = " << neq   << std::endl;
      size_t dofs = variableSet.degreesOfFreedom(0,nvars);
      std::cout << "number of degrees of freedom = " << dofs << std::endl;
    }
    
    //construct Galerkin representation
    Assembler assembler(gridManager,spaces);
    
    gridManager.enforceConcurrentReads(true);
    assembler.setNSimultaneousBlocks(blocks);
    assembler.setRowBlockFactor(rowBlockFactor);
    
    //keep mesh sizes
    double hs[maxRefSteps];
    //keep errors
    double ress[maxRefSteps];
    
    //calculate solution for different refinements
    int refSteps = 0;
    do
    {
      gridManager.globalRefine(1);
      boost::timer::cpu_timer assembTimer;
      VariableSet::VariableSet u(variableSet);

      size_t nnz = assembler.nnz(0,1,0,1,onlyLowerTriangle);
      if(verbosity > 0) {
        std::cout << "number of nonzero elements in the stiffness matrix: " << nnz << std::endl << std::endl;
      }
      
      CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,1>::init(spaces));
      solution = 0;
      
      assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
      CoefficientVectors rhs(assembler.rhs());
      AssGalOp A(assembler, onlyLowerTriangle);
      if(verbosity > 0) {
        std::cout << "assemble finished \n";
        std::cout << "computing time for assemble: " << (double)(assembTimer.elapsed().user)/1e9 << "s\n";
      }
     
      if (direct)
      {
        solverTypeII = directSolver;
        boost::timer::cpu_timer directTimer;
        directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
        u.data = solution.data;
        if(verbosity > 0) 
        {
          std::cout << "computing time for direct solve: " << boost::timer::format(directTimer.elapsed()) << "\n";
        }
      }
      else
      {
        solverTypeII = "CG, " + preconditioner;
        if(!quick) {
          std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl 
                    << "ATTENTION: Test will take a long time. Use --quick 1 for a fast test." << std::endl
                    << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
        boost::timer::cpu_timer iteTimer;
        int iteSteps = 1000;
        double iteEps = 1.0e-20;
        Dune::InverseOperatorResult res;
        const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
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
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl << "ATTENTION: Test does not work without preconditioner."
                      << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: NONE" << std::endl;
            }
            TrivialPreconditioner<AssGalOp > trivial;
            CG<LinearSpace,LinearSpace> cg(A,trivial,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ADDITIVESCHWARZ:
          {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl 
              << "ATTENTION: Test does not work with preconditioner PrecondType::ADDITIVESCHWARZ." << std::endl 
              << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ADDITIVESCHWARZ" << std::endl;
            }
            std::pair<size_t,size_t> idx = 
              temperatureSpace.mapper().globalIndexRange(gridManager.grid().leafIndexSet().geomTypes(dim)[0]);
            AdditiveSchwarzPreconditioner<AssGalOp > 
              addschwarz(A,idx.first,idx.second,verbosity);
            CG<LinearSpace,LinearSpace> cg(A,addschwarz,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ILUT:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ILUT" << std::endl;
              std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
              std::cout << "needs matrix.property = GENERAL" << std::endl;
            }
            int lfil = getParameter(pt, "solver.ILUT.lfil", 140);
            double dropTol = getParameter(pt, "solver.ILUT.dropTol", 0.01);
            ILUTPreconditioner<AssGalOp > ilut(A,lfil,dropTol,verbosity);
            Dune::BiCGSTABSolver<LinearSpace> cg(A,ilut,iteEps,iteSteps,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ILUK:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ILUK" << std::endl;
              std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
              std::cout << "needs matrix.property = GENERAL" << std::endl;
            }
            int fill_lev = getParameter(pt, "solver.ILUK.fill_lev", 3);
            ILUKPreconditioner<AssGalOp > iluk(A,fill_lev,verbosity);
            Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ICC:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ICC" << std::endl;
            }
            if (property != MatrixProperties::SYMMETRIC) 
            {
              std::cout << "PrecondType::ICC preconditioner of TAUCS lib has to be used with matrix.property==MatrixProperties::SYMMETRIC\n";
              std::cout << "i.e., call the executable with option --solver.property MatrixProperties::SYMMETRIC\n\n";
            }
            double dropTol = getParameter(pt, "solver.ICC.dropTol", 0.01);;
            ICCPreconditioner<AssGalOp > icc(A,dropTol);
            CG<LinearSpace,LinearSpace> cg(A,icc,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ICC0:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ICC0" << std::endl;
            }
            ICC_0Preconditioner<AssGalOp > icc0(A);
            CG<LinearSpace,LinearSpace> cg(A,icc0,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::HB:
          {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl << 
              "ATTENTION: Test does not work with preconditioner PrecondType::HB." << std::endl  << 
             "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: HB" << std::endl;
            }
            HierarchicalBasisPreconditioner<Grid,AssGalOp::range_type, 
              AssGalOp::range_type > hb(gridManager.grid());
            CG<LinearSpace,LinearSpace> cg(A,hb,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ARMS:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ARMS" << std::endl;
            }
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl << 
                         "ATTENTION: Test does not work with preconditioner PrecondType::ARMS." << std::endl <<
                         "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            int lfil = getParameter(pt, "solver.ARMS.lfil", 140);
            int lev_reord = getParameter(pt, "solver.ARMS.lev_reord", 1);
            double dropTol = getParameter(pt, "solver.ARMS.dropTol", 0.01);
            double tolind = getParameter(pt, "solver.ARMS.tolind", 0.2);
            ARMSPreconditioner<AssGalOp > iluk(A,lfil,dropTol,lev_reord,tolind,verbosity);
            CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::BOOMERAMG:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: BOOMERAMG" << std::endl;
            }
            int steps = getParameter(pt, "solver.BOOMERAMG.steps", iteSteps);
            int coarsentype = getParameter(pt, "solver.BOOMERAMG.coarsentype", 21);
            int interpoltype = getParameter(pt, "solver.BOOMERAMG.interpoltype", 0);
            int cycleType = getParameter(pt, "solver.BOOMERAMG.cycleType", 1);
            int relaxType = getParameter(pt, "solver.BOOMERAMG.relaxType", 3);
            int variant = getParameter(pt, "solver.BOOMERAMG.variant", 0);
            int overlap = getParameter(pt, "solver.BOOMERAMG.overlap", 1);
            double tol = getParameter(pt, "solver.BOOMERAMG.tol", iteEps);
            double strongThreshold = getParameter(pt, "solver.BOOMERAMG.strongThreshold", (dim==2)?0.25:0.6);
            BoomerAMG<AssGalOp >
                BoomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
                strongThreshold,variant,overlap,1,verbosity);
            Dune::LoopSolver<LinearSpace> cg(A,BoomerAMGPrecon,iteEps,iteSteps,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::EUCLID:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: EUCLID" << std::endl;
            }
            int level      = getParameter(pt, "solver.EUCLID.level",1);
            double droptol = getParameter(pt, "solver.EUCLID.droptol",0.01);
            int printlevel = 0;
            if (verbosity>2) printlevel=verbosity-2;
            printlevel = getParameter(pt,"solver.EUCLID.printlevel",printlevel);
            int bj = getParameter(pt, "solver.EUCLID.bj",0);
            Euclid<AssGalOp > EuclidPrecon(A,level,droptol,printlevel,bj,verbosity);
            CG<LinearSpace,LinearSpace> cg(A,EuclidPrecon,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::JACOBI:
          default:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: JACOBI" << std::endl;
            }
            JacobiPreconditioner<AssGalOp > jacobi(A,1.0);
            CG<LinearSpace,LinearSpace> cg(A,jacobi,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
        }
        solution *= -1.0;
        u.data = solution.data;
  
        if(verbosity > 0) {
          std::cout << "iterative solve eps= " << iteEps << ": " 
              << (res.converged?"converged":"failed") << " after "
              << res.iterations << " steps, rate="
              << res.conv_rate << ", computing time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
        }
            }
      
      //calculate error in l2norm
      VariableSet::VariableSet func( variableSet ) ;
      interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<0>(func.data),PoissonSolution());
      u -= func;
      L2Norm l2 ;
      double nrm2 = l2( boost::fusion::at_c<0>(u.data) ) ;
      if(verbosity > 0) {
        std::cout << "error in l2norm: " << nrm2 << "  has to be smaller than  " << errors[i][refSteps] << std::endl;
      }
      ress[refSteps] = nrm2;
      //calculate mesh size
      hs[refSteps] = 1.414213562/pow(2,refSteps+1+refinements);
      
      // query whether error is low enough
      if(!(nrm2 <= errors[i][refSteps])) {
      valid = false;
      message << "Test failed: The error after " << refSteps+2 << 
                 " refinements was too high at the test with ansatz functions of order " << order[i] << ".";
      }
      
      refSteps++;
    } while(refSteps < maxRefSteps);
        
     // calculate order of convergence by linear regression
    double a11=0, a12=0, a22=maxRefSteps, b1=0,b2=0;
    for(int i=0;i<maxRefSteps;i++) {
      hs[i] = std::log(hs[i]);
      ress[i] = std::log(ress[i]);
      a11 += hs[i]*hs[i];
      a12 += hs[i];
      b1 += hs[i]*ress[i];
      b2 += ress[i];
    }
    double det = 1.0/(a11*a22-a12*a12);
    double convOrd = det*(a22*b1-a12*b2);
    double logc = det*(a11*b2-a12*b1);
    if(verbosity > 0) {
      std::cout << "error = c*h^p with" << std::endl;
      std::cout << "p = " << convOrd << std::endl;
      std::cout << "log(c) = " << logc << std::endl;
    }
    
    if(convOrd < conOrder[i]) {
      valid = false;
      message << "Test failed: The order of convergence was too low at the test with ansatz functions of order "
              << order[i] << ".";
    }
  }

  std::cout << "total computing time: " << (double)(totalTimer.elapsed().user)/1e9 << "s\n";
  std::cout << message.str() << std::endl;
  if(result) {
    std::string description = "Test with simple poisson equation. Used " + solverTypeI + " solver " + solverTypeII + ":";
    std::ofstream outfile("../testResult.txt", std::ofstream::out | std::ofstream::app);
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }
}
