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
 * @file
 * @ingroup tests
 * @brief  Test with simple heat transfer equation, 3D.
 * 
 * Testprogram for Kaskade: tests wether the error between calculated and exact solution behaves as expected. 
 * This test was built by by using the stationary heat transfer example.
 * Boundary conditions and f are adapted in statHeatTransfer.hh.
 * 
 * The test does not work with the following preconditioners: PrecondType::NONE, PrecondType::ARMS, PrecondType::ADDITIVESCHWARZ, PrecondType::HB, PrecondType::ICC. 
 * The test works with: PrecondType::JACOBI, PrecondType::ILUT, PrecondType::ILUK,PrecondType::BOOMERAMG, PrecondType::EUCLID, PrecondType::ICC0 and with direct solvers DirectType::MUMPS, DirectType::UMFPACK and DirectType::UMFPACK3264. 
 * Uses ContinuousHierachicMapper.
 */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/istl/solvers.hh"

#include "fem/assemble.hh"
#include "fem/embedded_errorest.hh"
#include "fem/norms.hh"
//#include "fem/lagrangespace.hh"
#include "fem/hierarchicspace.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/partialDirectPreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/iccprecond.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG
#include "linalg/iluprecond.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "utilities/enums.hh"
#include "utilities/kaskopt.hh"


#include "cubus.hh"

using namespace Kaskade;
#include "statHeatTransfer.hh"
#undef max

/**
 * @brief  This represents the known exact solution. 
 * It is \f$ u = x(x-1)\exp (-(x-0.5)^2) \f$.
 */
struct StatHeatTransferSolution
{
   using Scalar = double;
   static int const components = 1;
   using ValueType = Dune::FieldVector<Scalar,components>;

   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   ValueType value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
     Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
     return x[0]*(x[0]-1) * std::exp(-(x[0]-0.5)*(x[0]-0.5));
   }
};


int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer test program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 0;
  bool dump = false; 
  constexpr int dim=3; 
  constexpr int heapSize=1024;
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  //using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = StatHeatTransferFunctional<double,VariableSet>;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,1>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,1>::type;

  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);
  
  constexpr int orderLength = 5;
  //ansatz orders to be tested
  int order[orderLength] = {1,2,3,4,5};
  int  maxRefSteps[orderLength] = {3,3,1,1,1};
  //error values for corresponding ansatz order and refinement
  double errorsI[2][4] = {{9e-3,3.3e-3,9.8e-4,2.5e-4}, //order 1
                          {3.2e-3,2.7e-4,2.4e-5,2.4e-6}}; //order 2
  double errorsII[3] = {4.5e-5, //order3
                        9.7e-7, //order4
                        5.9e-8}; //order5
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
  
  //  three-dimensional space: dim=3
  //  - more complex: 1 cube defined by 48 tetrahedra, provided in cubus.hh
    
  // definition of a more complex mesh using cubus.hh
  
 
  //iterate over the ansatz orders
  for(int i = 0;i < orderLength && valid;i++) {
    if(verbosity > 0) {
      std::cout << "original mesh shall be refined : " << maxRefSteps[i] << " times" << std::endl;
      std::cout << "discretization order           : " << order[i] << std::endl;
      std::cout << "output level (verbosity)       : " << verbosity << std::endl;
    }

     std::unique_ptr<Grid> grid( RefineGrid<Grid>(0, heapSize) );

    // some information on the refined mesh
    if(verbosity > 0) {
      std::cout << std::endl << "Grid: " << grid->size(0) << " tetrahedra, " << std::endl;
      std::cout << "      " << grid->size(1) << " triangles, " << std::endl;
      std::cout << "      " << grid->size(dim-1) << " edges, " << std::endl;
      std::cout << "      " << grid->size(dim) << " points" << std::endl;
    }
    // a gridmanager is constructed 
    // as connector between geometric and algebraic information
    GridManager<Grid> gridManager(std::move(grid));
      
      
    // construction of finite element space for the scalar solution T.
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
      
    double kappa = getParameter(pt, "kappa", 1.0);
    double q = getParameter(pt, "q", 1.0);
    Functional F(kappa,q);
    if(verbosity > 0) {
      std::cout << std::endl << "no of variables = " << nvars << std::endl;
      std::cout << "no of equations = " << neq   << std::endl;
      size_t dofs = variableSet.degreesOfFreedom(0,nvars);
      std::cout << "number of degrees of freedom = " << dofs   << std::endl;
    }

      
    //construct Galerkin representation
    Assembler assembler(gridManager,spaces);
    
    gridManager.enforceConcurrentReads(true);
    assembler.setNSimultaneousBlocks(blocks);
    assembler.setRowBlockFactor(rowBlockFactor);
    
    //keep mesh sizes
    //double hs[maxRefSteps+1];
    //keep errors
    //double ress[maxRefSteps+1];
    for(int refSteps = 0;refSteps<maxRefSteps[i]+1;refSteps++) {
      if(refSteps>0) {
        gridManager.globalRefine(1);
      } else {
        if(i>=2) continue;
      }
      
      boost::timer::cpu_timer assembTimer;
      VariableSet::VariableSet u(variableSet);
      
      size_t nnz = assembler.nnz(0,1,0,1,onlyLowerTriangle);
      if(verbosity > 0) {
        std::cout << "number of nonzero elements in the stiffness matrix: " << nnz << std::endl << std::endl;
      }
      
      CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,1>::init(spaces));
      solution = 0;
      
      assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
      if(verbosity > 0) {
        std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed()) << "\n";
      }
      
      CoefficientVectors rhs(assembler.rhs());
      AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, onlyLowerTriangle);
      
      if (direct)
      {
        solverTypeII = directSolver;
        boost::timer::cpu_timer directTimer;
        directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
        u.data = solution.data;
        if(verbosity > 0) {
          std::cout << "computing time for direct solve: " << boost::timer::format(directTimer.elapsed()) << "\n";
        }
      }
      else
      {
        solverTypeII = "CG, " + preconditioner;
        boost::timer::cpu_timer iteTimer;
        int iteSteps = getParameter(pt, "solver.iteMax", 1000);
        double iteEps = getParameter(pt, "solver.iteEps", 1.0e-10);
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
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl 
                      << "ATTENTION: Test does not work without preconditioner." << std::endl
                      << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: NONE" << std::endl;
            }
            TrivialPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > trivial;
            CG<LinearSpace,LinearSpace> cg(A,trivial,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ADDITIVESCHWARZ:
          {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl << "ATTENTION: Test does not work with preconditioner PrecondType::ADDITIVESCHWARZ." << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ADDITIVESCHWARZ" << std::endl;
            }
            std::pair<size_t,size_t> idx = temperatureSpace.mapper().globalIndexRange(gridManager.grid().leafIndexSet().geomTypes(dim)[0]);
            AdditiveSchwarzPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > addschwarz(A,idx.first,idx.second,verbosity);
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
            ILUTPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > ilut(A,lfil,dropTol,verbosity);
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
            ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > iluk(A,fill_lev,verbosity);
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
            ICCPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > icc(A,dropTol);
            CG<LinearSpace,LinearSpace> cg(A,icc,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ICC0:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ICC0" << std::endl;
            }
            ICC_0Preconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > icc0(A);
            CG<LinearSpace,LinearSpace> cg(A,icc0,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::HB:
          {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl << "ATTENTION: Test does not work with preconditioner PrecondType::HB." << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            if(verbosity > 0) {
              std::cout << "selected preconditioner: HB" << std::endl;
            }
            HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type, AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type > hb(gridManager.grid());
            CG<LinearSpace,LinearSpace> cg(A,hb,defaultScalarProduct,termination,verbosity);
            cg.apply(solution,rhs,res);
          }
          break;
          case PrecondType::ARMS:
          {
            if(verbosity > 0) {
              std::cout << "selected preconditioner: ARMS" << std::endl;
            }
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl 
              << "ATTENTION: Test does not work with preconditioner PrecondType::ARMS." << std::endl 
              << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            int lfil = getParameter(pt, "solver.ARMS.lfil", 140);
            int lev_reord = getParameter(pt, "solver.ARMS.lev_reord", 1);
            double dropTol = getParameter(pt, "solver.ARMS.dropTol", 0.01);
            double tolind = getParameter(pt, "solver.ARMS.tolind", 0.2);
            ARMSPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > iluk(A,lfil,dropTol,lev_reord,tolind,verbosity);
            CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
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
            //Dune::LoopSolver<LinearSpace> cg(A,BoomerAMGPrecon,iteEps,iteSteps,verbosity);
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
            Euclid<AssembledGalerkinOperator<Assembler,0,1,0,1> > EuclidPrecon(A,level,droptol,printlevel,bj,verbosity);
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
            JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > jacobi(A,1.0);
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
      interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<0>(func.data),StatHeatTransferSolution());
      u -= func;
      L2Norm l2 ;
      double nrm2 = l2( boost::fusion::at_c<0>(u.data) ) ;
      if(verbosity > 0) {
        std::cout << "error in l2norm: " << nrm2 << std::endl;
      }
      
      // query whether error is low enough
      if(i<2 && !(nrm2<=errorsI[i][refSteps])) {
        valid = false;
        message << "Test failed: The error after " << refSteps << " refinements was too high at the test with ansatz functions of order " << order[i] << ".";
      } else if(i>=2 && !(nrm2<=errorsII[i-2])) {
        valid = false;
        message << "Test failed: The error after 1 refinement was too high at the test with ansatz functions of order " << order[i] << ".";
      }
      //ress[refSteps] = nrm2;
      //calculate mesh size
      //hs[refSteps] = 1.414213562/pow(2,refSteps+1);
    }
  }
  
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << message.str() << std::endl;
  if(result) {
    std::string description = "Test with example stationary heat transfer, 3D. Used " + solverTypeI + " solver " + solverTypeII + ":";
    std::ofstream outfile("../testResult.txt", std::ofstream::out | std::ofstream::app);
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }
}
