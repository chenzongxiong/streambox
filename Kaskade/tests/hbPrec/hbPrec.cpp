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
 * @brief  Test for the PrecondType::HB preconditioner with simple heat transfer equation.
 * 
 * Testprogram for Kaskade: tests wether the error between calculated and exact solution behaves as expected. 
 * This test was built by changing the stationary heat transfer example. 
 * 
 * Variational functional, boundary conditions and f are adapted in hbPrec.hh. 
 * Only uses ansatz functions of order 1.
 */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"
#include "fem/norms.hh"
#include "linalg/direct.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "utilities/enums.hh"
#include "utilities/kaskopt.hh"
#include "mg/hb.hh"

#include "utilities/kaskopt.hh"

using namespace Kaskade;
#include "hbPrec.hh"

/**
 * @brief  This represents the known exact solution. 
 * It is \f$ u = \cos (2 \pi x) + \exp (4(y-0.5)^2) \f$.
 */
struct HBPrecSolution
{
   using Scalar = double;
   static int const components = 1;
   using ValueType = Dune::FieldVector<Scalar,components>;

   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   ValueType value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
     Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
     return std::cos(2*PI*x[0]) + std::exp(4*(x[1]-0.5)*(x[1]-0.5));
   }
};


int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer test program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 0;
  bool dump = false; 
  constexpr int dim=2; 
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  //using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = HBPrecFunctional<double,VariableSet>;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;

  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);
  
  //ansatz orders to be tested
  int order = 1;
  constexpr int maxRefSteps = 9;
  //error values for corresponding refinement
  double errors[maxRefSteps] = {0.65,0.13,0.034,0.0086,0.0021,5.4e-4,1.35e-4,3.4e-5,8.4e-6};
  bool valid = true;
  std::stringstream message("Test succeeded", std::stringstream::out);
  
  int verbosity   = getParameter(pt, "verbosity", 1);
  // if true, then the test result will be written in a file
  bool result = getParameter(pt, "result",0);
  int onlyLowerTriangle = false;
    
  MatrixProperties property = MatrixProperties::SYMMETRIC;
  
  int blocks = getParameter(pt,"blocks",40);
  int nthreads = getParameter(pt,"threads",4);
  double rowBlockFactor = getParameter(pt,"rowBlockFactor",2.0);

  property = MatrixProperties::SYMMETRIC;
  

  if(verbosity > 0) {
    std::cout << "original mesh shall be refined : " << maxRefSteps << " times" << std::endl;
    std::cout << "discretization order           : " << order << std::endl;
    std::cout << "output level (verbosity)       : " << verbosity << std::endl;
  }
      
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
  // the coarse grid will be refined three times
  grid->globalRefine(0);
  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(2) << " points" << std::endl;

  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
    
    
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
  Functional F;
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
  
  for(int refSteps = 0;refSteps<maxRefSteps;refSteps++) {
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
    if(verbosity > 0) {
      std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed()) << "\n";
    }
    
    CoefficientVectors rhs(assembler.rhs());
    AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, onlyLowerTriangle);
    
    boost::timer::cpu_timer iteTimer;
    int iteSteps = getParameter(pt, "solver.iteMax", 1000);
    double iteEps = getParameter(pt, "solver.iteEps", 1.0e-10);
    Dune::InverseOperatorResult res;
    const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
    StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
    int lookAhead = getParameter(pt, "solver.lookAhead", 50);
    termination.setLookAhead(lookAhead);

    HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type, AssembledGalerkinOperator<Assembler,0,1,0,1>::range_type > hb(gridManager.grid());
    CG<LinearSpace,LinearSpace> cg(A,hb,defaultScalarProduct,termination,verbosity);
    cg.apply(solution,rhs,res);
  
    solution *= -1.0;
    u.data = solution.data;
    
    if(verbosity > 0) {
      std::cout << "iterative solve eps= " << iteEps << ": " 
    << (res.converged?"converged":"failed") << " after "
    << res.iterations << " steps, rate="
    << res.conv_rate << ", computing time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
    }
    
    //calculate error in l2norm
    VariableSet::VariableSet func( variableSet ) ;
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<0>(func.data),HBPrecSolution());
    u -= func;
    L2Norm l2 ;
    double nrm2 = l2( boost::fusion::at_c<0>(u.data) ) ;
    if(verbosity > 0) {
      std::cout << "error in l2norm: " << nrm2 << std::endl;
    }
    // query whether error is low enough
    if(!(nrm2<=errors[refSteps])) {
      valid = false;
      message << "Test failed: The error after " << refSteps << " refinements was too high at the test with ansatz functions of order 1.";
    }
  }
       
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << message.str() << std::endl;
  if(result) {
    std::string description = "Test with example stationary heat transfer, 2D. Used iterate solver IterateType::CG, PrecondType::HB preconditioner:";
    std::ofstream outfile("../testResult.txt", std::ofstream::out | std::ofstream::app);
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }
}
