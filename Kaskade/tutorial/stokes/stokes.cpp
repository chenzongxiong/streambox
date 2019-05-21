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

//#define BOOST_DISABLE_ASSERTS

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/lagrangespace.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/uzawa.hh"
#include "linalg/direct.hh"
#include "linalg/cg.hh"
#include "io/vtk.hh"
#include "utilities/kaskopt.hh" // property_tree
#include "utilities/gridGeneration.hh"

using namespace Kaskade;
#include "stokes.hh"

int main(int argc, char *argv[])
{
  using namespace boost::fusion;


  int verbosity = 1;  // print to console if arguments are changed
  bool dump = false; // do not write properties into file
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosity, dump);

  std::cout << "Start stokes tutorial programm." << std::endl;

  boost::timer::cpu_timer totalTimer;

  constexpr int dim = 2;
  constexpr int uIdx = 0;
  constexpr int pIdx = 1;

  // command line parameters
  int refinements = getParameter(pt, "refinements", 5);
  int order = getParameter(pt, "order", 2); // order for velocity space, order for pressure space is (order-1)
  std::string empty;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  bool direct = getParameter(pt, s, 0);
    
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  DirectType directType = static_cast<DirectType>(getParameter(pt, s, 4));  // 4: DirectType::UMFPACK3264

  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "direct solver                  : " << directType << std::endl;

  boost::timer::cpu_timer gridTimer;
  // grid generation
  using Grid = Dune::UGGrid<dim>;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
  using Spaces = vector<H1Space const*,H1Space const*>;
  using VariableDescriptions = vector<Variable<SpaceIndex<1>,Components<2>,VariableId<uIdx> >,
                                      Variable<SpaceIndex<0>,Components<1>,VariableId<pIdx> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<>::type;
  using Functional = StokesFunctional<double,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;

  Dune::FieldVector<double,dim> x0(0.0), length(1.0);
  GridManager<Grid> gridManager( createRectangle<Grid>(x0,length,1.0));
  gridManager.globalRefine(refinements);
  std::cout << "computing time for generation of initial mesh: " << boost::timer::format(gridTimer.elapsed());


  // construct involved spaces.
  H1Space pressureSpace(gridManager,gridManager.grid().leafGridView(),order-1);
  H1Space velocitySpace(gridManager,gridManager.grid().leafGridView(),order);

  Spaces spaces(&pressureSpace,&velocitySpace);

  // construct variable list.
  std::string varNames[2] = { "u", "p" };

  VariableSet variableSet(spaces,varNames);

  // construct variational functional.
  Functional F;

  // construct Galerkin representation
  Assembler assembler(spaces);
  VariableSet::VariableSet x(variableSet);
  VariableSet::VariableSet dx(variableSet);

  size_t nnz = assembler.nnz(0,2,0,2,false);
  size_t dof = variableSet.degreesOfFreedom(0,2);
  std::cout << "overall degrees of freedom: " << dof << std::endl;
  std::cout << "(structurally) nonzero elements: " << nnz << std::endl;

  boost::timer::cpu_timer assembleTimer;
  assembler.assemble(linearization(F,x));
  std::cout << "computing time for assemble: " << boost::timer::format(assembleTimer.elapsed());

  if(direct)
  {
    CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<>::init(spaces));
    CoefficientVectors rhs(assembler.rhs());

    boost::timer::cpu_timer directTimer;
    // solve performing one Newton step
    directInverseOperator(AssembledGalerkinOperator<Assembler>(assembler),directType).applyscaleadd(-1.0,rhs,solution);
    std::cout << "computing time for directsolve: " << boost::timer::format(directTimer.elapsed());
    x.data = solution.data;
  }
  else // Uzawa Solver
  {
    using VectorOfU = VariableSet::CoefficientVectorRepresentation<uIdx,uIdx+1>::type;
    using VectorOfP = VariableSet::CoefficientVectorRepresentation<pIdx,pIdx+1>::type;
    using Assembler_UU = AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,uIdx,uIdx+1>;
    using Assembler_PU = AssembledGalerkinOperator<Assembler,pIdx,pIdx+1,uIdx,uIdx+1>;
    using Assembler_UP = AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,pIdx,pIdx+1>;
    using PreconAdapt = MatrixRepresentedOperator<MatrixAsTriplet<double>, VectorOfP, VectorOfP>;
    using UzSo = UzawaSolver<VectorOfU,VectorOfP>;

    Assembler_UU A(assembler);
    Assembler_PU B(assembler);
    Assembler_UP Bt(assembler);

    boost::timer::cpu_timer iterativeTimer;
    Dune::InverseOperatorResult res;
    const DefaultDualPairing<VectorOfU,VectorOfU> defaultScalarProduct{};
    StrakosTichyPTerminationCriterion<double> termination(1e-14,300);
    int lookAhead = getParameter(pt, "solver.lookAhead", 3);
    termination.setLookAhead(lookAhead);

    JacobiPreconditioner<Assembler_UU> jacobiPreconditioner(A);
    // inexact inner solver for upper left block
    //Dune::CGSolver<VectorOfU> cg(A,jacobiPreconditioner,1e-14,300,0);
    CG<VectorOfU,VectorOfU> cg(A,jacobiPreconditioner,defaultScalarProduct,termination,verbosity);
    TrivialPreconditioner<PreconAdapt> trivialPreconditioner;

    VectorOfU f(assembler.rhs<uIdx,uIdx+1>());
    VectorOfP g(assembler.rhs<pIdx,pIdx+1>());
    VectorOfU u(VariableSet::CoefficientVectorRepresentation<uIdx,uIdx+1>::init(spaces));
    VectorOfP p(VariableSet::CoefficientVectorRepresentation<pIdx,pIdx+1>::init(spaces));

    UzSo uzawa(A,cg,B,Bt,trivialPreconditioner,1e-4,100,2);

    UzSo::Domain solution(vector<VectorOfU,VectorOfP>(u,p));
    UzSo::Range rhs(vector<VectorOfU,VectorOfP>(f,g));
    rhs *= -1.0; // change sign as rhs is -'F while the assembler returns F'
    uzawa.apply(solution,rhs,res);

    std::cout << "computing time for iterative solve: " << boost::timer::format(iterativeTimer.elapsed());

    at_c<uIdx>(x.data).coefficients() = at_c<0>(at_c<uIdx>(solution.data).data);
    at_c<pIdx>(x.data).coefficients() = at_c<0>(at_c<pIdx>(solution.data).data);
  }

  boost::timer::cpu_timer outputTimer;
  writeVTKFile(x,"stokes", IoOptions().setOrder(order));
  std::cout << "graphical output finished, data in VTK format is written into file stokes.vtu \n";
  std::cout << "computing time for output: " << boost::timer::format(outputTimer.elapsed()) << "\n";

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End stokes tutorial program" << std::endl;
}
