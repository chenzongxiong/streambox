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
#include "dune/istl/solvers.hh"

#include "fem/assemble.hh"
#include "fem/lagrangespace.hh"
#include "fem/embedded_errorest.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/uzawa.hh"
#include "linalg/direct.hh"
#include "io/vtk.hh"
#include "utilities/kaskopt.hh" // property_tree
#include "utilities/gridGeneration.hh"

using namespace Kaskade;
#include "stokes-adaptive.hh"

int main(int argc, char *argv[])
{
  using namespace boost::fusion;


  int verbosity = 1;  // print to console if arguments are changed
  bool dump = false; // do not write properties into file
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosity, dump);

  std::cout << "Start Stokes tutorial programm" << std::endl;
  std::cout << "using addaptive refinement based on embedded error estimation." << std::endl;
  std::cout << "Note: use FixedFractionCriterion as refinement strategy!" << std::endl << std::endl;

  constexpr int dim = 2;
  constexpr int uIdx = 0;
  constexpr int pIdx = 1;

  // command line parameters
  int refinements = getParameter(pt, "refinements", 2);   // initial mesh is refined "refinements" times
  int order = getParameter(pt, "order", 3);               // order for velocity space, order for pressure space is (order-1)
  int steps = getParameter(pt, "steps", 5);               // adaptive refinement steps

  bool direct = getParameter(pt, "direct", true);
  DirectType directType = static_cast<DirectType>(getParameter(pt, "directSolver", 4)); // 4: DirectType::UMFPACK3264


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

  // accuracy in space
  std::vector<std::pair<double,double> > tol(variableSet.noOfVariables);
  tol[0] = std::make_pair(1e0,1e-2);   // u
  tol[1] = std::make_pair(1e0,1e-2);   // p
  bool accurate = 0;

  for (int i=0; i<steps; ++i)    // adaptive refinement loop
  {
    std::cout << "\nstep = " << i+1 << std::endl;
    size_t nnz = assembler.nnz(0,2,0,2,false);
    size_t dof = variableSet.degreesOfFreedom(0,2);
    std::cout << "overall degrees of freedom: " << dof << std::endl;
    //std::cout << "(structurally) nonzero elements: " << nnz << std::endl;

    x = 0;
    assembler.assemble(linearization(F,x));

    if(direct)
    {
      CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<>::init(spaces));
      CoefficientVectors rhs(assembler.rhs());

      // solve performing one Newton step
      directInverseOperator(AssembledGalerkinOperator<Assembler>(assembler),directType).applyscaleadd(-1.0,rhs,solution);
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

      JacobiPreconditioner<Assembler_UU> jacobiPreconditioner(A);
      // inexact inner solver for upper left block
      Dune::CGSolver<VectorOfU> cg(A,jacobiPreconditioner,1e-14,300,0);

      TrivialPreconditioner<PreconAdapt> trivialPreconditioner;

      VectorOfU f(assembler.rhs<uIdx,uIdx+1>());
      VectorOfP g(assembler.rhs<pIdx,pIdx+1>());
      VectorOfU u(VariableSet::CoefficientVectorRepresentation<uIdx,uIdx+1>::init(spaces));
      VectorOfP p(VariableSet::CoefficientVectorRepresentation<pIdx,pIdx+1>::init(spaces));

      UzSo uzawa(A,cg,B,Bt,trivialPreconditioner,1e-4,100,2);

      Dune::InverseOperatorResult res;
      UzSo::Domain solution(vector<VectorOfU,VectorOfP>(u,p));
      UzSo::Range rhs(vector<VectorOfU,VectorOfP>(f,g));
      rhs *= -1.0; // change sign as rhs is -'F while the assembler returns F'
      uzawa.apply(solution,rhs,res);

      at_c<uIdx>(x.data).coefficients() = at_c<0>(at_c<uIdx>(solution.data).data);
      at_c<pIdx>(x.data).coefficients() = at_c<0>(at_c<pIdx>(solution.data).data);
    }

    writeVTKFile(gridManager.grid().leafGridView(),x,"stokes-adaptive", IoOptions(), 1);

    VariableSet::VariableSet e = x;
    projectHierarchically(variableSet,e);
    e -= x;

    accurate = embeddedErrorEstimator(variableSet,e,x,IdentityScaling(),tol,gridManager,1);
    dof = variableSet.degreesOfFreedom(0,2);
    std::cout << "degrees of freedom after refinement = " << dof << std::endl;

  }

  return 0;
}
