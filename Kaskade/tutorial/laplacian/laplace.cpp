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

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/variables.hh"
#include "io/vtk.hh"
#include "linalg/direct.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;
#include "laplace.hh"

int main()
{
  std::cout << "Start Laplacian tutorial program" << std::endl;

  constexpr int dim = 2;
  int refinements = 5,
      order       = 2;

  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = HeatFunctional<double,VariableSetDesc>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using Operator = AssembledGalerkinOperator<Assembler>;
  using CoefficientVectors = VariableSetDesc::CoefficientVectorRepresentation<0,1>::type;


  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);

  // construction of finite element space for the scalar solution u.
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
  Spaces spaces(&temperatureSpace);
  VariableSetDesc variableSetDesc(spaces,{ "u" });

  Functional F;

  //construct Galerkin representation
  Assembler assembler(spaces);
  VariableSetDesc::VariableSet u(variableSetDesc);
  assembler.assemble(linearization(F,u));

  Operator A(assembler);
  CoefficientVectors solution(VariableSetDesc::CoefficientVectorRepresentation<>::init(spaces));
  CoefficientVectors rhs(assembler.rhs());

  directInverseOperator(A).applyscaleadd(-1.0,rhs,solution);
  component<0>(u) = component<0>(solution);

  writeVTKFile(u,"temperature",IoOptions().setOrder(order).setPrecision(7));

  std::cout << "graphical output finished, data in VTK format is written into file temperature.vtu \n";
  std::cout << "End Laplacian tutorial program" << std::endl;
}
