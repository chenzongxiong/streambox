/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "timestepping/limexWithoutJens.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;

#include "mygrid.hh"
#include "integrate.hh"
#include "movingsource.hh"

struct InitialValue 
{
  using Scalar = double;
  static int const components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  InitialValue(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::ctype,Cell::dimension> const& localCoordinate) const 
  {
    Dune::FieldVector<typename Cell::ctype,Cell::dimensionworld> x = cell.geometry().global(localCoordinate);
    x -= 0.5;
    return 0.0;
  }

private:
  int component;
};

int main(int argc, char *argv[])
  {
    int const dim = 2;
    int refinements = 5, order = 2, extrapolOrder = 2, maxSteps = 100,
      verbosity=1;
    double dt = 0.1, maxDT = 1.0, T = 1.0, rTolT = 1.0e-3, aTolT = 1.0e-3, rTolX = 2.0e-5, aTolX = 2.0e-5, writeInterval = 1.0;

    std::cout << "Start moving source tutorial program" << std::endl;

    using Grid = Dune::UGGrid<dim>;
    using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
    using Spaces = boost::fusion::vector<H1Space const*>;
    using VariableDescriptions = boost::fusion::vector<VariableDescription<0,1,0> >;
    using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
    using Equation = MovingSourceEquation<double,VariableSet>;

    GridManager<Grid> gridManager(RefineGrid<Grid>(refinements));
    std::cout << "Grid: " << gridManager.grid().size(0) << " " << gridManager.grid().size(1) << " " 
              << gridManager.grid().size(2) << std::endl;

  
// construct involved spaces.

  H1Space temperatureSpace(gridManager,gridManager.grid().leafView(),order);

    Spaces spaces(&temperatureSpace);

    std::string varNames[1] = { "u" };
  
    VariableSet variableSet(spaces,varNames);

    Equation Eq;

  std::vector<VariableSet::VariableSet> solutions;
  // std::vector<VariableSet::VariableSet> devnull;

  Eq.time(0);
  VariableSet::VariableSet x(variableSet);
  Eq.scaleInitialValue<0>(InitialValue(0),x);
  
  x = integrate(gridManager,Eq,variableSet,spaces,
          dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
          std::back_inserter(solutions),writeInterval,x,DirectType::SUPERLU,verbosity);

    std::cout << "End moving source tutorial program" << std::endl;

  }
