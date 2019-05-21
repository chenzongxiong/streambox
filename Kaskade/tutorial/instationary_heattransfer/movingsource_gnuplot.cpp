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

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "timestepping/limexWithoutJens.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;

#include "integrate_gnuplot.hh"
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
                  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
  {
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);

    double vx = x[0]-0.5;
    double vy = x[1]-0.75;

    return 0.8*exp(-80.0*(vx*vx + vy*vy));
  }

private:
  int component;
};

int main(int argc, char *argv[])
{
  int const dim = 2;
  int refinements = 4, order = 2, extrapolOrder = 2, maxSteps = 1000,
      verbosity=1;
  double dt = 0.1, maxDT = 1.0, T = 1.0, rTolT = 1.0e-3, aTolT = 1.0e-3, rTolX = 2.0e-5, aTolX = 2.0e-5, writeInterval = 1.0;
  if ( argc > 1 ) { sscanf(argv[1],"%i",&verbosity); };

  boost::timer::cpu_timer totalTimer;
  std::cout << "Start moving source tutorial program" << std::endl;

  using Grid = Dune::UGGrid<dim>;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,1,0> >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Equation = MovingSourceEquation<double,VariableSet>;
  
  GridManager<Grid> gridManager(createUnitSquare<Grid>(0.25));
  gridManager.globalRefine(refinements);
  std::cout << "Grid: " << gridManager.grid().size(0) << " " << gridManager.grid().size(1) << " " << gridManager.grid().size(2) << std::endl;
 
  gridManager.enforceConcurrentReads(true);
  
  // construct involved spaces.

  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);

  Spaces spaces(&temperatureSpace);

  std::string varNames[1] = { "u" };
  
  VariableSet variableSet(spaces,varNames);

  Equation Eq;

  std::vector<VariableSet::VariableSet> solutions;

  Eq.time(0);
  VariableSet::VariableSet x(variableSet);
  Eq.scaleInitialValue<0>(InitialValue(0),x);
    
  x = integrate(gridManager,Eq,variableSet,spaces,
                dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::SUPERLU,
                verbosity);

  std::cout << "End moving source tutorial program" << std::endl;
  std::cout << "used cpu time: " << (double)(totalTimer.elapsed().user)/1e9 << "s\n";
}
