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

#ifdef Cygwin
// availability of std::to_string
#define _GLIBCXX_USE_C99 1
#endif

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"
#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;

#include "integrate-navierStokes.hh"
//#include "integrate-Limex.hh"

#include "navierStokes.hh"


struct Initial_1Value 
{
  using Scalar = double;
  static int const components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  Initial_1Value(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
  {
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);

    if (component==0)               // vector u = (u_1,u_2)
    {
      std::cout << "Error: wrong call of Initial_1Value" << std::endl;
      return 0;
    }
    else if (component==1)          // scalar p
    {
      return 0;
    }
    else
      assert("wrong index!\n"==0);
      
    return 0;
  }

private:
  int component;
};

struct Initial_2Value 
{
  using Scalar = double;
  static int const components = 2;
  using ValueType = Dune::FieldVector<Scalar,components>;

  Initial_2Value(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
  {
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);

    ValueType result(0);
    if (component==0)                    // vector u = (u_1,u_2), concrete
    {
      return result;
    }
    else
      assert("wrong index!\n"==0);
      
    return result;
  }

private:
  int component;
};



int main(int argc, char *argv[])
{
  using namespace boost::fusion;


  int verbosityOpt = 1;  // print to console if arguments are changed
  bool dump = false; // do not write properties into file
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  std::cout << "Start Navier-Stokes tutorial programm." << std::endl;

  boost::timer::cpu_timer totalTimer;

  constexpr int dim = 2; 
  
  constexpr int uIdx  = 0;
  constexpr int pIdx  = 1;
   
  constexpr double rho    = 1.0;
  constexpr double lambda = 0.0;
  constexpr double mu     = 0.01;   // ----> Reynold number Re = 100

  using Grid = Dune::UGGrid<dim>;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;
  //using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,Grid::LeafGridView> >;
  using Spaces = vector<H1Space const*,H1Space const*>;
  using VariableDescriptions = vector<Variable<SpaceIndex<1>,Components<2>,VariableId<uIdx> >,
                                      Variable<SpaceIndex<0>,Components<1>,VariableId<pIdx> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<>::type;
  using Functional = NavierStokesFunctional<double,VariableSet>;

  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq   = Functional::TestVars::noOfVariables;
 
  std::cout << std::endl;
  std::cout << "density rho     = "  << rho    << std::endl;
  std::cout << "viscosity mu    = "  << mu     << std::endl;

  double dt, dtMax, tStart, tEnd, rTolT, aTolT, rTolX, aTolX, writeInterval;

  tStart = 0.0; 
  
  // command line parameters
  int verbosity = getParameter(pt, "verbosity", 1);
  int extrapolOrder = getParameter(pt, "extrapolOrder", 0);
  int maxTimeSteps  = getParameter(pt, "maxTimeSteps", 51);

  int refinements = getParameter(pt, "refinements", 6);
  int order = getParameter(pt, "order", 2); // order for velocity space, order for pressure space is (order-1)
  std::string empty;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  bool direct = getParameter(pt, s, 0);
    
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  DirectType directType = static_cast<DirectType>(getParameter(pt, s, 4));  // 4: DirectType::UMFPACK3264

  std::cout << std::endl;
  std::cout << "order of extrapolation in Limex : " << extrapolOrder << std::endl;
  std::cout << "original mesh shall be refined  : " << refinements << " times" << std::endl;
  std::cout << "discretization order in space   : " << order << std::endl;
  //std::cout << "direct solver                   : " << directType << std::endl;
  std::cout << std::endl;

  boost::timer::cpu_timer gridTimer;
  // grid generation
  Dune::FieldVector<double,dim> x0(0.0), length(1.0);
  //GridManager<Grid> gridManager( createRectangle<Grid>(x0,length,0.5));   // mesh 2
  GridManager<Grid> gridManager( createRectangle<Grid>(x0,length,1.0));     // mesh 1
  gridManager.globalRefine(refinements);

  std::cout << "Initial grid: trs=" << gridManager.grid().size(0) 
            << " eds=" << gridManager.grid().size(1) << " pts=" << gridManager.grid().size(2) 
            << std::endl << std::endl;
  //std::cout << "computing time for generation of initial mesh: " << boost::timer::format(gridTimer.elapsed());


  // construct involved spaces.

  H1Space pressureSpace(gridManager,gridManager.grid().leafGridView(),order-1);
  H1Space velocitySpace(gridManager,gridManager.grid().leafGridView(),order);


  Spaces spaces(&pressureSpace,&velocitySpace);

  // construct variable list.
  // concrete 0: u = (u_1,u_2), velocity 
  //          1: p              pressure
  std::string varNames[2] = { "u", "p" };

  VariableSet variableSet(spaces,varNames);


  // construct variational functional.
  Functional F(tStart,rho,lambda,mu);	

  VariableSet::VariableSet x(variableSet);
  VariableSet::VariableSet dx(variableSet);

  F.time(tStart);
  x = 0; 
  F.scaleInitialValue<uIdx>(Initial_2Value(uIdx),x);      // u
  F.scaleInitialValue<pIdx>(Initial_1Value(pIdx),x);      // p
  
  std::cout << "nvars = " << nvars << std::endl;
  std::cout << "  neq = " << neq   << std::endl;

  
  tStart = 0.0;
  tEnd   = 800.0;           	
  dt     = 2e-2; 		
  dtMax  = 2e-2; 
  
  std::cout << "start time = "   << tStart       << std::endl;
  std::cout << "final time = "   << tEnd         << std::endl;
  std::cout << "step size = "    << dt           << std::endl;
  std::cout << "maxTimeSteps = " << maxTimeSteps << std::endl;
  std::cout << std::endl;
  
  rTolT = 1.0e+30;		//2e-5;		//1e-5/8.0;
  aTolT = 1.0e+30;		//2e-5;		//1e-5/8.0;
  rTolX = 1.0e+12;		//1e-4;   1e-5*2.0;
  aTolX = 1.0e+12;		//1e-4;   1e-5*2.0;

  std::cout << "aTolT = " << aTolT << std::endl;
  std::cout << "rTolT = " << rTolT << std::endl;
  std::cout << "aTolX = " << aTolX << std::endl;
  std::cout << "rTolX = " << rTolX << std::endl;
  std::cout << std::endl;

  std::vector<VariableSet::VariableSet> solutions;
  
  x = integrate(gridManager,F,variableSet,spaces,
                dt,dtMax,tEnd,maxTimeSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::MUMPS,
                verbosity);


  //writeVTKFile(gridManager.grid().leafGridView(),x,"navierStokes", IoOptions(), order);
  //std::cout << "graphical output finished, data in VTK format is written into file stokes.vtu \n";

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End Navier-Stokes tutorial program" << std::endl;
}
