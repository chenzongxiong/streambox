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

#define FUSION_MAX_VECTOR_SIZE 25
#include <iostream>
//#include <fstream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/istl/solvers.hh"

#include "algorithm/giant_gbit.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "io/vtk.hh"
#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;

#include "sst.hh"

struct InitialValue 
{
  using Scalar = double;
  static constexpr int components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  InitialValue(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> const& localCoordinate) const
  {
  // Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
  if (component==0) 
    return 1.0e9;
  else if (component==1) 
    return 1.0e9;
  else if (component==2) 
    return 1.0e13;
  else if (component==3) 
    return 1.0e7;
  else
    assert("wrong index!\n"==0);
  return 0;
  
  }

private:
  int component;
};

int main(int argc, char *argv[])
{
  using Scalar = double;
  using namespace boost::fusion;

  std::cout << "Start sst transfer tutorial program with inexact damped Newton iteration " << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int refinements = getParameter(pt, "refinement", 5),
      order =  getParameter(pt, "order", 2);
//      verbosity   = getParameter(pt, "verbosity", 0);
  Scalar tol = getParameter(pt, "tolerance", 1.0e-10),
         rho = getParameter(pt, "safetyfactor",0.0625);
  //   IterateType iterateType = IterateType::CG;
  //   PrecondType precondType = PrecondType::NONE;
  std::string empty;

  std::cout << "refinements of original mesh   : " << refinements << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "tolerance for Newton iteration : " << tol << std::endl;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  //direct = getParameter(pt, s, 0);

  //   s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
  //   iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  //   s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  //   precondType = static_cast<PrecondType>(getParameter(pt, s, 0));

  constexpr int dim=2;    
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<Scalar,LeafView> >;
// using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<Scalar,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<1> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<2> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<3> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = SSTFunctional<Scalar,VariableSet>;
  using NonlinearSolver = Giant<Grid,Functional,VariableSet,Spaces>;

  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " triangles, " << std::endl;
  std::cout << "      " << gridManager.grid().size(1) << " edges, " << std::endl;
  std::cout << "      " << gridManager.grid().size(2) << " points" << std::endl;

  // construction of finite element space for the scalar solution T
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),
               order);
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  std::string varNames[4] = { "u0", "u1", "u2", "u3" };
  VariableSet variableSet(spaces,varNames);

  Functional F;
  VariableSet::VariableSet x(variableSet);
  
  F.scaleInitialValue<0>(InitialValue(0),x);
  F.scaleInitialValue<1>(InitialValue(1),x);
  F.scaleInitialValue<2>(InitialValue(2),x);
  F.scaleInitialValue<3>(InitialValue(3),x);
  
  LeafView leafGridView = gridManager.grid().leafGridView();
  writeVTKFile(x,"graph/sst_giant_start",IoOptions().setOrder(std::min(order,2)).setPrecision(7));
  gridManager.enforceConcurrentReads(false);
  
//   constexpr int nvars = Functional::AnsatzVars::noOfVariables;
//   constexpr int n     = variableSet.degreesOfFreedom(0,nvars);
//   std::vector<double> scale(n);
//   x.write(scale.begin());
//   for (int i=0;i<n;i++) scale[i]=2*scale[i];

//   std::fstream monitorStream;
//   monitorStream.open("sst_giant.mon",std::fstream::out);
//   if ( !monitorStream.is_open() )
//   {
//      std::cout << " Failed to open monitorStream file sst_giant.mon" << std::endl;
//      exit(-1);
//   };

  NonlinearSolver nonlinearSolver;
  
  s = "names.nonlinType." + getParameter(pt, "solver.nonlinType", empty);
  NonlinearSolver::NonlinProblemType 
    nonlinType=static_cast<NonlinearSolver::NonlinProblemType>(getParameter(pt,s,(int) NonlinearSolver::NonlinProblemType::highlyNonlinear));

  int preconFillLevel=getParameter(pt,"preconlevel",0);
  nonlinearSolver.setTolerance(tol);
//   nonlinearSolver.setMaximumNoIterations(50);
  nonlinearSolver.setPreconFillLevel(preconFillLevel);
  nonlinearSolver.setSafetyFactor(rho);
//   nonlinearSolver.setErrorLevel(NonlinearSolver::verbose);
  nonlinearSolver.setMonitorLevel(NonlinearSolver::verbose);
  nonlinearSolver.setDataLevel(NonlinearSolver::verbose);
//  nonlinearSolver.setErrorStream(std::cout);
//  nonlinearSolver.setMonitorStream(std::cout);
//   nonlinearSolver.setMonitorStream(monitorStream);
   nonlinearSolver.setNonlinProblemType( nonlinType );
//   nonlinearSolver.setRestricted(false);
//  nonlinearSolver.setScalingVector(&scale);
  nonlinearSolver.setOutFilePrefix((std::string) "graph/sst_giant_");

  boost::timer::cpu_timer nonlinTimer;
  nonlinearSolver.giantGbit(gridManager,F,variableSet,&x,spaces);
  std::cout << "computing time for nonlinear solver: " << boost::timer::format(nonlinTimer.elapsed()) << "\n";
  struct NonlinearSolver::NleqInfo info=nonlinearSolver.getInfo();
//  monitorStream << " The total number of iterative linear solver steps done is " << 
  std::cout << " The total number of iterative linear solver steps done is " << 
               info.noOrdLinIt+info.noSimLinIt << std::endl;
  std::string finalOutputFilename;
  if ( info.returnCode==0 ) 
    finalOutputFilename="graph/sst_giant_solution";
  else 
    finalOutputFilename="graph/sst_giant_final";
  writeVTKFile(x,finalOutputFilename,IoOptions().setOrder(std::min(order,2)).setPrecision(7));

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End sst transfer tutorial program" << std::endl;
}
