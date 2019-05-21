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

#include <complex>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <utility> // std::move

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"

#include "dune/common/stdstreams.hh"
#include "dune/grid/sgrid.hh"
#include "dune/grid/uggrid.hh"
#include "dune/grid/common/gridinfo.hh"

#include "fem/assemble.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/functional_aux.hh"
#include "fem/coarsening.hh"
#include "fem/norms.hh"

#include "linalg/umfpack_solve.hh"

#include "timestepping/limexWithoutJens.hh"
#include "timestepping/semieuler.hh"
#include "timestepping/extrapolation.hh"

#include "io/vtk.hh"

#include "integrate.hh"
#include "movingsource.hh"

//#include "io/lossystorageDUNE.hh"
#include "io/lossystorageWithoutTemporalPred.hh"
// #include "io/lossystorage.hh"

using namespace Kaskade;

template <class Grid>
std::unique_ptr<Grid> RefineGrid(int refinements, int heapSize=500)
{
  Grid::setDefaultHeapSize(heapSize);
  Dune::GridFactory<Grid> factory;
  Dune::FieldVector<double,2> v;
  // vertices
  v[0]=0; v[1]=0;
  factory.insertVertex(v);
  v[0]=1; v[1]=0;
  factory.insertVertex(v);
  v[0]=1; v[1]=1;
  factory.insertVertex(v);
  v[0]=0; v[1]=1;
  factory.insertVertex(v);
  // elements
  std::vector<unsigned int> vid(3);
  vid[0]=0; vid[1]=1; vid[2]=2;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=2; vid[2]=3;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  grid->globalRefine(refinements);  
  return grid;
}


struct InitialValue 
{
  typedef double Scalar;
  static int const components = 1;
  typedef Dune::FieldVector<Scalar,components> ValueType;

  InitialValue(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::ctype,Cell::dimension> const& localCoordinate) const 
  {
    Dune::FieldVector<typename Cell::ctype,Cell::dimensionworld> x = cell.geometry().global(localCoordinate);

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
    int refinements = 7, order = 1, extrapolOrder = 1, maxSteps = 10, heapSize = 1000 ;
    double dt = 0.1, maxDT = 1.0, T = 10.0, rTolT = 1.0e-2, aTolT = 1.0e-2, rTolX = 1.0e-4, aTolX = 1.0e-4, writeInterval = 1.0;

    typedef Dune::UGGrid<dim> Grid;
    std::unique_ptr<Grid> grid( RefineGrid<Grid>(refinements,heapSize) );
    
    std::cout << "Grid: " << grid->size(0) << " " << grid->size(1) << " " << grid->size(2) << std::endl << std::endl;

    GridManager<Grid> gridManager(std::move(grid));
  
    // construct involved spaces and define equation
    typedef Grid::LeafGridView LeafView;
    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> > H1Space;

    H1Space temperatureSpace(gridManager,gridManager.grid().leafView(),order);

    typedef boost::fusion::vector<H1Space const*> Spaces;
    Spaces spaces(&temperatureSpace);

    typedef boost::fusion::vector<VariableDescription<0,1,0> > VariableDescriptions;
    std::string varNames[1] = { "u" };
  
    typedef VariableSetDescription<Spaces,VariableDescriptions> VariableSet;
    VariableSet variableSet(spaces,varNames);

    typedef MovingSourceEquation<double,VariableSet> Equation;
    Equation Eq;

    std::vector<VariableSet::VariableSet> solutions;
    
    // prepare lossy storage
    std::vector<double> times ;
    
    int coarseLevel = 0 ;
    double qTol = 1e-5 ;   // quantization error tolerance
    
    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LevelGridView> > HierarchicH1Space;
    typedef UniformQuantizationPolicy<VariableSet,Grid> QuantizationPolicy ;
    typedef LossyStorage<Grid,VariableSet,HierarchicH1Space,QuantizationPolicy> Storage ;
    
    QuantizationPolicy quantizationPolicy ;
    Storage lossyStorage( gridManager, variableSet, coarseLevel, qTol, true, quantizationPolicy ) ;  

    // integrate equation
    Eq.time(0);
    VariableSet::VariableSet x(variableSet);
    Eq.scaleInitialValue<0>(InitialValue(0),x);
    
    x = integrate(gridManager,Eq,variableSet,spaces,gridManager.grid(),
		  dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
		  std::back_inserter(solutions),writeInterval,x,DirectType::UMFPACK,times,lossyStorage);

    // decode files & check error 
    // decode backwards in time (necessary if temporal prediction/differential
    // encoding is used)
    std::cout << "\n\nckecking error:\n" ;
    
    VariableSet::VariableSet state_data(variableSet) ;
    L2Norm l2 ; 
    
    std::ostringstream fn ; 
    fn.width(3);  fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);  
    for( int t = times.size()-1 ; t >= 0 ; t-- )
    {
      fn << "graph/quant" ;
      fn << t ;
      fn.flush();
      lossyStorage.decode( gridManager, state_data, fn.str() ) ;
      
      // write reconstructed states as vtu-files
//       writeVTKFile(gridManager.grid().leafView(), variableSet, state_data, fn.str());
      
      fn.clear() ; fn.str("");     
      
      // compute and print L^\infty error 
      state_data -= solutions[t] ;
      std::vector<double> foo( gridManager.grid().size(dim) ) ;
      state_data.write( foo.begin() ) ;
      double absQuantErr = fabs(*std::max_element(foo.begin(),foo.end(),abscompare)) ;
      std::cout << "t = " << times[t] << "\tL^\\infty = " << absQuantErr ;  
      absQuantErr = sqrt(l2.square( boost::fusion::at_c<0>(state_data.data) ));
      std::cout << "\tL^2 = " << absQuantErr << "\n" ;
    }
    
    return 0;
    
    
//     std::cout << times.size() << " timesteps.\n" ;
//     for( int i = 0 ; i < times.size() ; i++ ) std::cout << times[i] << "   " ;
//     std::cout << "\n\n" ;
//     std::cout << "Solution vector contains " << solutions.size() << " entries.\n" ;
//     
//     int windowSize = 2 ; 
//     
//     for( int i = 0 ; i < times.size()-windowSize ; i++ )
//     {
//     
//       VariableSet::VariableSet a = solutions[i+1] , b = solutions[i+0] ;
//       a -= solutions[i+0] ; a *= 1.0/(times[i+1]-times[i+0]) ;
// 
//       VariableSet::VariableSet pred = b ;
//       pred.axpy(times[i+2]-times[i+0],a);
//     
//       typedef Grid::LeafGridView LeafGridView;
//       LeafGridView leafView = gridManager.grid().leafView();
//       std::ostringstream fn ; 
//       fn.width(3);  fn.fill('0');
//       fn.setf(std::ios_base::right,std::ios_base::adjustfield);
//       fn << "graph/pred" ;
//       fn << i+2 ;
//       fn.flush();
//       writeVTKFile(leafView,variableSet,pred,fn.str() );
//       
//       lossyStorage.encode( pred, fn.str() ) ;
//     
//       fn.str("") ; 
//       fn << "graph/diff" ;
//       fn << i+2 ;
//       fn.flush();
//       pred -= solutions[i+2] ; 
//       writeVTKFile(leafView,variableSet,pred,fn.str());
//       
//       lossyStorage.encode( pred, fn.str() ) ;
//     }
  }
