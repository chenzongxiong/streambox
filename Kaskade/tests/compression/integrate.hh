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

#ifndef INTEGRATE_HH
#define INTEGRATE_HH

#include "utilities/enums.hh"

using namespace Kaskade;

template <class Grid, class Equation, class VariableSet, class Spaces, class OutIter, class Storage>
typename VariableSet::VariableSet  integrate(GridManager<Grid>& gridManager,
                                                Equation& eq, VariableSet const& variableSet, Spaces const& spaces,
                                                Grid const& grid, double dt, double dtMax,
                                                double T, int maxSteps, double rTolT, double aTolT,
                                                double rTolX, double aTolX, int extrapolOrder, 
                                                OutIter out, double outInterval,
						typename VariableSet::VariableSet x,
						DirectType directType,
						std::vector<double> &times,
						Storage& lossyStorage )
{
  // don't write initial position (initial value is known)
  //   *out = x;
  //   ++out;
  //   times.push_back(0) ;
  
  
  std::vector<std::pair<double,double> > tolX(variableSet.noOfVariables);
  std::vector<std::pair<double,double> > tolXC(variableSet.noOfVariables);

  for (int i=0; i<tolX.size(); ++i) {
    tolX[i] = std::make_pair(aTolX,rTolX);
    tolXC[i] = std::make_pair(aTolX/100,rTolX/100);
  }
  
  // for testing: no adaptivity in space
  tolX.clear() ; tolXC.clear() ; 
  
  std::vector<std::pair<double,double> > tolT(variableSet.noOfVariables);
  for (int i=0; i<tolT.size(); ++i)
    tolT[i] = std::make_pair(aTolT,rTolT);
  
  Limex<Equation> limex(gridManager,eq,variableSet,directType);

  Dune::FieldVector<double,2> samplePos(0);
  
  int steps;
  bool done = false;
  for (steps=0; !done && steps<maxSteps; ++steps) {
    if (eq.time()>T-1.1*dt) {
      dt = T-eq.time();
      done = true;
    }
    
    typename VariableSet::VariableSet dx(x);

    int redMax = 5;
    int red = 0;
    double factor = 1.0;
    double err;
    

    do {
      dt *= factor;
      dx = limex.step(x,dt,extrapolOrder,tolX);

      std::vector<std::pair<double,double> > errors = limex.estimateError(x,extrapolOrder,extrapolOrder-1);
      err = 0;
      for (int i=0; i<errors.size(); ++i)
        err = std::max(err, errors[i].first/(tolT[i].first+tolT[i].second*errors[i].second));
      
      factor = std::pow(0.5/err,1./(extrapolOrder+2));
      factor = std::max(0.5,std::min(factor,1.33));
      std::cout << eq.time() << ' ' << dt << ' ' << err << ' ' << red << ' ' << variableSet.degreesOfFreedom()
                << ' ' << boost::fusion::at_c<0>(x.data).value(samplePos)[0] << '\n';
      ++red;
    } while(err>1 && red<=redMax);
    
    std::cerr << "t= " << eq.time() << " dt = " << dt << " factor = " << factor << " red=" << red << '\n';
 
    
    // step ahead
    x += dx;
    eq.time(eq.time()+dt);
    dt = std::min(dt*factor,dtMax);

    *out = x ;
    ++out;
    times.push_back(eq.time());

    // store to disk
    std::ostringstream predfn ; 
    predfn.width(3);  predfn.fill('0');
    predfn.setf(std::ios_base::right,std::ios_base::adjustfield);
    predfn << "graph/quant" ;
    predfn << steps ;
    predfn.flush();
    
    lossyStorage.encode(x, predfn.str());
       
    if (1) 
    {
      // debugging output
      std::ostringstream fn;
      fn << "graph/outScaledGrid";
      fn.width(3);
      fn.fill('0');
      fn.setf(std::ios_base::right,std::ios_base::adjustfield);
      fn << steps ;
      fn.flush();
      typedef typename Grid::LeafGridView LeafGridView;
      LeafGridView leafView = gridManager.grid().leafView();
      writeVTKFile(leafView,variableSet,x,fn.str());
      std::cout << variableSet.degreesOfFreedom() << " values written to " << fn.str() << "\n\n";
    }    
  }
  
  std::cout << '\n';

  limex.reportTime(std::cerr);

  if (!done)
    std::cerr << "*** maxSteps reached ***\n";
  
  lossyStorage.flush();
  lossyStorage.finish();
  
  return x;
}

#endif
