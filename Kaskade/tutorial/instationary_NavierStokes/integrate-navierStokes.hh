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

//#include "fem/coarsening.hh"
#include "io/vtk.hh"
#include "timestepping/limexWithoutJens.hh"
#include "utilities/enums.hh"

template <class Grid, class Equation, class VariableSet, class Spaces, class OutIter>
typename VariableSet::VariableSet  integrate(GridManager<Grid>& gridManager,
                                             Equation& eq, VariableSet const& variableSet, Spaces const& spaces,
                                             double dt, double dtMax,
                                             double T, int maxSteps, double rTolT, double aTolT,
                                             double rTolX, double aTolX, int extrapolOrder, 
                                             OutIter out, double outInterval,
                                             typename VariableSet::VariableSet x,
                                             DirectType directType,
                                             int verbosity = 0)
{
  using LeafGridView = typename Grid::LeafGridView;
  
  // write initial position
  *out = x;
  ++out;
  double outTime = eq.time()+outInterval;
  
  std::vector<std::pair<double,double> > tolX(variableSet.noOfVariables);
  std::vector<std::pair<double,double> > tolXC(variableSet.noOfVariables);

  for (int i=0; i<tolX.size(); ++i) {
    tolX[i] = std::make_pair(aTolX,rTolX);
    tolXC[i] = std::make_pair(aTolX/100,rTolX/100);
  }
  
  std::vector<std::pair<double,double> > tolT(variableSet.noOfVariables);
  for (int i=0; i<tolT.size(); ++i)
    tolT[i] = std::make_pair(aTolT,rTolT);
  
  Limex<Equation> limex(gridManager,eq,variableSet,directType,PrecondType::ILUK,verbosity);

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
    

    dt *= factor;
    dx = limex.step(x,dt,extrapolOrder,tolX);

/*
    do {
      dt *= factor;
      dx = limex.step(x,dt,extrapolOrder,tolX);
      std::cout << "limex.step is finished " << std::endl;

      std::vector<std::pair<double,double> > errors(0);
      std::vector<std::pair<double,double> > errors = limex.estimateError(x,extrapolOrder,extrapolOrder-1);
      err = 0;
      for (int i=0; i<errors.size(); ++i)
        err = std::max(err, errors[i].first/(tolT[i].first+tolT[i].second*errors[i].second));
      
      factor = std::pow(0.5/err,1./(extrapolOrder+2));
      factor = std::max(0.5,std::min(factor,1.33));
      if ( verbosity > 0 ) 
      {
        std::cout << eq.time() << ' ' << dt << ' ' << err << ' ' << red << ' '
        << variableSet.degreesOfFreedom() << ' ' << std::scientific << std::setprecision(3) <<
        boost::fusion::at_c<0>(x.data).value(samplePos)[0]  << 
        std::resetiosflags(std::ios::scientific) << std::setprecision(6) << '\n';
      };
      ++red;
    } while(err>1 && red<=redMax);
*/

    if ( verbosity>0 ) 
    {
      std::cout.flush();
      std::cout << "t= " << eq.time() << " dt = " << dt << " factor = " << factor << " red=" << red << '\n';
    };

    
    // step ahead
    x += dx;
    eq.time(eq.time()+dt);
    dt = std::min(dt*factor,dtMax);

    //L2Norm l2 ;
    //double rCnorm = l2( boost::fusion::at_c<1>(x.data) );
    //std::cout << "l2norm(p)^2 = " << rCnorm*rCnorm << std::endl;


    // controlling of range of p
    //  std::cout << "   normalizing pressure p" << std::endl << std::flush;
      double pMin = 1e+20;
      for (int j=0; j<boost::fusion::at_c<1>(x.data).space().degreesOfFreedom(); ++j) 
      {
        double p = (boost::fusion::at_c<1>(x.data)).coefficients()[j];
        if ( p < pMin ) pMin = p;
      }
      //std::cout << " pMin = " << pMin << std::endl;
      for (int j=0; j<boost::fusion::at_c<1>(x.data).space().degreesOfFreedom(); ++j) 
      {
        (boost::fusion::at_c<1>(x.data)).coefficients()[j] -= pMin;
      }
   
    // end of control
    

    if (0) {
      // debugging output
      auto fname = "graph/outScaledGrid"+paddedString(2*steps);
      writeVTKFile(x,fname);
      if ( verbosity>0 ) 
      {
        std::cout << variableSet.degreesOfFreedom() << " values written to " << fname << std::endl;
      }
    } 
    
    // perform mesh coarsening
    std::vector<bool> dummy(0);
//    coarsening(variableSet,x,eq.scaling(),tolX,gridManager,dummy,verbosity);
    
    if (1) {
      // debugging output
      auto fname = "graph/navierStokes-order2-pEx1-tau1e-2-Re100-"+paddedString(steps);
      writeVTKFile(x,fname,IoOptions().setOrder(1));
      
      if ( verbosity>0 ) 
        std::cout << "   " << variableSet.degreesOfFreedom() << " values written to " << fname << std::endl;
    }
    
  }
  std::cout << '\n';

  limex.reportTime(std::cout);

  if (!done)
    std::cout << "*** maxSteps reached ***\n";
  
  
  return x;
}

#endif
