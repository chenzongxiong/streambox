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

#include "io/vtk.hh"
//#include "io/amira.hh"
#include "fem/coarsening.hh"
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
    

    do {
      dt *= factor;
      dx = limex.step(x,dt,extrapolOrder,tolX);

      std::vector<std::pair<double,double> > errors = limex.estimateError(x,extrapolOrder,extrapolOrder-1);
      err = 0;
      for (int i=0; i<errors.size(); ++i)
        err = std::max(err, errors[i].first/(tolT[i].first+tolT[i].second*errors[i].second));
      
      factor = std::pow(0.5/err,1./(extrapolOrder+2));
      factor = std::max(0.5,std::min(factor,1.33));
      if ( verbosity > 0 ) 
      {
        std::cout << "t=" << eq.time() << ", dt=" << dt << ", err=" << err << ", factor=" << factor <<
                     ", red=" << red << ", dof=" << variableSet.degreesOfFreedom() << '\n';
      };
      ++red;
    } while(err>1 && red<=redMax);
    
    if ( verbosity>0 ) 
    {
      std::cout.flush();
      std::cout << "t= " << eq.time() << ", dt = " << dt << ", factor = " << factor << ", red=" << red << '\n';
    };

    // write linearly interpolated equidistant output
    assert(eq.time()<=outTime);
    while (outTime<=eq.time()+dt) {
      typename VariableSet::VariableSet z(x);
      z.axpy((outTime-eq.time())/dt,dx);
      // *out = z;
      // ++out;
      outTime += outInterval;
    }


    
    
    // step ahead
    x += dx;
    eq.time(eq.time()+dt);
    dt = std::min(dt*factor,dtMax);

    if (0) {
    // debugging output
    std::ostringstream fn;
    fn << "graph/outScaledGrid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << 2*steps;
    fn.flush();
    writeVTKFile(x,fn.str(),IoOptions());
    if ( verbosity>0 ) 
    {
      std::cout << variableSet.degreesOfFreedom() << " values written to " << fn.str() << std::endl;
    }
    }
    
    // perform mesh coarsening
    std::vector<bool> dummy(0);
    coarsening(variableSet,x,eq.scaling(),tolX,gridManager,dummy,verbosity);
    
    if (1) {
    // debugging output
    std::ostringstream fn;
    fn << "graph/outScaledGrid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << 2*steps+1;
    fn.flush();
    
    writeVTKFile(x,fn.str(),IoOptions());

    // output of solution for Amira visualization,
    // the data are written in binary format into file temperature.am,
    // possible is also ascii
//  IoOptions options;
//  options.outputType = IoOptions::ascii;
//      writeAMIRAFile(leafGridView,variableSet,x,fn.str(),options);

    if ( verbosity>0 ) 
    {
      std::cout << variableSet.degreesOfFreedom() << " values written to " << fn.str() << std::endl;
    }
    }
    
  }
  std::cout << '\n';

  limex.reportTime(std::cout);

  if (!done)
    std::cout << "*** maxSteps reached ***\n";
  
  
  return x;
}

#endif
