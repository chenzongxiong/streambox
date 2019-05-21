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

template <class Grid, class Equation, class VariableSet, class Spaces, class OutIter>
typename VariableSet::VariableSet  integrate(GridManager<Grid>& gridManager,
                                                Equation& eq, VariableSet const& variableSet, Spaces const& spaces,
                                                Grid const& grid, double dt, double dtMax,
                                                double T, int maxSteps, double rTolT, double aTolT,
                                                double rTolX, double aTolX, int extrapolOrder, 
                                                OutIter out, double outInterval,
                                                typename VariableSet::VariableSet x,
                                                DirectType directType,
                                                int verbosity = 1)
{
  // write initial position
  *out = x;
  ++out;
  double outTime = eq.time()+outInterval;
  
  
  std::vector<std::pair<double,double> > tolX(variableSet.noOfVariables);
  std::vector<std::pair<double,double> > tolXC(variableSet.noOfVariables);
  std::vector<std::pair<double,double> > tolXdt(variableSet.noOfVariables);

  for (int i=0; i<tolX.size(); ++i) {
    tolX[i] = std::make_pair(aTolX,rTolX);
    tolXC[i] = std::make_pair(aTolX/100,rTolX/100);
  }

  std::cout << "aTolX = " << aTolX << ", rTolX = " << rTolX << std::endl;
  for (int i=0; i<tolX.size(); ++i)
      std::cout << "integrate: tolX = (" << tolX[i].first <<", " << tolX[i].second <<")" << std::endl;
  
  std::vector<std::pair<double,double> > tolT(variableSet.noOfVariables);
  for (int i=0; i<tolT.size(); ++i)
    tolT[i] = std::make_pair(aTolT,rTolT);
  
  Limex<Equation> limex(gridManager,eq,variableSet,directType,PrecondType::ILUK,verbosity);

	// heart Glenn
   Dune::FieldVector<double,3> samplePos(0.648);
   //samplePos[0] =  0.648;  
   //samplePos[1] =  0.521;  
   //samplePos[2] =  1.000;  
   samplePos[0] =  -0.9;  
   samplePos[1] =   3.0;  
   samplePos[2] =   1.1;  

	// heart smooth
//   Dune::FieldVector<double,3> samplePos(0.648);
//   samplePos[0] =  0.3;  
//   samplePos[1] =  -2.4;  
//   samplePos[2] =  182.5;  

  	 // slab
//   Dune::FieldVector<double,3> samplePos(1.9);
//   samplePos[0] =  1.9;  
//   samplePos[1] =  1.9;  
//   samplePos[2] =  1.1;  

  std::cout << "samplePos = " << samplePos[0] << "," << samplePos[1] << "," << samplePos[2] << std::endl;
  
  static int steps = -1;
  bool done = false;
  
  //for (steps=0; !done && steps<maxSteps; ++steps) {
  do {
  	steps++;

    if (eq.time()>T-1.1*dt) {
      dt = T-eq.time();
      done = true;
    }
    
    typename VariableSet::VariableSet dx(x);

    int redMax = 5;
    int red = 0;
    double factor = 1.0;
    double err;
    
	std::cout << "Start limex loop" << std::endl;
    do {
      dt *= factor;

	  // Kopplung der absoluten Ortsgenauigkeit an den Zeitschritt
  	  for (int i=0; i<tolX.size(); ++i) {
    		tolXdt[i] = std::make_pair(aTolX*dt*sqrt(0.3),rTolX);
  	  	}
      //dx = limex.step(x,dt,extrapolOrder,tolXdt);
      dx = limex.step(x,dt,extrapolOrder,tolX);
	  
      std::vector<std::pair<double,double> > errors = limex.estimateError(x,extrapolOrder,extrapolOrder-1);
      err = 0;
      for (int i=0; i<errors.size(); ++i)
        err = std::max(err, errors[i].first/(tolT[i].first+tolT[i].second*errors[i].second));
      
      factor = std::pow(0.5/err,1./(extrapolOrder+2));
      factor = std::max(0.5,std::min(factor,1.33));
      std::cout << "   limex loop: t=" << eq.time() << " dt=" << dt << " err=" << err << ' ' << red << " dof=" << variableSet.degreesOfFreedom()
                << " vertices=" << gridManager.grid().size(3) << " val(sample)=" << boost::fusion::at_c<0>(x.data).value(samplePos)[0] << '\n';
      ++red;
    } while(err>1 && red<=redMax);
    
    std::cout.flush();
	std::cout << "End limex loop" << std::endl;
    std::cout << "t= " << eq.time() << " dt = " << dt << " factor = " << factor << " red=" << red << '\n';

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
    
	double val_at_samplePos = boost::fusion::at_c<0>(x.data).value(samplePos)[0];
	double trans_val_at_samplePos = eq.theta(val_at_samplePos);
	std::cout << "val_at_samplePos: " << eq.time() << "   " << val_at_samplePos << std::endl;
	std::cout << "rescaled val_at_samplePos: " << eq.time() << "   " << trans_val_at_samplePos << std::endl;
    std::cout << "   step: t=" << eq.time() << " dt=" << dt << " err=" << err << ' '  << " dof=" << variableSet.degreesOfFreedom()
              << " vertices=" << gridManager.grid().size(3) << " rescaled_val(sample_t)=" << trans_val_at_samplePos << '\n';
    
    dt = std::min(dt*factor,dtMax);

    if (1) {
    // debugging output
    std::ostringstream fn;
    fn << "graph/outScaledGrid-1e-2-";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << steps;
    fn.flush();
    typedef typename Grid::LeafGridView LeafGridView;
    LeafGridView leafView = gridManager.grid().leafView();
    writeVTKFile(leafView,x,fn.str());
    if ( verbosity>0 ) 
    {
      std::cout << variableSet.degreesOfFreedom() << " values written to " << fn.str() << std::endl;
    }
    }
    
    // perform mesh coarsening
    std::vector<bool> dummy(0);
    coarsening(variableSet,x,eq.scaling(),tolX,gridManager,dummy,verbosity);
    
    
  }   while ( !done && (steps<(maxSteps-1)) );

  std::cout << '\n';

  limex.reportTime(std::cout);

  if (!done)
    std::cout << "*** maxSteps reached ***\n";
  
  
  return x;
}

#endif
