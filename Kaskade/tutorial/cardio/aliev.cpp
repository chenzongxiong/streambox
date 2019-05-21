#include <stdio.h>

#include <complex>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <utility>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"

#include "dune/common/stdstreams.hh"
#include "dune/grid/sgrid.hh"
#include "dune/grid/uggrid.hh"
#include "dune/grid/io/file/amirameshreader.hh"
#include "dune/grid/io/file/amirameshwriter.hh"
#include "dune/grid/common/gridinfo.hh"
#include "dune/istl/solvers.hh"

#include "utilities/enums.hh"
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
#include "io/amira.hh"

#include "utilities/kaskopt.hh"
#include "amiramesh.hh"

using namespace Kaskade;

#include "integrate.hh"
#include "aliev.hh"


boost::timer::cpu_timer partTimer;
double assTime, solveTime;


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
	double R=1.0;		// slab: R=0.09		heart: R=1.0
	double beta=1000.0;
	double f=0.0;
	double r;
	
	//r = sqrt( (x[0]-(0.3))*(x[0]-(0.3)) + (x[1]-(-2.4))*(x[1]-(-2.4)) + (x[2]-182.5)*(x[2]-182.5) );	// heart smooth
	r = sqrt( (x[0]-(-0.9))*(x[0]-(-0.9)) + (x[1]-1.7)*(x[1]-1.7) + (x[2]-1.1)*(x[2]-1.1) );			// heart glenn
	//r = sqrt( (x[0]-1.5)*(x[0]-1.5) + (x[1]-1.5)*(x[1]-1.5) + (x[2]-1.1)*(x[2]-1.1) );				// slab
	beta = -1e4*log(0.001);
	if (r < R) f = 1;
	else f = exp(-beta*(R-r)*(R-r));
	
    if (component==0) 
      return 0.0;		// 0.95*f;
    else if (component==1) 
      return 0.0;
    else
      assert("wrong index!\n"==0);
    return 0;
  }

private:
  int component;
};




template <class Eq>
typename Eq::AnsatzVars::VariableSet transformBack(Eq const& eq, typename Eq::AnsatzVars::VariableSet const& vars) 
{
  typename Eq::AnsatzVars::VariableSet  x(vars);
  
  for (int j=0; j<boost::fusion::at_c<0>(x.data).space().degreesOfFreedom(); ++j) {
    (boost::fusion::at_c<0>(x.data)).coefficients()[j] = eq.theta((boost::fusion::at_c<0>(x.data)).coefficients()[j]);
  }
  
  return x;
}

template <class Eq>
typename Eq::AnsatzVars::VariableSet transformBackDeriv(Eq const& eq, typename Eq::AnsatzVars::VariableSet const& vars) 
{
  typename Eq::AnsatzVars::VariableSet  x(vars);
  
  for (int j=0; j<boost::fusion::at_c<0>(x.data).space().degreesOfFreedom(); ++j) {
    (*boost::fusion::at_c<0>(x.data))[j] = eq.dtheta((*boost::fusion::at_c<0>(x.data))[j]);
  }
  
  return x;
}

template <class Eq>
typename Eq::AnsatzVars::VariableSet transformBackDerivMul(Eq const& eq, typename Eq::AnsatzVars::VariableSet const& vars1, typename Eq::AnsatzVars::VariableSet const& vars2) 
{
  typename Eq::AnsatzVars::VariableSet  x(vars1);
  
  for (int j=0; j<boost::fusion::at_c<0>(x.data).space().degreesOfFreedom(); ++j) {
    (*boost::fusion::at_c<0>(x.data))[j] = eq.dtheta((*boost::fusion::at_c<0>(x.data))[j]) * (*boost::fusion::at_c<0>(vars2.data))[j];
  }
  
  return x;
}

template <class Eq>
typename Eq::AnsatzVars::VariableSet transformForward(Eq const& eq, typename Eq::AnsatzVars::VariableSet const& vars) 
{
  typename Eq::AnsatzVars::VariableSet  x(vars);
  
  for (int j=0; j<boost::fusion::at_c<0>(x.data).space().degreesOfFreedom(); ++j) {
    (*boost::fusion::at_c<0>(x.data))[j] = eq.zeta((*boost::fusion::at_c<0>(x.data))[j]);
  }
  
  return x;
}

template <class Eq>
typename Eq::AnsatzVars::VariableSet transformInitial(Eq const& eq, typename Eq::AnsatzVars::VariableSet const& vars) 
{
  typename Eq::AnsatzVars::VariableSet  x(vars);
  
  for (int j=0; j<boost::fusion::at_c<0>(x.data).space().degreesOfFreedom(); ++j) {
    (*boost::fusion::at_c<0>(x.data))[j] = eq.zetaInitial((*boost::fusion::at_c<0>(x.data))[j]);
  }
  
  return x;
}



int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  std::cout << "Start cardiac simulation (Aliev-Panfilov)" << std::endl;

  fexcept_t flag;
  fegetexceptflag(&flag,FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
  fesetexceptflag(&flag,FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);  
  
  boost::timer::cpu_timer totalTimer;

  int const dim = 3;
  int refinements;
  int  order, extrapolOrder, maxSteps;
  double dt, maxDT, T, rTolT, aTolT, rTolX, aTolX, writeInterval;
  std::string logname = "cardio.log";
  bool refineTimeStep;
  std::string zetaFile, thetaFile;
  double xi, a, gs, ga, mu1, mu2, eps1, kappa, v_rest, v_peak;
  	
  std::string empty, geoFile("slab242.am");
  geoFile = getParameter(pt, "file", geoFile);


  refinements = 0;   // slab: 3,  Herz: 0
  order = 2;
  extrapolOrder = 2;
  T = 1.0;
  dt = 0.01;
  maxDT = 4.0;
  writeInterval = 1.0;
  refineTimeStep = true;
  rTolT = 1e-3;
  aTolT = 1e-3;
  rTolX = 1e-3;
  aTolX = 1e-3;
  maxSteps = 10000;
// 	thetaFile = "zetaLinearInv.gnu";
// 	zetaFile  = "zetaLinear.gnu";
// 	thetaFile = "zetaInv.gnu";     			// original nonlinear scaling
// 	zetaFile  = "zeta.gnu";        			// original nonlinear scaling
// 	thetaFile = "theta-tanh5-100000.gnu";   // tanh nonlinear scaling
// 	zetaFile  = "zeta-tanh5-100000.gnu";    // tanh nonlinear scaling
	thetaFile = "theta-tanh2-100000.gnu";   // tanh nonlinear scaling
	zetaFile  = "zeta-tanh2-100000.gnu";    // tanh nonlinear scaling
// 	thetaFile = "zetaInv-MW.gnu";     		// MW nonlinear scaling
// 	zetaFile  = "zeta-MW.gnu";        		// MW nonlinear scaling
// 	thetaFile = "theta-2d-1e-5.gnu";		// MW nonlinear scaling
// 	zetaFile  = "zeta-2d-1e-5.gnu";			// MW nonlinear scaling
//	thetaFile = "theta-2d-1e-7.gnu";		// MW nonlinear scaling
//	zetaFile  = "zeta-2d-1e-7.gnu";			// MW nonlinear scaling
// 	thetaFile = "zetaInvQuadratic.gnu";     // Quadratic nonlinear scaling
// 	zetaFile  = "zetaQuadratic.gnu";        // Quadratic nonlinear scaling
// 	thetaFile = "zetaIdentityInv.gnu";
// 	zetaFile  = "zetaIdentity.gnu";
// 	thetaFile = "";
// 	zetaFile  = "";
	
	
  xi   = 1.0;
  kappa = 0.001;
  a = 0.1;
  gs = 8.0;
  ga = 8.0;
  mu1 = 0.07;
  mu2 = 0.3;
  eps1 = 0.01;
	
  v_rest = 0.0;   //-85.0;
  v_peak = 1.0;   //35.0;
	

  int heapSize=16384;
  typedef Dune::UGGrid<dim> Grid;
  Grid::setDefaultHeapSize(heapSize);
  Dune::GridFactory<Grid> factory;

  DuneAmiraMesh<dim,dim,Grid> mesh(geoFile.c_str());
  mesh.InsertUGGrid(factory);
      
  std::unique_ptr<Grid> grid( factory.createGrid() );
  for (int k=0; k<refinements; k++)
  {
    grid->globalRefine(1);
  }

  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " tetrahedra, " << std::endl;
  std::cout << "      " << grid->size(1) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(dim-1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(dim) << " points" << std::endl;
  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
  gridManager.enforceConcurrentReads(true);
    
  typedef Grid::LeafGridView LeafView;
		
  	// construct involved spaces.
    //typedef ContinuousHierarchicMapper<double,Grid>  LMapper;
  typedef ContinuousLagrangeMapper<double,LeafView> LMapper;  
  typedef FEFunctionSpace<LMapper> H1Space;
	
  H1Space h1Space(gridManager,gridManager.grid().leafView(),order);
	
  	typedef boost::fusion::vector<H1Space const*> Spaces;
  	Spaces spaces(&h1Space);
	
  	// construct variable list.
  	typedef boost::fusion::vector<VariableDescription<0,1,0>,
    							  VariableDescription<0,1,1> > VariableDescriptions;
  	std::string varNames[2] = { "u", "v" };

  	typedef VariableSetDescription<Spaces,VariableDescriptions> VariableSet;
  	VariableSet variableSet(spaces,varNames);

  	typedef AlievPanfilovEquation<double,VariableSet> Equation;
  	Equation Eq(xi,kappa,a,gs,ga,mu1,mu2,eps1,v_rest,v_peak,zetaFile,thetaFile);
  	// Eq.check();
  

  std::vector<VariableSet::VariableSet> solutions;
  std::cout << "solutions.size = " << solutions.size() << std::endl;
  // std::vector<VariableSet::VariableSet> devnull;

  Eq.time(0);
  VariableSet::VariableSet x(variableSet);
  Eq.scaleInitialValue<0>(InitialValue(0),x);
  Eq.scaleInitialValue<1>(InitialValue(1),x);

  T     = 1.0;      	
  dt    = 1e-3;			//0.001;
  maxDT = 1.0;
  rTolT = 5e-3;		// slab: 1e-2; herz(Glenn): 5e-3;		//1e-5/8.0;
  aTolT = 5e-3;		// slab: 1e-2; herz(Glenn): 5e-3;		//1e-5/8.0;
  rTolX = 1e-2;		// slab: 1e-2; herz(Glenn): 5e-3;		//1e-5*2.0;
  aTolX = 1e-2;		// slab: 1e-2; herz(Glenn): 5e-3;		//1e-5*2.0;
  maxSteps = 10000;
  x = integrate(gridManager,Eq,variableSet,spaces,gridManager.grid(),
                dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::MUMPS);

  T     = 225.0;      	
  dt    = 5e-3;		//0.001;
  rTolT = 2e-3;		//1e-5/8.0;
  aTolT = 2e-3;		//1e-5/8.0;
  rTolX = 1.0e-2;   //1.5e-2;		//1e-5*2.0;
  aTolX = 1.0e-2;   //1.5e-2;		//1e-5*2.0;
  maxSteps = 10000;
  x = integrate(gridManager,Eq,variableSet,spaces,gridManager.grid(),
                dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::MUMPS);

  T     = 226.0;      	
  dt    = 1e-3;			//0.001;
  rTolT = 2e-3;		//1e-5/8.0;
  aTolT = 2e-3;		//1e-5/8.0;
  rTolX = 1.5e-2;		//1e-5*2.0;
  aTolX = 1.5e-2;		//1e-5*2.0;
  maxSteps = 10000;
  x = integrate(gridManager,Eq,variableSet,spaces,gridManager.grid(),
                dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::MUMPS);

  T     = 1000.0;      	
  dt    = 1e-3;		//0.001;
  rTolT = 2e-3;		//1e-5/8.0;
  aTolT = 2e-3;		//1e-5/8.0;
  rTolX = 1.5e-2;		//1e-5*2.0;
  aTolX = 1.5e-2;		//1e-5*2.0;
  maxSteps = 10000;
  x = integrate(gridManager,Eq,variableSet,spaces,gridManager.grid(),
                dt,maxDT,T,maxSteps,rTolT,aTolT,rTolX,aTolX,extrapolOrder,
                std::back_inserter(solutions),writeInterval,x,DirectType::MUMPS);

  std::cout << " End of integrate" << std::endl;

//   if (order>1)
//     gridManager.globalRefine(order);
  
//   for (int i=0; i<solutions.size(); ++i) {
//     //VariableSet::VariableSet x(transformBack(Eq,solutions[i]));
//     VariableSet::VariableSet xT(transformBack(Eq,x));
//     
//     std::ostringstream fn;
//     fn << "graph/out-rescaled";
//     fn.width(3);
//     fn.fill('0');
//     fn.setf(std::ios_base::right,std::ios_base::adjustfield);
//     fn << i;
//     fn.flush();
//     writeVTKFile(variableSet,xT,fn.str());
//   }
  
//   
//   for (int i=0; i<solutions.size(); ++i) {
//     std::ostringstream fn;
//     fn << "graph/outScaled";
//     fn.width(3);
//     fn.fill('0');
//     fn.setf(std::ios_base::right,std::ios_base::adjustfield);
//     fn << i;
//     fn.flush();
//     writeVTKFile(variableSet,solutions[i],fn.str());
//   }
  
  std::cout << "End Aliev-Panfilov integration" << std::endl;

  return 0;
}
