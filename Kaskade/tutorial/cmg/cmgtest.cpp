/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/istl/preconditioners.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/lagrangespace.hh"
#include "fem/hierarchicspace.hh"
#include "fem/embedded_errorest.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/iterate_grid.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/mumps_solve.hh"

#include "io/vtk.hh"
#include "io/amira.hh"

#include "mg/pcg.hh"
#include "linalg/apcg.hh"

#include "utilities/enums.hh"
#include "utilities/kaskopt.hh"
#include "utilities/geometric_sequence.hh"

using namespace Kaskade;
#include "poisson.hh"
#include "createGrid.hh"

int problemNo = 1;


bool compareAbs(const double x1, const double  x2)
  {
    return fabs(x1)>fabs(x2);
  }
  
// Implements a test problem for the cascadic multigrid algorithm as of 
// Deuflhard/Weiser chapter 7.

int main(int argc, char *argv[])
{
  int verbosity = 1;
  bool dump = false;
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosity, dump);
  
  std::cout << "Start cascadic multigrid test program" << std::endl;
  
  int  direct, refs, order, onlyLowerTriangle = false, maxAdaptSteps, gebiet;
  DirectType directType;
  IterateType iterateType = IterateType::PCG;
  MatrixProperties property;
  std::string empty, problem, geometry, functional;
  
  problem = getParameter(pt, "problem", empty);
  refs = getParameter(pt, problem+".refs", 0),
  order =  getParameter(pt, problem+".order", 1),
  maxAdaptSteps = getParameter(pt, problem+".maxAdaptSteps", 10);
  geometry = getParameter(pt, problem+".geometry", empty);
  functional = getParameter(pt, problem+".functional", empty);
  gebiet = getParameter(pt, "names.geometry."+geometry, 1);
  problemNo = getParameter(pt, "names.functional."+functional, 1);
  std::cerr << "selected problem " << problem << ", functional=" <<
  functional << "(" << problemNo <<  "), geometry=" <<
  geometry << "(" << gebiet << ")" << std::endl;
  
  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  std::cerr << "s = " << s << std::endl;
  direct = getParameter(pt, s, 0);
  std::cerr << "s = " << s << std::endl;
  
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  std::cerr << "s = " << s << std::endl;
  directType = static_cast<DirectType>(getParameter(pt, s, 2));
  std::cerr << "s = " << s << std::endl;
  
  s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
  iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  
  property = MatrixProperties::POSITIVEDEFINITE;
  
  if (((property == MatrixProperties::SYMMETRIC)||(property == MatrixProperties::POSITIVEDEFINITE))&&
    ((directType == DirectType::MUMPS)||(directType == DirectType::PARDISO)))
  {
    onlyLowerTriangle = true;
  }
  
  int const dim2=2;
  typedef Dune::UGGrid<dim2> Grid;
  Dune::GridFactory<Grid> factory;
  
  switch (gebiet)
  {
    case 0:
    {
      createUnitSquare<Grid,dim2>(factory);
    }
    break;
    case 1:
    {
      createUnitCrack<Grid,dim2>(factory);
    }
    break;
    case 2:
    {
      createOriginalBoDD<Grid,dim2>(factory);
    }
    break;
    case 3:
    {
      createAreaGrid<Grid,dim2>(factory);
    }
    break;
    default:
      std::cout << "Unknown gebiet" << std::endl;
      exit(2);
  }
  
  std::auto_ptr<Grid> grid( factory.createGrid() );
  // the coarse grid will be refined refs times
  grid->globalRefine(refs);
  
  // some information on the refined mesh
  std::cout << "Grid: " << grid->size(0) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(2) << " points" << std::endl;
  
  // a gridmanager is constructed 
  GridManager<Grid> gridManager(std::move(grid));    
  
  typedef Grid::LeafGridView LeafView;
  // construction of finite element space for the scalar solution T
  typedef FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> > H1Space;
  
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
  typedef boost::fusion::vector<H1Space const*> Spaces;
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  typedef boost::fusion::vector<VariableDescription<0,1,0> >
  VariableDescriptions;
  std::string varNames[1] = { "T" };
  typedef VariableSetDescription<Spaces,VariableDescriptions> VariableSet;
  VariableSet variableSet(spaces,varNames);
  
  typedef VariableSet::CoefficientVectorRepresentation<>::type CoefficientVector;
  
  // Define a higher order space for transfering the error estimate
  H1Space temperatureSpace2(gridManager,gridManager.grid().leafView(),order+1);
  
  // Define the variational functional
  typedef PoissonFunctional<double,VariableSet> Functional;
  Functional F;
  typedef VariationalFunctionalAssembler<LinearizationAt<Functional> > Assembler;
  typedef Assembler::RhsArray Rhs;
  Assembler assembler(gridManager,spaces);
  
  typedef FEFunctionSpace<ContinuousHierarchicExtensionMapper<double,LeafView> > H1ExSpace;
  H1ExSpace spaceEx(gridManager,gridManager.grid().leafView(), order+1);
  
  typedef boost::fusion::vector<H1Space const*,H1ExSpace const*> H1ExSpaces;
  H1ExSpaces exSpaces(&temperatureSpace,&spaceEx);
  
  typedef boost::fusion::vector<VariableDescription<1,1,0> > ExVariableDescriptions;
  typedef VariableSetDescription<H1ExSpaces,ExVariableDescriptions> ExVariableSet;
  std::string exVarNames[1] = { "e"};
  ExVariableSet exVariableSet(exSpaces, exVarNames);
  
  typedef ExVariableSet::CoefficientVectorRepresentation<>::type ExCoefficientVector;
  
  typedef HierarchicErrorEstimator<LinearizationAt<Functional>,ExVariableSet> ErrorEstimator;
  typedef VariationalFunctionalAssembler<ErrorEstimator> EstGOP;
  
  EstGOP estGop(gridManager,exSpaces);
  
  typedef VariableSet::Grid::Traits::LeafIndexSet IS ;
  IS const& is = gridManager.grid().leafIndexSet();
  
  VariableSet::VariableSet x(variableSet), dx(variableSet);
  
  double const TOL = getParameter(pt, problem+".TOL", 1.0-2),
  minRefine = getParameter(pt, problem+".minRefine", 0.0);
  std::cerr << "TOL = " << TOL << ", minRefine = " << minRefine << std::endl ;
  assert(TOL>0);
  
  IoOptions options;
  options.outputType = IoOptions::ascii;
  
  constexpr int neq = Functional::TestVars::noOfVariables;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  int refSteps = 0;
  bool accurate = false;
  int iteSteps = getParameter(pt, "solver.iteMax", 1000);
  int verbose = getParameter(pt, "solver.verbose", 1);
  double iteEps = getParameter(pt, "solver.iteEps", 1.0e-9);
  int lookahead = getParameter(pt,"solver.APCG.lookahead",6);
  double eefactor = getParameter(pt,"ccg.eefactor",0.0);
  Dune::InverseOperatorResult result;
  double errNorm = -1;
  
  size_t size = variableSet.degreesOfFreedom(0,1);
  
  double const gamma = iterateType==IterateType::APCG? 1.0: 0.5;   // Deuflhard/Weiser (7.60)
  int const d = 2;                                    // two dimensional problem
  assert(d*gamma>1);                                  // Deuflhard/Weiser Satz 7.36
  double const beta = 1.0/sqrt(d*gamma);              // Deuflhard/Weiser (7.67)
  double const alpha = (d*gamma-1.0)/(d*(1.0+gamma)); // Deuflhard/Weiser S. 303
  double q = 2;                                       // Deuflhard/Weiser (7.62)
  long N = size;
  double zk = 0.0;
  double const requested = sqrt(1-beta*beta)*TOL;     // Deuflhard/Weiser Algorithmus 7.37
  double yk = pow(N,alpha);
  double pcgerr = -1;
  double pcgerr_fact = -1;
  int required = 0;
  
  std::vector<VariableSet::VariableSet> solutions;

  do {
    // Zero coefficient vectors for initialization
    ExCoefficientVector exZero(ExVariableSet::CoefficientVectorRepresentation<>::init(exVariableSet)); exZero = 0;
    CoefficientVector zero(VariableSet::CoefficientVectorRepresentation<>::init(variableSet)); zero = 0;
    
    std::cerr << "----------------------------------------\nStarting round\n";
    CoefficientVector solution(zero), hilfe(zero);
    assembler.assemble(linearization(F,x));
    Rhs rhs(assembler.rhs());
    AssembledGalerkinOperator<Assembler,0,neq,0,nvars> A(assembler, onlyLowerTriangle);
    
    if (direct) 
      directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
    else  {
      switch (iterateType) {
	case IterateType::CG: {
	  JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobi(A,1.0);
	  Dune::CGSolver<CoefficientVector> cg(A,jacobi,iteEps,iteSteps,verbose);
	  cg.apply(hilfe,rhs,result);
	  break;
	}
	case IterateType::PCG: {
	  JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobi(A,1.0);
	  NMIIIPCGSolver<CoefficientVector> pcg(A,jacobi,iteEps,iteSteps,verbose);
	  pcg.apply(hilfe,rhs,result);
	  break;
	}
	case IterateType::APCG: {
	  JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobiPCG(A,1.0);
	  DefaultDualPairing<CoefficientVector,CoefficientVector> dp;
	  PCGEnergyErrorTerminationCriterion<double> terminate(iteEps,iteSteps);
	  terminate.lookahead(lookahead);
	  Pcg<CoefficientVector,CoefficientVector> apcg(A,jacobiPCG,dp,terminate,verbose);
	  apcg.apply(hilfe,rhs,result);
	  pcgerr = terminate.error();
	  
	  CoefficientVector hilfe2(zero);
	  PCGEnergyErrorTerminationCriterion<double> terminate2(iteEps/1000,iteSteps);
	  Pcg<CoefficientVector,CoefficientVector> apcg2(A,jacobiPCG,dp,terminate2,verbose);
          Dune::InverseOperatorResult result2;
	  apcg2.apply(hilfe2,rhs,result2);
    std::ostringstream fn;
    fn << "cmg-gamma-";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();
    std::ofstream out(fn.str().c_str());
    terminate.clear();
    required = 0;
    for (int i=0; i<terminate2.gamma2().size(); ++i) {
      terminate.step(terminate2.gamma2()[i]);
      out << terminate2.gamma2()[i] << ' ' << terminate.error() << ' ' << std::sqrt(std::accumulate(terminate2.gamma2().begin()+i,terminate2.gamma2().end(),0.0)) << '\n';
      if (std::sqrt(std::accumulate(terminate2.gamma2().begin()+i,terminate2.gamma2().end(),0.0)) >= iteEps)
	required = i;
    }    
	  
	  
// 	  PCGEnergyErrorTerminationCriterion<double> terminate3(0.0,required+1);
// 	  Pcg<CoefficientVector,CoefficientVector> apcg3(A,jacobiPCG,dp,terminate3,verbose);
// 	  apcg3.apply(hilfe,rhs,result);
// 	  pcgerr = terminate3.error();
// 	  //pcgstep = std::sqrt(std::accumulate(terminate3.gamma2().begin(),terminate3.gamma2().end(),0.0));
	  
	  CoefficientVector dhilfe(hilfe2); dhilfe -= hilfe;
	  CoefficientVector Adhilfe(zero);
	  A.apply(dhilfe,Adhilfe);
	  pcgerr_fact = std::sqrt(dhilfe*Adhilfe);
	  
	  break;
	}
// 	case IterateType::SGS: {
// 	  int addedIterations = 3;
// 	  typedef Dune::BlockVector<Dune::FieldVector<double,1> > NakedCoefficientVector;
// 	  Dune::SeqSSOR<MatrixAsTriplet<double>,NakedCoefficientVector,NakedCoefficientVector> sgs(A.getmat(),1,1.0);
// 	  sgs.pre(boost::fusion::at_c<0>(hilfe.data),boost::fusion::at_c<0>(rhs.data));
// 	  CoefficientVector y(zero), z(zero), r(rhs);
// 	  boost::circular_buffer<double> gamma(addedIterations);
// 	  for (int i=0; i<iteSteps; ++i) {
// 	    sgs.apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(r.data));
// 	    A.apply(y,z);
// 	    r -= z;
// 	    hilfe += y;
// 	    result.iterations = i+1;
// 	    gamma.push_back(std::sqrt(y*z));
// 	    if (i>=addedIterations && std::accumulate(gamma.begin(),gamma.end(),0.0)<iteEps)
// 	      break;
// 	  }
// 	  std::pair<double,double> geoest = estimateGeometricSequence(gamma.begin(),gamma.end());
// 	  std::cerr << "geometric estimate: c=" << geoest.first << " q=" << geoest.second << " delta=" << iteEps << '\n';
// 	  std::copy(gamma.begin(),gamma.end(),std::ostream_iterator<double>(std::cerr," ")); std::cerr << '\n';
// 	  break;
// 	}
	default:
	  std::cerr << "Solver " << iterateType << " not available" << std::endl;
	  throw -111;
      }
      solution.axpy(-1,hilfe);
    }
    dx.data = solution.data;
    
    
    // Do hierarchical error estimation. Remember to provide the very same underlying problem to the
    // error estimator functional as has been used to compute dx (do not modify x!).
    {

    
    estGop.assemble(ErrorEstimator(LinearizationAt<Functional>(F,x),dx));
      
      // iterative solution of error estimator      
      AssembledGalerkinOperator<EstGOP> E(estGop);
      Dune::InverseOperatorResult estRes;
      ExCoefficientVector const estRhside(estGop.rhs()) ;
      ExCoefficientVector estSol(exZero) ;
      JacobiPreconditioner<EstGOP> jprec(estGop, 1.0);
      jprec.apply(estSol,estRhside); //single Jacobi iteration

      // Represent error estimator as FE function
      H1ExSpace::Element<1>::type erroresttmp(spaceEx);
      erroresttmp = boost::fusion::at_c<0>(estSol.data);
      H1Space::Element<1>::type errorest(temperatureSpace2);
      errorest = erroresttmp;
      
      
      // Transfer error indicators to cells.
      std::vector<double> errorDistribution(is.size(0));
      typedef VariableSet::GridView::Codim<0>::Iterator CellIterator ;
      double maxErr = 0.0;
      for (CellIterator ci=variableSet.gridView.begin<0>(); ci!=variableSet.gridView.end<0>(); ++ci) {
	typedef H1ExSpace::Mapper::GlobalIndexRange GIR;
	double err = 0;
	GIR gix = spaceEx.mapper().globalIndices(*ci);
	for (GIR::iterator j=gix.begin(); j!=gix.end(); ++j) 
	  err += boost::fusion::at_c<0>(estSol.data)[*j] * boost::fusion::at_c<0>(estRhside.data)[*j]; // remember that what adds up here is positive (due to diagonal solve)
	assert(err>=0);
	errorDistribution[is.index(*ci)] = std::sqrt(err); // Deuflhard/Weiser (6.37)
	if (err>maxErr) maxErr = err;
      }

      // Select refinement threshold: Refine all elements which 
      // - have at least an error of half of the maximum error or
      // - are among the minRefine fraction with largest error
      double errLevel = maxErr/(order+1);
      if (minRefine>0.0) {
	std::vector<double> eSort(errorDistribution);
	std::sort(eSort.begin(),eSort.end(),std::greater<double>()); // sort decreasingly
	int minRefineIndex = minRefine*(eSort.size()-1);
	errLevel = std::min(errLevel,eSort[minRefineIndex]+1.0e-14);
      }

      // Compute the error estimator norm according to Deuflhard/Weiser (6.36)
      errNorm = 0;
      for (int i=0; i<boost::fusion::at_c<0>(estSol.data).N(); ++i)
	errNorm += boost::fusion::at_c<0>(estSol.data)[i] * boost::fusion::at_c<0>(estRhside.data)[i]; // remember that what adds up here is positive (due to diagonal solve)
      errNorm = std::sqrt(errNorm);
      
      // Termination check
      if (errNorm<requested) {
	accurate = true;
	
	for (int i=0; i<solutions.size()-1; ++i) {
          solutions[i] -= solutions.back();
	  
	  CoefficientVector du(zero), Adu(zero);
	  boost::fusion::at_c<0>(du.data) = boost::fusion::at_c<0>(solutions[i].data).coefficients();
	  A.apply(du,Adu);
	  std::cerr << "i= " << i << " err= " << std::sqrt(du*Adu) << '\n';
	}
      } else {
	// Refine mesh.
	int count = 0;
	for (CellIterator ci=variableSet.gridView.begin<0>(); ci!=variableSet.gridView.end<0>(); ++ci)
	  if (errorDistribution[is.index(*ci)] >= errLevel) {
	    gridManager.mark(1,*ci);
	    ++count;
	  }
	accurate = !gridManager.adaptAtOnce(); 
	
    
    
    std::ostringstream fn;
    fn << "graph/cmg-dx";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();
    writeVTKFile(gridManager.grid().leafView(),variableSet,dx,fn.str(),options,order);
    
      }
    // apply the Newton correction here
    x += dx;
    
    std::ostringstream fn;
    fn << "graph/cmg-sol";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();
    writeVTKFile(gridManager.grid().leafView(),variableSet,x,fn.str(),options,order);
    // add coarse grid error estimator to fine grid prolongation
    errorest *= eefactor;
    boost::fusion::at_c<0>(x.data) -= errorest;
    }
    
    // store solution for later checks
    solutions.push_back(x);
    
    
    std::cout << "step= " << refSteps << " [err]= " << errNorm << " zk= " << zk << " q= " << q << " alpha= " << alpha 
              << " N= " << N << " yk= " << yk << " iter= " << result.iterations << " delta= " << iteEps << " pcgerr= " << pcgerr << " pcgerr_fact= " << pcgerr_fact  
              << " required= " << required << '\n';
    std::cout.flush();

        refSteps++;

	
    // Modify tolerances according to Deuflhard/Weiser Algorithmus 7.37
    q = variableSet.degreesOfFreedom(0,1)/static_cast<double>(N); // DOF growth factor
    N = variableSet.degreesOfFreedom(0,1);                        // new DOF
    yk += pow(N,alpha);                                           // Deuflhard/Weiser S. 303
    zk = pow(N,alpha)*(pow(errNorm/requested,d*alpha)-pow(q,alpha))/(pow(q,alpha)-1.0); // Deuflhard/Weiser S. 304
    iteEps = beta*TOL*yk/(yk+zk);
    
    if (getParameter(pt, "strategy", 1) == 0)
      iteEps = beta*TOL;
    if (getParameter(pt, "strategy", 1) == 2)
      iteEps = errNorm*std::pow(q,1.0/d);
    
    
    
    if (refSteps>maxAdaptSteps) {
      std::cout << "number of steps " << refSteps << " > maxAdaptSteps (" << maxAdaptSteps << ")\n";
      break;
    }
    
    std::ostringstream fn;
    fn << "graph/cmg-grid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps-1;
    fn.flush();
    writeVTKFile(gridManager.grid().leafView(),variableSet,x,fn.str(),options,order);
    
  } while (!accurate);    
  
    

std::cout << "End cmgtest" << std::endl;
}
