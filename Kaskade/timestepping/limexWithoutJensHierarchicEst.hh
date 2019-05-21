#ifndef LIMEX_HH
#define LIMEX_HH

#include <fstream>
#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/istl/solvers.hh"

#include "fem/embedded_errorest.hh"
#include "fem/iterate_grid.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/hierarchicspace.hh"
#include "fem/istlinterface.hh"

#include "dune/istl/solvers.hh"
#include "dune/istl/preconditioners.hh"

#include "timestepping/extrapolation.hh"
#include "timestepping/semieuler.hh"

#include "linalg/factorization.hh"
#include "linalg/umfpack_solve.hh"
// #include "linalg/mumps_solve.hh"
// #include "linalg/superlu_solve.hh"

#include "linalg/trivialpreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/iluprecond.hh"
#include "linalg/additiveschwarz.hh"
// #include "linalg/hyprecond.hh"
#include "linalg/jacobiPreconditioner.hh"

#include "utilities/enums.hh"

#include "io/vtk.hh"

bool fabscompare( double a, double b ) { return fabs( a ) < fabs( b ) ; }

DirectType xyz = DirectType::UMFPACK;

namespace Kaskade
{
  
  
  template <class Functional>
  struct CardioD2
  {
    template <int row, int col>
    class D2: public Functional::template D2<row,col>
    {
      typedef typename Functional::template D2<row,col> d2;
      
      public:
	static bool const present = d2::present && (row==0) && (col==0);
	static bool const lumped = true;
    };
  };
  
/**
 * \ingroup timestepping
 * \brief Extrapolated linearly implicit Euler method.
 * 
 * This class implements the extrapolated linearly implicit Euler
 * method for integrating time-dependent evolution problems. The
 * implementation follows Deuflhard/Bornemann Chapter 6.4.3.
 */
template <class Eq>
class Limex 
{
public:
  typedef Eq EvolutionEquation;
  typedef typename EvolutionEquation::AnsatzVars::VariableSet State;

private:
  typedef SemiLinearizationAt<SemiImplicitEulerStep<EvolutionEquation> > Linearization;
  typedef VariationalFunctionalAssembler<Linearization> Assembler;
  
  // for hierarchic error estimator
  typedef FEFunctionSpace<ContinuousHierarchicExtensionMapper<double,typename Eq::AnsatzVars::Grid::LeafGridView> > SpaceEx;
  typedef boost::fusion::vector< typename SpaceType<typename Eq::AnsatzVars::Spaces,0>::type const*, SpaceEx const*> ExSpaces;
  
//   two components, 1 space -- how to generalize?
//   typedef boost::fusion::vector<VariableDescription<1,1,0>, VariableDescription<1,1,1> > ExVariableDescriptions;
  
  // 1 component, 1 space -- how to generalize?
  typedef boost::fusion::vector<VariableDescription<1,1,0> > ExVariableDescriptions;
  
  typedef VariableSetDescription<ExSpaces,ExVariableDescriptions> ExVariableSet;
  typedef HierarchicErrorEstimator<Linearization,ExVariableSet,ExVariableSet,CardioD2<Linearization> > ErrorEstimator;
  typedef VariationalFunctionalAssembler<ErrorEstimator> EstimatorAssembler;
  
//   typedef typename EstimatorAssembler::template AnsatzVariableRepresentation<> Ansatz;
//   typedef typename EstimatorAssembler::template TestVariableRepresentation<> Test;
/*  typedef typename ExVariableSet::template CoefficientVectorRepresentation<0,1>::type ExCoefficientVectors;*/
  
  typedef typename Eq::AnsatzVars::Grid::Traits::LeafIndexSet IS ;
  

public:
  /**
   * Constructs an ODE integrator. The arguments eq and ansatzVars
   * have to exist during the lifetime of the integrator.
   */
  Limex(GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager_,
        EvolutionEquation& eq_, typename EvolutionEquation::AnsatzVars const& ansatzVars_, DirectType st=DirectType::UMFPACK):
    gridManager(gridManager_), ansatzVars(ansatzVars_), eq(&eq_,0),
    assembler(gridManager,ansatzVars.spaces), 
    iteSteps(10000), iteEps(1e-6), extrap(0),
    rhsAssemblyTime(0.0), matrixAssemblyTime(0.0), factorizationTime(0.0), solutionTime(0.0), adaptivityTime(0.0),
    solverType(st)
    {}


  /**
   * In order to maintain compatibility with existing code,
   * this overload is needed, as non const reference parameters (refinements) can not
   * have default values. All ist does is call the original method with a temporary
   * parameter.
   */
  State const& step(State const& x, double dt, int order,
                    std::vector<std::pair<double,double> > const& tolX)
  {
    std::vector<std::vector<bool> > tmp ;
    return step(x,dt,order,tolX,tmp);
  }

  /**
   * Computes a state increment that advances the given state in
   * time. The time in the given evolution equation is increased by
   * dt.
   *
   * \param x the initial state to be evolved
   * \param dt the time step
   * \param order the extrapolation order >= 0 (0 corresponds to the linearly implicit Euler)
   * \param tolX 
   * \param refinements keep track of the cells marked for refinement
   * \return the state increment (references an internal variable that will be
   *         invalidated by a subsequent call of step)
   *
   * \todo (i) check for B constant, do not reassemble matrix in this case
   *       (ii) implement fixed point iteration instead of new factorization
   *            in case B is not constant
   */
  State const& step(State const& x, double dt, int order,
                    std::vector<std::pair<double,double> > const& tolX,
		    
		    std::vector< std::vector<bool> > &refinements )
  {
    boost::timer::cpu_timer timer;
    
    std::vector<double> stepFractions(order+1);
    for (int i=0; i<=order; ++i) stepFractions[i] = 1.0/(i+1); // harmonic sequence
    extrap.clear();

    typedef typename EvolutionEquation::AnsatzVars::Grid Grid;
    
    int const dim = EvolutionEquation::AnsatzVars::Grid::dimension;
    int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
    int const neq = EvolutionEquation::TestVars::noOfVariables;
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    size_t  size = ansatzVars.degreesOfFreedom(0,nvars);
    
    std::vector<double> rhs(size), sol(size);

    State dx(x), dxsum(x), tmp(x);
    double const t = eq.time();
    
    double estTime = 0 ;

    eq.temporalEvaluationRange(t,t+dt);
    
    typedef AssembledGalerkinOperator<Assembler,0,neq,0,neq> Op;
        
    bool iterative = true ;
    int fill_lev = /*0*/ /*1*/2;
        
    typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::type
      CoefficientVectors;
    typedef typename EvolutionEquation::TestVars::template CoefficientVectorRepresentation<0,neq>::type 
      LinearSpaceX;

    for (int i=0; i<=order; ++i) {
      double const tau = stepFractions[i]*dt;
      eq.setTau(tau);
      eq.time(t);

      // mesh adaptation loop. Just once if i>0.
      bool accurate = true;
      int refinement_count = 0 ;             
	
      
      do {	
        // Evaluate and factorize matrix B(t)-tau*J
        dx/* * */= 0;
        timer.start();
	
// 	assembler.assemble(Linearization(eq,x,x,dx)/*,(Assembler::VALUE|Assembler::RHS|Assembler::MATRIX),4*/);
	
// #ifndef KASKADE_SEQUENTIAL
// 	gridManager.enforceConcurrentReads(std::is_same<Grid,Dune::UGGrid<dim> >::value);
// 	assembler.setNSimultaneousBlocks(40);
// 	assembler.setRowBlockFactor(2.0);
// #endif
// 	std::cout << "assemble linear system ..." ; std::cout.flush();
	assembler.assemble(Linearization(eq,x,x,dx)/*,(Assembler::VALUE|Assembler::RHS|Assembler::MATRIX),4*/);
// 	std::cout << " done\n"; std::cout.flush();

	matrixAssemblyTime += (double)(timer.elapsed().user)/1e9;
	
	Op A(assembler);
// 	std::cout << "solve linear system\n"; std::cout.flush();
	  
	if (iterative) 
	{  	  
	  timer.start();

	  CoefficientVectors  solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars.spaces)
	    );
	  solution = 0;
	  CoefficientVectors rhs(assembler.rhs());
	  
	  ILUKPreconditioner<Op> p(A,fill_lev,0);
// 	  ILUTPreconditioner<Op> p(A,240,1e-2,0);
// 	  JacobiPreconditioner<Op> p(A,1.0);

	  Dune::BiCGSTABSolver<LinearSpaceX> cg(A,p,1e-7,2000,0); //verbosity: 0-1-2 (nothing-start/final-every it.)
	  Dune::InverseOperatorResult res;
	  cg.apply(solution,rhs,res);
	  if ( !(res.converged) || (res.iterations == 2001) ) {
	    std::cout << "   no of iterations in cg = " << res.iterations << std::endl;
	    std::cout << " convergence status of cg = " << res.converged  << std::endl;
	    assert(0);
	  }
	  dx.data = solution.data;
	  solutionTime += (double)(timer.elapsed().user)/1e9;
	}
	else
	{
	  // direct solution
	  timer.start();
// 	  assembler.toTriplet(0,neq,0,nvars,ridx.begin(),cidx.begin(),data.begin(),false);
	  MatrixAsTriplet<double> triplet = A.template get<MatrixAsTriplet<double> >();
	  
	  Factorization<double> *matrix = 0;
	  switch (solverType) {
	    case DirectType::UMFPACK:
	      matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
	      break;
	    case DirectType::UMFPACK3264:
	      matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
	      break;
// 	    case MUMPS:
// 	      matrix = new MUMPSFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
// 	      break;
// 	    case SUPERLU:
// 	      matrix = new SUPERLUFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
// 	      break;
	    default:
	      throw -321;
	      break;
	  }
	  factorizationTime += (double)(timer.elapsed().user)/1e9;
	  
	  // First right hand side (j=0) has been assembled together with matrix.
	  timer.start();
	  A.getAssembler().toSequence(0,neq,rhs.begin());
	  for (int k=0; k<rhs.size(); ++k) assert(std::isfinite(rhs[k]));
	  matrix->solve(rhs,sol);
	  delete matrix;
	  for (int k=0; k<sol.size(); ++k) assert(std::isfinite(sol[k]));
	  dx.read(sol.begin());
	  solutionTime += (double)(timer.elapsed().user)/1e9;
	  

// 	  Factorization<double> *matrix = 0;
// 	  switch (solverType)
// 	    {
// 	      case UMFPACK:
// 		matrix = new UMFFactorization<double>(size,0,ridx,cidx,data);
// 		break;
// 	      case UMFPACK3264:
// 		matrix = new UMFFactorization<double>(size,0,ridx,cidx,data);
// 		break;
// 	      case MUMPS:
// 		matrix = new MUMPSFactorization<double>(size,0,ridx,cidx,data);
// 		break;
// 	      case SUPERLU:
// 		matrix = new SUPERLUFactorization<double>(size,0,ridx,cidx,data);
// 		break;
// 	      default:
// 		throw -321;
// 		break;
// 	    }
// 	  factorizationTime += (double)(timer.elapsed().user)/1e9;
// 
// 	  // First right hand side (j=0) has been assembled together with matrix.
// 	  timer.start();
// 	  assembler.toSequence(0,neq,rhs.begin());
// 	  for (size_t k=0; k<rhs.size(); ++k) assert(finite(rhs[k]));
// 	  matrix->solve(rhs,sol);
// 	  delete matrix;
// 	  for (size_t k=0; k<sol.size(); ++k) assert(finite(sol[k]));
// 	  dx.read(sol.begin());
// 	  solutionTime += (double)(timer.elapsed().user)/1e9;
	  //end: direct solution
	}

        // Mesh adaptation only if requested and only for the full
        // implicit Euler step (first stage).
        
	typedef typename EvolutionEquation::AnsatzVars::GridView::template Codim<0>::Iterator CellIterator ;
        
	accurate = true ;
	timer.start();
	IS const& is = gridManager.grid().leafIndexSet();
	std::vector<double> errorDistribution(is.size(0),0.0);
	double maxErr = 0.0, errNorm = 0.0 ;

        if (!tolX.empty() && i==0) {

	  // embeddedErrorEstimator can not be used for linear finite elements
	  // use HierarchicErrorEstimator instead
	  // preparation for HierarchicErrorEstimator
	  if(true) // to avoid mesh adaptation for exVariables
	  {
	  //TODO: this should be using the gridviews/index sets of the original space, not the leaf!
	  SpaceEx  spaceEx(gridManager,gridManager.grid().leafGridView(), boost::fusion::at_c<0>(ansatzVars.spaces)->mapper().getOrder()+1);
	  typename SpaceType<typename Eq::AnsatzVars::Spaces,0>::type spaceH1 = *(boost::fusion::at_c<0>(ansatzVars.spaces)) ;
	  ExSpaces exSpaces(&spaceH1,&spaceEx);
	  std::string exVarNames[2] = { "ev", "ew" };
	  ExVariableSet exVariableSet(exSpaces, exVarNames);
	  EstimatorAssembler estAssembler(gridManager,exSpaces);   
	   
	  tmp/* * */= 0 ;
// 	  estAssembler.assemble(ErrorEstimator(Linearization(eq,x,x,tmp/*dx*/),dx));
	  estAssembler.assemble(ErrorEstimator(semiLinearization(eq,x,x,/*dx*/tmp),dx));

	  
	  int const estNvars = ErrorEstimator::AnsatzVars::noOfVariables;
// 	  std::cout << "estimator: nvars = " << estNvars << "\n";
	  int const estNeq = ErrorEstimator::TestVars::noOfVariables;
// 	  std::cout << "estimator:   neq = " << estNeq << "\n";
	  size_t  estNnz = estAssembler.nnz(0,estNeq,0,estNvars,false);
	  size_t  estSize = exVariableSet.degreesOfFreedom(0,estNvars);
	  
	  std::vector<int> estRidx(estNnz), estCidx(estNnz);
	  std::vector<double> estData(estNnz), estRhs(estSize), estSolVec(estSize);
	  
	  estAssembler.toSequence(0,estNeq,estRhs.begin());
	  
	  typedef typename ExVariableSet::template CoefficientVectorRepresentation<0,estNeq>::type ExCoefficientVectors;
	  
	  // iterative solution of error estimator
// 	  timer.start();
// 	  AssembledGalerkinOperator<EstimatorAssembler> E(estAssembler);
// 	  typename Dune::InverseOperatorResult estRes;
// 	  typename Test::type estRhside( Test::rhs(estAssembler) ) ;
// 	  typename Ansatz::type estSol( Ansatz::init(estAssembler) ) ;
// 	  estSol = 1.0 ;
// 	  JacobiPreconditioner<EstimatorAssembler> jprec(estAssembler, 1.0);
// 	  jprec.apply(estSol,estRhside); //single Jacobi iteration
// 	  estTime += (double)(timer.elapsed().user)/1e9;
// 	  estSol.write(estSolVec.begin());
	  
	  typedef AssembledGalerkinOperator<EstimatorAssembler> AssEstOperator;
	  AssEstOperator agro(estAssembler);
	  Dune::InverseOperatorResult estRes;
	  ExCoefficientVectors estRhside(estAssembler.rhs());
	  ExCoefficientVectors estSol(ExVariableSet::template CoefficientVectorRepresentation<0,estNeq>::init(exVariableSet.spaces));
	  estSol = 1.0 ;
	  JacobiPreconditioner<AssEstOperator> jprec(agro, 1.0);
	  jprec.apply(estSol,estRhside); //single Jacobi iteration
// 	  estTime += (double)(timer.elapsed().user)/1e9;
	  estSol.write(estSolVec.begin());
	  
	  //--
	  
	  // Transfer error indicators to cells.
	  CellIterator ciEnd = ansatzVars.gridView.template end<0>() ;
	  for (CellIterator ci=ansatzVars.gridView.template begin<0>(); ci!=ciEnd; ++ci) {
	    typedef typename SpaceEx::Mapper::GlobalIndexRange GIR;
	    double err = 0.0;
	    GIR gix = spaceEx.mapper().globalIndices(*ci);
	    for (typename GIR::iterator j=gix.begin(); j!=gix.end(); ++j)
	      err += fabs(boost::fusion::at_c<0>(estSol.data)[*j]);
	    // only for 1st component -- second in monodomain is ODE
	    errorDistribution[is.index(*ci)] = err;
	    if (fabs(err)>maxErr) maxErr = fabs(err);
	  }
// 	  std::cout << "maxErr: " << maxErr << "\n";
	  
// 	  for( size_t cno = 0 ; cno < is.size(0); cno++ )
// 	    std::cout << cno << "\t\t" << errorDistribution[cno] << "\n";
// 	  std::cout.flush();
	  
	  //TODO: estimate only error in PDE, not ODE!
	  // what about relative accuracy?
	  for (size_t k = 0; k < estRhs.size() ; k++ ) errNorm +=  fabs( estRhs[k] * estSolVec[k] ); 
	  // this is for all variables?!
// 	  errNorm = tau/*dt*/*sqrt(errNorm);
	  std::cout << "errNorm = " << errNorm << std::endl ;
	  std::cout.flush();

	  }
	
	  double overallTol = tolX[0].first ;

	  
// 	  for( int i = 0 ; i < 1 /*nvars*/ ; i++ )
// 	    overallTol += tolX[i].first*tolX[i].first ; 
// 	  overallTol = sqrt(overallTol);
	  
// 	  std::cout << "overallTol = " << overallTol << std::endl ;
	  
	  
// 	  alpha = 0 ; // for uniform refinement

	  int nRefMax = /*10*/7 ;/*7*//*5*//*3*//*15*/ /*4*/ /*5*/ /*10*/
	  double fractionOfCells = 0.05; /*0.1,0.2*/
	  unsigned long noToRefine = 0, noToCoarsen = 0;
	  double alpha = 1.0; /*0.5; //test april 13*/
	  
	  double errLevel = 0.5*maxErr ; //0.75
// 	  double minRefine = 0.05;
	  
// 	  if (minRefine>0.0)
// 	  {
// 	    std::vector<double> eSort(errorDistribution);
// 	    std::sort(eSort.begin(),eSort.end(),fabscompare);
// 	    int minRefineIndex = minRefine*(eSort.size()-1);
// 	    double minErrLevel = fabs(eSort[minRefineIndex])+1.0e-15;
// 	    if (minErrLevel<errLevel)
// 	      errLevel = minErrLevel;
// 	  }
	  
// 	  double minRefine = fractionOfCells; //* gridManager.grid().size(0);
// 	  if (minRefine >0.0)
// 	  {
// 	    std::vector<double> eSort(errorDistribution);
// 	    std::sort(eSort.begin(),eSort.end(),fabscompare);
// 	    int minRefineIndex = minRefine*(eSort.size()-1);
// 	    double minErrLevel = fabs(eSort[minRefineIndex])+1.0e-15;
// 	    if (minErrLevel<errLevel)
// 	      errLevel = minErrLevel;
// 	  }

	  size_t maxNoOfVertices =/* 500000*/ /*2000000*//*12000*/30000;
	  size_t maxNoOfElements = 2000000;
	  int maxLevel = /*20*/6 /*25*/ /*18*/; // additional refinement levels
	    
	  if( errNorm > overallTol && refinement_count < nRefMax 
	    && gridManager.grid().size(/*dim*/0) < /*maxNoOfVertices*/maxNoOfElements )
	    //current setup: for 2D BFGS; for 3D fibrillation: use maxNoOfElements
	  {    
// 	    accurate = false;
// 	    std::vector<bool> toRefine( is.size(0), false ) ; 
// 	    size_t noCells = gridManager.grid().size(0);
// // 	    std::cout << is.size(0) << "   " << noCells << "\n"; std::cout.flush();
// 	    for( size_t cno = 0 ; cno < noCells ; cno++ )
// 	    {
// 	      if (fabs(errorDistribution[cno]) >= alpha*errLevel)
// 	      {
// 		noToRefine++;
// 		toRefine[cno] = true ;
// 	      }
// 	    }
	    std::vector<bool> toRefine( is.size(0), false ) ; //for adaptivity in compression
	    
	    accurate = false ;
	    size_t noCells = gridManager.grid().size(0);
	    
// 	    // Nagaiah paper:
// 	    unsigned minLevel = 10 ; // do not coarsen below that level     
// 	    CellIterator ciEnd = ansatzVars.gridView.template end<0>() ;
// 	    for (CellIterator ci=ansatzVars.gridView.template begin<0>(); ci!=ciEnd; ++ci) 
// 	    {
// 	      double vol = sqrt( ci->geometry().volume() );
// 	      if (fabs(errorDistribution[is.index(*ci)])/vol >= 0.2 /*0.1*/ )
// 	      {
// 		if( ci->level() < maxLevel && gridManager.grid().size(dim) < maxNoOfVertices )
// 		{
// 		  gridManager.mark(1,*ci);
// 		  noToRefine++;
// 		}
// 	      }
// 	      else if (fabs(errorDistribution[is.index(*ci)])/vol < 0.1 /* /2 */ && ci->level() > minLevel ) 
// 	      {
// 		gridManager.mark(-1,*ci);
// 		noToCoarsen++;
// 	      }
// 	    }
// 	    std::cout << "refine: " << noToRefine << "    coarsen: " << noToCoarsen << "\n";
// 	    gridManager.countMarked();
	    
// 	    last implementation was the following:
	    unsigned long noToRefineOld = 0;
	    do{
	      noToRefineOld = noToRefine;
	      for( size_t cno = 0 ; cno < noCells ; cno++ )
	      {
		      if (fabs(errorDistribution[cno]) >= alpha*errLevel) 
		      { // TODO: check and bugfix this!
			      if( !toRefine[cno] ) noToRefine++;
			      toRefine[cno] = true ;
		      }
	      }
// 		std::cout << "alpha: " << alpha << "\trefine: " << noToRefine << "\n";
// 		std::cout.flush();
		alpha *= 0.75; // refine more aggressively: reduce alpha by larger amount
	    } while ( noToRefine < fractionOfCells * noCells && (!noToRefineOld == noToRefine) ) ;
//   	    std::cout << "to refine: " << noToRefine << "\n";
	    std::cout.flush();

	    
// 	    alpha = 0.5; //50
// 	    do {
// 	      for( size_t cno = 0 ; cno < noCells ; cno++ )
// 		if (fabs(errorDistribution[cno]) > alpha*maxErr/*overallTol/noCells*/) {
// 		    /*if( !toRefine[cno] )*/ noToRefine++;
// 		    toRefine[cno] = true;
// 		  }
// 	      
// 		  alpha *=0.9;
// 	    } while (noToRefine < fractionOfCells * noCells);
	    
// 	    long int nRefined = (long int) std::count( toRefine.begin(), toRefine.end(), true );
// 	    std::cout << "   noToRefine = " << nRefined << std::endl;
// 	    refinements.push_back(toRefine);

	    for (CellIterator ci=ansatzVars.gridView.template begin<0>(); ci!=ansatzVars.gridView.template end<0>(); ++ci)
	      if( toRefine[is.index(*ci)] && ci->level() < maxLevel )  gridManager.mark(1,*ci);
	    
	    bool refok = gridManager.adaptAtOnce();  // something has been refined?
	    if( !refok ) 
	    {
	      // nothing has been refined
	      std::cout << "nothing has been refined (probably due to maxLevel constraint), stopping\n";
	      accurate = true ;
	    }
	    
	    if (!accurate) {
	      refinement_count++ ;
	      nnz = assembler.nnz(0,neq,0,nvars,false);
	      size = ansatzVars.degreesOfFreedom(0,nvars);
	      rhs.resize(size);
	      sol.resize(size);
	    }
	  }
	  else
	  {
	   if( errNorm > overallTol && refinement_count > nRefMax ) 
	     std::cout << "max no of refinements exceeded\n";
// 	   if( errNorm > overallTol && gridManager.grid().size(2) > maxNoOfVertices ) 
// 	     std::cout << "max no of nodes exceeded\n";
	    
	   accurate = true ;
	  }
        }
      } while (!accurate); 
      adaptivityTime += (double)(timer.elapsed().user)/1e9;

//       std::cout << gridManager.grid().size(dim) << " nodes on " << gridManager.grid().maxLevel() << "levels\n";

      dxsum = dx;

      // propagate by linearly implicit Euler 
      for (int j=1; j<=i; ++j) {
	// Assemble new right hand side tau*f(x_j)
        eq.time(eq.time()+tau);
	tmp = x; tmp += dxsum;
        State zero(x); zero /* * */= 0;
        timer.start();
// 	gridManager.enforceConcurrentReads(true);
// 	assembler.setNSimultaneousBlocks(40);
// 	assembler.setRowBlockFactor(2.0);
        assembler.assemble(Linearization(eq,tmp,x,zero),Assembler::RHS|Assembler::MATRIX,4);
        rhsAssemblyTime += (double)(timer.elapsed().user)/1e9;
	
	Op A(assembler);
	
	if( !iterative )
	{
	  timer.start();
	  //direct solution
	  A.getAssembler().toSequence(0,neq,rhs.begin());
	  MatrixAsTriplet<double> triplet = A.template get<MatrixAsTriplet<double> >();
	  Factorization<double> *matrix = 0;
	  switch (solverType) {
	    case DirectType::UMFPACK:
	      matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
	      break;
	    case DirectType::UMFPACK3264:
	      matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
	      break;
// 	    case MUMPS:
// 	      matrix = new MUMPSFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
// 	      break;
// 	    case SUPERLU:
// 	      matrix = new SUPERLUFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
// 	      break;
	    default:
	      throw -321;
	      break;
	  }
	  factorizationTime += (double)(timer.elapsed().user)/1e9;
	  timer.start();
	  matrix->solve(rhs,sol);
	  delete matrix;
	  for (int k=0; k<rhs.size(); ++k) assert(std::isfinite(rhs[k]));
	  for (int k=0; k<sol.size(); ++k) assert(std::isfinite(sol[k]));
	  dx.read(sol.begin());
	  solutionTime += (double)(timer.elapsed().user)/1e9;
	}
	else
	{
	  timer.start();
	  CoefficientVectors  solution(EvolutionEquation::AnsatzVars::template
		    CoefficientVectorRepresentation<0,neq>::init(ansatzVars.spaces) );
	  solution = 0;
	  // 	  assembler.assemble(linearization(F,x));
	  CoefficientVectors rhs(assembler.rhs());
	  
	  ILUKPreconditioner<Op> p(A,fill_lev,0);
// 	  ILUTPreconditioner<Op> p(A,240,1e-2,0);
// 	  JacobiPreconditioner<Op> p(A,1.0);

	  
	  Dune::BiCGSTABSolver<LinearSpaceX> cg(A,p,1e-7,2000,0); //verbosity: 0-1-2 (nothing-start/final-every it.)
	  Dune::InverseOperatorResult res;
	  cg.apply(solution,rhs,res);
	  if ( !(res.converged) || (res.iterations == 2001) ) {
	    std::cout << "   no of iterations in cg = " << res.iterations << std::endl;
	    std::cout << " convergence status of cg = " << res.converged  << std::endl;
	    assert(0);
	  }
	  dx.data = solution.data;
// 	// iterative solution
//         timer.start();
// 	typedef typename Assembler::template TestVariableRepresentation<>::type Rhs;
// 	typedef typename Assembler::template AnsatzVariableRepresentation<>::type Sol;
// 	Sol solution(Assembler::template AnsatzVariableRepresentation<>::init(assembler));
// 	Rhs rhside(Assembler::template TestVariableRepresentation<>::rhs(assembler));
// 	AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, false);
// 	typename AssembledGalerkinOperator<Assembler,0,1,0,1>::matrix_type tri(A.getmat());
// 
// 	
//         typedef typename Assembler::template TestVariableRepresentation<>::type LinearSpace;
// 	Dune::InverseOperatorResult res;
// 	solution = 1.0;
//      	JacobiPreconditioner<Assembler,0,1,0> p(assembler,1.0);
// // 	TrivialPreconditioner<LinearSpace> trivial;
// 	Dune::CGSolver<LinearSpace> cg(A,/*trivial*/p,iteEps,iteSteps,0);
// 	cg.apply(solution,rhside,res);
// 	A.applyscaleadd(-1.0,rhside,solution);
// 	double slntime = (double)(timer.elapsed().user)/1e9;
// 	solutionTime += slntime ;
// 	dx.data = solution.data ;  
// 	// end: iterative solution
	}
	
        dxsum += dx;
        solutionTime += (double)(timer.elapsed().user)/1e9;
      }

      // insert into extrapolation tableau
      extrap.push_back(dxsum,stepFractions[i]);

      // restore initial time
      eq.time(t);
    }

    return extrap.back();
  }

  /**
   * Estimates the time discretization error of the previously
   * computed step by taking the difference between the diagonal and
   * subdiagonal extrapolation values of maximal order. This requires
   * that order>1 has been given for the last step.
   */
  std::vector<std::pair<double,double> > estimateError(State const& x,int i, int j) const 
  {
    assert(extrap.size()>1);

    std::vector<std::pair<double,double> > e(ansatzVars.noOfVariables);
    
    relativeError(typename EvolutionEquation::AnsatzVars::Variables(),extrap[i].data,
                  extrap[j].data,x.data,
                  ansatzVars.spaces,eq.scaling(),e.begin());
    
    return e;
  }
  
  template <class OutStream>
  void reportTime(OutStream& out) const {
    out << "Limex time: " << matrixAssemblyTime << "s matrix assembly\n"
        << "            " << rhsAssemblyTime << "s rhs assembly\n"
        << "            " << factorizationTime << "s factorization\n"
        << "            " << solutionTime << "s solution\n"
	<< "            " << adaptivityTime << "s adaptivity\n";
  }
  
  void advanceTime(double dt) { eq.time(eq.time()+dt); }
  
  
private:
  GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager;
  typename EvolutionEquation::AnsatzVars const& ansatzVars;
  SemiImplicitEulerStep<EvolutionEquation>      eq;
  Assembler                                     assembler;
  int iteSteps ;
  double iteEps ;
    
public:
  ExtrapolationTableau<State>  extrap;
  double rhsAssemblyTime, matrixAssemblyTime, factorizationTime, solutionTime, adaptivityTime;
  DirectType solverType;
};

} //namespace Kaskade
#endif
