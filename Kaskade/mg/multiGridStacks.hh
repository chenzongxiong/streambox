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

/**
 * \file
 * \author Felix Lehmann, Sebastian Götschel, Martin Weiser
 */

#ifndef MULTIGRIDSTACKS_HH
#define MULTIGRIDSTACKS_HH

#include <memory> // std::unique_ptr
#include <vector>
#include <map>
#include <time.h>
#include <math.h>
#include <fstream>

#include <iostream> // std::cout, std::endl
using std::cout;
using std::endl;


#include "dune/grid/common/grid.hh"
#include "dune/grid/common/entitypointer.hh"
#include "dune/grid/common/entity.hh"
#include "dune/grid/io/file/dgfparser/dgfparser.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/iteratorfacades.hh"
#include "dune/istl/matrixindexset.hh"

#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvers.hh"

// #include "fem/barycentric.hh"
#include "fem/fixdune.hh"
#include "fem/mllgeometry.hh"
#include "fem/fetransfer.hh"
#include "linalg/conjugation.hh"
#include "linalg/triplet.hh"

/**
 * Creates a vector of the prolongation matrices (base transformation..) and a 
 * vector of the restricted stiffness matrices for each gridlevel
 *
 *
 * Access to BCRS prolongation matrix from level k to k+1:
 * ProlongationStack<Dune::UGGrid<dim> > ps( gridManager.grid() );
 * BCRSMatrix<FieldMatrix<double,1,1> >  prolongation = ps.prolStack[k];
 *
 * Access to k times restricted BCRS stiffness matrix with assembled Matrix A:
 * ProlongationStack<Dune::UGGrid<dim> > ps( gridManager.grid() );
 * MlStack<Dune::UGGrid<dim> > myStack( ps , A );
 * BCRSMatrix<FieldMatrix<double,1,1> >  levelStiff = myStack.levelMatrixStack[k];
 *
 *
 * We have three classes:  ParentalNodes
 *			   ProlongationStack
 *			   MlStack
 */

void bcrsPrint( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& bcrsMatrix  );
double energyNorm( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& matrix_, 
		   const Dune::BlockVector<Dune::FieldVector<double,1> >& vec_ );

//--------------------------------------------------------------------------------------------------------------------------


/**
 * \brief Finds the parent nodes and their interpolation weight for each node in the grid.
 * usage eg.:  ParentalNodes<Grid> parentalNodes( gridManager.grid() );
 *
 * creates a "level-nested" vector "levelParentsVector" whose components are vectors, containing
 * the indices of vertices and barycentric coordinates corresponding to the parent cells

 * creates also the "level-nested" Vector "levelDisappear" whose entries will be set 1 if 
 * corresponding node disappears in current level and 0 else
 */
template < class Grd >
class ParentalNodes
{
  typedef typename Grd::LevelGridView::IndexSet::IndexType 	IndexType;
  typedef typename std::pair<IndexType,float> 			FatherIndexValue;
  typedef typename std::vector<FatherIndexValue> 		FatherSimplexPairs;
  
  private:
    void findParents( const Grd& grd, int levelNr );
  public:
    std::vector< std::vector<FatherSimplexPairs> >  	levelParentsVector;
    std::vector< std::vector<int> > 			levelDisappear;
    
    ParentalNodes( Grd const& grd ): levelParentsVector(grd.maxLevel()), levelDisappear(grd.maxLevel())
    {
      for(int levelNr=1; levelNr<=grd.maxLevel(); levelNr++)   
        findParents(grd,levelNr);
    }
};

template < class Grd >
void ParentalNodes<Grd>::findParents( const Grd& grd, int levelNr )
{
  // First we iterate through all elements in
  // the current level. In each element we iterate through each corner node. First we check if we
  // visited it before (vector "processed"). If not we collect the barycentric coordinates in 
  // the father element and the corresponding level indices of the father corner vertex in the
  // level above ( using vector<pair<int index,double barycentricCoordinate> > >, each vector 
  // consisting of maximal four pairs ). We also collect the nodes, which disappear from one to
  // another level, we need them (!) to build up consistent prolongations.
  
  typedef typename Grd::template Codim<0>::LevelIterator   		ElementLevelIterator;
  typedef typename Grd::LevelGridView::IndexSet::IndexType 		IndexType;
  
  std::vector<FatherSimplexPairs> levelParents;
  levelParents.resize( grd.size(levelNr,Grd::dimension) );
    
  // helpvector to avoid multiple processing of identical vertices
  std::vector<bool> processed(grd.size(Grd::dimension),false);
  
  // disappearance monitor: necessary to keep track of disappearing nodes
  std::vector<bool> disappears( grd.size(levelNr-1,Grd::dimension) , true );

  typename Grd::LevelGridView::IndexSet const& indexSet = grd.levelGridView(levelNr).indexSet() ;
  typename Grd::LevelGridView::IndexSet const& fatherIndexSet = grd.levelGridView(levelNr-1).indexSet();
  
  // iterate through all cells on current level
  ElementLevelIterator itEnd = grd.template lend<0>(levelNr) ;
  for (ElementLevelIterator it = grd.template lbegin<0>(levelNr) ; it != itEnd; ++it)
  {
    // for storing the father level indices of the father's corner points
    std::vector<IndexType> fatherIndices;
    
    // VertexPointer to father Corners
    typename Grd::template Codim<0>::EntityPointer itFather = (*it).father();
    
    for(int i=0 ; i<=Grd::dimension ; i++ ) 
    {
      fatherIndices.push_back( fatherIndexSet.subIndex(*itFather, i, Grd::dimension) );
      // corners of the father cell won't disappear on coarser level
      disappears[ indexSet.subIndex( *itFather,i,Grd::dimension ) ] = false;
    }
    
    // now regard corners of subentity
    for( int i=0 ; i<=Grd::dimension ; i++ )
    {
      IndexType idx_child = indexSet.subIndex(*it, i, Grd::dimension);
      if (!processed[idx_child])
      { 
	processed[idx_child] = true;
	
        std::vector<float> relativeCoord;
        relativeCoord.push_back(1);
	for( int component = 0 ; component < Grd::dimension ; component++ ) 
	{
	  relativeCoord.push_back( it->geometryInFather().corner(i)[component] );
	  relativeCoord.front() -=  relativeCoord.back();
	}
                
	for( int component = 0 ; component <= Grd::dimension ; component++ ) 
	{
          FatherIndexValue couple;
          couple.first =  fatherIndices[component];
          couple.second = relativeCoord[component];
          levelParents[idx_child].push_back(couple);
	}
	
	// TODO: instead use barycentric() from fem/barycentric.hh, something like this:
	// Dune::FieldVector<typename Grd::ctype,Grd::dimension+1> relativeCoord 
        //                    = barycentric(it->geometryInFather().corner(i));
        // ...for() {...
	//       couple.second = relativeCoord[(component+1)%(Grd::dimension+1)];
	// ... }
	//
        // somehow barycentric coordinates are in different order than mine
	
      }
    }
  }
  
  // store which nodes will not turn up on the coarser level
  for (int k = 0 ; k < disappears.size(); k++ )
    if(disappears[k])
      levelDisappear[levelNr-1].push_back(k);
  
  // move data to the global store.
  levelParentsVector[levelNr-1].swap(levelParents);
}

//--------------------------------------------------------------------------------------------------------------------------

/**
 * (II) ProlongationStack:
 *
 * Now we use the resulting vectors "levelParentsVector" and "levelDisappear" to create the 
 * prolongation matrices. UGGrid creates the leaf indices as follows: The indices are consecutive,
 * starting at zero. They are numerated exactly in the order, in which the nodes disappeared level 
 * after level. Finally the indices of the maximal grid level are attached in natural order.
 * Following this idea, we build up a sparse BCRS matrix, row-wise:
 * 
 * The first (m) rows enclose a simple Identity, needed for carrying the (m) old nodes, which 
 * disappeared the level before the previous one.
 * 
 * The next (n) rows take care of the (n) nodes which disappeared right before. Each of these 
 * rows also consists of exactly one 1, at the position, which ist the index of the 
 * corresponding node in previous level (where it disappeared) plus (m), the number of 
 * already carried with nodes.
 * 
 * The last rows provide information care of the still present nodes.
 * Each row consists of four entries ( 3-DIM, for 2-DIM there are 3).
 * Each value ist a barycentric coordinate and its position ist the level index of the 
 * corresponding father node plus (m).
 * 
 * I called the numbers (m) and (m)+(n) of disappeared nodes disappearSumBefore and disappearSumAfter
 */

// example usage:  ProlongationStack<Grid> prolstack( gridManager.grid() );

// creates Vector prolStack of prolongationmatrices
template < class Grd >
class ProlongationStack
{
  typedef Dune::FieldMatrix<double,1,1> M;
  
  public:
    std::vector<Dune::BCRSMatrix<M> > prolStack;
    ProlongationStack( const Grd& grd );
};


template < class Grd >
ProlongationStack<Grd>::ProlongationStack( const Grd& grd )
{
  int maxiLevel = grd.maxLevel();
  //  creates levelParentsVector
  ParentalNodes<Grd> parentalNodes( grd ); 
  
  std::vector<std::vector<int> > disappear(parentalNodes.levelDisappear);
  
  int disappearSumBefore = 0;
  int disappearSumAfter  = 0;
  
  for(int levelNr = 1 ; levelNr <= maxiLevel ; levelNr++)
  {
    disappearSumBefore = disappearSumAfter;
    disappearSumAfter += disappear[levelNr-1].size(); 
    
    int prolCols = grd.size( levelNr-1 , Grd::dimension ) + disappearSumBefore ;
    int prolRows = grd.size( levelNr , Grd::dimension ) + disappearSumAfter ;
    
    // In order to create a BCRS Prolongation matrix we first compute an upper bound for the 
    // number of nonzeros. The nodes on lower levels just retain their values and are not
    // interpolated, hence there is exactly one nonzero entry in the row. 
    // Each child node's value is determined by the dim+1 values of its father simplex's corners, such that
    // each row has at most dim+1 nonzeros. Most refinement schemes lead to even sparser prolongation matrices,
    // as not all father's corners contribute. We start with the worst case here.
    int nz = disappearSumAfter + grd.size(levelNr,Grd::dimension)*(Grd::dimension+1);
    Dune::BCRSMatrix<M> prolMatrix(prolRows,prolCols,nz,Dune::BCRSMatrix<M>::row_wise);
    
    // Now we create the sparsity pattern. Here we can be more precise and allocate only those entries
    // which are actually nonzero. Remember that "very close to zero" barycentric coordinates have been
    // set to exactly zero.
    typedef Dune::BCRSMatrix<M>::CreateIterator Iter;
    
    for (Iter row=prolMatrix.createbegin(); row!=prolMatrix.createend(); ++row)
    {
      if (row.index() < disappearSumBefore)
	row.insert( row.index() );
      else if (row.index() >= disappearSumBefore && row.index() < disappearSumAfter  )
	row.insert(  disappearSumBefore + disappear[levelNr-1][ row.index() - disappearSumBefore ] );
      else
	for( int k = 0 ; k <= Grd::dimension ; k++ ) 
          if (parentalNodes.levelParentsVector[levelNr-1][row.index()-disappearSumAfter][k].second > 0)
            row.insert( parentalNodes.levelParentsVector[levelNr-1][row.index()-disappearSumAfter][k].first + disappearSumBefore );
    }
    
    // now enter the values
    for(int row = 0 ; row < disappearSumBefore ; row++)
      prolMatrix[row][row] = 1;
    for(int indx = 0 ; indx < disappear[levelNr-1].size(); indx++)
      prolMatrix[ disappearSumBefore + indx ][ disappear[levelNr-1][indx] + disappearSumBefore ] = 1;
    for(int row = disappearSumAfter ; row < prolRows ; row++)
    {
      for( int k = 0 ; k <= Grd::dimension ; k++ ) 
	if (parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter /*row.index()*/][k].second > 0)
	{
	  prolMatrix[row][ parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter][k].first + disappearSumBefore ]
	    = parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter][k].second;
	}
    }
    
    prolStack.push_back( prolMatrix );
  }
}

// ------------------------------------------------------------------------------------------------------

/**
 * (III) MlStack:
 *
 * The obtained prolongationmatrices are used to build up a vector of n times restricted 
 * stiffness matrices A[k]
 * I used conjugation.hh to multiply P^T A P recursively.
 *
 * In conjugation.hh parameter onlyLowerTriangle can be set true, so that only the lower triangle of A will be visited
 * and only lower triangle of the level stiffness matrices will be created.
 *
 */

// example usage: ProlongationStack<Grid> prolStack( gridManager.grid() );  MlStack<Grid> myStack( prolStack , bcrsMatrix );

// creates vector levelMatrixStack of level stiffness matrices
template < class Grd >
class MlStack
{
  typedef Dune::FieldMatrix<double,1,1> M;
  public:
    // sets the member levelMatrixStack, a vector of level stiffness matrices
    MlStack( ProlongationStack<Grd> const& ps , Dune::BCRSMatrix<M> const& bcrsA , bool onlyLowerTriangle = false);
    std::vector<Dune::BCRSMatrix<M> > levelMatrixStack;
};


template < class Grd >
MlStack<Grd>::MlStack(ProlongationStack<Grd> const& ps , Dune::BCRSMatrix<M> const& bcrsA , bool onlyLowerTriangle)
{ 
  // number of prolongations ( = # gridlevels - 1 )
  int prols = ps.prolStack.size();
  // sets the first entry of levelMatrixStack to A (stiffness matrix on leaflevel)
  levelMatrixStack.push_back(bcrsA);

  for(int i = 0 ; i < prols ; i++)
  { 
    // format needed for toolbox conjugation.hh --> C=C+P^T*A*P
    std::unique_ptr<Dune::MatrixIndexSet> C = conjugationPattern<double,M>( ps.prolStack[prols-1-i] , levelMatrixStack[ i ] , onlyLowerTriangle);
    Dune::BCRSMatrix<M> Cnew;
    // exportIdx from Dune::MatrixIndexSet
    C->exportIdx(Cnew);
    // initialize Cnew
    Cnew *= 0.0;
    // now Cnew = Cnew + P^T * A * P , where P and A are level corresponding 
    // prolongation- and stiffness matrices
    conjugation( Cnew , ps.prolStack[prols-1-i] , levelMatrixStack[ i ] , onlyLowerTriangle );
    // append new stiffness matrix
    levelMatrixStack.push_back( Cnew );
  }
}

// ------------------------------------------------------------------------------------------------------

/**
 * \ingroup linalgsolution
 * MultigridSolver realizes a classical v-cycle multigrid solver
 * 
 * yet truncation criterion is relative defect tolerance ||Ax-b||_2/||b||_2 < maxTol
 *
 * example usage:
 *
 *   ProlongationStack<Dune::UGGrid<dim> > ps( gridManager.grid() );
 *   MlStack<Dune::UGGrid<dim> > myStack( ps , bcrsMatrix );
 *   MultigridSolver<Grid> MGSolver( ps , myStack , JacobiRelaxFactor , maxIter , mgSmoothings, maxTol  );
 *   MGSolver.mgSolve( solution , rhs );
 */

template < class Grd , class domain_type >
class MultigridSolver
{
  typedef Dune::BlockVector<Dune::FieldVector<double,1> > CoeffVector ;
  typedef Dune::FieldMatrix<double,1,1> M;
  typedef Dune::SeqJac<Dune::BCRSMatrix<M>,CoeffVector,CoeffVector> JacobiSmoother;
  private:
    int    mNumOfLevels;
    const  ProlongationStack<Grd> &mProlongations;
    const  MlStack<Grd> &mLevelMatrices;
    double mJacRelaxFactor, mMgIter, mMgSmoothings, mMaxTol, mTimeDirectSolve , mExactSolEnergyNorm ;
  public:
    MultigridSolver( const ProlongationStack<Grd>& prolongations,
                     const MlStack<Grd>& levelMatrices,
                     double relax, 
                     double mgIter,
                     double mgSmoothings,
                     double maxTol );
    void mgSolve( domain_type& solution , domain_type& rightHand );
  private:
    void   mgSolveRecursive( CoeffVector& solution , CoeffVector& rightHand , int Level );
    void   jacobiSmooth( CoeffVector& solution , CoeffVector& rightHand , int Level );
};

// constructor
template < class Grd , class domain_type >  
MultigridSolver<Grd,domain_type>::MultigridSolver( const  ProlongationStack<Grd>& prolongations, 
                                                   const  MlStack<Grd>& levelMatrices, 
                                                   double relax, 
                                                   double mgIter,
                                                   double mgSmoothings,
                                                   double maxTol
                                                 ): mProlongations( prolongations ), 
                                                    mLevelMatrices( levelMatrices )
{
  mNumOfLevels    = levelMatrices.levelMatrixStack.size();
  mJacRelaxFactor = relax;
  mMgIter         = mgIter;
  mMgSmoothings   = mgSmoothings;
  mMaxTol         = maxTol;
}

//member mgSolve()
template < class Grd , class domain_type >
void MultigridSolver<Grd,domain_type>::mgSolve( domain_type& domain_type_solution , domain_type& domain_type_rightHand )
{
  int Level = mNumOfLevels - 1;
  double relDefect;
  
  // references Dune::BlockVector
  CoeffVector &solution  = boost::fusion::at_c<0>(domain_type_solution.data);
  CoeffVector &rightHand = boost::fusion::at_c<0>(domain_type_rightHand.data);
  
  for(int i=0;i<mMgIter;i++)
  {
    mgSolveRecursive( solution , rightHand , Level );
    // calculate relative Defect
    CoeffVector defect(rightHand);
    mLevelMatrices.levelMatrixStack.front().mmv( solution , defect );
    relDefect = defect.two_norm()/rightHand.two_norm();
    if( relDefect < mMaxTol)
    { 
      cout << "MG converged: " << i << " iterations.   rel.Defect = "<< relDefect<< "\n";
      break; 
    }
    else if( i==mMgIter-1 ) 
      cout << "MG not converged: " << mMgIter << " iterations. rel. Defect = "<< relDefect<< '\n';
  }
}

//member mgSolveRecursive()
template < class Grd, class domain_type >
void MultigridSolver<Grd,domain_type>::mgSolveRecursive( CoeffVector& solution , CoeffVector& rightHand , int Level )
{
  if(Level==0)
  {
    // direct solve on Level 0
    Dune::MatrixAdapter<Dune::BCRSMatrix<M>, CoeffVector, CoeffVector>  matrixAdapter( mLevelMatrices.levelMatrixStack.back() ); 
    TrivialPreconditioner<CoeffVector> trivialPrecond;
    Dune::InverseOperatorResult resultInfo;
    CoeffVector rhsDummy(rightHand);
    Dune::CGSolver<CoeffVector> cg( matrixAdapter , trivialPrecond , 1.0e-14 , 1000 , 0 );
    cg.apply( solution , rhsDummy , resultInfo );
    mTimeDirectSolve += resultInfo.elapsed;
        
  }
  else
  {
    // pre-smoothing
    jacobiSmooth( solution , rightHand , Level );
    CoeffVector levelRightHand( mProlongations.prolStack[Level-1].M() ), tmpRhs(rightHand);
    mLevelMatrices.levelMatrixStack[mNumOfLevels-Level-1].mmv( solution , tmpRhs );// mmv(const X &x, Y &y) ... y -= A x
    // restriction residual:
    mProlongations.prolStack[Level-1].mtv( tmpRhs , levelRightHand );// mtv(const X &x, Y &y) ... y = A^T x 
    // initialize level-error
    CoeffVector levelError(levelRightHand);
    levelError *= 0;    
    mgSolveRecursive( levelError , levelRightHand , Level - 1 );
    // prolongate correction
    mProlongations.prolStack[Level-1].umv( levelError , solution );// umv(const X &x, Y &y) ... y += A x 
    // post-smoothing
    jacobiSmooth( solution , rightHand , Level );
  }
}


// member jacobiSmooth()
template < class Grd ,class domain_type>
void MultigridSolver<Grd,domain_type>::jacobiSmooth( CoeffVector& solution , CoeffVector& rightHand , int Level )
{
  JacobiSmoother jacobi( mLevelMatrices.levelMatrixStack[ mNumOfLevels-Level-1 ] , mMgSmoothings , mJacRelaxFactor );
  CoeffVector rhsDummy(rightHand);
  jacobi.pre( solution , rhsDummy );
  jacobi.apply( solution , rhsDummy );
  jacobi.post( solution );
}

// ------------------------------------------------------------------------------------------------------

// for output in console
void bcrsPrint( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& bcrsMatrix  )
{
  for(int i=0;i<bcrsMatrix.N();i++)
  {
    for(int k =0;k<bcrsMatrix.M();k++)
    {
      if(bcrsMatrix.exists(i,k)) cout << bcrsMatrix[i][k] << " ";
      else cout << "0 ";
    }
    cout << endl;
  }
  cout<<endl;
}

// ------------------------------------------------------------------------------------------------------


double energyNorm( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& matrix_, const Dune::BlockVector<Dune::FieldVector<double,1> >& vec_ )
{
  Dune::BlockVector<Dune::FieldVector<double,1> > Avec_( vec_.size() );
  Avec_=0;
  matrix_.mv(vec_,Avec_);// mv (const X &x, Y &y) ... y = A x 
  return sqrt(Avec_*vec_);
}

#endif
