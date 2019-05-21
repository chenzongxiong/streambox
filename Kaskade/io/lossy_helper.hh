#ifndef LOSSY_HELPER
#define LOSSY_HELPER

#include <map>

#include "dune/grid/common/grid.hh"
#include "dune/grid/common/entitypointer.hh"
#include "dune/grid/common/entity.hh"
#include "dune/grid/io/file/dgfparser/dgfparser.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/iteratorfacades.hh"
#include "dune/istl/matrixindexset.hh"

#include "linalg/conjugation.hh"


namespace Lossy_Detail {

  void bcrsPrint( const Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& bcrsMatrix  )
  {
    std::cout << std::endl;
    for(int i=0;i<bcrsMatrix.N();i++)
    {
      for(int k =0;k<bcrsMatrix.M();k++)
      {
	if(bcrsMatrix.exists(i,k)) 
	/*if( bcrsMatrix[i][k] > 0 )*/ std::cout << "(" << i << "," << k << ") : " << bcrsMatrix[i][k] << "\n";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }


  template < class Grid >
  class ParentalNodes
  {
      static const int dim = Grid::dimension;
      typedef typename Grid::LevelGridView::IndexSet::IndexType IndexType;
      typedef typename std::pair<IndexType,double> 	       FatherIndexValue;
      typedef typename std::vector<FatherIndexValue> 	       FatherSimplexPairs;
	  
    private:
      void findParents( const Grid& grid, int level );
      int maxLevel;
      
    public:
      std::vector< std::vector<FatherSimplexPairs> >  	levelParentsVector;
      // levelParentsVector[level][vertexIdx at level][fatherVertexNo] = pair(fatherVertexIdx, barycentricCoor)
      // level > 0 (level 0 vertices have no father), fatherVertexNo < dim+1 (in 2D for splitting edges at midpoint/red
      // refinement there should be only 2 father nodes?!)

      ParentalNodes( Grid const& grid ) : maxLevel(grid.maxLevel()), levelParentsVector(maxLevel)
      {
	for(int level=1; level<=maxLevel; level++)   findParents( grid, level );
      }
	  
  };
  
  template < class Grid >
  void ParentalNodes<Grid>::findParents( const Grid& grid, int level )
  {
    typedef typename Grid::template Codim<0>::LevelIterator   ElementLevelIterator;
    typedef typename Grid::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename Grid::template Codim<0>::EntityPointer   ElementPointer;   
    typedef typename Grid::LevelGridView::IndexSet::IndexType IndexType;
	  
    std::vector<FatherSimplexPairs> levelParents( grid.size(level,dim), FatherSimplexPairs(dim+1) );
	  
    
    // helpvector to avoid multiple calls of identical vertices
    std::vector<bool> toProcess(grid.size(dim), true);
    
    typename Grid::LevelGridView::IndexSet const& indexSet = grid.levelIndexSet(level) ;
    typename Grid::LevelGridView::IndexSet const& fatherIndexSet = grid.levelIndexSet(level-1);
//          typename Grid::LevelGridView::IndexSet const& indexSet = grid.levelView(level).indexSet() ;
//          typename Grid::LevelGridView::IndexSet const& fatherIndexSet = grid.levelView(level-1).indexSet();
          
    Dune::FieldVector<double, dim> geomInFatherCorner ;
    
    
    // iterate through all cells on current level
    ElementLevelIterator itEnd = grid.template lend<0>(level) ;
    for(ElementLevelIterator it = grid.template lbegin<0>(level) ; it != itEnd; ++it)
    {     
      std::vector<IndexType> fatherIndices(dim+1);
      ElementPointer itFather = it->father();
      for(int i=0 ; i<=dim ; i++) 
      {
	fatherIndices[i] = fatherIndexSet/*grid.levelView(level-1).indexSet()*/.subIndex(*itFather, i, dim);
      }
      
      // now regard corners of subentity
      for( int i=0 ; i<=dim ; i++ )
      {
	IndexType idx_child = indexSet/*grid.levelView(level).indexSet()*/.subIndex(*it, i, dim);
	
	if(toProcess[idx_child])
	{ 
	  toProcess[idx_child]=false;
	  
	  double relativeCoor = 1 ;
	  geomInFatherCorner = it->geometryInFather().corner(i);
	  
	  for( int component = dim/*-1*/ ; component >/*=*/ 0; component-- ) 
// 	  for( int component = 1 ; component <= dim; component++ ) 
	  {
	    double val = geomInFatherCorner[component-1];
	    relativeCoor -= val;
	    levelParents[idx_child][component/*+1*/] = std::make_pair(fatherIndices[component/*+1*/],val);
	  }
	  levelParents[idx_child][0] = std::make_pair(fatherIndices[0],relativeCoor);
	}
      }
    }
    
//     std::cout << "Element loop: " << time << "s\n";
//     levelParentsVector.push_back( levelParents );
    levelParentsVector[level-1] = levelParents;
  }
	
	
  template<class Grid>
  class Prolongation
  {
    typedef typename Grid::template Codim<Grid::dimension>::LevelIterator VertexLevelIterator;
    typedef typename Grid::LevelGridView::IndexSet::IndexType 	IndexType;
    typedef typename std::pair<IndexType,double> 		FatherIndexValue;
    typedef Dune::FieldMatrix<double,1,1> M;
	    
    public:
      Prolongation(Grid const& grid_) : parentalNodes(grid_) 
      {
	// test: compute prolongation matrix, dropping entries < 1 
	auto maxLevel = grid_.maxLevel();
	for(int level = 1 ; level <= maxLevel ; level++)
	{
	  int prolCols = grid_.size( level-1, Grid::dimension ) ;
	  int prolRows = grid_.size( level,   Grid::dimension ) ;
      
	  // we now create a BCRS Prolongation matrix
	  // we need three (triangle) or four (tetraeder) nnz's per row, because barycentric coordinates
	  // therefore:  prolRows*((Grid::dimension)+1)
	  Dune::BCRSMatrix<M> prolMatrix(prolRows,prolCols, prolRows*((Grid::dimension)+1),Dune::BCRSMatrix<M>::row_wise);
	  typedef Dune::BCRSMatrix<M>::CreateIterator Iter;
	  // now set sparsity pattern
	  for(Iter row=prolMatrix.createbegin(); row!=prolMatrix.createend(); ++row)
	  {
	    for( int k = 0 ; k <= Grid::dimension ; k++ ) 
	    {
	      row.insert( parentalNodes.levelParentsVector[level-1][row.index()][k].first );
	    }
	  }
	  // now set values
	  for(int row = 0 ; row < prolRows ; row++)
	  {
	    for( int k = 0 ; k <= Grid::dimension ; k++ ) 
	    {
	      if( parentalNodes.levelParentsVector[level-1][row][k].second > 0.9 ) // consider only 1s
		prolMatrix[row][ parentalNodes.levelParentsVector[level-1][row][k].first ]
		      = parentalNodes.levelParentsVector[level-1][row][k].second;
	    }
	  }
	  prolStack.push_back( prolMatrix );
// 	  std::cout << "modified prolongation matrix to level " << level << "\n";
// 	  bcrsPrint(prolMatrix);
	}
// 	std::cerr << "prolStack computed\n";
	
	// build level index mapper: levelIndexMapper[level l][idx on l][i] = idx on level l+i+1, i =0,..,maxLevel-l
	// not only on maxLevel, but for all levels from l+1 up to maxLevel, e.g. for use with solutions defined on coarser levels
	// as im MLSDC/PFASST
	indexOnNextLevel.resize(maxLevel, std::vector<size_t>());
        for( int l = 0; l < maxLevel; l++ )
	{
	  auto bcrsMatrix = prolStack[l]; // prolongation to level l+1 
	  indexOnNextLevel[l].resize(bcrsMatrix.M());
	  for(size_t k =0; k < bcrsMatrix.M(); k++) // iterate through columns
	  {
	    for( size_t i = 0; i < bcrsMatrix.N(); i++ ) // find existing row entry
	    {
	      if(bcrsMatrix.exists(i,k) && bcrsMatrix[i][k] > 0.9 ) // there are zeros and ones only in the matrix
	      {
		indexOnNextLevel[l][k] = i;
		break; // leave for i loop
	      }
	    }
	  }
	}
// 	std::cerr << "indexOnNextLevel computed\n";
	
	levelIndexMapper.resize(maxLevel);
	for( int l = 0 ; l < maxLevel; l++ )
	{
// 	  std::cout << "\nLEVEL " << l << "\n";
	  for( size_t k = 0; k < indexOnNextLevel[l].size(); k++ )
	  {
// 	    std::cout << "idx " << k << " -> ";
	    levelIndexMapper[l].resize(indexOnNextLevel[l].size());
	    size_t tmpIdx = k ;
	    for( int j = 0; j < maxLevel-l; j++ )
	    {
	      // determine index on level l+j+1
	      tmpIdx = indexOnNextLevel[l+j][tmpIdx];
	      levelIndexMapper[l][k].push_back(tmpIdx);
// 	      std::cout << tmpIdx << ", " ;
	    }
// 	    std::cout << "\n";
	  }
// 	  std::cout << "\n";
	}
// 	std::cerr << "levelIndexMapper computed\n";
      }
      
      ~Prolongation() {}
      
      
      // transfer source from level-1 to level
      // target is resized to have the correct size
      template<class CoeffVector> 
      void mv( unsigned level, CoeffVector const& source, CoeffVector& target)
      {
	unsigned int maxIndex = parentalNodes.levelParentsVector[level].size();
	target.resize(maxIndex);
	for( unsigned int i = 0 ; i <maxIndex ; i++ )
	{
	  target[i] = 0 ;
	  unsigned int maxIndexOld = parentalNodes.levelParentsVector[level][i].size();
	  for( unsigned int j = 0 ; j < maxIndexOld; j++ )
	  {
	    FatherIndexValue tmp = parentalNodes.levelParentsVector[level][i][j];
	    target[i] += source[ tmp.first ] * tmp.second ;
	  }
	}
      }
    
      // leaf has to be handeled separately!
      size_t getIndexOnLevel(int currentLevel, size_t currentIndex, int targetLevel)
      {
	if( targetLevel <= currentLevel ) return currentIndex;
	return levelIndexMapper[currentLevel][currentIndex][(targetLevel-currentLevel-1)];
      }
      
    private:
      Prolongation( Prolongation const & ) {} 
      
      ParentalNodes<Grid> parentalNodes;
      
      std::vector<Dune::BCRSMatrix<M>> prolStack;
      std::vector<std::vector<size_t>> indexOnNextLevel;
      std::vector<std::vector<std::vector<size_t>>> levelIndexMapper;
  };

  //------------------------------------------------------------------------------------------------------------------------

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
  template < class Grid >
  class ProlongationStack
  {
      typedef Dune::FieldMatrix<double,1,1> M;
  
    public:	
	std::vector<Dune::BCRSMatrix<M> > prolStack;
	ProlongationStack( const Grid& grid );
  };
  
  
  template < class Grid >
  ProlongationStack<Grid>::ProlongationStack( const Grid& grid )
  {
    using namespace Dune; 
    
    int maxLevel = grid.maxLevel();
        
    ParentalNodes<Grid> parentalNodes( grid ); 
    
    for(int level = 1 ; level <= maxLevel ; level++)
    {
      int prolCols = grid.size( level-1, Grid::dimension ) ;
      int prolRows = grid.size( level,   Grid::dimension ) ;
      
      // we now create a BCRS Prolongation matrix
      // we need three (triangle) or four (tetraeder) nnz's per row, because barycentric coordinates
      // therefore:  prolRows*((Grid::dimension)+1)
      BCRSMatrix<M> prolMatrix(prolRows,prolCols, prolRows*((Grid::dimension)+1),BCRSMatrix<M>::row_wise);
      typedef BCRSMatrix<M>::CreateIterator Iter;
      // now set sparsity pattern
      for(Iter row=prolMatrix.createbegin(); row!=prolMatrix.createend(); ++row)
      {
	for( int k = 0 ; k <= Grid::dimension ; k++ ) 
	{
	  row.insert( parentalNodes.levelParentsVector[level-1][row.index()][k].first );
	}
      }
      // now set values
      for(int row = 0 ; row < prolRows ; row++)
      {
	for( int k = 0 ; k <= Grid::dimension ; k++ ) 
	{
	  prolMatrix[row][ parentalNodes.levelParentsVector[level-1][row][k].first ]
		= parentalNodes.levelParentsVector[level-1][row][k].second;
	}
      }
      prolStack.push_back( prolMatrix );
    }
  }
  
  
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
    using namespace Kaskade;
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
  
  bool abscompare( double a, double b ) { return fabs( a ) < fabs( b ) ; }
  template <class T> bool greaterZero(T v) {return v > 0 ;}
  double ld( double val ) { return log(val)/log(2.0) ; }
	
} //namespace Lossy_Detail

#endif
