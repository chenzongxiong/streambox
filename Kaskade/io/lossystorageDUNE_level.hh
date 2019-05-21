#ifndef LOSSYSTORAGEDUNE_HH
#define LOSSYSTORAGEDUNE_HH

#include <map>
#include <cassert>

#include "dune/grid/common/grid.hh"
#include "dune/grid/common/entitypointer.hh"
#include "dune/grid/common/entity.hh"
#include "dune/grid/io/file/dgfparser/dgfparser.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/iteratorfacades.hh"
#include "dune/istl/matrixindexset.hh"

#include "rangecoder.hh"

#include "lossy_helper.hh"

// #define USEBOUNDS
#define USEFREQ
// #define LEVELWISE

template <class Grid, class CoeffVector>
class LossyStorage 
{
public:
  static const int dim = Grid::dimension ;
  // degree of interpolation polynomials; currently only linear interpolation is supported
  static const int order = 1 ; 
  
  // some typedefs used throughout the class
  typedef typename Grid::template Codim<dim>::LevelIterator VertexLevelIterator ;
  typedef typename Grid::template Codim<dim>::LeafIterator  VertexLeafIterator ;
  typedef typename Grid::LeafIndexSet IndexSet;
  typedef typename Grid::LevelGridView::IndexSet::IndexType IndexType;
  typedef typename Grid::LeafGridView::IndexSet::IndexType  LeafIndexType;
  
  typedef typename Grid::Traits::GlobalIdSet::IdType IdType;
  
  typedef unsigned long int ULong; 
	
  LossyStorage( int coarseLevel_, double aTol_) :  ps(NULL), coarseLevel(coarseLevel_), aTol(aTol_),
						   accumulatedEntropy(0), accumulatedOverhead(0),
						   accumulatedUncompressed(0), accumulatedCompressed(0)
  {   
/*#ifdef USEFREQ
    std::cout << "USES PRECOMPUTED FREQUENCIES FOR ENCODING A SINGLE FUNCTION!\nBEWARE OF THE OVERHEAD!\n";
    std::cout << "Due to implementation, not using pre-computed frequencies fails, as the vectors used for\n"
	      << "Alphabet are too large.\n";
#else
    std::cout << "Uses dictionary for encoding. Beware of the overhead!\n";
    std::cout << "Due to implementation, not using the dictionary fails, as the vectors used for\n"
	      << "Alphabet would be too large for fine quantization.\n";
#endif*/
  }

  ~LossyStorage()  { if( ps != NULL ) delete ps; }
  
  // returns the entropy (summed since the last call to resetEntropy/resetSizes)
  // of the data (= average symbol size in byte needed to store the file)
  double reportEntropy()
  {
    return accumulatedEntropy;
  }
  
  
  void resetEntropy()
  {
    accumulatedEntropy = 0;
  }
  
  
  // reports the overhead (summed since last call to resetOverhead/resetSizes)
  // i.e. interval bounds, frequencies, etc.
  double reportOverhead()
  {
    return accumulatedOverhead;
  }
  
  void resetOverhead()
  {
    accumulatedOverhead = 0;
  }
  
  // returns the compressed size 
  // which is (in the best case) the size of the encoded, compressed file
  // without the overhead from storing intervals etc.
  double reportCompressedSize()
  {
    return accumulatedCompressed;
  }
  
  void resetCompressedSize()
  {
    accumulatedCompressed = 0 ;
  }
  
  // returns the compressed size + overhead, which is (in the best case/in the limit)
  // the size of the compressed file (including all side information)
  double reportOverallSize()
  {
    return accumulatedCompressed + accumulatedOverhead;
  }
  
  // reset everything
  void resetSizes()
  {
    accumulatedEntropy = 0;
    accumulatedOverhead = 0;
    accumulatedCompressed = 0;
    accumulatedUncompressed = 0;
  }
  
  
  // returns the size needed for uncompressed storage of the data
  double reportUncompressedSize()
  {
    return accumulatedUncompressed;
  }
  
  void resetUncompressedSize()
  {
    accumulatedUncompressed = 0 ;
  }
  
  
  // returns the compression factor (=uncompressed size/compressed size)
  double reportRatio()
  {
    if( accumulatedCompressed + accumulatedOverhead > 0 )
      return accumulatedUncompressed / (accumulatedCompressed + accumulatedOverhead );
    
    return -1;
  }


  // TODO: has to be faster for adaptive meshes!
  void setup( Grid const& grid )
  {
//     boost::timer::cpu_timer timer;
    if( ps != NULL ) delete ps;
    ps = new Lossy_Detail::Prolongation<Grid>(grid);
//     std::cerr << "time for ps: " << timer.format() << "\n";
    
//     auto foo = Lossy_Detail::ProlongationStack<Grid>(grid); // for debugging, should use ps = new Prolongation<Grid> instead!
   
    typename Grid::GlobalIdSet const& idSet = grid.globalIdSet();  
    
    levelInfo.clear();
    for( int level = grid.maxLevel(); level >= 0; --level )
    {
      VertexLevelIterator lEnd =  grid.template lend <dim>( level );
      for( VertexLevelIterator it = grid.template lbegin <dim>( level ); it != lEnd; ++it )
      {
	levelInfo[idSet.id( *it )] = level ;
      }
    }
    
  }
  
  

  /**
  *  Encode a given state, e.g. the difference between to timesteps.
  *  Quantized values are written to a file specified by fn.
  */
  void encode( Grid const& grid, CoeffVector const& sol, std::string fn, double aTol_ = 0, int maxLevel_ = -1 )
  {
     std::ofstream out( fn.c_str(), std::ios::binary ) ;  
     encode( grid, sol, out, aTol_, maxLevel_);
     out.close();
  }

//   void encode( Grid const& grid, CoeffVector const& sol, std::ostream& out, double aTol_ = 0, int maxLevel_ = -1)
//   {         
//     // use maxLevel_ instead of maxLevel member
//     int maxLevel = maxLevel_;
//     if( maxLevel < 0 || maxLevel >= grid.maxLevel() ) 
//     {
//       maxLevel = grid.maxLevel() ;
//       encode( grid, sol, out, aTol_, maxLevel, grid.leafIndexSet());
//     }
//     else
//     {
//       encode( grid, sol, out, aTol_, maxLevel, grid.levelIndexSet(maxLevel));
//     }
//   }
//   
//   
//   template<class SolutionIndexSet>  
  void encode( Grid const& grid, CoeffVector const& sol, std::ostream& out, double aTol_, int maxLevel_/*, 
	       SolutionIndexSet const& solutionIndexSet*/)
  {         
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;
    
    // use maxLevel_ instead of maxLevel member
    int maxLevel = maxLevel_;
    if( maxLevel < 0 ) 
    {
      maxLevel = grid.maxLevel() ;
    }
        
    double overhead = 0, entropy = 0, compressed = 0 ;
//     size_t nNodes = solutionIndexSet.size(dim);
    size_t nNodes = grid.size(maxLevel, dim);
       
    
    
    typename Grid::GlobalIdSet const& idSet = grid.globalIdSet();
//     typename Grid::LeafGridView::IndexSet const& leafIndexSet = grid.leafIndexSet();
    
    
    std::vector<long int> indices ;
    std::vector<double> predictionError ;
    
    std::vector<std::vector<long int> > levelIndices;

//    std::cout << "\nLEVEL 0\n";
    
    VertexLevelIterator itEnd = grid.template lend<dim>(coarseLevel);
//     for( VertexLevelIterator it = grid.template lbegin<dim>( coarseLevel ) ;  it != itEnd; ++it )
    for( size_t i = 0; i < grid.size(coarseLevel, dim); i++ )
    {
//       predictionError.push_back( -sol[/*leafIndexSet*/solutionIndexSet.index(*it)] ) ;
//      std::cout << -sol[leafIndexSet.index(*it)] << " (" <<  leafIndexSet.index(*it) << ") ";
      predictionError.push_back( -sol[ ps->getIndexOnLevel(coarseLevel, i, maxLevel) ] ); 
      // if sol lives on leaf and leaf is not maxlevel (i.e. adaptive mesh) this is not covered yet!
    }    

    quantize( predictionError, indices) ;    
    reconstruct( predictionError, indices ) ;
    
    levelIndices.push_back(indices);
  
    CoeffVector reconstruction(grid.size(coarseLevel, dim) ) ;
    
    // assign reconstructed coarse grid values    
    size_t nEntries = sol.size();  
    std::vector<long int> intervalIndicesTmp( nEntries );
    long int minIndex = 0;
    
    size_t vertexCount = 0 ;
//     for( VertexLevelIterator it = grid.template lbegin<dim>(coarseLevel); it != itEnd; ++it)
    for( size_t i = 0; i < grid.size(coarseLevel, dim); i++ )
    {
      reconstruction[vertexCount] = -predictionError[vertexCount] ;
      
      long int tmpIndex = indices[vertexCount];
      intervalIndicesTmp[/*solutionIndexSet.index(*it)*/ps->getIndexOnLevel(coarseLevel, i, maxLevel)] = tmpIndex;
      if( tmpIndex < minIndex ) minIndex = tmpIndex;
            
      vertexCount++ ;
    }
      
    CoeffVector prediction;
        
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
      ps->mv(l, reconstruction, prediction) ;
//      std::cout << "\n\nTO LEVEL " << l+1 << "\n";
            
      std::vector<std::vector<size_t> > vertexInfo( prediction.size(), std::vector<size_t>(3) );
      
      predictionError.clear();
      predictionError.reserve( prediction.size() ); // avoid reallocations
      
      vertexCount = 0 ;
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = grid.levelGridView(l+1).indexSet();

      itEnd = grid.template lend<dim>(l+1);
      for( VertexLevelIterator it = grid.template lbegin<dim>(l+1);  it != itEnd; ++it)
      {
        unsigned char vertexLevel = levelInfo[idSet.id( *it )] ;
	auto levelIdx = levelIndexSet.index(*it);
	if( vertexLevel == l+1 )  
        {	  
// 	  auto solIdx =  /*leafIndexSet*/solutionIndexSet.index(*it);
// 	  if( ! solutionIndexSet.contains(*it) ) 
// 	    std::cerr << "LEVEL " << l+1 << " -- Entity idx " << levelIndexSet.index(*it) 
// 		      << " not contained in solutionIndexSet (determined idx = " << leafIdx << ")!\n" ;
	  auto solIdx = ps->getIndexOnLevel(l+1, levelIdx, maxLevel);
	  predictionError.push_back(prediction[/*levelIndexSet.index(*it)*/levelIdx] - sol[/*solutionIndexSet.index(*it)*/solIdx] );
// 	  std::cout << prediction[levelIndexSet.index(*it)] -sol[leafIndexSet.index(*it)] << " " ;
        }
        vertexInfo[vertexCount][0] = levelIdx;
	vertexInfo[vertexCount][1] = /*solutionIndexSet.index(*it); */ps->getIndexOnLevel(l+1, levelIdx, maxLevel);
	vertexInfo[vertexCount][2] = vertexLevel ;
	
//	std::cout << sol[leafIndexSet.index(*it)] << " (" << leafIndexSet.index(*it) << ") ";
	
        vertexCount++;
      }
           
      quantize( predictionError, indices) ;
      reconstruct( predictionError, indices) ;
    
      levelIndices.push_back(indices);
      
      if( (l+1) < maxLevel )   
      {
        // prepare prediction on next level -- use reconstructed values
        reconstruction.resize( prediction.size() );
        vertexCount = 0 ;

	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
        {
          // correction only for the nodes on level l+1
	  unsigned char vertexLevel = vertexInfo[ii][2]; 
	  IndexType levelIndex = vertexInfo[ii][0];
          if( vertexLevel < l+1 ) 
          { 
            reconstruction[levelIndex] =  prediction[levelIndex];	    
          }
          else
          {
	    long int tmpIndex = indices[vertexCount];
	    intervalIndicesTmp[vertexInfo[ii][1] ] = tmpIndex;
	    if( tmpIndex < minIndex ) minIndex = tmpIndex;
	    reconstruction[levelIndex] = prediction[levelIndex]-predictionError[vertexCount] ;
	    
            vertexCount++;
          }
        }
      }
      else
      {
	vertexCount = 0 ;
	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
	{
	  unsigned char vertexLevel = vertexInfo[ii][2];//levelInfo[idSet.id( *it )] ;
	  if( vertexLevel < l+1 )  continue; 
	  
	  long int tmpIndex = indices[vertexCount];
	  intervalIndicesTmp[vertexInfo[ii][1] ] = tmpIndex;
	  if( tmpIndex < minIndex ) minIndex = tmpIndex;
	  vertexCount++ ;
	}
      }
    }
    
      
    std::vector<ULong> intervalIndices( nEntries ) ;
    for( size_t i = 0; i <  nEntries; i++ ) intervalIndices[i] = intervalIndicesTmp[i] - minIndex ;
    
    out.write(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;
    overhead += sizeof(long int);
//     std::cout << "minIndex = " << minIndex << "\n"; std::cout.flush();
    
    std::map<unsigned long, unsigned long> freqMap;	  
    unsigned int nnz = 0;
    for( size_t i = 0 ; i < nEntries ; i++ ) 
    {
      if( freqMap.find(intervalIndices[i] ) != freqMap.end() ) // already there, increment
      {
        freqMap[intervalIndices[i]]++ ;
      }
      else // new nonzero, add to map
      {
        nnz++;
        freqMap[intervalIndices[i]] = 1;
      }
    }
    
    out.write(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;
    overhead += sizeof(unsigned int);
//     std::cout << "nnz = " << nnz << "\n"; std::cout.flush();

    std::vector<unsigned long> frequenciesForRangeCoder(nnz);
  
    std::map<unsigned long, unsigned long> dictMap;
    unsigned long curNo = 0 ;
	  
    std::map<unsigned long,unsigned long>::iterator mapIt;
    std::map<unsigned long,unsigned long>::iterator mapEnd = freqMap.end();
    for( mapIt = freqMap.begin(); mapIt != mapEnd; ++mapIt )
    {
      unsigned long lv = mapIt->first;
      frequenciesForRangeCoder[curNo] = mapIt->second;
      dictMap.insert( dictMap.begin(), std::pair<unsigned long, unsigned long>(lv, curNo) ) ;
      out.write(reinterpret_cast<char*>( &lv ), sizeof(unsigned long) ) ;
      overhead += sizeof(unsigned long);
      curNo++;
    }


    for( unsigned int i = 0 ; i < nnz ; i++ ) 
    {
#ifdef USEFREQ
      out.write(reinterpret_cast<char*>( &frequenciesForRangeCoder[i] ), sizeof(unsigned long) ) ; 
      overhead += sizeof(unsigned long);
#endif
      
      double tmp =  -ld(frequenciesForRangeCoder[i]/((double)nNodes)) *
			frequenciesForRangeCoder[i]/((double)nNodes);
      entropy += tmp;
      compressed += tmp;
    }
    accumulatedCompressed += compressed/8*nNodes;
    compressed = 0;

#ifndef USEFREQ
    frequenciesForRangeCoder.clear();
    freqMap.clear();
#endif
    
    std::cout.flush();
#ifdef USEFREQ
    Alphabet<unsigned long> alphabet( frequenciesForRangeCoder.begin(), frequenciesForRangeCoder.end() ) ;
    RangeEncoder<unsigned long> encoder( out ) ;

    for( size_t i = 0 ; i < nEntries; i++ ) 
    {
      encodeSymbol( encoder, alphabet, dictMap.find(intervalIndices[i])->second/* dictMap[intervalIndices[i]]*/ );
    }

#else

    
    size_t symbolCounter = 0 ;                
    std::vector<unsigned long> count(nnz, 1) ; 
    Alphabet<unsigned long> alphabet( count.begin(), count.end() ) ;
    
    RangeEncoder<unsigned long> encoder(out) ;
    
    for( size_t i = 0 ; i < nEntries; i++ ) 
    {
      encodeSymbol( encoder, alphabet, dictMap.find(intervalIndices[i])->second);
      ++count[dictMap[intervalIndices[i]]];
      ++symbolCounter ;
      if (symbolCounter>0.1*nEntries)
      {
	alphabet.update(count.begin(),count.end());
	symbolCounter=0;
      }
    }   
#endif

    accumulatedEntropy += entropy/8 ;
    accumulatedOverhead += overhead; 
    accumulatedUncompressed += 8*nNodes;
    
    aTol = aTolSave;
  }
  
  /**
  *  Decode quantized values and store in a VariableSet::Representation.
  *  The file to be read is specified by fn.
  */
  void decode( Grid const& gridState, CoeffVector& sol, std::string fn, double aTol_ = 0, int maxLevel_ = -1 ) 
  {
    std::ifstream in( fn.c_str(), std::ios::binary ) ;  
    decode(gridState, sol, in, aTol_, maxLevel_);
    in.close() ;
  }
  
//   void decode( Grid const& gridState, CoeffVector& sol, std::istream& in, double aTol_ = 0, int maxLevel_ = -1 ) 
//   {
//     // use maxLevel_ instead of maxLevel member
//     int maxLevel = maxLevel_;
//     if( maxLevel < 0 || maxLevel >= gridState.maxLevel() ) 
//     {
//       maxLevel = gridState.maxLevel() ;
//       decode( gridState, sol, in, aTol_, maxLevel, gridState.leafIndexSet());
//     }
//     else
//     {
//       decode( gridState, sol, in, aTol_, maxLevel, gridState.levelIndexSet(maxLevel));
//     }
//   }
// 
//   template<class SolutionIndexSet>
  void decode( Grid const& gridState, CoeffVector& sol, std::istream& in, double aTol_, int maxLevel_/*, 
	       SolutionIndexSet const& solutionIndexSet */) 
  {
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;
    
    int maxLevel = maxLevel_;
    if( maxLevel < 0 ) maxLevel = gridState.maxLevel();
            
    // prepare prediction
    typename Grid::GlobalIdSet const& idSet = gridState.globalIdSet();
//     IndexSet const& indexSet = gridState.leafIndexSet();
    
    std::vector<double> values ;
    std::vector<long int> intervalIndices( gridState.size(maxLevel, dim)/*solutionIndexSet.size(dim) */);
    size_t nEntries = intervalIndices.size();
    
    VertexLevelIterator itEnd = gridState.template lend<dim>(coarseLevel);
    
    // read in indices from file   
    long int minIndex ;
    in.read(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;   
//     std::cout << "minIndex = " << minIndex << "\n"; std::cout.flush();

    unsigned int nnz ;
    in.read(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;   // # non-empty intervals
//     std::cout << "nnz = " << nnz << "\n"; std::cout.flush();
    
    std::vector<long int> dictionary( nnz ) ;
    for( int i = 0 ; i < nnz ; i++ )
    {
      in.read(reinterpret_cast<char*>( &dictionary[i] ), sizeof(unsigned long) ) ; // existing intervals (dictionary)
    }
#ifdef USEFREQ     
    std::vector<unsigned long> frequencies( nnz, 0 ) ;
    for( int i = 0 ; i < nnz ; i++ )
    {
      in.read(reinterpret_cast<char*>( &frequencies[i] ), sizeof(unsigned long) ) ; // frequencies
    } 
      
    Alphabet<unsigned long> alphabet( frequencies.begin(), frequencies.end() ) ;    
    try
    {
      RangeDecoder<unsigned long> decoder( in ) ;
      for( int i = 0 ; i < intervalIndices.size() ; i++ ) 
      {
	unsigned long s = decodeSymbol(decoder,alphabet) ;
	intervalIndices[i] = dictionary[s] + minIndex ;
      }
    }
    catch( std::ios_base::failure& e )
    {
      if (in.rdstate() & std::ifstream::eofbit) 
      {
	std::cout << "EOF reached.\n";
      }
      else
      {
	std::cout  << " Decoding error\n" << e.what() << "\n"; 
      }
    }
#else   
    size_t symbolCounter = 0 ;

    std::vector<unsigned long> symcount(nnz, 1) ;
    Alphabet<unsigned long> alphabet( symcount.begin(), symcount.end() ) ;
    
    try
    {
      RangeDecoder<unsigned long> decoder( in ) ;
      for( size_t i = 0 ; i < nEntries ; i++ ) 
      {
	unsigned long s = decodeSymbol(decoder,alphabet) ;
	intervalIndices[i] = dictionary[s] + minIndex ;
	
	++symcount[s];
	++symbolCounter ;
	
	if (symbolCounter>0.1*nEntries) 
	{
	  alphabet.update(symcount.begin(),symcount.end());
	  symbolCounter=0;
	}
      }
    }
    catch( std::ios_base::failure& e )
    {
      if (in.rdstate() & std::ifstream::eofbit) 
      {
	std::cout << "EOF reached.\n";
      }
      else
      {
	std::cout  << " Decoding error\n" << e.what() << "\n"; 
      }
    }
#endif

//     in.close() ;
   
    // start reconstruction
    int vertexCount = 0;
    
    CoeffVector reconstruction( gridState.size(coarseLevel, dim) );
    sol.resize( gridState.size(maxLevel, dim) ); // if on leaf this has to be handeled differently

//     auto const& leafIndexSet = gridState.leafGridView().indexSet();

    
//    std::cout << "\n\nLEVEL 0\n";
//     for( VertexLevelIterator it = gridState.template lbegin<dim>( coarseLevel ); it != itEnd; ++it)
    for( size_t i = 0; i < gridState.size(coarseLevel, dim); i++ )
    {

      double recVal = reconstruct( intervalIndices[/*solutionIndexSet.index(*it)*/ps->getIndexOnLevel(coarseLevel, i, maxLevel)]);
      reconstruction[/*gridState.levelGridView(coarseLevel).indexSet().index(*it)*/i] = -recVal ; 
            
      // store predicted values for the coarse grid in solution vector
      sol[/*solutionIndexSet.index(*it)*/ ps->getIndexOnLevel(coarseLevel, i, maxLevel)] = -recVal ;      
//      std::cout << -recVal << " " << " (" << leafIndexSet.index(*it) << ") ";
      vertexCount++ ;
    }
    
    
   
    CoeffVector prediction;
    // perform prediction and correction
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
      ps->mv( l, reconstruction, prediction ) ;
      
//      std::cout << "\n\nDECODE TO LEVEL " << l+1 << "\n";
      
      reconstruction.resize( prediction.size() );
      
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = gridState.levelGridView(l+1).indexSet();
            
      vertexCount = 0 ;

      itEnd = gridState.template lend<dim>(l+1);
      for( VertexLevelIterator it = gridState.template lbegin<dim>(l+1); it != itEnd ; ++it)
      {
	IndexType levelIndex = levelIndexSet.index(*it);
	if(levelInfo[gridState.globalIdSet().id(*it)] == l+1) //correct only vertices on level l+1
	{
	  reconstruction[levelIndex] = prediction[levelIndex] 
	  - reconstruct( intervalIndices[/*solutionIndexSet.index(*it)*/ps->getIndexOnLevel(l+1, levelIndex, maxLevel)]);
// 	  std::cout << prediction[levelIndex] - reconstruction[levelIndex] << " ";
	}
	else
	{
	  reconstruction[levelIndex] = prediction[levelIndex];
	}
	sol[/*solutionIndexSet.index(*it)*/ps->getIndexOnLevel(l+1, levelIndex, maxLevel)] = reconstruction[ levelIndex] ;
//	std::cout << sol[leafIndexSet.index(*it)] << " " << " (" << leafIndexSet.index(*it) << ") ";

	vertexCount++ ;
      }       
    }
    
    //exit(1);
    aTol = aTolSave ;
  }

  
  
private:
  LossyStorage( LossyStorage const& );
  LossyStorage& operator=( LossyStorage const& ) ;
    
  
  /** Helper method to perform the actual quantization for a whole level. */
  void quantize( std::vector<double> const& values, std::vector<long int>& indices)
  {           
    indices.clear() ;
    indices.resize( values.size(), 0 ) ;
    
    for( size_t i = 0 ; i < values.size() ; i++ )
    {
      indices[i] = static_cast<long int>( floor( values[i] / (2*aTol) + 0.5 ) ); 
    }
  }
  
  
    /** Helper method to perform the actual reconstruction of quantized values without lb, ub. */
  void reconstruct( std::vector<double>& values, std::vector<long int> const& indices)
  {   
    values.clear() ;
    values.resize( indices.size() ) ;
    for( size_t i = 0 ; i < indices.size() ; i++ )
    {
      values[i] = indices[i] * 2* aTol ;
    }
  }
  

  
  /** Helper method to perform the actual reconstruction of quantized values without lb, ub. */
  double reconstruct( long int const& index )
  {
    return index*2*aTol;
  }
  
  /** Helper method returning the base 2-logarithm. */
  double ld( double val ) { return log(val)/log(2.0) ; }
   
   Lossy_Detail::Prolongation<Grid> * ps;
   std::map<IdType, unsigned char> levelInfo ;
   int coarseLevel ;
   double aTol ;
   double accumulatedEntropy, accumulatedOverhead, accumulatedUncompressed, accumulatedCompressed;
   
   std::vector<std::vector<int>> idToLevelIdx;
} ;

#endif
