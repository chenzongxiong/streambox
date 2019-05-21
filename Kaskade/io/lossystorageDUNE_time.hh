#ifndef LOSSYSTORAGEDUNE_HH
#define LOSSYSTORAGEDUNE_HH

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

#include "rangecoder.hh"

#include "lossy_helper.hh"

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
  typedef typename Grid::LevelGridView::IndexSet::IndexType 	IndexType;
  typedef typename Grid::LeafGridView::IndexSet::IndexType 	LeafIndexType;
  
  
  typedef typename Grid::Traits::GlobalIdSet::IdType IdType;
  
  typedef unsigned long int ULong; 
	
  LossyStorage( int coarseLevel_, double aTol_ ) :  ps(NULL), coarseLevel(coarseLevel_), aTol(aTol_),
						    accumulatedEntropy(0), accumulatedOverhead(0),
						    accumulatedUncompressed(0), accumulatedCompressed(0)
  {
//     std::cout << "BEWARE: USES LB, UB FOR ENCODING TWO FUNCTIONS. USING ENCODE(X) NEEDS TO BE ADAPTED!!!\n";
//     std::cerr << "BEWARE: USES LB, UB FOR ENCODING TWO FUNCTIONS. USING ENCODE(X) NEEDS TO BE ADAPTED!!!\n";
  }
		  
  ~LossyStorage()  { if( ps != NULL ) delete ps; }
  
  
  // returns the entropy (summed since the last call to resetEntropy/resetSizes)
  // of the data (= average symbol size needed to store the file)
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
    boost::timer::cpu_timer timer;
    if( ps != NULL ) delete ps;
    ps = new Lossy_Detail::Prolongation<Grid>(grid);
//     ps = new Lossy_Detail::ProlongationStack<Grid>(grid); // for debugging, should use ps = new Prolongation<Grid> instead!

//     std::cout << "ps: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  

    timer.start();    
    maxLevel = grid.maxLevel();
    
    typename Grid::GlobalIdSet const& idSet = grid.globalIdSet();  
    

    levelInfo.clear();
//     IdType ID ;
    for( int level = maxLevel; level >= 0; --level )
    {
      VertexLevelIterator lEnd =  grid.template lend <dim>( level );
      for( VertexLevelIterator it = grid.template lbegin <dim>( level ); it != lEnd; ++it )
      {
// 	ID = idSet.id( *it ) ;
	levelInfo[/*ID*/idSet.id( *it )] = level ;
      }
    }
//     std::cout << "levelInfo: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  

  }
  
  

  /**
  *  Encode a given state, e.g. the difference between to timesteps.
  *  Quantized values are written to a file specified by the member variables 
  *  filePrefix and count, keeping track of the already encoded timesteps.
  */
  
  /** REVISION AUG. 13, 2013: DON'T USE LB, UB ANYMORE 
   *           AUG  23, 2013: DON'T USE PRECOMPUTED FREQUENCIES ANYMORE (in encoding/decoding a single function!)
   */
  
  /** TODO: pull out range encoding, to be able to handle temporal prediction more easily
   *        -> this nethod should just preform transform+quantization, encoding should be done
   *           separately! (how/where to write ub/lb etc?) 
   **/
  void encode( Grid const& grid, CoeffVector const& sol, std::string fn, double aTol_ = 0 )
  {
    boost::timer::cpu_timer timer;
    boost::timer::cpu_timer fooTimer;
     
     
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;
    
    double overhead = 0, entropy = 0, compressed = 0 ;
    size_t nNodes = grid.size(dim);
    
    std::ofstream out( fn.c_str(), std::ios::binary ) ;  

//     int maxLevel = grid.maxLevel() ;
    typename Grid::GlobalIdSet const& idSet = grid.globalIdSet();
    typename Grid::LeafGridView::IndexSet const& leafIndexSet = grid.leafIndexSet();
    
    // for adaptivity in space, prolongation matrices have to be recomputed for the current grid
    // --> call setup(grid) first!
//      setup(grid);
    
//      std::cout << "init: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  
     timer.start();
  
    std::vector</*ULong*/long int> indices ;
    std::vector<double> predictionError ;

    VertexLevelIterator itEnd = grid.template lend<dim>(coarseLevel);
    for( VertexLevelIterator it = grid.template lbegin<dim>( coarseLevel ) ;  it != itEnd; ++it )
    {
      predictionError.push_back( -sol[leafIndexSet.index(*it)] ) ;
    }    

//     double ub, lb;
//     ub = *std::max_element( predictionError.begin(), predictionError.end(), Lossy_Detail::abscompare ) ;
//     if( ub < 0 ) { lb = ub ; ub = -lb  ; }
//     else         { lb = -ub ; }
     
// //     out.write(reinterpret_cast<char*>( &lb ), sizeof(double) ) ; 
//     out.write(reinterpret_cast<char*>( &ub ), sizeof(double) ) ; 
//     overhead += sizeof(double);

    quantize( predictionError, indices/*, lb, ub*/ ) ;    
    reconstruct( predictionError, indices/*, lb, ub*/ ) ;
    
//     std::cout << "\nreconstructed prediction error for level 0:\n";
//     for( size_t i = 0; i < predictionError.size(); i++ )
//       std::cout << predictionError[i] << "\n";
//     std::cout << "\n";
    
  
//     CoeffVector reconstruction(ps->prolStack[coarseLevel].M()) ;
    CoeffVector reconstruction(grid.size(coarseLevel, dim) ) ;
    
    // assign reconstructed coarse grid values    
    size_t nEntries = sol.size();  
    std::vector<long int/*ULong*/> intervalIndicesTmp( nEntries );
    long int minIndex = 0;

    size_t vertexCount = 0 ;
//     itEnd = grid.template lend<dim>(coarseLevel);
    for( VertexLevelIterator it = grid.template lbegin<dim>(coarseLevel); it != itEnd; ++it)
    {
      reconstruction[vertexCount] = -predictionError[vertexCount] ;
      
      long int tmpIndex = indices[vertexCount];
      intervalIndicesTmp[leafIndexSet.index(*it)] = tmpIndex;//indices[vertexCount];
      if( tmpIndex < minIndex ) minIndex = tmpIndex;
      
      vertexCount++ ;
    }
      
//      std::cout << "coarse grid: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  
     timer.start();

    CoeffVector prediction;
    
    double quantTime = 0, mvTime = 0, prepTime = 0, predErrTime = 0;
    
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
//       CoeffVector prediction( ps->prolStack[l].N() );
//       CoeffVector prediction( grid.size(l+1,dim) );
//       ps->prolStack[l].mv( reconstruction, prediction ) ;

      fooTimer.start();
      ps->mv(l, reconstruction, prediction) ;
      
//       std::cout << "prediction vs solution vs error for level " << l+1 << "\n";
//       for( VertexLevelIterator it = grid.template lbegin<dim>(l+1);  it != itEnd; ++it)
//       {
// 	std::cout << prediction[grid.levelView(l+1).indexSet().index(*it)] << "\t\t" <<  sol[leafIndexSet.index(*it)] 
// 		  << "\t\t" <<  prediction[grid.levelView(l+1).indexSet().index(*it)] - sol[leafIndexSet.index(*it)]<< "\n";
//       }
//       std::cout << "\n";
      
//       prediction.resize( ps->prolStack[l].N() );
//       prediction = 0;
//       ps->prolStack[l].mv( reconstruction, prediction ) ;

      mvTime += (double)(fooTimer.elapsed().user)/1e9 ;
      
      std::vector<std::vector<size_t> > vertexInfo( prediction.size(), std::vector<size_t>(3) );
      
      fooTimer.start();
      predictionError.clear();
      predictionError.reserve( prediction.size() ); // avoid reallocations
      
      vertexCount = 0 ;
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = grid.levelView(l+1).indexSet();
		
//       ub = 0 ;
//       double tmpval;
      itEnd = grid.template lend<dim>(l+1);
      for( VertexLevelIterator it = grid.template lbegin<dim>(l+1);  it != itEnd; ++it)
      {
        unsigned char vertexLevel = levelInfo[idSet.id( *it )] ;
	if( vertexLevel == l+1 )  
        {
	  //         tmpval = prediction[levelIndexSet.index(*it)] - sol[leafIndexSet.index(*it)] ;
	  // 	if( fabs(tmpval) > fabs(ub) ) ub = tmpval;
	  predictionError.push_back(/*tmpval*/prediction[levelIndexSet.index(*it)] - sol[leafIndexSet.index(*it)] );
        }
        vertexInfo[vertexCount][0] = levelIndexSet.index(*it) ;
	vertexInfo[vertexCount][1] = leafIndexSet.index(*it) ;
	vertexInfo[vertexCount][2] = vertexLevel ;
	
        vertexCount++;
      }
      

//       // symmetric quantization
//       ub = *std::max_element( predictionError.begin(), predictionError.end(), Lossy_Detail::abscompare ) ;
      
      predErrTime += (double)(fooTimer.elapsed().user)/1e9;
      
//       if( ub < 0 ) { lb = ub ; ub = -lb  ; }
//       else         { lb = -ub ; }
//       overhead += sizeof(double);
     
// //       out.write(reinterpret_cast<char*>( &lb ), sizeof(double) ) ; 
//       out.write(reinterpret_cast<char*>( &ub ), sizeof(double) ) ; 

      fooTimer.start();
      quantize( predictionError, indices/*, lb, ub*/) ;
      reconstruct( predictionError, indices/*, lb, ub*/ ) ;
      quantTime += (double)(fooTimer.elapsed().user)/1e9;
    
      fooTimer.start();
      if( (l+1) < maxLevel )   
      {
        // prepare prediction on next level -- use reconstructed values
        reconstruction.resize( prediction.size() );
        vertexCount = 0 ;
//         itEnd = grid.template lend<dim>(l+1);
//         for( VertexLevelIterator it = grid.template lbegin<dim>(l+1); it != itEnd ; ++it)
	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
        {
          // correction only for the nodes on level l+1
	  unsigned char vertexLevel = vertexInfo[ii][2]; //levelInfo[idSet.id( *it )] ;
	  IndexType levelIndex = vertexInfo[ii][0];//levelIndexSet.index(*it) ;
          if( vertexLevel < l+1 ) 
          { 
            reconstruction[levelIndex] =  prediction[levelIndex];	    
          }
          else
          {
	    long int/*ULong*/ tmpIndex = indices[vertexCount];
	    intervalIndicesTmp[/*leafIndexSet.index(*it)*/vertexInfo[ii][1] ] = tmpIndex;//indices[vertexCount];
	    if( tmpIndex < minIndex ) minIndex = tmpIndex;
   
	    reconstruction[levelIndex] = prediction[levelIndex]-predictionError[vertexCount] ;
	    
            vertexCount++;
          }
        }
      }
      else
      {
	vertexCount = 0 ;
// 	itEnd = grid.template lend<dim>(l+1);
// 	for( VertexLevelIterator it = grid.template lbegin<dim>(l+1); it != itEnd; ++it)
	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
	{
	  unsigned char vertexLevel = vertexInfo[ii][2];//levelInfo[idSet.id( *it )] ;
	  if( vertexLevel < l+1 )  continue; 
	  
	  long int/*ULong*/ tmpIndex = indices[vertexCount];
	  intervalIndicesTmp[/*leafIndexSet.index(*it)*/vertexInfo[ii][1] ] = tmpIndex;//indices[vertexCount];
	  if( tmpIndex < minIndex ) minIndex = tmpIndex;
	  vertexCount++ ;
	}
      }
      prepTime += (double)(fooTimer.elapsed().user)/1e9;
    }
  
//      std::cout << "prediction: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  
//      std::cout << "          mv: " << mvTime << "s\n"
// 	       << "     predErr: " << predErrTime << "s\n"
// 	       << "       quant: " << quantTime << "s\n"
// 	       << "        prep: " << prepTime << "s\n";
    
    // In order to maximize efficiency of the range coder, precalculate the
    // frequency of each interval, befor encoding.
    
  
    //TODO: check size/data type for nnz, maxIndex etc!

    timer.start();
    fooTimer.start();
//     long int minIndex = *std::min_element( intervalIndicesTmp.begin(), intervalIndicesTmp.end() ) ;
      
    std::vector<ULong> intervalIndices( nEntries ) ;
    for( size_t i = 0; i <  nEntries; i++ ) intervalIndices[i] = intervalIndicesTmp[i] - minIndex ;
    
    out.write(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;
    overhead += sizeof(long int);
//     std::cout << "    minIndex: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";

    fooTimer.start();

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

//     out.write(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;
//     overhead += sizeof(unsigned int);

    std::vector<unsigned long> frequenciesForRangeCoder(nnz);
//     std::cout << "    freq: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";
  
    fooTimer.start();

    std::map<unsigned long, unsigned long> dictMap;
    unsigned long curNo = 0 ;
	  
    std::map<unsigned long,unsigned long>::iterator mapIt;
    std::map<unsigned long,unsigned long>::iterator mapEnd = freqMap.end();
    for( mapIt = freqMap.begin(); mapIt != mapEnd; ++mapIt )
    {
      unsigned long lv = mapIt->first;
      frequenciesForRangeCoder[curNo] = mapIt->second;
      dictMap.insert( dictMap.begin(), std::pair<unsigned long, unsigned long>(lv, curNo) ) ;
      curNo++;
//       out.write(reinterpret_cast<char*>( &lv ), sizeof(unsigned long) ) ;
//       overhead += sizeof(unsigned long);
    }
//     std::cout << "    dict: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";

	  
    fooTimer.start();

    for( unsigned int i = 0 ; i < nnz ; i++ ) 
    {
//       out.write(reinterpret_cast<char*>( &frequenciesForRangeCoder[i] ), sizeof(unsigned long) ) ; 
//       overhead += sizeof(unsigned long);
	
      double tmp =  -ld(frequenciesForRangeCoder[i]/((double)nNodes)) *
			frequenciesForRangeCoder[i]/((double)nNodes);
      entropy += tmp;
      compressed += tmp;
    }
    accumulatedCompressed += compressed/8*nNodes;
    compressed = 0;

//     std::cout << "    write freq: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";

//      std::cout << "preparation for rc: " << (double)(timer.elapsed().user)/1e9 << "s\n";
     timer.start();
    
    std::cout << "\nEncoding " << intervalIndices.size() << " values, " 
	      << frequenciesForRangeCoder.size() <<  " different values\n" ;  
    std::cout << "Overhead: " << overhead << " B\n";
//      
//     Alphabet<unsigned long> alphabet( frequenciesForRangeCoder.begin(), frequenciesForRangeCoder.end() ) ;
//     RangeEncoder<unsigned long> encoder( out ) ;
// 
//     for( size_t i = 0 ; i < nEntries; i++ ) 
//     {
//       encodeSymbol( encoder, alphabet, dictMap.find(intervalIndices[i])->second/* dictMap[intervalIndices[i]]*/ );
//     }
// //     std::cout << "entropy coding: " << (double)(timer.elapsed().user)/1e9 << "s\n";

    size_t maxIntervals = *std::max_element(intervalIndices.begin(), intervalIndices.end()) + 1;
    out.write(reinterpret_cast<char*>( &maxIntervals ), sizeof(size_t) ) ;
    overhead += sizeof(size_t);
    
    size_t symbolCounter = 0 ;                
    std::vector<unsigned long> count(maxIntervals, 1) ; 
    Alphabet<unsigned long> alphabet( count.begin(), count.end() ) ; 
    RangeEncoder<unsigned long> encoder(out) ;
    
    for( size_t i = 0 ; i < nEntries; i++ ) 
    {
      encodeSymbol( encoder, alphabet, intervalIndices[i] ); 
      ++count[intervalIndices[i]];
      ++symbolCounter ;
      if (symbolCounter>0.1*alphabet.size())
      {
	alphabet.update(count.begin(),count.end());
	symbolCounter=0;
      }
    }   
    
    accumulatedEntropy += entropy/8 ;
    accumulatedOverhead += overhead; 
    accumulatedUncompressed += 8*nNodes;
    
    aTol = aTolSave;
  }
  

	
  void flush(Grid const& /*grid*/)
  {}

  
  /**
  *  Decode quantized values and store in a VariableSet::VariableSet.
  *  The file to be read is specified by the member variables 
  *  filePrefix and count, keeping track of the already decoded timesteps.
  */
  void decode( Grid const& gridState, CoeffVector& sol, std::string fn, double aTol_ = 0 ) 
  {
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;
    
    std::ifstream in( fn.c_str(), std::ios::binary ) ;  
    
    // for adaptivity in space, prolongation matrices have to be recomputed for the current grid
    // --> call setup(gridState)
//     setup(gridState);
    
    // prepare prediction
    IndexSet const& indexSet = gridState.leafIndexSet();
    
    std::vector<double> values ;
    std::vector<long int> intervalIndices( gridState.size(dim) );
    
    // read lb, ub from file
//     std::vector<double> lb(maxLevel-coarseLevel+1), ub(maxLevel-coarseLevel+1) ;
//     for( int i = 0 ; i <= maxLevel-coarseLevel ; i++ )
//     {
// //       in.read( reinterpret_cast<char*>( &lb[i] ), sizeof(double) ) ; 
//       in.read( reinterpret_cast<char*>( &ub[i] ), sizeof(double) ) ; 
//       lb[i] = -ub[i];
//     }
    
    // read in indices from file   
    long int minIndex ;
    in.read(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;   
    
//     unsigned int nnz ;
//     in.read(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;   // # non-empty intervals
//     std::vector<long int> dictionary( nnz ) ;
//     for( int i = 0 ; i < nnz ; i++ )
//     {
//       in.read(reinterpret_cast<char*>( &dictionary[i] ), sizeof(unsigned long) ) ; // existing intervals (dictionary)
//     }
//     std::vector</*ULong*/unsigned long> frequencies( nnz, 0 ) ;
//     for( int i = 0 ; i < nnz ; i++ )
//     {
//       in.read(reinterpret_cast<char*>( &frequencies[i] ), sizeof(unsigned long) ) ; // frequencies
//     } 
//       
//     Alphabet</*ULong*/unsigned long> alphabet( frequencies.begin(), frequencies.end() ) ;    
//     try
//     {
//       RangeDecoder</*ULong*/unsigned long> decoder( in ) ;
//       for( int i = 0 ; i < intervalIndices.size() ; i++ ) 
//       {
// 	/*ULong*/unsigned long s = decodeSymbol(decoder,alphabet) ;
// 	intervalIndices[i] = dictionary[s] + minIndex ;
//       }
//     }
//     catch( std::ios_base::failure& e )
//     {
//       if (in.rdstate() & std::ifstream::eofbit) 
//       {
// 	std::cout << "EOF reached.\n";
//       }
//       else
//       {
// 	std::cout  << " Decoding error\n" << e.what() << "\n"; 
//       }
//     }
    size_t symbolCounter = 0 ;
    size_t maxIntervals;
    in.read(reinterpret_cast<char*>( &maxIntervals ), sizeof(size_t) ) ;   

    std::vector<unsigned long> symcount(maxIntervals, 1) ;
    Alphabet<unsigned long> alphabet( symcount.begin(), symcount.end() ) ;
    try
    {
      RangeDecoder<unsigned long> decoder( in ) ;
      for( size_t i = 0 ; i < intervalIndices.size() ; i++ ) 
      {
	unsigned long s = decodeSymbol(decoder,alphabet) ;
	intervalIndices[i] = s + minIndex ;
	++symcount[s];
	++symbolCounter ;
	if (symbolCounter>0.1*alphabet.size()) 
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
      
    in.close() ;
      
    VertexLevelIterator itEnd = gridState.template lend<dim>(coarseLevel);
//     std::vector<std::vector<ULong> > allIndices( maxLevel+1-coarseLevel ) ;
    
//     for( int l = coarseLevel ; l <= maxLevel ; l++ )
//     {
//       itEnd = gridState.template lend<dim>(l);
//       for ( VertexLevelIterator it = gridState.template lbegin<dim>(l); it!= itEnd ; ++it)
//       {
// 	allIndices[l-coarseLevel].push_back(intervalIndices[indexSet.index(*it)] ) ;           
//       }
//     }
   
    // start reconstruction
    int vertexCount = 0;
//     reconstruct( values, allIndices[0], lb[0], ub[0]);
    
//     CoeffVector reconstruction( ps->prolStack[coarseLevel].M() );
    CoeffVector reconstruction( gridState.size(coarseLevel, dim) );
    sol.resize( gridState.size(dim) );
    
//     itEnd = gridState.template lend<dim>(coarseLevel);
    
//     typename Grid::LevelGridView::IndexSet const& coarseIndexSet = gridState.levelIndexSet(coarseLevel);
    
    
    for( VertexLevelIterator it = gridState.template lbegin<dim>( coarseLevel ); it != itEnd; ++it)
    {
      double recVal = reconstruct( intervalIndices[indexSet.index(*it)]/*, lb[0], ub[0]*/);
      reconstruction[gridState.levelView(coarseLevel).indexSet().index(*it)] = -recVal ; 
      
//       reconstruction[/*coarseIndexSet*/gridState.levelView(coarseLevel).indexSet().index(*it)] = -values[vertexCount] ; 
      
      // store predicted values for the coarse grid in solution vector
      sol[ indexSet.index(*it) ] = -recVal ;
//       sol[ indexSet.index(*it) ] = -values[vertexCount] ;
      
      vertexCount++ ;
    }
    
   
    CoeffVector prediction;
    // perform prediction and correction
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
//       prediction.resize( ps->prolStack[l].N(), false );
//       prediction = 0;
//       CoeffVector prediction( gridState.size(l+1,dim) );
//       ps->prolStack[l].mv( reconstruction, prediction ) ;
      ps->mv( l, reconstruction, prediction ) ;
      
//       reconstruct( values, allIndices[l-coarseLevel+1], lb[l+1-coarseLevel], ub[l+1-coarseLevel] );     
      reconstruction.resize( prediction.size() );
      
//       typename Grid::LevelGridView::IndexSet const& levelIndexSet = gridState.levelIndexSet(l+1);
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = gridState.levelView(l+1).indexSet();
            
      vertexCount = 0 ;

      itEnd = gridState.template lend<dim>(l+1);
      for( VertexLevelIterator it = gridState.template lbegin<dim>(l+1); it != itEnd ; ++it)
      {
	IndexType levelIndex = /*gridState.levelView(l+1).indexSet()*/levelIndexSet.index(*it);
	if(levelInfo[gridState.globalIdSet().id(*it)] == l+1) //correct only vertices on level l+1
	{
	  reconstruction[levelIndex] = prediction[levelIndex] - //values[vertexCount] ;
		reconstruct( intervalIndices[indexSet.index(*it)]/*, lb[l+1-coarseLevel], ub[l+1-coarseLevel]*/);
	}
	else
	{
	  reconstruction[levelIndex] = prediction[levelIndex];
	}
	sol[indexSet.index(*it)] = reconstruction[ levelIndex] ;
	vertexCount++ ;
      }       
    }
    aTol = aTolSave ;
  }

  
  /** encode two FE functions at once for cardiac optimization;
      just encoding, not writing to disk
  */
  
  void encode( Grid const& grid, CoeffVector const& v, CoeffVector const& w,
	       std::vector<std::vector<long int> > & intervalIndicesTmp,
	       /*std::string fn_v, std::string fn_w, */ double aTol_ = 0 )
  {
    boost::timer::cpu_timer timer;
    boost::timer::cpu_timer fooTimer;
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;
    
    double overhead = 0, entropy = 0, compressed = 0 ;
    size_t nNodes = grid.size(dim);

    
//     std::ofstream out[2];
//     out[0].open(fn_v.c_str(), std::ios::binary ) ;
//     out[1].open(fn_w.c_str(), std::ios::binary ) ;

    typename Grid::GlobalIdSet const& idSet = grid.globalIdSet();
    typename Grid::LeafGridView::IndexSet const& leafIndexSet = grid.leafIndexSet();
    
    timer.start();
    
    std::vector<long int> indices_v ;
    std::vector<long int> indices_w ;
    std::vector<double> predictionError_v ;
    std::vector<double> predictionError_w ;
    
    VertexLevelIterator itEnd = grid.template lend<dim>(coarseLevel);
    for( VertexLevelIterator it = grid.template lbegin<dim>( coarseLevel ) ;  it != itEnd; ++it )
    {
      predictionError_v.push_back( -v[leafIndexSet.index(*it)] ) ;
      predictionError_w.push_back( -w[leafIndexSet.index(*it)] ) ;
    }    
       
    quantize( predictionError_v, indices_v) ;    
    reconstruct( predictionError_v, indices_v) ;
    
    quantize( predictionError_w, indices_w) ;    
    reconstruct( predictionError_w, indices_w) ;
    
    CoeffVector reconstruction_v( grid.size(coarseLevel, dim) ) ;
    CoeffVector reconstruction_w( reconstruction_v.size() ) ;
    
    // assign reconstructed coarse grid values    
    size_t nEntries = v.size();  
//     std::vector<std::vector<long int> > intervalIndicesTmp( 2, std::vector<long int>(nEntries) );
    intervalIndicesTmp.resize( 2, std::vector<long int>(nEntries) );
//     long int minIndex[2] = {0,0}; 
    
    size_t vertexCount = 0 ;
    for( VertexLevelIterator it = grid.template lbegin<dim>(coarseLevel); it != itEnd; ++it)
    {
      reconstruction_v[vertexCount] = -predictionError_v[vertexCount] ;
      reconstruction_w[vertexCount] = -predictionError_w[vertexCount] ;
      
      
      long int tmpIndex = indices_v[vertexCount];
      intervalIndicesTmp[0][leafIndexSet.index(*it)] = tmpIndex;
//       if( tmpIndex < minIndex[0] ) minIndex[0] = tmpIndex;
      
      tmpIndex = indices_w[vertexCount];
      intervalIndicesTmp[1][leafIndexSet.index(*it)] = tmpIndex;
//       if( tmpIndex < minIndex[1] ) minIndex[1] = tmpIndex;
      
      vertexCount++ ;
    }
    
    timer.start();
    
    CoeffVector prediction_v, prediction_w;
    
    double quantTime = 0, mvTime = 0, prepTime = 0, predErrTime = 0, boundsTime = 0, allocTime = 0,
	   vertexInfoTime = 0, levelIndexTime = 0, leafIndexTime = 0, levelInfoTime = 0;
    
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
      fooTimer.start();
      ps->mv(l, reconstruction_v, prediction_v) ;
      ps->mv(l, reconstruction_w, prediction_w) ;
      mvTime += (double)(fooTimer.elapsed().user)/1e9 ;
      
      fooTimer.start();
      size_t levelSize = prediction_v.size() ; 
      std::vector<std::vector<size_t> > vertexInfo( levelSize, std::vector<size_t>(3) );
      predictionError_v.clear();
      predictionError_v.reserve( levelSize); // avoid reallocations
      predictionError_w.clear();
      predictionError_w.reserve( levelSize); // avoid reallocations

      vertexCount = 0 ;
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = grid.levelView(l+1).indexSet();
      
      itEnd = grid.template lend<dim>(l+1);
      
      allocTime += (double)(fooTimer.elapsed().user)/1e9;
      fooTimer.start();
	
      for( VertexLevelIterator it = grid.template lbegin<dim>(l+1);  it != itEnd; ++it)
      {
	predErrTime += (double)(fooTimer.elapsed().user)/1e9;
	
	fooTimer.start();
	vertexInfo[vertexCount][0] = levelIndexSet.index(*it) ;
	levelIndexTime += (double)(fooTimer.elapsed().user)/1e9;
	
	fooTimer.start();
	vertexInfo[vertexCount][2] = levelInfo[idSet.id( *it )] ;
	levelInfoTime += (double)(fooTimer.elapsed().user)/1e9;
	
	fooTimer.start();
	if( vertexInfo[vertexCount][2] == l+1 )  
	{
	  fooTimer.start();
	  auto tmpLeafIxd = leafIndexSet.index(*it) ;
	  vertexInfo[vertexCount][1] = tmpLeafIxd;
	  leafIndexTime += (double)(fooTimer.elapsed().user)/1e9;
	  
	  predictionError_v.push_back(prediction_v[vertexInfo[vertexCount][0]] - v[tmpLeafIxd/*vertexInfo[vertexCount][1]*/] );
	  predictionError_w.push_back(prediction_w[vertexInfo[vertexCount][0]] - w[tmpLeafIxd/*vertexInfo[vertexCount][1]*/] ); 
	}
	vertexCount++;
      }
      predErrTime += (double)(fooTimer.elapsed().user)/1e9;
      
      fooTimer.start();
      quantize( predictionError_v, indices_v) ;
      reconstruct( predictionError_v, indices_v) ;
      quantize( predictionError_w, indices_w) ;
      reconstruct( predictionError_w, indices_w) ;
      quantTime += (double)(fooTimer.elapsed().user)/1e9;
      
      
      fooTimer.start();
      if( (l+1) < maxLevel )   
      {
	// prepare prediction on next level -- use reconstructed values
	reconstruction_v.resize( prediction_v.size() );
	reconstruction_w.resize( prediction_w.size() );
	
	vertexCount = 0 ;
	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
	{
	  // correction only for the nodes on level l+1
	  unsigned char vertexLevel = vertexInfo[ii][2]; 
	  IndexType levelIndex = vertexInfo[ii][0];
	  if( vertexLevel < l+1 ) 
	  { 
	    reconstruction_v[levelIndex] =  prediction_v[levelIndex];
	    reconstruction_w[levelIndex] =  prediction_w[levelIndex];	    
	  }
	  else
	  {
	    long int tmpIndex = indices_v[vertexCount];
	    intervalIndicesTmp[0][vertexInfo[ii][1] ] = tmpIndex;
// 	    if( tmpIndex < minIndex[0] ) minIndex[0] = tmpIndex;
	    
	    tmpIndex = indices_w[vertexCount];
	    intervalIndicesTmp[1][vertexInfo[ii][1] ] = tmpIndex;
// 	    if( tmpIndex < minIndex[1] ) minIndex[1] = tmpIndex;
	    
	    reconstruction_v[levelIndex] = prediction_v[levelIndex]-predictionError_v[vertexCount] ;
	    reconstruction_w[levelIndex] = prediction_w[levelIndex]-predictionError_w[vertexCount] ;
	    
	    vertexCount++;
	  }
	}
      }
      else
      {
	vertexCount = 0 ;
	for( size_t ii = 0 ; ii < vertexInfo.size(); ii++ ) 
	{
	  unsigned char vertexLevel = vertexInfo[ii][2];
	  if( vertexLevel < l+1 )  continue; 
	  
	  long int/*ULong*/ tmpIndex = indices_v[vertexCount];
	  intervalIndicesTmp[0][vertexInfo[ii][1] ] = tmpIndex;
// 	  if( tmpIndex < minIndex[0] ) minIndex[0] = tmpIndex;
	  
	  tmpIndex = indices_w[vertexCount];
	  intervalIndicesTmp[1][vertexInfo[ii][1] ] = tmpIndex;
// 	  if( tmpIndex < minIndex[1] ) minIndex[1] = tmpIndex;
	  
	  vertexCount++ ;
	}
      }
      prepTime += (double)(fooTimer.elapsed().user)/1e9;
    }

    
//          std::cout << "prediction: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  
//          std::cout << "          mv: " << mvTime << "s\n"
// 	           << "       alloc: " << allocTime << "s\n"
// 		   << "  levelIndex: " << levelIndexTime << "s\n"
// 		   << "   leafIndex: " << leafIndexTime << "s\n"
// 		   << "   levelInfo: " << levelInfoTime << "s\n"
// // 		   << "  vertexInfo: " << vertexInfoTime << "s\n"
//     	           << "     predErr: " << predErrTime << "s\n"
// 		   << "      bounds: " << boundsTime << "s\n"
//     	           << "       quant: " << quantTime << "s\n"
//     	           << "        prep: " << prepTime << "s\n";

    aTol = aTolSave;

    // now build difference to previous indices in integrate()
  }

  /** write quantized indices to files using range encoding */
  void write(std::vector<std::vector<long int> > const & intervalIndicesTmp,
	     std::string fn_v, std::string fn_w )
	     
  {
    std::ofstream out[2];
    out[0].open(fn_v.c_str(), std::ios::binary ) ;
    out[1].open(fn_w.c_str(), std::ios::binary ) ; 
    
    size_t nFunctions = intervalIndices.size();
    assert(nFunctions == 2);
    
    size_t nEntries = intervalIndicesTmp[0].size();
    
    boost::timer::cpu_timer timer;   
    boost::timer::cpu_timer fooTimer;

    long int minIndex[2] = {0,0}; 
    minIndex[0] = *std::min_element( intervalIndicesTmp[0].begin(), intervalIndicesTmp[0].end() ) ;
    minIndex[1] = *std::min_element( intervalIndicesTmp[1].begin(), intervalIndicesTmp[1].end() ) ;
    
    std::vector<ULong> intervalIndices( nEntries ) ;
    
    // first v, then w
    for( int function = 0 ; function < 2 ; function++ )
    {
    
      for( size_t i = 0; i <  nEntries; i++ ) 
	intervalIndices[i] = intervalIndicesTmp[function][i] - minIndex[function] ;
    
      out[function].write(reinterpret_cast<char*>( &minIndex[function] ), sizeof(/*ULong*/long int) ) ;
      overhead += sizeof(long int);
//         std::cout << "    minIndex: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";
    
      fooTimer.start();
    
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
      
      out[function].write(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;
      overhead += sizeof(unsigned int);
//       std::cout << "nonzero frequencies: " << nnz << "\n";
    
      std::vector<unsigned long> frequenciesForRangeCoder(nnz);
//         std::cout << "    freq: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";
    
      fooTimer.start();
    
      std::map<unsigned long, unsigned long> dictMap;
      unsigned long curNo = 0 ;
    
      std::map<unsigned long,unsigned long>::iterator mapIt;
      std::map<unsigned long,unsigned long>::iterator mapEnd = freqMap.end();
      for( mapIt = freqMap.begin(); mapIt != mapEnd; ++mapIt )
      {
	unsigned long lv = mapIt->first;
	frequenciesForRangeCoder[curNo] = mapIt->second;
	dictMap.insert( dictMap.begin(), std::pair<unsigned long, unsigned long>(lv, curNo) ) ;
	curNo++;
	out[function].write(reinterpret_cast<char*>( &lv ), sizeof(unsigned long) ) ;
	overhead += sizeof(unsigned long);
      }
//     std::cout << "    dict: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";
    
      fooTimer.start();        
      for( unsigned int i = 0 ; i < nnz ; i++ ) 
      {
	out[function].write(reinterpret_cast<char*>( &frequenciesForRangeCoder[i] ), sizeof(unsigned long) ) ; 
	overhead += sizeof(unsigned long);
	
	double tmp =  -ld(frequenciesForRangeCoder[i]/((double)nNodes)) *
			  frequenciesForRangeCoder[i]/((double)nNodes);
	entropy += tmp;
	compressed += tmp;/*/8 * frequenciesForRangeCoder[i];*/
      }
      accumulatedCompressed += compressed/8*nNodes;
      compressed = 0;
//       std::cout << "    write freq: " << (double)(fooTimer.elapsed().user)/1e9 << "s\n";
//       std::cout << "preparation for rc: " << (double)(timer.elapsed().user)/1e9 << "s\n";
      
      timer.start();
    
      Alphabet<unsigned long> alphabet( frequenciesForRangeCoder.begin(), frequenciesForRangeCoder.end() ) ;
      RangeEncoder<unsigned long> encoder( out[function] ) ;
      
      for( size_t i = 0 ; i < /*2**/nEntries; i++ ) 
      {
	  encodeSymbol( encoder, alphabet, dictMap.find(intervalIndices[i])->second/*[intervalIndices[i]]*/ );
      }  
      
//       std::cout << "range coder: " << (double)(timer.elapsed().user)/1e9 << "s\n";	  
    } // for( int function ...     
    
    accumulatedEntropy += entropy/8 ;
    accumulatedOverhead += overhead; 
    accumulatedUncompressed += 2*8*nNodes;    
  }  
  
  
  /** read quantized indices from files using range encoding */
  void read( size_t nEntries, // degrees of freedom to read (given by the grid)
	     std::vector<std::vector<long int> > & intervalIndices,
	     std::string fn_v, std::string fn_w )
  {
    std::ifstream in[2];
    in[0].open( fn_v.c_str(), std::ios::binary ) ;
    in[1].open( fn_w.c_str(), std::ios::binary ) ;  
    
    intervalIndices.resize( 2, std::vector<long int>(nEntries) );
  
    long int minIndex[2];
    for( int function = 0 ; function < 2 ; function++)
    {
      // read in indices from file   
      in[function].read(reinterpret_cast<char*>( &minIndex[function] ), sizeof(long int) ) ;  
    
      unsigned int nnz ;
      in[function].read(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;   // # non-empty intervals
      std::vector<ULong> dictionary( nnz ) ;
      for( int i = 0 ; i < nnz ; i++ )
      {
	in[function].read(reinterpret_cast<char*>( &dictionary[i] ), sizeof(unsigned long) ) ;
      }
      std::vector<ULong> frequencies( nnz, 0 ) ;
      for( int i = 0 ; i < nnz ; i++ )
      {
	in[function].read(reinterpret_cast<char*>( &frequencies[i] ), sizeof(/*ULong*/unsigned long) ) ; // frequencies
      } 
           
      Alphabet<unsigned long> alphabet( frequencies.begin(), frequencies.end() ) ;    
      
      try
      {
	RangeDecoder<unsigned long> decoder( in[function] ) ;
	for( int i = 0 ; i < nEntries ; i++ ) 
	{
	  unsigned long s = decodeSymbol(decoder,alphabet) ;
	  intervalIndices[function][i] = dictionary[s] + minIndex[function] ;
	}
      }
      catch( std::ios_base::failure& e )
      {
	if (in[function].rdstate() & std::ifstream::eofbit) 
	{
	  std::cout << "EOF reached.\n";
	}
	else
	{
	  std::cout  << " Decoding error\n" << e.what() << "\n"; 
	}
      }
    in[function].close() ;
    } // for( int function = ...    
  }
  
  /** decode two FE functions at once for cardiac optimization;
      reads each function from a separate file */
  
  void decode( Grid const& gridState, CoeffVector& v,  CoeffVector& w,
	       std::vector<std::vector<long int> const& intervalIndices,
	       /*std::string fn_v, std::string fn_w,*/ double aTol_ = 0) 
  {
    double aTolSave = aTol ;
    if( aTol_ > 0 ) aTol = aTol_ ;

    assert( intervalIndices.size() == 2 );
    size_t nEntries = intervalIndices[0].size();
    assert( nEntries = gridState.size(dim);
    
    // prepare prediction
    IndexSet const& indexSet = gridState.leafIndexSet();
    std::vector<double> values_v, values_w ;

    VertexLevelIterator itEnd = gridState.template lend<dim>(coarseLevel);
   
    // start reconstruction
    int vertexCount = 0;

    CoeffVector reconstruction_v( gridState.size(coarseLevel, dim) );
    v.resize( gridState.size(dim) );
    CoeffVector reconstruction_w( gridState.size(coarseLevel, dim) );
    w.resize( gridState.size(dim) );

    
    for( VertexLevelIterator it = gridState.template lbegin<dim>( coarseLevel ); it != itEnd; ++it)
    {
      double recVal = reconstruct( intervalIndices[0][indexSet.index(*it)]/*, lb[0][0], ub[0][0]*/);
      reconstruction_v[gridState.levelView(coarseLevel).indexSet().index(*it)] = -recVal ;      
      v[ indexSet.index(*it) ] = -recVal ;
      
      recVal = reconstruct( intervalIndices[1][indexSet.index(*it)]/*, lb[1][0], ub[1][0]*/);
      reconstruction_w[gridState.levelView(coarseLevel).indexSet().index(*it)] = -recVal ;      
      w[ indexSet.index(*it) ] = -recVal ;
      
      vertexCount++ ;
    }
    
    
    CoeffVector prediction_v, prediction_w;
    // perform prediction and correction
    for( int l = coarseLevel ; l < maxLevel ; l++ )
    {
      ps->mv( l, reconstruction_v, prediction_v ) ;
      ps->mv( l, reconstruction_w, prediction_w ) ;
      
      reconstruction_v.resize( prediction_v.size() );
      reconstruction_w.resize( prediction_w.size() );
      
      //       typename Grid::LevelGridView::IndexSet const& levelIndexSet = gridState.levelIndexSet(l+1);
      typename Grid::LevelGridView::IndexSet const& levelIndexSet = gridState.levelView(l+1).indexSet();
      
      vertexCount = 0 ;
      
      itEnd = gridState.template lend<dim>(l+1);
      for( VertexLevelIterator it = gridState.template lbegin<dim>(l+1); it != itEnd ; ++it)
      {
	IndexType levelIndex = levelIndexSet.index(*it);
	if(levelInfo[gridState.globalIdSet().id(*it)] == l+1) //correct only vertices on level l+1
	{
	  reconstruction_v[levelIndex] = prediction_v[levelIndex] - 
		reconstruct( intervalIndices[0][indexSet.index(*it)]/*, lb[0][l+1-coarseLevel], ub[0][l+1-coarseLevel]*/);
	  reconstruction_w[levelIndex] = prediction_w[levelIndex] - 
		reconstruct( intervalIndices[1][indexSet.index(*it)]/*, lb[1][l+1-coarseLevel], ub[1][l+1-coarseLevel]*/);
		
	}
	else
	{
	  reconstruction_v[levelIndex] = prediction_v[levelIndex];
	  reconstruction_w[levelIndex] = prediction_w[levelIndex];
	}
	v[indexSet.index(*it)] = reconstruction_v[ levelIndex] ;
	w[indexSet.index(*it)] = reconstruction_w[ levelIndex] ;
	vertexCount++ ;
      }       
    }
    aTol = aTolSave ;
  }
  
  
  
private:
  LossyStorage( LossyStorage const& );
  LossyStorage& operator=( LossyStorage const& ) ;
    
  /** Helper method to perform the actual quantization for a whole level. */
  void quantize( std::vector<double> const& values, std::vector</*ULong*/long int>& indices, double lb, double ub )
  {           
    double range = ub-lb;
    long int nIntervals = (long int) ceil( range / aTol / 2.0 ) ; 
    if( nIntervals % 2 == 0 ) nIntervals += 1 ;
//     
// //     std::cout << "[" << lb << ", " << ub << "] -> " << nIntervals << " intervals,   " <<  values.size() << "nodes\n";
//     
    indices.clear() ;
    indices.resize( values.size(), 0 ) ;
//     
//     // handle range == 0 !
    if( range > 0 )
    {
      for( size_t i = 0 ; i < values.size() ; i++ )
      {
// 	indices[i] = /*(ULong)*/ floor( values[i] / (2*aTol) + 0.5 ) ; // without lb/ub
	indices[i] = ((values[i] - lb) * nIntervals / range ) + 1 ; 
      }
    }
  }
  
  /** Helper method to perform the actual quantization for a whole level without lb, ub. */
  void quantize( std::vector<double> const& values, std::vector</*ULong*/long int>& indices)
  {           
    indices.clear() ;
    indices.resize( values.size(), 0 ) ;
    
    for( size_t i = 0 ; i < values.size() ; i++ )
    {
	indices[i] = /*(ULong)*/ floor( values[i] / (2*aTol) + 0.5 ) ; // without lb/ub
    }
  }
  
  /** Helper method to perform the actual reconstruction of quantized values. */
  void reconstruct( std::vector<double>& values, std::vector</*ULong*/long int> const& indices, double lb, double ub )
  {
    double range = ub-lb;
    long int nIntervals = (long int) ceil( range / aTol / 2.0 ) ; 
    if( nIntervals % 2 == 0 ) nIntervals += 1 ;
    
    double delta;
    if( nIntervals <= 0 ) delta = 0 ;
    else                  delta = range / nIntervals ;
    
    values.clear() ;
    values.resize( indices.size() ) ;
    for( size_t i = 0 ; i < indices.size() ; i++ )
    {
// 	values[i] = indices[i] * 2* aTol ;
      values[i] = lb + (/*2**/indices[i]-0.5/*1*/) * delta/*/2 */;
    }
  }
  
    /** Helper method to perform the actual reconstruction of quantized values without lb, ub. */
  void reconstruct( std::vector<double>& values, std::vector</*ULong*/long int> const& indices)
  {   
    values.clear() ;
    values.resize( indices.size() ) ;
    for( size_t i = 0 ; i < indices.size() ; i++ )
    {
      values[i] = indices[i] * 2* aTol ;
    }
  }
  
  /** Helper method to perform the actual reconstruction of quantized values. */
  double reconstruct( /*ULong*/long int const& index, double lb, double ub)
  {
    double range = ub-lb;
    long int nIntervals = (long int) ceil( range / aTol / 2.0 ) ; 
    if( nIntervals % 2 == 0 ) nIntervals += 1 ;
    
    double delta;
    if( nIntervals <= 0 ) delta = 0 ;
    else                  delta = range / nIntervals ;
    
    return lb + (index-0.5) * delta;
//     return index*2*aTol;
  }
  
  /** Helper method to perform the actual reconstruction of quantized values without lb, ub. */
  double reconstruct( /*ULong*/long int const& index )
  {
    return index*2*aTol;
  }
  
  /** Helper method returning the base 2-logarithm. */
  double ld( double val ) { return log(val)/log(2.0) ; }
   
//    Lossy_Detail::ProlongationStack<Grid> * ps; // for debugging
   Lossy_Detail::Prolongation<Grid> * ps;
   std::map<IdType, unsigned char> levelInfo ;
   int coarseLevel, maxLevel ;
   double aTol ;
   double accumulatedEntropy, accumulatedOverhead, accumulatedUncompressed, accumulatedCompressed;
} ;

#endif