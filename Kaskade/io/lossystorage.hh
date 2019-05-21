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


/**
*
* Use this file for lossy compression on uniformly refined, fixed grids.
* Temporal prediction/differential encoding is performed.
*
*/


#ifndef LOSSYSTORAGE_HH
#define LOSSYSTORAGE_HH


#include <limits>

#include "io/vtk.hh"
#include "io/rangecoder.hh"
#include "fem/mgtools.hh"

#include <boost/unordered_map.hpp>

typedef unsigned long ULong ;

bool abscompare( double a, double b ) { return fabs( a ) < fabs( b ) ; }
template<class T> bool greaterZero(T t) {return t > 0 ;}


// Use policies for determining the quantization intervals
template <class VariableSet, class Grid, int D=2>
struct UniformQuantizationPolicy 
{
   typedef Dune::FieldVector<double,D> Coor ;
   
   long int getIntervals(double range, double aTol, Coor const& coor = Coor(0), bool uniform = true ) 
   { 
//      return (long int) ceil( range / aTol / 2.0 ) ; 
      // symmetric mid-tread quantization : need odd number of intervals
      long int nIntervals = (long int) ceil( range / aTol / 2.0 ) ; 
      if( nIntervals % 2 == 0 ) return nIntervals + 1 ;
      return nIntervals ;
   }
} ;



template <class Grid, class VariableSet, class Space, class QuantizationPolicy=UniformQuantizationPolicy<VariableSet,Grid,Grid::dimension> >
class LossyStorage 
{   
  public:
    static const int dim = Grid::dimension ;

    // some typedefs used throughout the class
    typedef Dune::FieldVector<double,1> StorageValueType;
    typedef boost::fusion::vector<Space const*> Spaces ;
    typedef boost::fusion::vector<VariableDescription<0,1,0> > VariableDescriptions;
    typedef VariableSetDescription<Spaces,VariableDescriptions> PredVariableSet;
    typedef typename Grid::template Codim<dim>::LevelIterator VertexLevelIterator ;
    typedef typename Grid::template Codim<dim>::LeafIterator  VertexLeafIterator ;
    typedef typename Grid::template Codim<dim>::LeafIndexSet IndexSet ;
    typedef typename Grid::LevelGridView LevelView ;

  public:

    LossyStorage( GridManager<Grid>& gridManager_, VariableSet const& varSet_, int coarseLevel_, double aTol_ , bool uniform_ = false, 
		  QuantizationPolicy& quantizationPolicy_ = QuantizationPolicy()  ) 
      : gridManager(gridManager_), varSet(varSet_), coarseLevel(coarseLevel_), aTol(aTol_), 
        count(0), firstCall(true), uniform(uniform_), order(1), quantizationPolicy(quantizationPolicy_)
    {
      mlTransfer = new MultilevelTransfer<Space,Grid>( gridManager, order, coarseLevel ) ;
      ioOptions.outputType = IoOptions::ascii ;
      previousIndices = new typename VariableSet::VariableSet(varSet) ;
      
      // build map: vertex id - level the vertex is created
      unsigned long int ID ;
      for( int level = 0 ; level <= gridManager.grid().maxLevel() ; level++ )
      {
	for ( VertexLevelIterator it = gridManager.grid().template lbegin <dim>( level ); 
	      it != gridManager.grid().template lend <dim>( level ); ++ it)
	{
	  ID = gridManager.grid().globalIdSet().id( *it ) ;
	  if( levelInfo.find(ID) != levelInfo.end() ) continue ;
	  levelInfo[ID] = level ;
	}
      }
    }


    ~LossyStorage()
    {
      delete mlTransfer ;
      delete previousIndices ;
    }

    /**
     *  Encode a given state, e.g. the difference between to timesteps.
     *  Quantized values are written to a file specified by the member variables 
     *  filePrefix and count, keeping track of the already encoded timesteps.
     */
    typename VariableSet::VariableSet encode( typename VariableSet::VariableSet const& sol, std::string fn )
    {
      fn += ".dat" ;

      std::ostringstream debugfn ;

      IndexSet const& indexSet = varSet.indexSet ;
      typename VariableSet::VariableSet predictionError(varSet) ;
      predictionError *= 0 ; 

      // create variable set over coarse grid function space
      LevelView gv = gridManager.grid().levelView( coarseLevel ) ;
      Space space( gridManager, gv, order ) ;
      Spaces spaces( &space ) ;
      std::string varNames[1] = { "pred" };
      PredVariableSet predVarSet( spaces, varNames ) ;
      typename PredVariableSet::VariableSet x( predVarSet ) ;
      x *= 0 ; 

      // calculate quantization (for coarseLevel)
      std::vector<std::vector<ULong> > allIndices ;
      std::vector<ULong> indices ;
      std::vector<double> nodalValues ;
      
      VertexLevelIterator itEnd = gridManager.grid().template lend<dim>( coarseLevel );
      for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( coarseLevel ) ;
	    it != itEnd; ++it )
      {
	nodalValues.push_back( -(*boost::fusion::at_c<0>(sol.data))[indexSet.index(*it)] ) ;
      }    

      double ub = *std::max_element( nodalValues.begin(), nodalValues.end() ) ;
      double lb = *std::min_element( nodalValues.begin(), nodalValues.end() ) ;

      // symmetric quantization
      ub = *std::max_element( nodalValues.begin(), nodalValues.end(), abscompare ) ;
      if( ub < 0 ) { lb = ub ; ub = -lb  ; }
      else         { lb = -ub ; }
      // --
      
      double range = ub - lb ;
     
      double maxBits ;
      if( range > 0 ) maxBits = -ld(aTol) ;
      else            maxBits = 0 ;
      if( maxBits < 0 ) maxBits = 0 ; 
      maxIntervals = (long int)ceil(pow( 2, maxBits ));

      std::ofstream out( fn.c_str(), std::ios::binary ) ;  
      out.write(reinterpret_cast<char*>( &lb ), sizeof(double) ) ; 
      out.write(reinterpret_cast<char*>( &ub ), sizeof(double) ) ; 
      
      int vertexCount = 0 ;
      if( uniform )
      {
	quantize( nodalValues, indices, coarseLevel, lb, ub ) ;    
	// collect indices on each level in a common vector, and after finishing the
	// multilevel prediction, calculate occuring frequencies of interval numbers, and encode
	// using range coder
	allIndices.push_back( indices ) ;
	reconstruct( nodalValues, indices, coarseLevel, lb, ub ) ;
      }
      else
      {
	indices.clear() ; indices.resize(nodalValues.size(), 0);
	
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( coarseLevel ) ;
	  it != gridManager.grid().template lend<dim>( coarseLevel ); ++it )
	{ 
	  indices[vertexCount] = quantize( nodalValues[vertexCount], lb, ub, it->geometry().corner(0) );
	  nodalValues[vertexCount] = reconstruct( indices[vertexCount], lb, ub, it->geometry().corner(0) );
	  vertexCount++ ;
	}    
	allIndices.push_back( indices ) ;
      }
      
      
      // assign reconstructed coarse grid values of solution to hierarchic variable
      // to avoid error accumulation
      vertexCount = 0 ;
      for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( coarseLevel ); 
	      it != itEnd; ++it)
      {
	(*boost::fusion::at_c<0>(x.data))[gv.indexSet().index(*it)] -= nodalValues[vertexCount] ;
	vertexCount++ ;
      }

      int maxLevel = gridManager.grid().maxLevel() ;

      for( int l = coarseLevel ; l < maxLevel ; l++ )
      {
	itEnd = gridManager.grid().template lend<dim>( l+1 );
	
	gv = gridManager.grid().levelView(l+1) ;
	space.setGridView(gv);

	*boost::fusion::at_c<0>(x.data) = *(mlTransfer->apply(l, *boost::fusion::at_c<0>(x.data))) ;
  
	nodalValues.clear() ;
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l+1 ); 
		it != itEnd; ++it)
	{
	  double val = (*boost::fusion::at_c<0>(x.data))[gv.indexSet().index(*it)] 
		      -(*boost::fusion::at_c<0>(sol.data))[gridManager.grid().leafIndexSet().index(*it)] ;

	  unsigned ID = gridManager.grid().globalIdSet().id( *it ) ;
	  unsigned char vertexLevel = levelInfo[ID] ;
	  if( vertexLevel == l+1 ) nodalValues.push_back(val);
	}

	ub = *std::max_element( nodalValues.begin(), nodalValues.end() ) ;
	lb = *std::min_element( nodalValues.begin(), nodalValues.end() ) ;
	
	// for symmetric quantization (mid tread)
	ub = *std::max_element( nodalValues.begin(), nodalValues.end(), abscompare ) ;
	if( ub < 0 ) { lb = ub ; ub = -lb  ; }
	else         { lb = -ub ; }
	
	out.write(reinterpret_cast<char*>( &lb ), sizeof(double) ) ; 
	out.write(reinterpret_cast<char*>( &ub ), sizeof(double) ) ; 
	
	if( uniform ) 
	{
	  quantize( nodalValues, indices, l+1, lb, ub) ;
	  allIndices.push_back( indices ) ;
	  reconstruct( nodalValues, indices, l+1, lb, ub ) ;
	}
	else
	{
	  vertexCount = 0 ;
	  indices.clear() ; indices.resize(nodalValues.size(), 0);
	  	    
	  for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l+1 ) ;
	    it != itEnd; ++it )
	  {
	    
	    indices[vertexCount] = quantize( nodalValues[vertexCount], lb, ub, it->geometry().corner(0) );
	    nodalValues[vertexCount] = reconstruct( indices[vertexCount], lb, ub, it->geometry().corner(0) );
	    vertexCount++ ;
	  }    
	  allIndices.push_back( indices ) ;
	}
	
	// prepare prediction on next level -- use reconstructed values
	vertexCount = 0 ;
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l+1 ); 
		it != itEnd; ++it)
	{
	  // correction only for the nodes on level l+1
	  unsigned ID = gridManager.grid().globalIdSet().id( *it ) ;
	  unsigned char vertexLevel = levelInfo[ID] ;
	  if( vertexLevel < l+1 ) { continue ; }

	  int idx = gv.indexSet().index(*it) ;
	  (*boost::fusion::at_c<0>(x.data))[idx] -= nodalValues[vertexCount] ;
	  vertexCount++ ;
	}
      }
      
      // for debugging
      typename VariableSet::VariableSet pred(varSet) ;
      for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>(maxLevel); 
	      it != gridManager.grid().template lend<dim>(maxLevel); ++it)
      {
	(*boost::fusion::at_c<0>(pred.data))[ gridManager.grid().leafIndexSet().index(*it) ] = 
	  (*boost::fusion::at_c<0>(x.data))[ gv.indexSet().index(*it) ] ;
      }
      
      pred -= sol ;
//       std::vector<double> foo( gridManager.grid().size(dim) ) ;
//       pred.write( foo.begin() ) ;
//       double absQuantErr = fabs(*std::max_element(foo.begin(),foo.end(),abscompare)) ;
//       sol.write( foo.begin() ) ;
//       double relQuantErr = absQuantErr / fabs(*std::max_element(foo.begin(),foo.end(),abscompare)) ;
//       std::cout << "Linf quant err: abs = " << absQuantErr << "\trel = " << relQuantErr << "\n" ;
//       
//       L2Norm l2 ; 
//       double l2q2 = l2.square( boost::fusion::at_c<0>(pred.vars) ) ;
//       double l2y2 = l2.square( boost::fusion::at_c<0>(sol.vars) ) ;
//       relQuantErr = sqrt(l2q2) / sqrt(l2y2) ;
//       absQuantErr = sqrt(l2q2) ;
//       std::cout << "  L2 quant err: abs = " << absQuantErr << "\trel = " << relQuantErr << "\n" ;

      count++ ;
      
      // Now allIndices contains the collected interval numbers of all levels.
      // In order to maximize efficiency of the range coder, precalculate the
      // frequency of each interval, befor encoding.
      typename VariableSet::VariableSet intervals(varSet) ;
      for( int l = maxLevel ; l >= coarseLevel ; l-- )
      {
	itEnd = gridManager.grid().template lend<dim>( l );
	vertexCount = 0 ;
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l ); 
		it != itEnd; ++it)
	{
	  unsigned ID = gridManager.grid().globalIdSet().id(*it);
	  unsigned char vertexLevel = levelInfo[ID] ;
	  if( vertexLevel != l ) continue ;
	  
	  (*boost::fusion::at_c<0>(intervals.data))[indexSet.index(*it)] = allIndices[l-coarseLevel][vertexCount] ;
	  vertexCount++ ;
	}
      }

      out.close() ; // until now, only the interval bounds were written to "out"

      if( firstCall ) // first call of the method, must wait for next call to build difference
      { 
	prevFilename = fn ;
	*previousIndices *= 0 ; *previousIndices += intervals ;
	firstCall = false ;
	return pred ;
      }
      
      // write difference to previous file
      out.open( prevFilename.c_str(), std::ios::binary | std::ios::app ) ;
      
      double byteswritten = 0 ; // count overhead bytes

      // generate difference between previous indices and current indices
      typename VariableSet::VariableSet difference(varSet) ;
      difference *= 0 ;
      difference += intervals ;  difference -= *previousIndices ; 

      // prepare next write
      *previousIndices *= 0 ; *previousIndices += intervals ;  prevFilename = fn ;

      // and output the differences
      std::vector<long int> intervalIndicesTmp( gridManager.grid().size(dim) ) ;
      difference.write( intervalIndicesTmp.begin() ) ;

      long int minIndex = *std::min_element( intervalIndicesTmp.begin(), intervalIndicesTmp.end() ) ;
      std::vector<ULong> intervalIndices( gridManager.grid().size(dim) ) ;
      for( int i = 0 ; i < intervalIndicesTmp.size() ; i++ ) intervalIndices[i] = intervalIndicesTmp[i] - minIndex ;
      out.write(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;
      byteswritten += sizeof(long int) ;

      long int maxIndex = *std::max_element( intervalIndices.begin(), intervalIndices.end() ) + 1 ;
      
      std::vector<ULong> frequencies( maxIndex, 0 ) ;
      for( int i = 0 ; i < intervalIndices.size() ; i++ ) frequencies[ intervalIndices[i] ]++ ;

      unsigned int nnz = (unsigned int) std::count_if( frequencies.begin(), frequencies.end(), greaterZero<ULong> ) ;
      std::vector<int> dictionary( frequencies.size(), -1 ) ;
      int cur_no = 0 ;
      for( int i = 0 ; i < frequencies.size() ; i++ )
      {
	if( frequencies[i] > 0 ) 
	{ 
	  dictionary[i] = cur_no ; cur_no++ ; 
	}
      }
      out.write(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;
      byteswritten += sizeof(unsigned int) ;
      
      for( int i = 0 ; i < dictionary.size() ; i++ )
      {
	if( dictionary[i] >= 0 ) { byteswritten += sizeof(int) ; out.write(reinterpret_cast<char*>( &i ), sizeof(int) )  ; }
      }
      frequencies.erase( std::remove( frequencies.begin(), frequencies.end(), 0 ), frequencies.end() ) ;
      for( int i = 0 ; i < frequencies.size() ; i++ ) 
      {
	out.write(reinterpret_cast<char*>( &frequencies[i] ), sizeof(ULong) ) ; 
	byteswritten += sizeof(ULong) ;
      }
    
      byteswritten += maxLevel * sizeof(double) * 2 ;
      
//       std::cout << "Overhead written: " << byteswritten << " byte.\n" ;
    
      Alphabet<ULong> alphabet( frequencies.begin(), frequencies.end() ) ; 
      RangeEncoder<ULong> encoder( out ) ;      
      for( int i = 0 ; i < intervalIndices.size(); i++ ) 
      {
	encodeSymbol( encoder, alphabet, dictionary[intervalIndices[i]] );
      }
      
      // don't use precomputed frequencies for range encoding
//       int symbolCounter = 0 ;
//       std::vector<ULong> count(maxIntervals, 1) ;
//       Alphabet<ULong> alphabet( count.begin(), count.end() ) ;
//       RangeEncoder<ULong> encoder( out ) ;      
//       for( int i = 0 ; i < intervalIndices.size(); i++ ) 
//       {
// 	encodeSymbol( encoder, alphabet, intervalIndices[i] );
// 	++count[intervalIndices[i]];
// 	++symbolCounter ;
// 	if (symbolCounter>2*alphabet.size())
// 	{
// 	  alphabet.update(count.begin(),count.end());
// 	  symbolCounter=0;
// 	}
//       }

      return pred ; // return quantization error
    }

    /** In the end, write remaining indices to disc. */
    void flush()
    {
      std::ofstream out( prevFilename.c_str(), std::ios::binary | std::ios::app ) ;
      
      double byteswritten = 0 ;

      std::vector<long int> intervalIndicesTmp( gridManager.grid().size(dim) ) ;
      previousIndices->write( intervalIndicesTmp.begin() ) ;

      long int minIndex = *std::min_element( intervalIndicesTmp.begin(), intervalIndicesTmp.end() ) ;
      std::vector<ULong> intervalIndices( gridManager.grid().size(dim) ) ;
      for( int i = 0 ; i < intervalIndicesTmp.size() ; i++ ) intervalIndices[i] = intervalIndicesTmp[i] - minIndex ;
      out.write(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;

      long int maxIndex = *std::max_element( intervalIndices.begin(), intervalIndices.end() ) + 1 ;
      
      std::vector<ULong> frequencies( maxIndex, 0 ) ;
      for( int i = 0 ; i < intervalIndices.size() ; i++ ) frequencies[ intervalIndices[i] ]++ ;

      unsigned int nnz = (unsigned int) std::count_if( frequencies.begin(), frequencies.end(), greaterZero<ULong> ) ;
      std::vector<int> dictionary( frequencies.size(), -1 ) ;
      int cur_no = 0 ;
      for( int i = 0 ; i < frequencies.size() ; i++ )
      {
	if( frequencies[i] > 0 ) 
	{ 
	  dictionary[i] = cur_no ; cur_no++ ; 
	}
      }
      out.write(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;
      byteswritten += sizeof(unsigned int);     
      
      for( int i = 0 ; i < dictionary.size() ; i++ )
      {
	if( dictionary[i] >= 0 ) 
	{ 
	  out.write(reinterpret_cast<char*>( &i ), sizeof(int) ) ;
	  byteswritten += sizeof(int) ; 
	}
      }
      
      frequencies.erase( std::remove(frequencies.begin(), frequencies.end(), 0), frequencies.end() ) ;
      for( int i = 0 ; i < frequencies.size() ; i++ ) 
      {
	out.write(reinterpret_cast<char*>( &frequencies[i] ), sizeof(ULong) ) ; 
	byteswritten += sizeof(ULong);
      }
        
      byteswritten += gridManager.grid().maxLevel() * 2 * sizeof(double) ;
      
//       std::cout << "Overhead written: " << byteswritten << " byte.\n" ;
    
      Alphabet<ULong> alphabet( frequencies.begin(), frequencies.end() ) ; 
      RangeEncoder<ULong> encoder( out ) ;
      for( int i = 0 ; i < intervalIndices.size(); i++ ) 
      {
 	encodeSymbol( encoder, alphabet, dictionary[intervalIndices[i]] );
      }
      
      // don't use precomputed frequencies for range encoding
//       int symbolCounter = 0 ;
//       std::vector<ULong> count(maxIntervals, 1) ; 
//       Alphabet<ULong> alphabet( count.begin(), count.end() ) ;
//       
//       RangeEncoder<ULong> encoder( out ) ;
//       for( int i = 0 ; i < intervalIndices.size(); i++ ) 
//       {
// 	encodeSymbol( encoder, alphabet, intervalIndices[i] );
// 	++count[intervalIndices[i]];
// 	++symbolCounter ;
// 	if (symbolCounter>2*alphabet.size())
// 	{
// 	  alphabet.update(count.begin(),count.end());
// 	  symbolCounter=0;
// 	}
//       }
    }

    void finish() { firstCall = true ; } ;

    /**
     *  Decode quantized values and store in a VariableSet::VariableSet.
     *  The file to be read is specified by the member variables 
     *  filePrefix and count, keeping track of the already decoded timesteps.
     */
    void decode( GridManager<Grid>& gridManagerState, typename VariableSet::VariableSet& sol, std::string fn ) 
    {
      fn += ".dat" ;

      // build hierarchic variable set -- could probably be moved to constructor
      LevelView gv = gridManager.grid().levelView( coarseLevel ) ;
      Space space( gridManager, gv, order ) ;
      Spaces spaces( &space ) ;
      std::string varNames[1] = { "pred" };
      PredVariableSet predVarSet( spaces, varNames ) ;
      typename PredVariableSet::VariableSet x( predVarSet ) ;
      x *= 0 ; 
      
      // prepare prediction
      int maxLevel = gridManager.grid().maxLevel() ;
      IndexSet const& indexSet = varSet.indexSet ;

      std::vector<double> values ;
      
      std::vector<double> lb(maxLevel-coarseLevel+1), ub(maxLevel-coarseLevel+1) ;
      // read lb, ub from file
      std::ifstream in( fn.c_str(), std::ios::binary ) ;  
      for( int i = 0 ; i <= maxLevel-coarseLevel ; i++ )
      {
	in.read( reinterpret_cast<char*>( &lb[i] ), sizeof(double) ) ; 
	in.read( reinterpret_cast<char*>( &ub[i] ), sizeof(double) ) ; 
      }

      // read in indices from file
      typename VariableSet::VariableSet corr(varSet) ;
      typename VariableSet::VariableSet difference(varSet) ;

      std::vector<long int> intervalIndices( gridManager.grid().size(dim) ) ;

      long int minIndex ;
      in.read(reinterpret_cast<char*>( &minIndex ), sizeof(long int) ) ;   

      unsigned int nnz ;
      in.read(reinterpret_cast<char*>( &nnz ), sizeof(unsigned int) ) ;   // # non-empty intervals
      std::vector<int> dictionary( nnz ) ;
      for( int i = 0 ; i < nnz ; i++ )
      {
	in.read(reinterpret_cast<char*>( &dictionary[i] ), sizeof(int) ) ; // existing intervals (dictionary)
      }
      std::vector<ULong> frequencies( nnz, 0 ) ;
      for( int i = 0 ; i < nnz ; i++ )
	in.read(reinterpret_cast<char*>( &frequencies[i] ), sizeof(ULong) ) ; // frequencies

      Alphabet<ULong> alphabet( frequencies.begin(), frequencies.end() ) ; 
      
      try
      {
	RangeDecoder<ULong> decoder( in ) ;
	for( int i = 0 ; i < intervalIndices.size() ; i++ ) 
	{
	  ULong s = decodeSymbol(decoder,alphabet) ;
	  intervalIndices[i] = dictionary[s] + minIndex ;
	}
      }
      catch( std::ios_base::failure& e ) { std::cout << e.what() << " Decoding error\n" ; } 

      // don't use precomputed frequencies for range encoding
//       int symbolCounter = 0 ;
//       std::vector<ULong> symcount(maxIntervals, 1) ;
//       Alphabet<ULong> alphabet( symcount.begin(), symcount.end() ) ;
//       try
//       {
// 	RangeDecoder<ULong> decoder( in ) ;
// 	for( int i = 0 ; i < intervalIndices.size() ; i++ ) 
// 	{
// 	  ULong s = decodeSymbol(decoder,alphabet) ;
// 	  intervalIndices[i] = s + minIndex ;
// 	  ++symcount[s];
// 	  ++symbolCounter ;
// 	  if (symbolCounter>2*alphabet.size()) 
// 	  {
// 	    alphabet.update(symcount.begin(),symcount.end());
// 	    symbolCounter=0;
// 	  }
// 	}
//       }
//       catch( std::ios_base::failure& e ) { std::cout << e.what() << " Decoding error\n" ; } 


      difference.read( intervalIndices.begin() ) ;

      if( firstCall )
      {
	// on first call to decode, the first file contains the actual interval indices of the
	// prediction errors quantization
	*previousIndices *= 0 ;
	*previousIndices += difference ;
	corr *= 0 ;
	corr += difference ; 
	firstCall = false ;
      }
      else
      {
	// the file from which is read contains the difference between the previous interval
	// indices and the current indices. The interval indices of the current timestep are
	// previousIndices - difference.
	corr *= 0 ;
	corr += *previousIndices ; corr -= difference ;
	// prepare next call
	*previousIndices *= 0 ;  *previousIndices += corr ; 
      }
      
      
      // assign the overall vector to the single levels using levelInfo
      std::vector<std::vector<ULong> > allIndices( maxLevel+1-coarseLevel ) ;
      for( int l = coarseLevel ; l <= maxLevel ; l++ )
      {
	VertexLevelIterator itEnd = gridManager.grid().template lend<dim>( l );
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l ); 
	      it != itEnd; ++it)
	{
	  unsigned int ID = gridManager.grid().globalIdSet().id(*it);
	  unsigned char vertexLevel = levelInfo[ID];
	  if( vertexLevel == l )
	    allIndices[l-coarseLevel].push_back( (*boost::fusion::at_c<0>(corr.data))[indexSet.index(*it)] ) ;
	}
      }

      // start reconstruction
      int vertexCount = 0;
      VertexLevelIterator itEnd = gridManager.grid().template lend<dim>( coarseLevel );
      
      if( uniform )
      {
	reconstruct( values, allIndices[0], coarseLevel, lb[0], ub[0] );
      }
      else
      {
	values.clear()  ; values.resize( allIndices[0].size() ) ;
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( coarseLevel ) ;
	  it != itEnd; ++it )
	{
	  values[vertexCount] = reconstruct( allIndices[0][vertexCount], lb[0], ub[0], it->geometry().corner(0) );
	  vertexCount++ ;
	}    
      }
      
      vertexCount = 0 ;
      for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( coarseLevel ); 
	      it != itEnd; ++it)
      {
	(*boost::fusion::at_c<0>(x.data))[gv.indexSet().index(*it)] = -values[vertexCount] ; 
	vertexCount++ ;
      }
  
      // perform prediction and correction
      for( int l = coarseLevel ; l < maxLevel ; l++ )
      {
	itEnd = gridManager.grid().template lend<dim>( l+1 );
	
	gv = gridManager.grid().levelView(l+1) ;
	space.setGridView(gv);
	
	*boost::fusion::at_c<0>(x.data) = *(mlTransfer->apply(l, *boost::fusion::at_c<0>(x.data))) ;

	if( uniform )
	{
	  reconstruct( values, allIndices[l-coarseLevel+1], l+1, lb[l+1-coarseLevel], ub[l+1-coarseLevel] );
	}
	else
	{
	  vertexCount = 0;
	  values.clear()  ; values.resize( allIndices[l-coarseLevel+1].size() ) ;
	  for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l+1 ) ;
	    it != itEnd;  ++it )
	  {	
	    values[vertexCount] = reconstruct( allIndices[l-coarseLevel+1][vertexCount], lb[l-coarseLevel+1], 
					       ub[l-coarseLevel+1], it->geometry().corner(0));
	    vertexCount++ ;
	  }    
	}

	vertexCount = 0 ;
	for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>( l+1 ); 
		it != itEnd ; ++it)
	{
	  //correct only vertices on level l+1
	  if( levelInfo[gridManager.grid().globalIdSet().id(*it)] == l+1 ) 
	  {
	    (*boost::fusion::at_c<0>(x.data))[gv.indexSet().index(*it)] -= values[vertexCount] ;
	    vertexCount++ ;
	  }
	}  
      }
  
      itEnd = gridManager.grid().template lend<dim>(maxLevel);
      for ( VertexLevelIterator it = gridManager.grid().template lbegin<dim>(maxLevel); 
	      it != itEnd; ++it)
      {
	(*boost::fusion::at_c<0>(sol.data))[ gridManager.grid().leafIndexSet().index(*it) ] = 
	  (*boost::fusion::at_c<0>(x.data))[ gv.indexSet().index(*it) ] ;
      }
      
      std::ostringstream debugfn ;
      count-- ;
      in.close() ;
    }


  private:
    LossyStorage( LossyStorage const& );
    LossyStorage& operator=( LossyStorage const& ) ;

    /** Helper method to perform the actual quantization for a whole level. */
    void quantize( std::vector<double> const& values, std::vector<ULong>& indices, int level, double lb, double ub )
    {
      double range = ub - lb ;
     
      long int nIntervals = quantizationPolicy.getIntervals( range, aTol );
      
      if( nIntervals > maxIntervals ) maxIntervals = nIntervals ; // for rangecoder without freq.

//      double bits = 0 ;
//      if( nIntervals>0) bits = int(ld(nIntervals)*100)/100.0 ;
//      std::cout << "  Level:  " << level << "   range: " << range << "  Intervals: " << nIntervals 
//                 << "   Bits:  "  << ( (nIntervals > 0) ? (int(ld(nIntervals)*100))/100.0 : 0 ) 
// 		<< "   nodes: " << values.size() << "\n"
// 		<< "   bits*nodes: " << bits*values.size() << "\n" ;

      indices.clear() ;
      indices.resize( values.size(), 0 ) ;
      if( range > 0 ) // if range is empty, indices[i]=0 for all values
      {
	ULong idx ;
	for( int i = 0 ; i < values.size() ; i++ )
	{
	  idx = (ULong)( (values[i] - lb) * nIntervals / range ) + 1 ; 
	  indices[i] = idx ;
	}
      }
    }
    
    /** Helper method to perform the actual quantization for a single node. */
    ULong quantize( double value, double lb, double ub, Dune::FieldVector<double,dim> const& coor, bool uniform = false )
    {
      double range = ub - lb ;
      
      if( range == 0 ) return 0 ;
      
      long int nIntervals = quantizationPolicy.getIntervals( range, aTol, coor, uniform ) ;
      return (ULong)( (value - lb) * nIntervals / range ) + 1 ; 
    }
    

    /** Helper method to perform the actual reconstruction of quantized values. */
    void reconstruct( std::vector<double>& values, std::vector<ULong> const& indices, int level, double lb, double ub )
    {
      double range = ub - lb ; 

      long int nIntervals = quantizationPolicy.getIntervals( range, aTol );

      double delta ;
      if( nIntervals <= 0 ) delta = 0 ;
      else                  delta = range / nIntervals ;

      values.clear() ;
      values.resize( indices.size() ) ;
      for( int i = 0 ; i < indices.size() ; i++ )
      {
	if( indices[i] == 0 ) 		  values[i] = lb ;
	else	               	          values[i] = lb + (2*indices[i]-1) * delta/2 ;
      }
    }

    /** Helper method to perform the actual reconstruction for a single node, used for adaptive quantization. */
    double reconstruct( ULong index, double lb, double ub, Dune::FieldVector<double,dim> const& coor, bool uniform = false )
    {
      double range = ub - lb ;
      long int nIntervals = quantizationPolicy.getIntervals( range, aTol, coor, uniform ) ;
      double delta ;
      if( nIntervals <= 0 ) delta = 0 ;
      else                  delta = range / nIntervals ;
      if( index == 0 ) return lb ;
      else             return lb + (2*index-1) * delta/2 ;
    }


    /** Helper method returning the base 2-logarithm. */
    double ld( double val ) { return log(val)/log(2.0) ; }
  
  private:
    GridManager<Grid>& gridManager ;
    VariableSet varSet ;
    MultilevelTransfer<Space,Grid>* mlTransfer ;

    int coarseLevel ;
    double aTol ;
    std::string filePrefix ;

    int count ; // count timesteps
    IoOptions ioOptions;

    typename VariableSet::VariableSet* previousIndices ;
    bool firstCall ;
    std::string prevFilename ;
    long int maxIntervals ;
    bool uniform ;
    boost::unordered_map<unsigned long int, unsigned char> levelInfo ;
    
    int order ;
    
  public:
    QuantizationPolicy quantizationPolicy ;
} ;

#endif
