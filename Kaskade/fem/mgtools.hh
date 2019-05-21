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

#ifndef MGTOOLS_HH
#define MGTOOLS_HH

#include <limits>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <iterator>
#include <boost/signal2.hpp>

#include "dune/grid/common/grid.hh"
#include "dune/grid/common/entitypointer.hh"
#include "dune/grid/common/entity.hh"
#include "dune/grid/io/file/dgfparser/dgfparser.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/iteratorfacades.hh"

#include "fem/fixdune.hh"
#include "fem/mllgeometry.hh"
#include "fem/fetransfer.hh"

/**
 * TODO: THIS CLASS IS NOT NEEDED ANYMORE!
 * Class providing a forward iterator for traversing the whole grid
 * hierarchy (only entities with codim=0, i.e. cells). Needed for
 * constructing FE function spaces with an HierarchicIndexSet.
 */

template< class GridImp, class ValueType >
class FullHierarchicIterator 
:  public Dune::ForwardIteratorFacade< FullHierarchicIterator<GridImp,ValueType>, ValueType  >

{
public:
  static const int codimension = 0 ;

  typedef typename GridImp::template Codim<codimension>::LevelIterator LevelIterator ;
  typedef typename GridImp::template Codim<codimension>::HierarchicIterator HierarchicIterator ;
  typedef ValueType Entity ;

  FullHierarchicIterator( GridImp const& grid_, int maxlevel_, bool end = false ) :
    grid( &grid_ ), coarseGridIt( grid_.template lbegin<0>(0) ),
    coarseGridEnd( grid_.template lend<0>(0) ), hierIt( coarseGridIt->hbegin(maxlevel_) ),
    hierEnd( coarseGridIt->hend(maxlevel_) ), maxlevel( maxlevel_ ), coarse(true)
  {
    if( end ) coarseGridIt = coarseGridEnd ;
  }

  FullHierarchicIterator( FullHierarchicIterator const& other ) :
    grid(other.grid), coarseGridIt(other.coarseGridIt), coarseGridEnd(other.coarseGridEnd), hierIt(other.hierIt),
    hierEnd(other.hierEnd), maxlevel(other.maxlevel), coarse(other.coarse)
  {}

  void increment()
  {
    if( coarse ) // last increment: ++coarseGridIt, now initialize new hierarchic iterator
    {
      if( maxlevel == 0 ) // special case, don't use hierIt
      {
        ++coarseGridIt ;
      }
      else
      {
        hierIt = coarseGridIt->hbegin(maxlevel) ;
        hierEnd = coarseGridIt->hend(maxlevel) ;
        coarse = false ;
      }
    }
    else
    {
      ++hierIt ;
      if( hierIt == hierEnd )
      {
        ++coarseGridIt ;
        coarse = true ;
      }
    }
  }

  ValueType& dereference() const
  {
    if( coarse || maxlevel==0 ) return *coarseGridIt ;
    return *hierIt ;
  }

  int level()
  {
    if( coarse || maxlevel==0 ) return 0 ;
    return hierIt.level() ;
  }

  bool equals( FullHierarchicIterator const& rhs ) const
  {
    if( rhs.maxlevel != maxlevel || rhs.grid != grid ) return false ;
    if( rhs.coarse == true && coarse == true ) return ( rhs.coarseGridIt == coarseGridIt ) ;
    return ( rhs.coarseGridIt == coarseGridIt && rhs.hierIt == hierIt ) ;
  }

  // cast to EntityPointer -- needed e.g. for TransferData
  operator typename GridImp::template Codim<0>::EntityPointer ()
        {
    if( coarse || maxlevel==0 ) return typename GridImp::template Codim<0>::EntityPointer( coarseGridIt ) ;
    return typename GridImp::template Codim<0>::EntityPointer( hierIt ) ;
        }

private:
  const GridImp* grid ;
  LevelIterator coarseGridIt, coarseGridEnd ;
  HierarchicIterator  hierIt, hierEnd ;
  int maxlevel ;
  bool coarse ;
};


/**
 * TODO: THIS CLASS IS NOT NEEDED ANYMORE!
 * Class providing an index set containing the whole grid hierarchy, 
 * i.e. every entity on every grid level has a unique index.
 * By default, the index set is connected to the grid signals, and
 * updates on grid refinement.
 * 
 */

template< class Grid > 
class HierarchicIndexSet
{
private:
  static const int dim = Grid::dimension;
  typedef typename Grid::template Codim<0>::LevelIterator ElementLevelIterator ;
  typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator ;
  typedef typename Grid::Traits::LevelIndexSet LevelIndexSet ;
  typedef typename Grid::template Codim<0>::Entity Entity;

public:
  template <int cd>
  struct Codim
  {
    template <Dune::PartitionIteratorType pitype>
    struct Partition
    {
      typedef FullHierarchicIterator<Grid, Entity> Iterator;
    };
  };

  typedef int IndexType;
  typedef HierarchicIndexSet<Grid> Self ;

  HierarchicIndexSet( GridManager<Grid>& gridMan_, bool connectMe_ = true ) :
    gridMan( gridMan_ ), grid(gridMan_.grid() ), maxlevel( gridMan_.maxLevel() ),connectMe( connectMe_ ), rtf(*this)
  {
    if( connectMe ) refConnection = gridMan.signals.informAboutRefinement.connect(0, rtf);
    init() ;
    update() ;
  }

  HierarchicIndexSet( GridManager<Grid>& gridMan_, int maxlevel_, bool connectMe_ = true ) :
    gridMan( gridMan_ ), grid(gridMan_.grid() ), maxlevel( maxlevel_ ), connectMe(connectMe_), rtf(*this)
  {
    if( connectMe ) refConnection = gridMan.signals.informAboutRefinement.connect(0, rtf);
    init();
    update() ;
  }

  HierarchicIndexSet(const HierarchicIndexSet& other ) :
    gridMan( other.gridMan ), grid(other.grid ), maxlevel( other.maxlevel ), connectMe(other.connectMe), rtf(*this)
  {
    if( connectMe ) refConnection = gridMan.signals.informAboutRefinement.connect(0, rtf);
    init() ;
    update() ;
  }

  HierarchicIndexSet& operator=(const HierarchicIndexSet& other)
  {
    gridMan = other.gridMan ;
    grid = other.grid ;
    maxlevel = other.maxlevel ;
    connectMe = other.connectMe ;
    if( connectMe ) refConnection = gridMan.signals.informAboutRefinement.connect(0, rtf);
    init();
    update() ;
  }

  ~HierarchicIndexSet()
  {
    refConnection.disconnect() ;
  }

public:

  void init()
  {
    std::vector< std::vector<int> > tmp( maxlevel+1, std::vector<int>(1, -1) ) ;

    if( indices.size() > 0 ) indices.clear() ;
    if( numEntities.size() > 0 ) numEntities.clear() ;

    indices.resize( 4, tmp ) ;
    numEntities.resize( 4, 0 ) ;
  }

  /// update hierarchic index set, e.g. after grid change
  void update( int maxlevel_ )
  {
    maxlevel = maxlevel_ ;
    for( int i = 0 ; i < 4 ; i++ ) indices[i].resize( maxlevel+1, std::vector<int>(1,-1) ) ;
    update() ;
  }

  void update( )
  {
    LevelIndexSet const& coarseLevelIndexSet = grid.levelIndexSet(0) ;

    for( ElementLevelIterator it = grid.template lbegin<0>( 0 );
        it != grid.template lend<0>( 0 ); ++it )
    {
      // add father to index set, including sub-entities of all available codimensions
      for( int codim = 0 ; codim <= dim ; codim++ )
      {
        int nSubentities = FixDuneDetail::GetNumberOfSubEntities<dim>::value( *it, codim ) ;
        for( int k = 0 ; k < nSubentities ; k++ )
        {
          int subindex = FixDuneDetail::GetIndexOfSubEntity<dim>::value(*it, codim, k, coarseLevelIndexSet ) ;
          if( subindex >= indices[codim][0].size() ) indices[codim][0].resize( subindex+1, -1 ) ;
          if( indices[codim][0][subindex] == -1 )
          {
            indices[codim][0][subindex] = numEntities[codim] ;
            numEntities[codim]++ ;
          }
        }
      }

      // traverse children with level <= maxlevel and add indices for all subentities
      HierarchicIterator end = it->hend(maxlevel);
      for( HierarchicIterator hi = it->hbegin(maxlevel); hi != end; ++hi)
      {
        int level = hi.level() ;
        LevelIndexSet const& levelIndexSet = grid.levelIndexSet( level ) ;

        for( int codim = 0 ; codim <= dim ; codim++ )
        {
          int nSubentities = FixDuneDetail::GetNumberOfSubEntities<dim>::value( *hi, codim ) ;
          for( int k = 0 ; k < nSubentities ; k++ )
          {
            int subindex = FixDuneDetail::GetIndexOfSubEntity<dim>::value( *hi, codim, k, levelIndexSet ) ;
            if( subindex >= indices[codim][level].size() ) indices[codim][level].resize( subindex+1, -1 ) ;
            if( indices[codim][level][subindex] == -1 )
            {
              indices[codim][level][subindex] = numEntities[codim] ;
              numEntities[codim]++ ;
            }
          }
        }
      }
    }

    // build vector of geometry types
    if( geometryTypes.size() > 0 ) geometryTypes.clear() ;
    geometryTypes.resize( numEntities.size(), std::vector<Dune::GeometryType>(0) ) ;
    std::vector<Dune::GeometryType> geotmp ;
    for( int codim = 0; codim <= dim ; codim++ )
    {
      for( int level = 0 ; level <= maxlevel ; level++ )
      {
        geotmp = grid.levelIndexSet(level).geomTypes(codim) ;
        geometryTypes[codim].insert( geometryTypes[codim].end(), geotmp.begin(), geotmp.end() ) ;
      }
      std::sort( geometryTypes[codim].begin(), geometryTypes[codim].end() ) ;
      std::vector<Dune::GeometryType>::iterator it
      = std::unique( geometryTypes[codim].begin(), geometryTypes[codim].end() ) ;
      geometryTypes[codim].resize( it-geometryTypes[codim].begin() ) ;
    }
  }


  template<int cc>
  int index( const typename Grid::Traits::template Codim<cc>::Entity &e ) const
  {
    int idx = grid.levelIndexSet(e.level()).index(e) ;
    int level = e.level() ;
    assert( level < indices[cc].size() ) ;
    assert( idx < indices[cc][level].size() );
    return indices[cc][level][idx] ;
  }

  template<class EntityType >
  int index( const EntityType &e ) const
  {
    int idx = grid.levelIndexSet(e.level()).index(e) ;
    int level = e.level() ;
    assert( level < indices[EntityType::codimension].size() ) ;
    assert( idx < indices[EntityType::codimension][level].size() );
    return indices[EntityType::codimension][level][idx] ;
  }

  template<int cc>
  int subIndex( const typename Grid::Traits::template Codim<0>::Entity &e, int i ) const
  {
    int sidx = grid.levelIndexSet(e.level()).template subIndex<cc>(e,i) ;
    int level = e.level() ;
    assert( level < indices[cc].size() ) ;
    assert( sidx < indices[cc][level].size() ) ;
    return indices[cc][level][sidx] ;
  }

  const std::vector<Dune::GeometryType>& geomTypes( int codim ) const
        {
    assert( codim < geometryTypes.size() ) ;
    return geometryTypes[codim] ;
        }

  int size( Dune::GeometryType type ) const
  {
    int s=0 ;
    for( int level = 0 ; level <= maxlevel ; level++ )
    {
      s += grid.levelIndexSet(level).size(type) ;
    }
    return s ;
  }

  int size( int codim ) const
  {
    assert( codim < numEntities.size() ) ;
    return numEntities[codim] ;
  }

  template<class EntityType>
  bool contains( const EntityType &e ) const
  {
    int level = e.level() ;
    if( level > maxlevel ) return false ;
    int lindex = grid.levelIndexSet(level).index(e) ;
    int codim = EntityType::codimension ;
    if( lindex > indices[codim][level].size() ) return false ;
    if( indices[codim][level][lindex] > -1 ) return true ;
    return false ;
  }

  // iterators
  template<int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return FullHierarchicIterator<Grid, Entity>( grid, maxlevel ) ;
  }

  template<int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return FullHierarchicIterator<Grid, Entity>( grid, maxlevel, true ) ;
  }

  // debug
  void out(std::ostream& o) const
  {
    for( int i = 0 ; i < indices.size() ; i++ )
    {
      o << "\nCODIM " << i << "\n" ;
      for( int j = 0 ; j < indices[i].size() ; j++ )
      {
        for( int k = 0 ; k < indices[i][j].size() ; k++ )
          o << indices[i][j][k] << "\t" ;
        o << "\n" ;
      }
      o << "\n" ;
    }
  }

private:

  struct ReactToRefinement
  {
    ReactToRefinement(Self &is_) : is(is_) {}

    void operator()(typename GridSignals::Status const ref)
    {
      if(ref==GridSignals::AfterRefinement) is.update( is.gridMan.grid().maxLevel() ) ;
    }

  private:
    Self& is ;
  } ;

  std::vector< std::vector< std::vector<int> > > indices ; //indices[codim][level][level_index] = hierarchic_index ;
  std::vector< int > numEntities ;
  GridManager<Grid>& gridMan ;
  Grid const& grid ;

  int maxlevel ;
  bool connectMe ;
  ReactToRefinement rtf;
  boost::signals::connection refConnection;

  std::vector< std::vector<Dune::GeometryType> > geometryTypes ;
};





/**
 *  Class that computes and stores prolongation matrices between grid levels.
 *
 *  Template parameter "Space" has to be a function space over a LevelGridView, e.g. 
 *  typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LevelGridView> > H1LevelSpace ;
 *
 *  Example usage (see io/lossystorage.hh):
\verbatim
      LevelGridView gv = gridManager.grid().levelView( coarseLevel ) ;
      Space space( gridManager, gv, order ) ;
      Spaces spaces( &space ) ;
      VariableSet varSet( spaces, "y" ) ;
      typename VariableSet::VariableSet x( varSet ) ;

      // [e.g. set values for x on the coarse level]

      MultilevelTransfer<H1LevelSpace,Grid>  mlTransfer( gridManager, order, coarseLevel ) ;

      // prolongation from coarse grid to levels 1, ..., maxLevel :

      for( int l = coarseLevel ; l < maxLevel ; l++ )  {
	gv = gridManager.grid().levelView(l+1) ;
	space.setGridView(gv) ;
	std::unique_ptr<Dune::BlockVector<StorageValueType> > newCoeff=mlTransfer.apply(l, *boost::fusion::at_c<0>(x.data)) ;
 *boost::fusion::at_c<0>(x.data) = *newCoeff ;
      }
 \endverbatim
 *
 *  TODO: connect to refinement signal for automatic re-calculation in case of 
 *  grid changes
 */
template<class Space, class Grid> 
class MultilevelTransfer
{
private:
  MultilevelTransfer( MultilevelTransfer const& ) { }

public:
  MultilevelTransfer( GridManager<Grid>& gridMan_, int order_ = 1, int startlevel_ = 0 )
  : gridMan( gridMan_ ), order(order_), startlevel(startlevel_)
  {
    MultilevelCoarseningPolicy multilevelCoarseningPolicy(startlevel) ;
    typename Grid::LevelGridView gv = gridMan.grid().levelView(startlevel);

    for( int level=startlevel ; level < gridMan.grid().maxLevel() ; level++ )
    {
      Space fatherSpace( gridMan, gridMan.grid().levelView(level+1), order ) ;

      TransferData<Space,MultilevelCoarseningPolicy> transferData( fatherSpace, multilevelCoarseningPolicy );
      matrices.push_back( transferData.transferMatrix().release() );
      multilevelCoarseningPolicy.update(level+1) ;
    }
  }

  ~MultilevelTransfer()
  {
    for( int i = 0 ; i < matrices.size() ; i++ ) delete matrices[i] ;
  }

  void out(std::ostream& o) const
  {
    for( int level=startlevel ; level < gridMan.grid().maxLevel() ; level++ )
    {
      o << "prolongation " << level << "->" << level+1 << "\n" ;
      matrices[level-startlevel]->out(o) ;
    }
  }

  /// Access transfer matrix
  typename TransferData<Space,MultilevelCoarseningPolicy>::TransferMatrix* operator[](int level) const
  {
    return matrices[level-startlevel] ;
  }

  /// Apply transfer matrix to coefficients
  template <class StorageValue>
  std::unique_ptr<Dune::BlockVector<StorageValue> > apply(int level, Dune::BlockVector<StorageValue> const& oldCoeff) const
  {
    std::unique_ptr<Dune::BlockVector<StorageValue> > newCoeff = matrices[level-startlevel]->apply(oldCoeff) ;
    return newCoeff ;
  }

private:
  GridManager<Grid>& gridMan ;
  int order ;
  int startlevel ;
  std::vector<typename TransferData<Space,MultilevelCoarseningPolicy>::TransferMatrix* > matrices ;
} ;




#endif
