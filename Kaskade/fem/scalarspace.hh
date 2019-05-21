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

#ifndef SCALARSPACE_HH
#define SCALARSPACE_HH

#include <algorithm>
#include <functional>
#include <map>
#include <vector>

#include <boost/compressed_pair.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "fem/firstless.hh"
#include "fem/gridcombinatorics.hh"
#include "fem/views.hh"

namespace Kaskade
{
  // forward declarations
  template <class,class> class ContinuousLagrangeBoundaryMapper;

  namespace ScalarSpaceDetail
  {
    struct Empty {};

    /**
     * \brief A functor for extracting the first component of a boost compressed pair.
     */
    template <class Pair>
    struct CompressedFirst
    {
      typedef typename Pair::first_const_reference result_type;

      typename Pair::first_const_reference operator()(Pair const& pair) const
      {
        return pair.first();
      }
    };

    // Overload here as compressed pair's entries are accessed by method call. 
    struct FirstLess
    {
      template <class Data>
      bool operator()(boost::compressed_pair<size_t,Data> const& a, boost::compressed_pair<size_t,Data> const& b)
      {
        return a.first() < b.first();
      }
    };

    /// Default policy, in general use this.
    /**
     * \tparam Implementation ...MapperImplementation
     * \tparam optional additional data stored for each shapefunction (see hierarchical and boundary spaces for examples)
     */
    template <class Implementation,class SFData=Empty>
    class MapperPolicy
    {
      typedef typename Implementation::Grid              Grid;
      typedef typename Implementation::GridView          GridView;
      typedef typename GridView::IndexSet                IndexSet;
      typedef typename Grid::template Codim<0>::Entity   Cell;
      typedef typename Implementation::ShapeFunctionSet  ShapeFunctionSet;
      typedef boost::compressed_pair<size_t,SFData>      Data;

    protected:
      static constexpr bool boundaryPolicy = false;

      explicit MapperPolicy(Implementation const& impl) : implementation(impl)
      {}

      /**
       * \brief Computes the number of dofs associated to each geometry type in the grid and 
       *        initializes the startIndex indices accordingly, and resizes globIdx and sortedIdx appropriately.
       */
      template <class GlobalIndices, class SortedIndices>
      void initIndices(std::map<Dune::GeometryType,size_t>& startIndex, GlobalIndices& globIdx, SortedIndices& sortedIdx)
      {
        globIdx.resize(implementation.indexSet().size(0));
        sortedIdx.resize(implementation.indexSet().size(0));
        size_t n = 0;
        for (int codim=0; codim<=Grid::dimension; ++codim)
          for(auto const& geoType: implementation.indexSet().geomTypes(codim))
          {
            startIndex[geoType] = n;
            n += implementation.indexSet().size(geoType) * implementation.dofOnEntity(geoType);
          }
      }

      /**
       * \brief Store the global indices of the ansatz functions on cell into globIdx and a pair of global and local indices in sortedIdx.
       * 
       * Sorting of sortedIdx will be performed in UniformScalarMapper::update() and thus must not be done here.
       *
       * \param startIndex holds offsets for different geometry types, useful if global indices are ordered according to geometry types, see below
       */
      template <class GlobalIndices, class SortedIndices>
      void computeIndices(Cell const& cell, std::map<Dune::GeometryType,size_t>& startIndex, GlobalIndices& globIdx, SortedIndices& sortedIdx, size_t cellIndex) const
      {
        ShapeFunctionSet const& sf = implementation.shapeFunctions(cell);
        size_t localNumberOfShapeFunctions = sf.size();
        globIdx[cellIndex].resize(localNumberOfShapeFunctions);
        sortedIdx[cellIndex].resize(localNumberOfShapeFunctions);

        // step through all shape functions.
        // TODO: pull common parts out of this loop
        for(size_t i=0; i<localNumberOfShapeFunctions; ++i)
        {
          Dune::GeometryType gt;
          int subentity, codim, indexInSubentity;
          SFData sfData;
          implementation.entityIndex(cell,sf[i],i,gt,subentity,codim,indexInSubentity,sfData);

          // compute the global index of the subentity to which the shape function is associated
          int gt_idx = subIndex(implementation.indexSet(),cell,codim,subentity);

          auto mi = startIndex.find(gt);
          assert(mi!=startIndex.end());

          globIdx[cellIndex][i] =  Data(mi->second+gt_idx*implementation.dofOnEntity(gt)+indexInSubentity,sfData);
          sortedIdx[cellIndex][i] = std::make_pair(globIdx[cellIndex][i].first(),i);
        }
      }

      template <class Container, class ShapeFunctionSet>
      void initShapeFunctionSet(Container& sfs, ShapeFunctionSet const& sf, size_t cellIndex) const
      {
        sfs[cellIndex] = &sf; // cache shape function set
      }

      Implementation implementation;
    };

    /**
     * \cond internal
     */
    inline bool onFace_3D(int codim, int id, int localFaceId)
    {
      // edges
      if(codim==2)
      {
        if(localFaceId==0) return (id>=0 && id<3);
        if(localFaceId==1) return (id==0 || id==3 || id==4);
        if(localFaceId==2) return (id==1 || id==3 || id==5);
        if(localFaceId==3) return (id==2 || id==4 || id==5);
      }
      // vertices
      if(codim==3)
      {
        if(localFaceId==0) return (id>=0 && id<3);
        if(localFaceId==1) return (id==0 || id==1 || id==3);
        if(localFaceId==2) return (id==0 || id==2 || id==3);
        if(localFaceId==3) return (id==1 || id==2 || id==3);
      }
      return false;
    }
    /**
     * \endcond 
     */

    /// check if subentity is subentity of face
    /**
     * \param dim space dimension
     * \param codim codimension of subentity
     * \param id index within the subentities with codimension codim
     * \param localFaceId local index of face in cell
     */
   inline bool onFace(int dim, int codim, int id, int localFaceId)
    {
      if(codim==0) return false;
      if(codim==1) return id==localFaceId;

      if(dim==2)
      {
        if(localFaceId==0) return (id==0 || id==1);
        if(localFaceId==1) return (id==0 || id==2);
        if(localFaceId==2) return (id==2 || id==1);
      }
      if(dim==3) return onFace_3D(codim,id,localFaceId);

      return false;
    }

    static constexpr int defaultIndex = -94279;

    template <class Face>
    int getBoundaryId(Face const& face, std::vector<int> const& boundaryIds)
    {
      size_t boundarySegmentIndex = face.boundarySegmentIndex();
      if(boundarySegmentIndex < boundaryIds.size()) return boundaryIds[boundarySegmentIndex];
      return defaultIndex;
    }

    inline bool usedId(int id, std::vector<int> const& usedIds)
    {
      if (id==defaultIndex) return false;
      return std::find(usedIds.begin(), usedIds.end(), id) != usedIds.end();
    }

    template <class Face>
    inline bool considerFace(Face const& face, std::vector<int> const& boundaryIds, std::vector<int> const& usedIds)
    {
      if(boundaryIds.empty()) return face.boundary();
      return (face.boundary() && usedId(getBoundaryId(face,boundaryIds),usedIds));
    }

    ///
    struct RestrictToBoundary
    {
      //static constexpr int defaultIndex = -99999;

      RestrictToBoundary() : boundaryIds(), usedIds() {}
      RestrictToBoundary(RestrictToBoundary const& other) : boundaryIds(other.boundaryIds), usedIds(other.usedIds) {}
      RestrictToBoundary(std::vector<int> const& boundaryIds_, std::vector<int> const& usedIds_) : boundaryIds(boundaryIds_), usedIds(usedIds_)
      {}

      template <class Data, class GridView, class Cell>
      void treatBoundary(Data& data, GridView const& gridView, Cell const& cell, int codim, int subentity) const
      {
        // ignore shape functions that are associated with codim-0 entities
        if(codim > 0) {
          for(auto const& face : intersections(gridView,cell)) {
            // only consider boundary faces
            if(considerFace(face,boundaryIds,usedIds)) {
              if(onFace(GridView::dimension,codim,subentity,face.indexInInside())) {
                data.first() = true;
              }
            }
          }
        }
      }

      std::vector<int> const boundaryIds;
      std::vector<int> const usedIds;
    };

    template <class Policy>
    constexpr bool isRestriction()
    {
      return std::is_same<Policy,RestrictToBoundary>::value;
    }

    ///
    struct AllShapeFunctions
    {
      template <class Data, class GridView, class Cell> void treatBoundary(Data&, GridView const&, Cell const&, int, int) const {}
    };
  } // End of namespace ScalarSpaceDetail
  
  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------

  /**
   * \brief A traits class defining the type of argument that is provided by
   * UniformScalarMapper:combiner on the call to the Combiner constructor.
   */
  template <class SFData>
  struct UniformScalarMapperCombinerArgument
  {
    // nice try :DD
    // private:
    typedef typename std::vector<boost::compressed_pair<size_t,SFData> >::const_iterator Iterator;
    /**
     * \brief The sequence type that is provided to the call of the Combiner constructor. 
     * 
     * The type is a lightweight view type and should be kept by value, not by reference.
     */
    typedef RangeView<Iterator> type;
  };

  /**
   * \ingroup fem
   * \brief Base class for uniform scalar local to global mappers.
   *
   * It manages degrees of freedom for ansatz spaces where to each type of entity
   * of the grid the same number of global degrees of freedom is
   * associated and on each cell live the same number of shape
   * functions.
   *
   * We call a degree of freedom associated to an entity e of the grid
   * - if the support of the ansatz function is contained within the union
   *   of the cells (codim 0 entitites) which are incident to e and
   * - the ansatz function vanishes on all subentities of incident cells
   *   with a codimension not less than e.
   * This includes the usual polynomial FE spaces of arbitrary but fixed order.
   *
   * With each shape function on each cell, we associate its global
   * degree of freedom and additional data the type of is specified as
   * template parameter \a SFData.
   * 
   * \tparam Implementation
   * \tparam SFData type of (optional) additional data associated to global degrees of freedom
   *
   */
  template <class Implementation, class SFData = ScalarSpaceDetail::Empty>
  class UniformScalarMapper : public ScalarSpaceDetail::MapperPolicy<Implementation,SFData>
  {
  public:
    typedef typename Implementation::Grid              Grid;
    typedef typename Grid::template Codim<0>::Entity   Cell;
    typedef typename Implementation::ShapeFunctionSet  ShapeFunctionSet;
    typedef typename Implementation::Converter         Converter;
    typedef typename Implementation::Combiner          Combiner;
    typedef typename Implementation::Scalar            Scalar;
    typedef typename Implementation::GridView          GridView;
    typedef typename GridView::IndexSet                IndexSet;
    typedef std::pair<size_t,int>                      IndexPair;

  private:
    typedef boost::compressed_pair<size_t,SFData> Data;

    typedef ScalarSpaceDetail::CompressedFirst<Data> First;

    typedef boost::transform_iterator<First,typename std::vector<Data>::const_iterator> GlobalIndexIterator;
    typedef std::vector<IndexPair>::const_iterator                                      SortedIndexIterator;

    static constexpr int dim = Grid::dimension;

  public:
    typedef RangeView<GlobalIndexIterator> GlobalIndexRange;
    typedef RangeView<SortedIndexIterator> SortedIndexRange;

    /**
     * \brief Whether the ansatz functions have global support (i.e. lead to dense matrices).
     */
    static bool const globalSupport = false;


    UniformScalarMapper(Implementation const& impl) : ScalarSpaceDetail::MapperPolicy<Implementation,SFData>(impl)
    {
      update();
    }

    /**
     * \brief Returns the maximal polynomial order of shape functions encountered in any cell.
     */
    int maxOrder() const
    {
      return order;
    }
    
    
    /**
     * \brief Returns an empty range just for initialization purposes, since RangeView is not default constructible.
     */
    GlobalIndexRange initGlobalIndexRange() const
    {
      return GlobalIndexRange(GlobalIndexIterator(globIdx[0].begin(),First()),
                              GlobalIndexIterator(globIdx[0].begin(),First()));
    }


    /**
     * \brief Returns an immutable sequence containing the global indices of the shape functions associated to this cell. 
     * 
     * Global indices start at 0 and are consecutive - in the range returned here, an unordered subset
     * is contained.
     */
    GlobalIndexRange globalIndices(Cell const& cell) const
    {
      return globalIndices(implementation.indexSet().index(cell));
    }

    /**
     * \brief Returns an immutable sequence containing the global indices of the shape functions associated to this cell. 
     * 
     * Global indices start at 0 and are consecutive - in the range returned here, an unordered subset
     * is contained.
     */
    GlobalIndexRange globalIndices(size_t n) const
    {
      return GlobalIndexRange(GlobalIndexIterator(globIdx[n].begin(),First()),
                              GlobalIndexIterator(globIdx[n].end(),  First()));
    }

    /**
     * \brief Returns an empty range just for initialization, since RangeView is not default constructible.
     */
    static SortedIndexRange initSortedIndexRange()
    {
      static std::vector<IndexPair> dummy; // empty
      return SortedIndexRange(dummy.begin(),dummy.end());
    }

    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(Cell const& cell) const
    {
      return sortedIndices(implementation.indexSet().index(cell));
    }

    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(size_t n) const
    {
      return SortedIndexRange(sortedIdx[n].begin(),sortedIdx[n].end());
    }

    /**
     * \brief Returns the number of global degrees of freedom managed. 
     * 
     * Note that this does not correspond directly to the number of
     * coefficients in a FE function (if the FE function has more than
     * one component).
     */
    size_t size() const { return n; }

    /**
     * \brief DEPRECATED. Use maxOrder instead. This method will be removed after 2016-12-31.
     */
    int getOrder() const { return order; }

    /**
     * \brief Returns the set of shape functions defined on this cell.
     * \param cell the codim 0 entity of the grid for wich the shape functions are to be retrieved
     * \param contained if true, the method may assume that the cell is contained in the index set of the space.
     *                  (The other case occurs during interpolation between different grids).
     */
    ShapeFunctionSet const& shapefunctions(Cell const& cell, bool contained=false) const
    {
      if (contained || implementation.indexSet().contains(cell))
        return shapefunctions(implementation.indexSet().index(cell));
      else
        return implementation.shapeFunctions(cell);
      // This is the previous implementation (as of Kaskade 7.1). I'm not sure the
      // condition of the if makes sense at all. First, cell has to be contained in
      // the index set, otherwise the index and hence globalIndices(cell) is undefined.
      // However, if the cell is contained in the index set, the globalIndices are not
      // empty (compare the construction in update()). Ws-2012-06-07.
      //
      //     if ( cell.isLeaf() && globalIndices(cell).empty() )
      //       return implementation.emptyShapeFunctionSet();
      //     else
      //       return implementation.shapeFunctions(cell);
    }

    ShapeFunctionSet& shapefunctions_non_const(Cell const& cell)
    {
      return shapefunctions_non_const(implementation.indexSet().index(cell));
    }

    /**
     * \brief Returns the set of shape functions defined on this cell.
     */
    ShapeFunctionSet const& shapefunctions(size_t n) const
    {
      return *sfs[n];    
    }

    ShapeFunctionSet& shapefunctions_non_const(size_t n)
    {
      return *sfs[n];
    }

    ShapeFunctionSet const& lowerShapeFunctions(Cell const& cell) const
    {
      if (globalIndices(cell).empty())
        return implementation.emptyShapeFunctionSet();
      else
        return implementation.lowerShapeFunctions(cell);
    }

    /**
     * \brief Returns a combiner for the given leaf cell.
     * \param cell the grid cell for which the combiner is requested
     * \param index the index of the cell
     */
    Combiner combiner(Cell const& cell, size_t index) const {
      assert(implementation.indexSet().index(cell)==index);
      return Combiner(rangeView(globIdx[index].begin(),globIdx[index].end()));
    }

    /**
     * \brief (Re)computes the internal data.
     *
     * This has to be called after grid modifications and on construction.
     */
    void update()
    {
      GridView const& gridView = implementation.gridView();
      IndexSet const& indexSet = implementation.indexSet();

      // For each codimension (i.e. type of subentity) compute the
      // number of global ansatz functions as well as an accumulated
      // index into an array of all ansatz functions.
      startIndex.clear();
      // compute global degrees of freedom and store the offsets for
      // different geometry types in startIndex
      // see MapperPolicy
      // Precompute and cache all the global indices. First allocate the
      // memory needed to prevent frequent reallocation.
      this->initIndices(startIndex,globIdx,sortedIdx);

      sfs.resize(indexSet.size(0));
      order = 0;

      // Step through all cells and compute the global ansatz function
      // indices of all shape functions on that cell.
      for(auto const& element : elements(gridView))
      {
        size_t const cellIndex = indexSet.index(element);
        ShapeFunctionSet const& sf = implementation.shapeFunctions(element);
        this->initShapeFunctionSet(sfs,sf,cellIndex);                           // initShapeFunctionSet is defined in template base class MapperPolicy.
        order = std::max(order,sf.order());
        this->computeIndices(element, startIndex, globIdx, sortedIdx, cellIndex);   // computeIndices is defined in template base class MapperPolicy.


        // sort the index pairs according to the global index
        std::sort(sortedIdx[cellIndex].begin(),sortedIdx[cellIndex].end(),FirstLess());

        // make sure that any assigned shape function forms at most one ansatz function.
        // This is a functionality restriction, but seems quite reasonable as a debugging check.
        // Note that not all shape functions need not be assigned. If a shape function is mapped     //Jakob: Are negative indices for unused shape function really allowed? Are'nt they simply not present in globIdx, such that the check idx[j]<0 below is nonsense?
        // to negative global index, it takes not part at all (e.g. in hp-methods).
#ifndef NDEBUG
        std::vector<int> idx(globIdx[cellIndex].size());
        for (int j=0; j<idx.size(); ++j) idx[j] = globIdx[cellIndex][j].first();
        std::sort(idx.begin(),idx.end());
        for (int j=1; j<idx.size(); ++j) assert(idx[j]>idx[j-1] || idx[j]<0);
#endif
      }

      // compute overall number of degrees of freedom
      n = 0;
      for(size_t i=0; i<globIdx.size(); ++i)
        if(globIdx[i].size() > 0)
          n = std::max(n,std::max_element(globIdx[i].begin(),globIdx[i].end(),ScalarSpaceDetail::FirstLess())->first());
      ++n;
    }

    /**
     * \brief Returns a half-open range of global indices that are associated
     * to the entities with given geometry type.
     *
     * This is useful to partition stiffness matrices when using higher
     * order hierarcical ansatz functions in order to use a direct
     * solver on the low order ansatz functions as a preconditioner.
     */
    std::pair<size_t,size_t> globalIndexRange(Dune::GeometryType gt) const
    {
      assert(startIndex.find(gt)!=startIndex.end());

      size_t first = startIndex.find(gt)->second;
      size_t last = size();

      for (auto const& si: startIndex)
        if (si.second>first && si.second<last)
          last = si.second;

      return std::make_pair(first,last);
    }

  private:
    // In startIndex, for any geometry type that occurs in the indexSet
    // there is an index stored, such that the dofs associated with
    // nodes on subentities of this type start at this index. We request
    // that there is a fixed number of dofs associated with each
    // subentity type, such that equal length index blocks are assigned
    // to each subentity of a given type. The layout of dofs is thus as
    // follows (example continuous 2D simplex mesh with two triangles,
    // order 4):
    //
    // [Interior Triangle Nodes, Interior Edge Nodes, Vertex Nodes] with substructuring
    // [ [[t1a,t1b,t1c],[t2a,t2b,t2c]], [[e1a,e1b],[e2a,e2b],[e3a,e3b],[e4a,e4b],[e5a,e5b]], [v1,v2,v3,v4] ]
    //     | start index triangles        | start index edges                                 | start index vertices
    // localDofs(triangles)=3, localDofs(edges)=2, localDofs(vertices)=1
    size_t                              n;
    std::map<Dune::GeometryType,size_t> startIndex;

    // In globIdx, for each cell there are the global ansatz function indices of the
    // shape functions on that cell.
    std::vector<std::vector<Data>>      globIdx;

    // In sortedIdx, for each cell there is a vector of (global index, local index) pairs, sorted ascendingly
    // by the global index. This could be computed outside on demand, but the sorting takes time and may be
    // amortized over multiple assembly passes/matrix blocks if cached here.
    std::vector<std::vector<IndexPair> > sortedIdx;

    // In sfs, for each cell there is a pointer to the shape function set on this cell cached.
    std::vector<typename std::conditional<ScalarSpaceDetail::MapperPolicy<Implementation,SFData>::boundaryPolicy,ShapeFunctionSet,ShapeFunctionSet const>::type*> sfs;
    //std::vector<ShapeFunctionSet const*> sfs;

    int                                 order;

  protected:
    using ScalarSpaceDetail::MapperPolicy<Implementation,SFData>::implementation; // direct access to base class member variable
  };
} // end of namespace Kaskade

#endif
