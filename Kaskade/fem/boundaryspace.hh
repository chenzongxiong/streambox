/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BOUNDARY_SPACE_HH
#define BOUNDARY_SPACE_HH

#include <map>
#include <utility>

#include <dune/common/fvector.hh>

#include "fem/scalarspace.hh"

/**
 * @file
 * @brief  Trace mapper
 * @author Lars Lubkoll
 *
 * This file contains a boundary mapper that restricts existing mappers to (possibly parts of) the domain's boundary.
 *
 */


namespace Kaskade
{
  // forward declarations
  /// \internal
  template <class,class> class ContinuousLagrangeMapper;
  template <class,class,class> class ContinuousLagrangeMapperImplementation;
  template <class,class> class DiscontinuousLagrangeMapper;
  template <class,class> class ContinuousHierarchicMapper;
  template <class,class> class DiscontinuousHierarchicMapper;
  template <class,class,class> class ContinuousHierarchicMapperImplementation;
  template <class,class> class ContinuousHierarchicExtensionMapper;
  template <class,class,class> class ContinuousHierarchicExtensionMapperImplementation;

  template <class,int,class> class LagrangeSimplexShapeFunctionSet;
  template <class,int,class> class HierarchicSimplexShapeFunctionSet;
  template <class,int,class> class HierarchicExtensionSimplexShapeFunctionSet;
  /// \endinternal 

  namespace BoundarySpace_Detail
  {
    using namespace ScalarSpaceDetail;

    /// Choose implementation for mapper. Not the best design, but works...
    template <class MapperImplementation> struct ChooseMapper;

    template <class Scalar, class GridView>
    struct ChooseMapper<ContinuousLagrangeMapper<Scalar,GridView> >
    {
      typedef ContinuousLagrangeMapperImplementation<Scalar,GridView,RestrictToBoundary> type;
    };

    struct BOUNDARY_MAPPERS_ARE_NOT_YET_IMPLEMENTED_FOR_DISCONTINUOUS_SPACES;

    template <class Scalar, class GridView>
    struct ChooseMapper<DiscontinuousLagrangeMapper<Scalar,GridView> >
    {
      typedef BOUNDARY_MAPPERS_ARE_NOT_YET_IMPLEMENTED_FOR_DISCONTINUOUS_SPACES type;
    };

    template <class Scalar, class GridView>
    struct ChooseMapper<ContinuousHierarchicMapper<Scalar,GridView> >
    {
      typedef ContinuousHierarchicMapperImplementation<Scalar,GridView,RestrictToBoundary> type;
    };

    template <class Scalar, class GridView>
    struct ChooseMapper<ContinuousHierarchicExtensionMapper<Scalar,GridView> >
    {
      typedef ContinuousHierarchicExtensionMapperImplementation<Scalar,GridView,RestrictToBoundary> type;
    };

    template <class Scalar, class GridView>
    struct ChooseMapper<DiscontinuousHierarchicMapper<Scalar,GridView> >
    {
      typedef BOUNDARY_MAPPERS_ARE_NOT_YET_IMPLEMENTED_FOR_DISCONTINUOUS_SPACES type;
    };

    template <class Mapper> struct GetSFDataType;

    template <class Scalar, class GridView>
    struct GetSFDataType<ContinuousLagrangeMapper<Scalar,GridView> >
    {
      typedef boost::compressed_pair<bool,Empty> type;
    };

    template <class Scalar, class GridView>
    struct GetSFDataType<ContinuousHierarchicMapper<Scalar,GridView> >
    {
      typedef boost::compressed_pair<bool,int> type;
    };

    template <class Scalar, class GridView>
    struct GetSFDataType<ContinuousHierarchicExtensionMapper<Scalar,GridView> >
    {
      typedef boost::compressed_pair<bool,int> type;
    };

    /// get restricted shape function set for simplicial grid
    template <class Mapper> struct ShapeFunctionSetRestriction;

    template <class Scalar, class GridView, class Data>
    struct ShapeFunctionSetRestriction<ContinuousLagrangeMapperImplementation<Scalar,GridView,Data> >
    {
      typedef RestrictedShapeFunctionSet<LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> > type;
    };

    template <class Scalar, class GridView, class Data>
    struct ShapeFunctionSetRestriction<ContinuousHierarchicMapperImplementation<Scalar,GridView,Data> >
    {
      typedef RestrictedShapeFunctionSet<HierarchicSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> > type;
    };

    template <class Scalar, class GridView, class Data>
    struct ShapeFunctionSetRestriction<ContinuousHierarchicExtensionMapperImplementation<Scalar,GridView,Data> >
    {
      typedef RestrictedShapeFunctionSet<HierarchicExtensionSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> > type;
    };

  }

  /**
   * \ingroup fem
   * \brief A local to global mapper implementation for boundary spaces, with functions defined on the domain boundary (or parts of it) only.
   *
   * This works by restriction of the shape and ansatz functions. From the regular shapefunctions those are removed which are zero on boundary faces.
   * (For cells which have no boundary intersection all shape functions are removed.) Therefore, those boundary functions can be also seen as FE functions over the
   * entire domain whose support is restricted to a neighbourhood of the boundary.
   *
   * Unfortunately even continuous boundary functions will be discontinuous if considered as functions over the entire domain. See the following
   *
   *                              f   2    f
   *                          1    f     f    3
   *                                f  f
   *                    bbbbbbbbbbbb * bbbbbbbbbbbbbbbbbb
   *
   * Since cell 2 touches the boundary only at vertex * all shapefunctions of cell 2 will be removed. Therefore, the ansatz function belonging to node * will be discontinuous
   * (its support will be cell 1 and cell 3).
   *
   * This means that at the moment boundary functions should not be passed as second argument to interpolateGlobally, since there the function
   * values on vertices are averaged leading to results which are probably unexpected (depending on what you expect).
   * TODO: Fix this. (If you have done you can also slightly simplify writeVTK() and the assignment operator for FunctionSpaceElement.)
   *
   * \tparam Implementation_: Underlying FE mapper (e.g. ContinuousLagrangeMapper or ContinuousHierarchicMapper)
   * \tparam Scalar_: Type of scalars.
   * \tparam GridView_: Type of the corresponding grid view.
   */
  template <template <class,class> class Implementation_, class Scalar_, class GridView_>
  class BoundaryMapper : public UniformScalarMapper<typename BoundarySpace_Detail::ChooseMapper<Implementation_<Scalar_,GridView_>>::type, 
                                                    typename BoundarySpace_Detail::GetSFDataType<Implementation_<Scalar_,GridView_>>::type>
  {
    typedef typename BoundarySpace_Detail::ChooseMapper<Implementation_<Scalar_,GridView_>>::type Implementation;
    typedef typename BoundarySpace_Detail::GetSFDataType<Implementation_<Scalar_,GridView_>>::type DataType;
    typedef UniformScalarMapper<Implementation,DataType> Base;

  public:
    typedef Scalar_ Scalar;
    typedef GridView_ GridView;
    typedef typename BoundarySpace_Detail::ShapeFunctionSetRestriction<Implementation>::type ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = 0;
    static constexpr int dim = GridView::dimension;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     *
     * \tparam m number of components
     */
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<BoundaryMapper>,m> type;
    };

    /**
     * \ingroup fem
     * @brief BoundaryMapper: Create a boundary mapper over the entire boundary.
     * @param gridView: The corresponding grid view.
     * @param order: Order of the underlying FE mapper.
     */
    BoundaryMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {
      assert(order >= 1);
    }

    /**
     * \ingroup fem
     * @brief BoundaryMapper: Create a boundary mapper over parts of the boundary.
     * @param gridView: The corresponding grid view.
     * @param order: Order of the underlying FE mapper.
     * @param boundaryIds_: Container which maps each boundary segment index, which occurs in the grid (consecutive from zero to #(boundary segments in level 0 grid view)-1) to an int property value.
     * @param usedIds_: States by int property value which boundary faces belong to the desired part of the boundary.
     */
    BoundaryMapper(GridView const& gridView, int order, std::vector<int> const& boundaryIds_, std::vector<int> const& usedIds_)
     : Base(Implementation(gridView,order,ScalarSpaceDetail::RestrictToBoundary(boundaryIds_,usedIds_)))
    {
      assert(order>=1);
    }

    template <class T, class Functor>
    BoundaryMapper(GridView const& gridView, int order, std::vector<T> const& boundaryIds_, std::vector<T> const& usedIds_, Functor boundaryIdsToInts)
     : Base(Implementation(gridView,order,ScalarSpaceDetail::RestrictToBoundary(boundaryIdsToInts(boundaryIds_),boundaryIdsToInts(usedIds_))))
    {
      assert(order>=1);
    }

    bool inDomain(typename GridView::Intersection const& face) const
    {
      return ScalarSpaceDetail::considerFace(face,this->implementation.boundaryIds,this->implementation.usedIds);
    }

    /**
     * \ingroup fem
     * @brief numFaces: Counts number of faces on which this boundaryMapper is defined. This costs some time so use with care.
     * @return number of faces
     */
    int numFaces() {
      int nFaces = 0;
      forEachBoundaryFace(this->implementation.gridView(),[&](typename GridView::Intersection const& intersection) {
        if(inDomain(intersection)) ++nFaces;
      });
      return nFaces;
    }
  };


  /**
   * \cond internal
   */
  namespace ScalarSpaceDetail
  {
    /// Own implementation for boundary spaces
    template <template <class,class,class> class MapperImplementation, class Scalar_, class GridView_, class SFData>
    class MapperPolicy<MapperImplementation<Scalar_,GridView_,RestrictToBoundary>, boost::compressed_pair<bool,SFData> >
    {
      typedef Scalar_ Scalar;
      typedef GridView_ GridView;
      typedef MapperImplementation<Scalar,GridView,RestrictToBoundary>                     Implementation;
      typedef typename Implementation::Grid                                                  Grid;
      typedef typename GridView::IndexSet                                                    IndexSet;
      typedef typename Grid::template Codim<0>::Entity                                      Cell;
      typedef typename Implementation::ShapeFunctionSet                                      ShapeFunctionSet;
      typedef boost::compressed_pair<bool,SFData>                                            Data;
      typedef typename BoundarySpace_Detail::ShapeFunctionSetRestriction<Implementation>::type  ShapeFunctionSet_Restricted;

    protected:
      static constexpr bool boundaryPolicy = true;

      explicit MapperPolicy(Implementation const& impl) : implementation(impl),
      insertionIndex(0), consideredCells(implementation.indexSet().size(0))
      {}

      template <class GlobalIndices, class SortedIndices>
      void initIndices(std::map<Dune::GeometryType,size_t>& startIndex, GlobalIndices& globIdx, SortedIndices& sortedIdx)
      {
        insertionIndex = 0;
        globIdx.clear();
        globIdx.resize(implementation.indexSet().size(0));
        sortedIdx.clear();
        sortedIdx.resize(implementation.indexSet().size(0));

        consideredCells.clear();
        consideredCells.resize(implementation.indexSet().size(0),false);
        shapeFunctionSets.clear();
        shapeFunctionSets.resize(implementation.indexSet().size(0));

        size_t n = 0;
        for (int codim=0; codim<=Grid::dimension; ++codim)
          for(auto const& geoType: implementation.indexSet().geomTypes(codim))
          {
            startIndex[geoType] = n;
            n += implementation.indexSet().size(geoType) * implementation.dofOnEntity(geoType);
          }

        indexMapper.clear();
        indexMapper.resize(n,std::make_pair(false,0));

        restrictCells();
      }

      template <class GlobalIndices, class SortedIndices>
      void computeIndices(Cell const& cell, std::map<Dune::GeometryType,size_t>& startIndex, GlobalIndices& globIdx, SortedIndices& sortedIdx, size_t cellIndex)
      {
        // No evaluation will be done on faces that do not intersect the boundary.
        // If evaluation should be allowed make sure to set a, possibly empty, vector using setRestriction()
        std::vector<int> idRestriction;
        ShapeFunctionSet& sf = static_cast<UniformScalarMapper<Implementation,Data>&>(*this).shapefunctions_non_const(cell);
        if(!consideredCells[cellIndex])
        {
          if(cell.type().isSimplex()) static_cast<ShapeFunctionSet_Restricted&>(sf).setRestriction(idRestriction);
          return;
        }

        int nominalOrder; // unused

        // step through all shape functions.
        Dune::GeometryType gt;
        int subentity, codim, indexInSubentity;
        size_t realLocalId = 0;
        for(size_t i=0; i<sf.size(); ++i)
        {
          std::tie(nominalOrder,codim,subentity,indexInSubentity) = sf[i].location();

          if(codim == 0) continue;

          Data data(false);
          assert(cell.type().isSimplex());
          implementation.entityIndex(cell,sf[i],i,gt,subentity,codim,indexInSubentity,data);
          if(data.first())
          {
            std::pair<bool,size_t> contains = std::make_pair(false,0);
            if(codim > 0)
            {
              // compute the global index of the subentity to which the shape function is associated
              int gt_idx = subIndex(implementation.indexSet(),cell,codim,subentity);

              std::map<Dune::GeometryType,size_t>::const_iterator mi = startIndex.find(gt);
              assert(mi!=startIndex.end());

              // On the boundary codim-0-entities do vanish
              // Shape functions associated with codim-1-entities are unique
              // thus we only have to check shape functions associated with entities
              // of codim > 1
              size_t unrestrictedIndex = mi->second+gt_idx*implementation.dofOnEntity(gt)+indexInSubentity;
              if(codim > 1) {
                contains = indexMapper[unrestrictedIndex];
                if(!contains.first) indexMapper[unrestrictedIndex] = std::make_pair(true,insertionIndex);
              } else {
                // for shapefunctions on codim 1 entities also store the restricted index (originally this map was only provided and used internally for codim >1 shapefunctions),
                indexMapper[unrestrictedIndex] = std::make_pair(true,insertionIndex);
              }
            }

            if(!contains.first) globIdx[cellIndex].push_back( boost::compressed_pair<size_t,Data>(insertionIndex,data) );
            else globIdx[cellIndex].push_back( boost::compressed_pair<size_t,Data>(contains.second,data) );

            sortedIdx[cellIndex].push_back( std::make_pair(globIdx[cellIndex].back().first(),realLocalId) );

            idRestriction.push_back(i);

            if(!contains.first) ++insertionIndex;
            ++realLocalId;
          }
        }

        if(cell.type().isSimplex()) static_cast<ShapeFunctionSet_Restricted&>(sf).setRestriction(idRestriction);
      }

      template <class Container, class ShapeFunctionSet>
      void initShapeFunctionSet(Container& sfs, ShapeFunctionSet const& sf, size_t cellIndex)
      {
        // TODO: probably it is not the most efficient idea to create for each cell an own restricted shape function set
        // the next line was replaced since it caused a memory leak
        //sfs[cellIndex] = new ShapeFunctionSet_Restricted(static_cast<ShapeFunctionSet_Restricted const&>(sf),nullptr); // cache shape function set
        assert(cellIndex < shapeFunctionSets.size());
        shapeFunctionSets[cellIndex] = std::make_unique<ShapeFunctionSet_Restricted>(static_cast<ShapeFunctionSet_Restricted const&>(sf),nullptr);
        sfs[cellIndex] = shapeFunctionSets[cellIndex].get();
      }

      Implementation implementation;

    public:
      /**
       * @brief unrestrictedToRestrictedIndex: This map provides an easy way to identify the functions from the restricted ansatzfunctionset.
       *
       * It relies on the assumption that here the unrestricted index is computed in the same way as in scalarspace.hh.
       *
       * @param unrestrictedIndex: Unique index of ansatz function in the FE-space over the whole domain.
       * @return pair of values: The first states whether the (restriction of the) corresponding ansatz function is also contained in the FE-space over the boundary,
       *                         the second gives the index of ansatz function in the basis of the boundary FE-space, if it is contained.
       */
      std::pair<bool,size_t> const& unrestrictedToRestrictedIndex(std::vector<std::pair<bool,size_t>>::size_type unrestrictedIndex) const {
        assert(unrestrictedIndex < indexMapper.size());
        return indexMapper[unrestrictedIndex];
      }

    private:
      void restrictCells()
      {
        auto cend = implementation.gridView().template end<0>();
        for(auto cell = implementation.gridView().template begin<0>(); cell!=cend; ++cell)
        {
          size_t cellIndex = implementation.indexSet().index(*cell);
          auto fend = implementation.gridView().iend(*cell);
          for(auto face = implementation.gridView().ibegin(*cell); face!=fend; ++face)
            if(ScalarSpaceDetail::considerFace(*face,implementation.boundaryIds,implementation.usedIds))
              consideredCells[cellIndex] = true;
        }
      }

      size_t insertionIndex;
      std::vector<bool> consideredCells;
      // use the global numbering of standard spaces for fast identification of shared degrees of freedoms
      std::vector<std::pair<bool,size_t> > indexMapper;
      std::vector<std::unique_ptr<ShapeFunctionSet>> shapeFunctionSets; //store non-const copies of restricted shape function sets
    };

  } // end namespace ScalarSpaceDetail
  /**
   * \endcond
   */
} // end namespace Kaskade

#endif
