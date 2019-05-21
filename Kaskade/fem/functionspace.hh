/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FUNCTIONSPACE_HH
#define FUNCTIONSPACE_HH

/**
 * @file
 * @brief  FEFunctionSpace and FunctionSpaceElement and the like
 * @author Martin Weiser
 *
 *
 * This file contains a number of core classes of Kaskade, needed to define a problem. The important ones are:
 *
 * - FEFunctionSpace: information about the finite element space used (reference element, grid,...)
 * - FunctionSpaceElement: FEFunction that lives on a FEFunctionSpace
 *
 * The other classes are auxilliary
 */

#include <cassert>
#include <map>
#include <sstream>

#include <boost/mpl/int.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/multi_array.hpp>
#include <boost/signals2.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>


#include "dune/geometry/quadraturerules.hh"
#include "dune/istl/bvector.hh"

#include "fem/gridmanager.hh"
#include "fem/pshapefunctions.hh"
#include "fem/shapefunctioncache.hh"
#include "linalg/crsutil.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  // forward declarations
  template <class LocalToGlobalMapper> class FEFunctionSpace;

  template <template <class, class> class, class, class> class BoundaryMapper;

  /// \cond internals
  namespace FunctionSpace_Detail
  {
    template <class Mapper>
    struct IsBoundaryMapperStruct : std::false_type {};
    template<template <class, class> class DomainMapper, class Scalar, class GridView>
    struct IsBoundaryMapperStruct<BoundaryMapper<DomainMapper,Scalar,GridView>> : std::true_type {};
    // if Mapper is a BoundaryMapper this returns true, otherwise false
    template <class Mapper>
    constexpr bool isBoundaryMapper() {
      return IsBoundaryMapperStruct<Mapper>::value;
    }

    template <class Mapper>
    struct ChooseDomainMapperStruct {
      using type = void;
    };
    template<template <class, class> class DomainMapper, class Scalar, class GridView>
    struct ChooseDomainMapperStruct<BoundaryMapper<DomainMapper,Scalar,GridView>> {
      using type = DomainMapper<Scalar,GridView>;
    };
    // if Mapper is a BoundaryMapper this gives the underlying FE mapper, otherwise it gives the void type
    template<class Mapper>
    using ChooseDomainMapper = typename ChooseDomainMapperStruct<Mapper>::type;

    // checks whether two LeafGridViews are the same (i.e. whether they belong to the same grid (tested by memory address of grid))
    template<class Grid>
    bool gridViewsAreSame(typename Grid::LeafGridView gv1, typename Grid::LeafGridView gv2) {
      return &(gv1.grid()) == &(gv2.grid());
    }
    // checks whether two LevelGridViews are the same (i.e. whether they belong to the same grid and are of the same level)
    template<class Grid>
    bool gridViewsAreSame(typename Grid::LevelGridView gv1, typename Grid::LevelGridView gv2) {
      return (gv1.template begin<0>()->level() == gv2.template begin<0>()->level()) && (&(gv1.grid()) == &(gv2.grid()));
    }
  }
  /// \endcond

  /** 
   * \ingroup functional
   * \brief A type for describing the nature of a PDE problem.
   *
   * Possibilities:
   * - VariationalFunctional: problem that comes from a minimization problem. This means that test space and ansatz space are identical
   * - WeakFormulation: this means that test space and ansatz space differ
   */
  enum ProblemType { VariationalFunctional, WeakFormulation };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \ingroup fem
   * \brief Extracts the type of FE space with given index from set of spaces.
   *
   * In PDE problems, pointers to the involved spaces are bundled into
   * static vectors to be fed into the assembler. This trait class
   * extracts the space type at given index position.
   *
   * \tparam Spaces a boost::fusion sequence of pointers to FE spaces.
   * \tparam Idx the index of the requested space.
   */
  template <class Spaces, int Idx>
  struct SpaceType
  {
    using type = std::remove_pointer_t<typename boost::fusion::result_of::value_at_c<Spaces,Idx>::type>;
  };


  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief Returns an entity pointer pointing to a cell of the given grid containing the specified global coordinate. 
   * 
   * This method makes use of the grid hierarchy, and is (on refined meshes) significantly faster than a linear search. 
   * If no cell containing the given point is found, a \ref LookupException is thrown.
   *
   * Uses now checkInside from fixdune.hh and a (fixed) tolerance.
   * 
   * \param gv the grid view on which to operate (TODO: does this make sense? We only use the grid itself!)
   * \param global the global position for which a cell shall be found
   */
  template <class GridView>
  typename GridView::template Codim<0>::EntityPointer
  findCell(GridView const& gv, Dune::FieldVector<typename GridView::Grid::ctype,GridView::Grid::dimensionworld> global)
  {
    const double tol = 1e-12 ;

    using namespace Dune;
    typedef typename GridView::Grid::template Codim<0>::Entity          Cell;

    // Do linear search on coarse grid (level 0).
    auto coarseIterator = std::find_if(gv.grid().template lbegin<0>(0), gv.grid().template lend<0>(0),
                                       [&](Cell const& cell){ return checkInside(cell.type(),cell.geometry().local(global)) < tol; });

    if (coarseIterator == gv.grid().template lend<0>(0)) // not found
    {
      std::ostringstream msg;
      msg << "findCell(): global point << " << global << " is not contained in any level 0 cell in the grid!";
      msg.flush();
      throw LookupException(msg.str(),__FILE__,__LINE__);
    }

    // Do a hierarchical search for a leaf cell containing the point.
    typename GridView::Grid::template Codim<0>::EntityPointer ci(coarseIterator);

    int lastFound = 0 ;
    
    for( int level = 0 ; level <= gv.grid().maxLevel() ; level++ )
    {
      if( ci->isLeaf() ) 
        return ci;
      
      for (auto hi = ci->hbegin(level+1); hi != ci->hend(level+1) ; ++hi)
        if( checkInside(hi->type(), hi->geometry().local(global)) < tol )
        {
          ci = typename GridView::Grid::template Codim<0>::EntityPointer(hi);
          lastFound = level;
          break ;
        }
    }

    // We should never get here... For a position outside the domain, no coarse grid cell should have been
    // found, leadin to an exception being raised above. If, on the other hand, a coarse grid cell has been 
    // found, there should be a leaf cell containing the given point as well.
    std::cerr << "findCell() at " << __FILE__ << ':' << __LINE__ << '\n'
        << "Global point << " << global
        << " is not contained in any cell in the grid!\n"
        << " last found on " << lastFound << " of " << gv.grid().maxLevel() << "\n";
    abort();
    return typename GridView::template Codim<0>::EntityPointer(gv.template begin<0>()); // never get here! (just for formal correctness return something)
  }


  //---------------------------------------------------------------------


  /**
   * \ingroup fem
   * \brief A class for representing finite element functions.
   * 
   * A class for representing finite element function space elements
   * (also known as finite element functions).  Finite element functions
   * are piecewisely smooth functions (on each codim 0 entity,
   * i.e. cell, of an associated grid) that may or may not be globally
   * continuous or differentiable.  FE functions support pointwise
   * evaluation and differentiation, but on sets of measure zero (more
   * precisely, on the codim 1 entities, i.e. faces, of the grid), these
   * values are undefined, and calling the corresponding methods may
   * yield arbitrary values.
   *
   * \tparam FunctionSpace Each FunctionSpaceElement is associated with an
   * FEFunctionSpace (it "lives" in this space).  This function space
   * is defined in the first template parameter 
   *	
   * \tparam m FunctionSpaceElements may have multiple components. The number of components 
   * attached to each shape function is given by this template parameter. Note that even 
   * with m=1, the function can be vector-valued if the shape functions themselves are 
   * vector-valued.
   *
   * Usually, FSElements
   * are created as part of a VariableSet, which in turn is created via
   * a VariableSetDescription, containing information about variables
   * and their spaces (VariableDescription s). FSElements contain the
   * data for the variables.
   *
   * 
   * If a GridManager is used, FunctionSpaceElements are automatically
   * prolongated, if the grid is refined, and approximated, if the grid is
   * coarsened.
   *
   * If a VariableSet vs is at hand, the FunctionSpaceElements can be
   * accessed as public data-members in a \ref fusion
   * "boost::fusion::vector" vs.data.  This means: \ref fusion
   * "boost::fusion::at_c<N>(vs.data)" gives a reference to the N'th
   * FSElement in vs.
   *
   * Simple assignments are performed by the assignment operator. For
   * more complex specification of FE function values see \ref spaceTransfer,
   * \ref interpolateGlobally, and \ref interpolateGloballyWeak.
   * 
   */
  template <class FunctionSpace, int m>
  class FunctionSpaceElement
  {
  public:
    /// FEFunctionSpace, this element lives on
    typedef FunctionSpace      Space;
    /// own type
    typedef FunctionSpaceElement <Space, m>  Self;
    /// scalar field type
    typedef typename Space::Scalar Scalar;
    
    /**
     * \brief scalar field type
     * 
     * This is provided for consistency with Dune::BlockVector.
     */
    using field_type = typename Space::field_type;
    
    /// type of the elements of the data vector
    typedef Dune::FieldVector<Scalar,m>                 StorageValueType;

    /**
     * \brief type of coefficient vector entries
     * 
     * This is provided for compatibility with Dune::BlockVector.
     */
    using block_type = StorageValueType;

    /// type of the data vector
    typedef Dune::BlockVector<StorageValueType>         StorageType;

    /// components at each point of evaluation
    static int const components = Space::sfComponents * m;

    /// return type of the function value(...)
    typedef Dune::FieldVector<Scalar,components>        ValueType;
    
    /**
     * \brief return type of the function derivative (...)
     * 
     * The derivative is always computed wrt the global world coordinate system, even for functions defined on
     * a lower dimensional manifold (e.g. in shell or plate problems). In that case, the function is extended 
     * to a neighborhood of the manifold in the canonical way - constant extension normal to the manifold. The 
     * derivative of that extension wrt the global world coordinate system is then returned. Hence, the number 
     * of entries in the derivative is always the world dimension.
     * 
     * Design rationale: The alternative would be to consider the derivative wrt a coordinate system of the 
     * tangent space. The choice of such a coordinate system is, however, not unique.
     */
    typedef Dune::FieldMatrix<Scalar,components,Space::dimworld> DerivativeType;

    /// deprecated - use DerivativeType instead
    typedef DerivativeType GradientType [[deprecated]];
    
    /// return type of the function hessian(...)
    typedef Tensor3<Scalar,components,Space::dimworld,Space::dimworld> HessianType;


    /**
     * \name Constructors & Destructor
     * @{
     */
    /**
     * \brief Constructor.
     * 
     * The function is initialized to zero. Space is a FEFunctionSpace.
     */
    FunctionSpaceElement(Space const& fs):
      data(fs.degreesOfFreedom()), sp(&fs), transferMe(*this), blk(transferConnection)
    {
      blk.unblock();
      data = 0;
      transferConnection = sp->requestToRefine.connect(transferMe);
    }

    /// Copy constructor
    FunctionSpaceElement(Self const& fse):
      data(fse.data), sp(fse.sp), transferMe(*this), blk(fse.blk)
    {
      transferConnection = sp->requestToRefine.connect(transferMe);
    }

    ~FunctionSpaceElement() 
    {
      transferConnection.disconnect();
      sp = 0;
    }
    
    /**
     * @}
     */
    
    /**
     * \name Assignment
     * @{
     */

    /**
     * \brief Copy assignment.
     * 
     * Note that even if the type of function (and hence of the function space) is the same, the actual function spaces
     * may be different (e.g. by different polynomial order). If different spaces are involved, assignment is just interpolation.
     */
    Self& operator=(Self const& fse) 
    {
      if(this != &fse) { // nothing to do on self-assignment
        if (sp == fse.sp)  // the same function space
          data = fse.data;
        else  // different function spaces
          interpolateGlobally<PlainAverage>(*this,fse);
      }
      return *this;
    }

    /**
     * \brief Move assignment.
     * 
     * Note that even if the type of function (and hence of the function space) is the same, the actual function spaces
     * may be different (e.g. by different polynomial order). If different spaces are involved, assignment is just interpolation.
     */
    Self& operator=(Self&& fse) 
    {
      if(this != &fse) { // nothing to do on self-assignment
        if (sp == fse.sp)  // the same function space
        {
          data = std::move(fse.data);
          fse.transferConnection.disconnect();
        }
        else  // different function spaces
          interpolateGlobally<PlainAverage>(*this,fse);
      }
      return *this;
    }

    /**
     * \brief Assignment from functions belonging to other function spaces.
     * 
     * This assignment is defined by interpolation and hence results in pointwise equality only if the
     * right hand side f lies inside the function space represented by our space. In cases where the
     * coefficients of the involved spaces are just subsets (either way round), the interpolation approach
     * is probably of suboptimal performance.
     */
    template <class OtherSpace>
    Self& operator=(FunctionSpaceElement<OtherSpace,m> const& f) { interpolateGlobally<PlainAverage>(*this,f); return *this; }

    /**
     * @brief Assignment from an FE function which is restricted to the boundary to a suitable FE function.
     *
     * If the underlying mappers are the same for both FE function, assignment can be done fast by simply copying the indices.
     * Otherwise interpolation will be used. (TODO: Currently interpolateGlobally with a boundary FE function as second parameter does not work as it should. Therefore, before doing the
     * interpolation the boundary FE function is assigned to a suitable FE function, which is then used for interpolation.)
     */
    template <template <class, class> class DomainMapper>
    Self& operator=(FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper,Scalar,typename Space::GridView>>,m> const& bf) {
      bool hasEquivalentUnderlyingMapper = !FunctionSpace_Detail::isBoundaryMapper<typename Space::Mapper>();
      hasEquivalentUnderlyingMapper = (hasEquivalentUnderlyingMapper && std::is_same<DomainMapper<Scalar,typename Space::GridView>, typename Space::Mapper>::value);
      hasEquivalentUnderlyingMapper = (hasEquivalentUnderlyingMapper && space().mapper().maxOrder()==bf.space().mapper().maxOrder());
      hasEquivalentUnderlyingMapper = (hasEquivalentUnderlyingMapper && FunctionSpace_Detail::gridViewsAreSame<typename Space::GridView::Grid>(space().gridView(),bf.space().gridView()));
      if(hasEquivalentUnderlyingMapper) {
        BoundaryMapper<DomainMapper, Scalar, typename Space::GridView> const& boundaryMapper = bf.space().mapper();
        *this = 0.0;
        for(typename StorageType::size_type vi=0;vi<size();++vi) {
          std::pair<bool,typename StorageType::size_type> restrictedIndex = boundaryMapper.unrestrictedToRestrictedIndex(vi);
          if(restrictedIndex.first) {
            data[vi] = bf[restrictedIndex.second];
          }
        }
      } else {
        // TODO: When interpolateGlobally works correctly (better said: as one would expect) for boundary functions as second argument remove the following three lines
        FEFunctionSpace<DomainMapper<Scalar,typename Space::GridView>> domainSpace(bf.space().gridManager(),bf.space().gridView(),bf.space().mapper().maxOrder());
        typename FEFunctionSpace<DomainMapper<Scalar,typename Space::GridView>>::template Element_t<m> df(domainSpace);
        df = bf;
        interpolateGlobally<PlainAverage>(*this,df);
      }
      return *this;
    }

    /**
     * @brief Assignment from a suitable FE function to an FE function which is restricted to the boundary.
     *
     * If the underlying mappers are the same for both FE function, assignment can be done fast by simply copying the indices.
     * Otherwise interpolation will be used.
     * This overload will only be available if this FE function belongs to a boundary space (otherwise ChooseDomainMapper
     * yields void which is not the same as Mapper and type substitution will fail).
     */
    template <class Mapper>
    std::enable_if_t<std::is_same<Mapper,FunctionSpace_Detail::ChooseDomainMapper<typename Space::Mapper>>::value,Self&>
    operator=(FunctionSpaceElement<FEFunctionSpace<Mapper>,m> const& f) {
      bool hasEquivalentUnderlyingMapper = FunctionSpace_Detail::isBoundaryMapper<typename Space::Mapper>();
      hasEquivalentUnderlyingMapper = (hasEquivalentUnderlyingMapper && space().mapper().maxOrder()==f.space().mapper().maxOrder());
      hasEquivalentUnderlyingMapper = (hasEquivalentUnderlyingMapper && FunctionSpace_Detail::gridViewsAreSame<typename Space::GridView::Grid>(space().gridView(),f.space().gridView()));
      if(hasEquivalentUnderlyingMapper) {
        typename Space::Mapper const& boundaryMapper = space().mapper();
        for(typename StorageType::size_type vi=0;vi<f.size();++vi) {
          std::pair<bool,typename StorageType::size_type> restrictedIndex = boundaryMapper.unrestrictedToRestrictedIndex(vi);
          if(restrictedIndex.first) {
            data[restrictedIndex.second] = f[vi];
          }
        }
      } else {
        interpolateGlobally<PlainAverage>(*this,f);
      }
      return *this;
    }

    /// Assignment from vector of coefficients.
    Self& operator=(StorageType const& v)
    {
      assert(data.size()==v.size());
      data = v;
      return *this;
    }

    /**
     * \brief Sets all coefficients to the given value. As StorageValueType is a vector one value for each component can be specified.
     *
     * Note that this does in general \em not result in a constant function (a notable exception
     * being Lagrange elements, where this property indeed holds).
     */
    Self& operator=(StorageValueType const& a)  { for (int i=0; i<data.N(); ++i) data[i] = a; return *this; }

    /**
     * \brief Sets all coefficients to the given value.
     *
     * Note that this does in general \em not result in a constant function (a notable exception
     * being Lagrange elements, where this property indeed holds).
     */
    Self& operator=(Scalar val) { for(size_t i=0; i<data.N(); ++i) data[i] = val; return *this; }
    
    /**
     * @}
     */
    
    /**
     * \name Basic linear algebra
     */

    /**
     * \brief Addition of a function from the same (type of) space.
     * Note that despite the same type the function spaces may actually be different (e.g. different polynomial order).
     * In this case, interpolation is employed.
     */
    Self& operator+=(Self const& f) 
    {
      if (sp==f.sp) // same function space
        coefficients() += f.coefficients();
      else {
        Self tmp(space()); tmp = f; *this += tmp;
      }
      return *this;
    }

    /**
     * \brief Addition of a function from a different function space, but with the same number of components.
     * 
     * Here, interpolation is used.
     */
    template <class OtherSpace>
    Self& operator+=(FunctionSpaceElement<OtherSpace,m> const& f) 
    {
      Self tmp(space()); tmp = f; *this += tmp;
      return *this;
    }

    /**
     * \brief Addition of coefficient vectors
     */
    Self& operator+=(StorageType const& v) { assert(data.size()==v.size()); data += v; return *this; }

    /**
     * \brief Subtraction of a function from the same (type of) space.
     * 
     * Note that despite the same type the function spaces may actually be different (e.g. different polynomial order).
     * In this case, interpolation is employed.
     */
    Self& operator-=(Self const& f)
    {
      if (sp==f.sp) // same function space
        coefficients() -= f.coefficients();
      else {
        Self tmp(space()); tmp = f; *this -= tmp;
      }
      return *this;
    }

    /**
     * \brief Subtraction of a function from a different function space, but with the same number of components.
     * 
     * Here, interpolation of the given function to our function space is used.
     */
    template <class OtherSpace>
    Self& operator-=(FunctionSpaceElement<OtherSpace,m> const& f) 
    {
      Self tmp(space()); tmp = f; *this -= tmp;
      return *this;
    }

    /**
     * \brief Subtraction of coefficient vectors
     */
    Self& operator-=(StorageType const& v) { assert(data.size()==v.size()); data -= v; return *this; }

//    /**
//     * \brief Computes \f$ y \leftarrow ax + y \f$.
//     */
//    Self& axpy(Scalar a, Self const& x) { assert(sp==x.sp); data.axpy(a,x.data); return *this; }

    /**
     * \brief Computes \f$ y \leftarrow ax + y \f$ for a function x from the same (type of) space.
     * Note that despite the same type the function spaces may actually be different (e.g. different polynomial order).
     * In this case, interpolation is employed.
     */
    Self& axpy(Scalar a, Self const& f)
    {
      if (sp==f.sp) // same function space
        coefficients().axpy(a,f.coefficients());
      else {
        Self tmp(space()); tmp = f; this->axpy(a,tmp);
      }
      return *this;
    }

    /**
     * \brief Computes \f$ y \leftarrow ax + y \f$ for a function x from a different function space, but with the same number of components.
     *
     * Here, interpolation is used.
     */
    template <class OtherSpace>
    Self& axpy(Scalar a, FunctionSpaceElement<OtherSpace,m> const& f)
    {
      Self tmp(space()); tmp = f; this->axpy(a,tmp);
      return *this;
    }

    /**
     * \brief Computes \f$ y \leftarrow av + y \f$ for coefficient vectors
     */
    Self& axpy(Scalar a, StorageType const& v) { assert(data.size()==v.size()); data.axpy(a,v); return *this; }

    /// Scaling
    Self& operator*=(Scalar a) { data *= a; return *this; }
    
    /**
     * @}
     */

    /**
     * \name Function and derivative evaluation
     * @{
     */
    
    /**
     * \brief Evaluates the FE function at the position used to construct the evaluator.  
     * 
     * On codimension 1 entitites, the value is undefined if
     * Space::continuity<0. For vectorial shape functions and m>1, the
     * components are grouped by shape function component first. Example
     * with shape functions with 3 components and m==2:
     * [sf0m0,sf1m0,sf2m0,sf0m1,sf1m1,sf2m1]
     */
    ValueType value(typename Space::Evaluator const& evaluator) const
    {
      ValueType y(0);

      size_t esize = evaluator.size();

      for (size_t i=0; i<esize; ++i) {
        StorageValueType const& coeff = data[evaluator.globalIndices()[i]];
        auto const& sfVal = evaluator.evalData[i].value;
        for (int n=0; n<m; ++n)
          for (int l=0; l<Space::sfComponents; ++l)
            y[n*Space::sfComponents+l] += coeff[n]*sfVal[l];
      }

      return y;
    }

    /**
     * \brief Evaluates the FE function at the specified local position in the given cell.
     */
    ValueType value(typename Space::Grid::template Codim<0>::Entity const& cell,
                    Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimension> const& local) const
    {
      typename Space::Evaluator eval(space());
      eval.moveTo(cell);
      eval.evaluateAt(local);

      return value(eval);
    }

    /**
     * \brief Evaluates the FE function at the specified global position. 
     * 
     * \WARNING This method is comparatively DEAD SLOW.  
     * The components are sorted as described for value().
     */
    ValueType value(Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimensionworld> const& global) const
    {
      auto ci = findCell(space().gridView(),global) ;
      return value(*ci,ci->geometry().local(global));
    }


    /**
     * \brief Evaluates the derivative of the FE function wrt global coordinates at the position used to construct the evaluator.  
     * 
     * On codimension 1 entitites, the value is undefined if Space::continuity<1. The components are sorted as
     * described for value().
     */
    DerivativeType derivative(typename FunctionSpace::Evaluator const& evaluator) const
    {
      DerivativeType dy(0);

      size_t esize = evaluator.size();
      for (size_t i=0; i<esize; ++i)                              // step through all restricted ansatz functions
      {
        auto const& coeff = data[evaluator.globalIndices()[i]];
        auto const& sfGrad = evaluator.evalData[i].derivative;
        for (int n=0; n<m; ++n)                                   // step through all FE function components
          for (int l=0; l<Space::sfComponents; ++l)               // step through all shape function components
            dy[n*Space::sfComponents+l] += coeff[n] * sfGrad[l];
      }

      return dy;
    }
    
    /// deprecated - use derivative instead (will be removed on 2015-12-31)
    DerivativeType gradient(typename FunctionSpace::Evaluator const& evaluator) const { return derivative(evaluator); }
    

    /**
     * \brief Evaluates the derivative of the FE function at the given local position inside the given cell.  
     * 
     * On codimension 1 entitites, the value is undefined if Space::continuity<1. The components are sorted as
     * described for value().
     */
    DerivativeType derivative(typename Space::Grid::template Codim<0>::Entity const& cell,
                              Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimension> const& local) const
    {
      typename Space::Evaluator eval(space());
      eval.moveTo(cell);
      eval.evaluateAt(local);

      return derivative(eval);
    }

    /// deprecated - use derivative instead (will be removed on 2015-12-31)
    GradientType gradient(typename Space::Grid::template Codim<0>::Entity const& cell,
                          Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimension> const& local) const { return derivative(cell,local); }

    /**
     * \brief Evaluates the FE function's gradient at the specified global position. 
     * 
     *  \warning This method is comparatively DEAD SLOW (it has to locate the cell in which the position is).  
     */
    DerivativeType derivative(Dune::FieldVector<typename Space::Grid::ctype,Space::Grid::dimensionworld> const& global) const
    {
      auto ci = findCell(space().gridView(),global) ;

      return derivative(*ci,ci->geometry().local(global));
    }

    /// deprecated - use derivative instead (will be removed on 2015-12-31)
    GradientType gradient(Dune::FieldVector<typename Space::Grid::ctype,Space::Grid::dimensionworld> const& global) const { return derivative(global); }

    /**
     * \brief Evaluates the Hessian of the FE function at the position used to construct the evaluator.  
     * 
     * On codimension 1 entitites, the value is undefined if Space::continuity<2. The components are sorted as
     * described for value().
     */
    HessianType hessian(typename FunctionSpace::Evaluator const& evaluator) const
    {
      HessianType ddy(0);

      assert(evaluator.derivatives() >= 2);

      size_t esize = evaluator.size();
      for (size_t i=0; i<esize; ++i) {
        StorageValueType const& coeff = data[evaluator.globalIndices()[i]];
        auto const& sfHess = evaluator.evalData[i].hessian;
        for (int n=0; n<m; ++n)
          for (int l=0; l<Space::sfComponents; ++l) {
            auto hess = sfHess[l];
            hess *= coeff[n];
            ddy[n*Space::sfComponents+l] += hess;
          }
      }

      return ddy;
    }

    /**
     * \brief Evaluates the Hessian of the FE function at the given local position inside the given cell.  
     * 
     * On codimension 1 entitites, the value is undefined if Space::continuity<2. The components are sorted as
     * described for value().
     */
    HessianType hessian(typename Space::Grid::template Codim<0>::Entity const& cell,
                        Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimension> const& local) const
    {
      typename Space::Evaluator eval(space(),nullptr,2);
      eval.moveTo(cell);
      eval.evaluateAt(local);

      return hessian(eval);
    }

    /**
     * \brief Evaluates the FE function's Hessian at the specified global position. 
     * 
     *  \warning This method is comparatively DEAD SLOW.  
     */
    HessianType hessian(Dune::FieldVector<typename Space::Grid::ctype, Space::Grid::dimensionworld> const& global) const
    {
      auto ci = findCell(space().gridView(),global);
      return hessian(*ci,ci->geometry().local(global));
    }
    
    /**
     * @}
     */


    /**
     * \name Access to coefficient vector
     * @{
     */
    /**
     * \brief Provides read-only access to the coefficients of the finite element basis functions.
     */
    StorageType const& coefficients() const { return data; }

    /**
     * \brief Provides access to the coefficients of the finite element basis functions.
     */
    StorageType& coefficients() { return data; }
    
    /**
     * \brief EXPERIMENTAL Provides access to the coefficients of the finite element basis functions.
     */
    StorageValueType& operator[](size_t i) { return data[i]; }

    /**
     * \brief EXPERIMENTAL Provides read-only access to the coefficients of the finite element basis functions.
     */
    StorageValueType const& operator[](size_t i) const { return data[i]; }
    
    /**
     * @}
     */
    
    /**
     * \name Miscellaneous
     * @{
     */

    /// Number of scalar degrees of freedom.
    int dim() const { return data.dim(); }

    /**
     * \brief Number of coefficient blocks.
     *
     * This is dim()/m. Method is added for interface consistency with
     * Dune::BlockVector.
     */
    int size() const { assert(data.size()*m==dim()); return data.size(); }


    /**
     * \brief Provides access to the FEFunctionSpace to which this function belongs.
     */
    Space const& space() const { return *sp; }

    //// Order of the shape functions
    template<class SFS>
    int order(SFS const& sfs) const
    {
      return sfs.order();
    }

    /**
     * Group of members that are necessary for automatic grid transfer
     * after refinement. The member transferMe is connected to a boost::signal in
     * function space automatically at construction
     */
    void transfer(typename TransferData<Space>::TransferMatrix const& transferMatrix)
    {
      std::unique_ptr<Dune::BlockVector<StorageValueType> > newCoeff = transferMatrix.apply(data);
      data = *newCoeff;
    }

    /// Block grid transfer, i.e., if the grid changes, this FSElement is not transferred
    // why would we need this functionality? It doesn't appear to be used anywhere. Should be
    // removed after 2015-12-31.
    void blockGridTransfer()
    {
      blk.block(); //.block();
    }
    
    /**
     * @}
     */

  private:

    struct TransferMe
    {
      TransferMe(Self& fe) : myFSE(fe) {}

      void operator()(typename TransferData<Space>::TransferMatrix const& transfer, GridManagerBase<typename Space::Grid> const& gm)
      {
        gm.transferFSE(transfer,myFSE);
      }

    private:
      Self &myFSE;
    };



  protected:
    Dune::BlockVector<StorageValueType> data;
    Space const* sp;

  private:
    TransferMe transferMe;
    boost::signals2::connection transferConnection;
    boost::signals2::shared_connection_block blk;
  };
  
  /**
   * \ingroup fem
   * \brief Writes the coefficients into a flat sequence.
   * \related FunctionSpaceElement<Space,m>
   */
  template <class Space, int m, class OutIter>
  OutIter vectorToSequence(FunctionSpaceElement<Space,m> const& v, OutIter& iter) 
  {
    return vectorToSequence(v.coefficients(),iter);
  }
  
  template <class Space, int m, class InIter>
  InIter vectorFromSequence(FunctionSpaceElement<Space,m>& v, InIter in) 
  {
    return vectorFromSequence(v.coefficients(),in);
  }

  namespace FunctionSpace_Detail {
  // if an FEFunction is restricted to the boundary transform it in to a new FEFunction over the whole domain, otherwise leave it unchanged. Used by writeVTK.
    template <class FEFunction>
    class ToDomainRepresentation {
    public:
      ToDomainRepresentation(FEFunction const& f_) : f(f_){}

      FEFunction const& get() const {
        return f;
      }

    private:
      FEFunction const& f;
    };

    template <template <class, class> class DomainMapper, typename Scalar, typename GridView, int m>
    class ToDomainRepresentation<FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper, Scalar, GridView>>, m>> {

      using BoundaryFEFunction = FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper, Scalar, GridView>>, m>;
      using DomainFEFunction = FunctionSpaceElement<FEFunctionSpace<DomainMapper<Scalar, GridView>>, m>;
      using DomainFESpace = typename DomainFEFunction::Space;

    public:
      ToDomainRepresentation(BoundaryFEFunction const& f_) : space(f_.space().gridManager(),f_.space().gridView(),f_.space().mapper().maxOrder()), f(space) {
        f = f_;
      }

      DomainFEFunction const& get() const {
        return f;
      }

    private:
      DomainFESpace space;
      DomainFEFunction f;
    };
  }

  //---------------------------------------------------------------------

  /// Represents a finite element function space, based on a \ref LocalToGlobalMapperConcept "mapper".
  /**
   * \ingroup fem
   * \brief A representation of a finite element function space defined over a domain covered by a grid.  
   * 
   * Elements of the space are defined in
   * terms of shape functions defined over cells (codim 0 entities) of
   * the grid and a LocalToglobalMapper (also known as node manager)
   * combining local shape function coefficients on different cells to
   * one global degree of freedom. Depending on the LocalToGlobalMapper,
   * the finite element functions may be globally continuous or
   * discontinuous, and may support pointwise evaluation of gradients
   * besides the pointwise evaluation of function values.
   *
   * The FEFunctionSpace defines an element type (a class that holds a
   * certain function from the space) as well as evaluation cache
   * classes that store global indices of local shape functions and
   * values and possibly derivatives of shape functions at certain
   * integration points.
   *
   */
  template <class LocalToGlobalMapper>
  class FEFunctionSpace
  {
    typedef FEFunctionSpace<LocalToGlobalMapper> Self;

  public:
    /// Type of grid view
    using GridView = typename LocalToGlobalMapper::GridView;
    /// Type of grid
    typedef typename GridView::Grid                    Grid;
    /// Type of IndexSet
    typedef  typename GridView::IndexSet               IndexSet;

    /** The scalar field type on which the linear (function) space is based. */
    typedef typename LocalToGlobalMapper::Scalar    Scalar;
    
    /**
     * \brief The scalar field type on which the linear (function) space is based.
     * This is provided for consistency with Dune.
     */
    using field_type = Scalar;

    /// Coordinate Type
    typedef typename Grid::ctype ctype;
    /// corresponding \ref LocalToGlobalMapperConcept "mapper"
    typedef LocalToGlobalMapper Mapper;
    /// Grid dimension
    static int const dim = Grid::dimension;
    /// spatial dimension
    static int const dimworld = Grid::dimensionworld;

    /**
     * \brief The number of components of shape functions.
     */
    static int const sfComponents = Mapper::ShapeFunctionSet::comps;

    /**
     * \brief Tells how often elements are globally continuously differentiable.
     * 
     * E.g., 0 means a globally continuous, but not
     * differentiable function, -1 means a discontinuous function (with
     * finite values), -2 may contain Dirac functionals etc.
     */
    static int const continuity = Mapper::continuity;

    /**
     * \brief DEPRECATED use Element_t<m> instead
     */
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<Self,m> type;
    };
    
    /**
     * \brief Defines the element function type.  
     * \tparam  m  number of "components" 
     * A "component" is understood as a full shape function value. This can be a
     * multi-component value as well. With vectorial shape functions
     * (e.g. Nedelec edge elements), m==2 means you end up with 2*dim
     * entries of type Scalar in the function's values.
     */
    template <int m>
    using Element_t = FunctionSpaceElement<Self,m>;

    /**
     * \brief  A helper class that stores all informations necessary to evaluate FE functions at a certain point.
     * 
     * The evaluator caches all the information necessary for evaluations on one cell/face. Information that
     * is independent of the function space is computed once and used for all spaces. Then evaluation essentially
     * reduces to look-ups, thus allowing efficient evaluations.
     * 
     * \related moveEvaluatorsToCell moveEvaluatorsToIntegrationPoint
     */
    struct Evaluator
    {
      typedef typename IndexSet::IndexType IndexType;

      /// The type of function space we belong to.
      typedef Self Space;
      
      /**
       * \brief The type of cells in the grid.
       */
      using Cell = typename Grid::template Codim<0>::Entity;

      /**
       * \brief The number of components of the shape functions.
       */
      static int const sfComponents = Mapper::ShapeFunctionSet::comps;

      /**
       * \brief The type of variational test/ansatz arguments storing value and derivatives/hessians of restricted ansatz functions.
       */
      typedef VariationalArg<Scalar,dimworld,sfComponents> VarArg;

      /**
       * \brief The type of ranges of global indices.
       */
      typedef typename Mapper::GlobalIndexRange GlobalIndexRange;
    
      /**
       * \brief A range of sorted (global index, local index) pairs.
       */
      typedef typename Mapper::SortedIndexRange SortedIndexRange;

      /**
       * \brief Construct a shape function evaluator for a given space, possibly using the given shape function cache for efficiency.
       * 
       * \param deriv how many derivatives of the shape function shall be evaluated (in [0,2]).
       * 
       * Note that calling moveTo() is necessary for the globalIndices
       * to have meaningful values and calling evaluateAt() is necessary
       * for evalData to have meaningful values.
       */
      Evaluator(Self const& space_, ShapeFunctionCache<Grid,Scalar>* cache_=nullptr, int deriv_ = 1) :
        space(space_), sf(nullptr), cache(cache_), sfValues(nullptr), combiner_(nullptr),
        globIdx(space.mapper().initGlobalIndexRange()), sortedIdx(space.mapper().initSortedIndexRange()),
        deriv(deriv_), currentCell(nullptr)
      { }
 
      /**
       * \brief Copy constructor.
       */
      Evaluator(Evaluator const& eval):
        evalData(eval.evalData),
        space(eval.space), sf(eval.sf), sfSize(eval.sfSize), cache(eval.cache), sfValues(eval.sfValues), lastQr(eval.lastQr), lastSubId(eval.lastSubId), combiner_(nullptr),
        globIdx(eval.globIdx), sortedIdx(eval.sortedIdx), deriv(eval.deriv), currentXloc(eval.currentXloc), currentCell(nullptr)
      {
        if (eval.combiner_) {
          combiner_ = std::make_unique<typename Mapper::Combiner>(*eval.combiner_);
        }
        if (eval.currentCell) {
          currentCell = std::make_unique<Cell>(*eval.currentCell);
        }
      }

      /**
       * \brief Assignment.
       */
      Evaluator& operator=(Evaluator const& eval)
      {
        assert(&space==&eval.space);
        sf = eval.sf;
        sfSize = eval.sfSize();
        cache = eval.cache;
        sfValues = eval.sfValues;
        lastQr = nullptr;
        lastSubId = -1;
        converter = eval.converter;
        if (eval.combiner_) {
          combiner_ = std::make_unique<typename Mapper::Combiner>(*eval.combiner_);
        } else {
          combiner_.reset();
        }
        globIdx = eval.globIdx;
        sortedIdx = eval.sortedIdx;
        evalData = eval.evalData;
        deriv = eval.deriv;
        currentXloc = eval.currentXloc;
        if (eval.currentCell) {
          currentCell = std::make_unique<Cell>(*eval.currentCell);
        } else {
          currentCell.reset();
        }
      
        return *this;
      }

      /**
       * \brief Tells the evaluator that subsequent calls refer to the given cell.
       * 
       * Obtains the global indices of the shape functions associated to
       * the given cell and a shape function value converter for
       * transforming values and gradients between global world
       * coordinates and local reference element coordinates.
       *
       * Warning: Never pass a temporary object as cell (as for example returned by intersection.inside() or dereferencing EntityPointer), since the evaluator stores
       * a raw pointer to the cell. Therefore, the cell object must be still available (not destroyed) when later evaluateAt() is called.
       * 
       * \param cell  a reference to a grid cell
       * \param index the index of the cell in the grid view index set. If not provided, it is looked up (which can be somewhat expensive).
       */
      void moveTo(Cell const& cell, IndexType const index)
      {
        sf = &space.mapper().shapefunctions(index);
        sfSize = sf->size();
        assert(sfSize>=0);

        currentCell = std::make_unique<Cell>(cell);

        converter.moveTo(*currentCell);

        combiner_ = std::make_unique<typename Mapper::Combiner>(space.mapper().combiner(*currentCell,index));

        globIdx = space.mapper().globalIndices(index);
        sortedIdx = space.mapper().sortedIndices(index);

        //this->moveToCell(cell);
      }
      
      /**
       * \brief A convenience overload, looking up the cell index itself.
       * 
       * Note that index lookups are not particularly cheap. Provide the index explicitly for performance reasons if already available.
       *
       * Warning: Never pass a temporary object as cell (as for example returned by intersection.inside() or dereferencing EntityPointer), since the evaluator stores
       * a raw pointer to the cell. Therefore, the cell object must be still available (not destroyed) when later evluateAt() is called.
       */
      void moveTo(Cell const& cell)
      {
        moveTo(cell,space.indexSet().index(cell));
      }
    
      /**
       * \brief Tells the evaluator that subsequent calls to evaluateAt(x,qr,ip,subId) refer to the given quadrature rule and subentity.
       * 
       * Providing this information to the evaluator allows to reuse data without further checking and hence improves performance.
       * 
       * \param qr the quadrature rule used
       * \param subId the subindex of the boundary entity (if integrating over a boundary face), or 0 for the cell
       */
      template <int subDim>
      void useQuadratureRule(Dune::QuadratureRule<typename Grid::ctype,subDim> const& qr, int subId)
      {
        // reusing data does not work for boundary mapper in its current (non optimal) implementation,
        // where every cell has its own restricted shape function set.
        // TODO: When the implementation of boundary mapper is improved (i.e. all cells use the same shape function set and
        // restriction is handled by the Combiner) remove the following line.
        if(FunctionSpace_Detail::isBoundaryMapper<Mapper>()) return;
        // check whether we have the data corresponding to the lookup key (qr,subId) already available
        if (sfValues && lastQr==&qr && lastSubId==subId) 
          return;
        
        if (cache) {
          sfValues = &cache->evaluate(*sf,qr,subId); // extract data from cache
          lastQr = &qr;                              // remember extraction key
          lastSubId = subId;
        }
      }

      /**
       * \brief Evaluates the shape functions at the given local coordinate.
       */
      void evaluateAt(Dune::FieldVector<ctype,dim> const& x)
      {
        currentXloc = x;
        evalData.clear();

        converter.setLocalPosition(x);

        evalData.reserve(sfSize);
        for (int i=0; i<sfSize; ++i) 
        {
          typedef VariationalArg<Scalar,dim,sfComponents> SfArg;
          
          typename Mapper::ShapeFunctionSet::value_type const& sfun = (*sf)[i];
          evalData.push_back( derivatives()<=1?
                                converter.global(SfArg(sfun.evaluateFunction(x),sfun.evaluateDerivative(x)),1) :
                                converter.global(SfArg(sfun.evaluateFunction(x),sfun.evaluateDerivative(x),sfun.evaluate2ndDerivative(x)),2) );
        }
        combiner_->rightTransform(evalData);
      }

      /**
       * \brief Evaluates the shape functions at the given integration point of the given quadrature rule. 
       * If possible, an evaluation cache stored in the shape function set is used.
       * \param x local position of integration point
       * \param qr the quadrature rule used
       * \param ip the number of the integration point of the quadrature rule
       * \param subId the subindex of the boundary entity (if integrating over a boundary face), or 0 for the cell
       */
      template <int subDim>
      void evaluateAt(Dune::FieldVector<ctype,dim> const& x,
                      Dune::QuadratureRule<typename Grid::ctype,subDim> const& qr, int ip, int subId)
      {
        if (cache) 
        {
          currentXloc = x;
          evalData.clear();
          converter.setLocalPosition(x);
          
          if (sfValues) 
            for (int i=0; i<sfSize; ++i)
              evalData.push_back( converter.global((*sfValues)[ip][i],derivatives()) );
          else
          {
            auto const& localData = cache->evaluate(*sf,qr,ip,subId);
            
            for (int i=0; i<sfSize; ++i)
              evalData.push_back( converter.global(localData[i],derivatives()) );
          }

          combiner_->rightTransform(evalData);
        } 
        else 
          evaluateAt(x);
      }

      /// Provides access to the global indices.
      GlobalIndexRange globalIndices() const { return globIdx; }
    
      /**
       * \brief Provides access to (global,local) index pairs sorted ascendingly by global index.
       */
      SortedIndexRange sortedIndices() const { return sortedIdx; }

      /// Returns the number of global ansatz functions. This is a
      /// convenience method equivalent to globalIndices().size(). Note
      /// that this may differ from the number of shape functions living
      /// on this cell.
      int size() const 
      { 
        return globalIndices().size(); 
      }


      /// Returns the polynomial order of the shape functions.
      int order() const { return sf->order(); }

      /// Returns the shape function set.
      typename Mapper::ShapeFunctionSet const& shapeFunctions() const { return *sf; }

      /// Returns the algebraic combiner associated to the cell.
      typename Mapper::Combiner const& combiner() const { return *combiner_; }

      /**
       * \brief The container holding values and gradients of the global ansatz functions at the evaluation point.
       * 
       * The entries are sorted according to the local ordering of restricted ansatz functions.
       */
      std::vector<VarArg> evalData;

      /**
       * \brief The grid view on which the space is defined.
       */
      GridView const& gridView() const
      {
        return space.gridView();
      }

      LocalToGlobalMapper const& mapper() const
      {
        return space.mapper();
      }
      
      /**
       * \brief Returns the cell containing the current evaluation point.
       */
      Cell const& cell() const { return *currentCell; }
      
      /**
       * \brief Returns the local coordinates of the current evaluation point.
       */
      Dune::FieldVector<ctype,dim> const& xloc() const { return currentXloc; }
      
      /**
       * \brief The number of derivatives that can be evaluated.
       */
      int derivatives() const { return deriv; }

    private:
      Self const&                                    space;
      typename Mapper::ShapeFunctionSet const*       sf;
      int                                            sfSize;
      ShapeFunctionCache<Grid,Scalar>*               cache;
      typename ShapeFunctionCache<Grid,Scalar>::template DataType<sfComponents>::type const* sfValues;
      void const*                                    lastQr;    // lookup key for cached shape function values
      int                                            lastSubId; // lookup key for cached shape function values
      typename Mapper::Converter                     converter;
      std::unique_ptr<typename Mapper::Combiner>     combiner_; // TODO: why held by pointer? Wouldn't a member variable make more sense and simplify things? (For this, something like a default constructor would be needed)
      GlobalIndexRange                               globIdx;
      SortedIndexRange                               sortedIdx;
      int                                            deriv; // number of derivatives to evaluate
      Dune::FieldVector<ctype,dim>                   currentXloc;      // local coordinates of current evaluation point
      std::unique_ptr<Cell>                          currentCell;   // cell containing current evaluation point (want Evaluator to store its own Entity object (previously it only stored a raw pointer),
                                                                    // unfortunately not all grid implementations (especially AluGrid) provide a default constructor for Entity, so we cannot hold it by value.
      
    }; // end of class Evaluator

    /**
     * \brief Constructs a finite element space from a GridManager, a grid view, and possibly an extra argument passed to the node manager. 
     * 
     * Note that the grid has to exist during the whole lifetime of the function
     * space.
     *
     * The construction from a GridManager, rather than a Grid makes
     * sure, that the FEFunctionSpace is informed about grid
     * modifications, such that data can be transfered from the original
     * grid to the modified grid.
     *
     * \param args a parameter pack that is forwarded to the local to global mapper.
     * If LocalToGlobalMapper is a LagrangeMapper or a HierarchicMapper,
     * then the last argument arg denotes the polynomial order of the
     * finite element space.
     */
    template <typename... Args>
    FEFunctionSpace(GridManagerBase<Grid>& gridMan,
                    GridView const& gridView_,
                    Args... args):
          gridman(gridMan),
          m_gridView(gridView_),
          is(m_gridView.indexSet()),
          mapper_(m_gridView,args...),
          rtf(*this,gridman)
    {
      refConnection = gridman.signals.informAboutRefinement.connect(int(GridSignals::allocResources),rtf);
    }

    /**
     * \brief Copy constructor.
     */
    FEFunctionSpace(const FEFunctionSpace<LocalToGlobalMapper> & fs) :
      gridman(fs.gridman),
      m_gridView(fs.m_gridView),
      is(m_gridView.indexSet()),
      mapper_(fs.mapper()),
      rtf(*this,gridman)
    {
      refConnection = gridman.signals.informAboutRefinement.connect(int(GridSignals::allocResources),rtf);
    }

    ~FEFunctionSpace() 
    {
      refConnection.disconnect();
    }

    /**
     * \brief Returns the dimension of the function space, i.e. the number of
     * degrees of freedom of a one-component function from this space.
     */
    size_t degreesOfFreedom() const { return mapper_.size(); }


    /**
     * \brief Provides access to the grid on which this finite element space is defined.
     */
    Grid const& grid() const { return gridman.grid(); }

    /**
     * \brief Access to the GridManager
     * \todo why is this a const function returning a non-const grid manager reference? 
     */
    GridManagerBase<Grid>&  gridManager() const { return gridman; }

    /**
     * \brief Provides access to the index set defining the part of the grid on which the space is defined.
     */
    IndexSet const& indexSet() const { return is; }

    /**
     * \brief Provides access to the index set defining the potential support of the space.
     */
    GridView const& gridView() const { return m_gridView; }

    /**
     * \brief Provides read access to the node manager.
     */
    LocalToGlobalMapper const & mapper() const { return mapper_; }
    
    /**
     * \brief Provides read-write access to the node manager.
     */
    LocalToGlobalMapper /* const */ & mapper() /* const */ { return mapper_; }

    /**
     * \brief Creates a new evaluator for this space.
     * \param cache an (optional) shape function cache to be used
     * \param deriv number of highest requested derivatives
     */
    Evaluator evaluator(ShapeFunctionCache<Grid,Scalar>* cache=nullptr, int deriv = 1) const { return Evaluator(*this,cache,deriv); }
    
    /**
     * \brief Creates a new function in this function space.
     * \tparam m the number of components
     */
    template <int m>
    Element_t<m> element() const { return Element_t<m>(*this); }

  private:

    // The slot class that connects to the grid manager signal to be informed about the 
    // state of the refinement process.
    struct ReactToRefinement
    {
      ReactToRefinement(Self & space, GridManagerBase<Grid>& gm_) : mySpace(space),gm(gm_)  {};

      void operator()(typename GridSignals::Status const ref)
      {
        switch (ref)
        {
          case GridSignals::BeforeRefinement:
            if(!mySpace.requestToRefine.empty())
              gm.getTransferData(transferData,mySpace);
            break;
          case GridSignals::AfterRefinement:
            mySpace.mapper().update();
            if(!mySpace.requestToRefine.empty())
            {
              assert(transferData!=0);
              gm.getTransferMatrix(transferData,transferMatrix,mySpace);
            }
            break;
          case GridSignals::TransferCompleted:
            if(!mySpace.requestToRefine.empty())
            {
              assert(transferMatrix!=0);
              assert(transferData!=0);
              delete transferMatrix;
              delete transferData;  //reset
            }
        }
      }
      
    private:
      Self&                   mySpace;
      GridManagerBase<Grid>&  gm;
      TransferData<Self>*     transferData;
      typename TransferData<Self>::TransferMatrix* transferMatrix;
    };

    GridManagerBase<Grid>& gridman;
    GridView               m_gridView;
    IndexSet const&        is;
    LocalToGlobalMapper    mapper_;
    ReactToRefinement      rtf;
    boost::signals2::connection refConnection;

  public:
    /**
     * \brief Signal for refinement.
     * 
     * This signal is emitted after mesh refinement/coarsening. Clients are informed with a transfer matrix that allows to
     * compute the new FE coefficients from the old ones. This is used by FEFunctions belonging to this space.
     */
    mutable typename boost::signals2::signal<void (typename TransferData<Self>::TransferMatrix &, GridManagerBase<Grid>const &) > requestToRefine;
  };



  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  
  /**
   * \brief A boost::fusion functor for creating an evaluator from the space,
   * using a shape function cache that is given on construction of the
   * functor.
   * 
   * \tparam Cache the type of  shape function cache to use for construction of evaluators
   * 
   * \param deriv a map from space address to the highest requested derivative.
   */
  template <class Cache>
  class GetEvaluators
  {
  public:
    GetEvaluators(Cache* cache_, std::map<void const*,int> const& deriv_ = std::map<void const*,int>()):
      cache(cache_), deriv(deriv_)
    {}

    template <class Space>
    typename Space::Evaluator operator()(Space const* space) const
    {
      auto it = deriv.find(space);
      int derivatives = it==deriv.end()? 1: it->second; // defaults to first derivatives
      return typename Space::Evaluator(*space,cache,derivatives);
    }

  private:
    Cache* cache;
    std::map<void const*,int> const& deriv;
  };

  /**
   * \brief A boost::fusion functor for extracting the evaluator type from a
   * pointer to a FE space. This is more or less redundant to
   * GetEvaluators, but no shape function cache type needs to be
   * provided.
   */
  class GetEvaluatorTypes
  {
  public:
    template <class Space>
    typename Space::Evaluator operator()(Space const*);
  };

  /**
   * \ingroup fem
   * \brief returns a heterogeneous sequence of evaluators for the given spaces
   * \tparam Spaces a heterogeneous sequence of FEFunctionSpace s
   * \tparam ShapeFunctionCache a type for caching shape function values
   * 
   * \param cache pointer to a shape function cache - if null, no cache is used
   * \param deriv a map from space addresses to highest derivatives requested
   */
  template <class Spaces, class ShapeFunctionCache>
  auto getEvaluators(Spaces const& spaces, ShapeFunctionCache* cache, std::map<void const*,int> const& deriv = std::map<void const*,int>())
  {
//     return boost::fusion::as_vector(boost::fusion::transform(spaces,GetEvaluators<ShapeFunctionCache>(cache,deriv)));
    return boost::fusion::as_vector(boost::fusion::transform(spaces,[&](auto const* space) 
    { 
      auto it = deriv.find(space);
      int derivatives = it==deriv.end()? 1: it->second; // defaults to first derivatives
      return space->evaluator(cache,derivatives); 
    }));
  }
  
  /**
   * \ingroup fem
   * \brief the heterogeneous sequence type of Evaluators for the given spaces
   * \tparam Spaces a heterogeneous sequence type of FEFunctionSpace s
   */
  template <class Spaces>
  using Evaluators = typename boost::fusion::result_of::as_vector<
                          typename boost::fusion::result_of::transform<Spaces,GetEvaluatorTypes>::type
                       >::type;

  //---------------------------------------------------------------------
  /**
   * \brief A boost::fusion functor for creating a variational arg from the space,
   * using a shape function cache that is given on construction of the
   * functor.
   * 
   * WARNING this appears to be unused - can we remove this (2015-05-12)?
   */
  template<class Evaluators>
  class [[deprecated]] GetVarArgs
  {
  public:
    template <class T> struct result {};

    template <class VD>
    struct result<GetVarArgs(VD)>
    {
      static int const sidx = VD::spaceIndex;
      typedef typename boost::fusion::result_of::value_at_c<Evaluators, sidx>::type::VarArg type;
    };
  };

  //---------------------------------------------------------------------
  
  /**
   * \brief Moves all provided evaluators to the given cell.
   */
  template <class Evaluators, class Cell>
  void moveEvaluatorsToCell(Evaluators& evals, Cell const& cell)
  {
    boost::fusion::for_each(evals,[&cell](auto& eval) { eval.moveTo(cell); });
  }

  /**
   * \brief Moves all provided evaluators to the given cell with provided index.
   * 
   * As index lookups may be relatively expensive, providing the index here once for all spaces/evaluators
   * may improve performance.
   */
  template <class Evaluators, class Cell, class Index>
  void moveEvaluatorsToCell(Evaluators& evals, Cell const& cell, Index idx)
  {
    boost::fusion::for_each(evals,[&cell,idx](auto& eval) { eval.moveTo(cell,idx); });
  }
  
//---------------------------------------------------------------------

  /**
   * \ingroup fem
   * \brief Tells all evaluators to use the given quadrature rule on the given subentity.
   */
  template <class Evaluator, class QuadratureRule>
  void useQuadratureRuleInEvaluators(Evaluator& evals, QuadratureRule const& qr, int subId)
  {
    boost::fusion::for_each(evals,[&qr,subId](auto& eval) { eval.useQuadratureRule(qr,subId); });
  }

  //---------------------------------------------------------------------

  /**
   * \ingroup fem
   * \brief Moves all provided evaluators to the given integration point, evaluating the shape functions there.
   * \param evals the evaluators to be moved
   * \param x the evaluation point (local coordinates)
   * \param qr the quadrature rule defining the integration point x
   * \param ip the number of the integration point in the quadrature rule qr
   * \param subId the sub-id of the entity on which the quadrature rule lives, w.r.t. the current cell
   */
  template <class Evaluators, class CoordType, int dim, int subDim>
  void moveEvaluatorsToIntegrationPoint(Evaluators& evals, Dune::FieldVector<CoordType,dim> const& x,
                                        Dune::QuadratureRule<CoordType,subDim> const& qr, int ip, int subId)
  {
    boost::fusion::for_each(evals,[&x,&qr,ip,subId](auto& eval) { eval.evaluateAt(x,qr,ip,subId); });
  }

  /**
   * \ingroup fem
   * \brief Moves all provided evaluators to the given integration point, evaluating the shape functions there.
   * \param evals the evaluators to be moved
   * \param x the evaluation point (local coordinates)
   */
  template <class Evaluators, class CoordType, int dim>
  void moveEvaluatorsToIntegrationPoint(Evaluators& evals, Dune::FieldVector<CoordType,dim> const& x)
  {
    boost::fusion::for_each(evals,[&x](auto& eval) { eval.evaluateAt(x); });
  }

  //---------------------------------------------------------------------
  
  /**
   * \ingroup fem
   * \brief Computes the maximum ansatz order used in a collection of evaluators.
   */
  template <class Evaluators>
  int maxOrder(Evaluators const& evals)
  {
    auto getOrder = [](auto const& eval) -> int { return eval.order(); };
    auto max = [](auto a, auto b) -> int { return std::max(a,b); };
    return boost::fusion::fold(boost::fusion::transform(evals,getOrder),0,max);
  }

} /* end of namespace Kaskade */

#endif
