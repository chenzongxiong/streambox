/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LAGRANGESPACE_HH
#define LAGRANGESPACE_HH

#include <cmath>
#include <numeric>
#include <tuple>

#include "boost/multi_array.hpp"

#include "dune/grid/common/capabilities.hh"

#include "fem/converter.hh"
#include "fem/fixdune.hh"
#include "fem/gridcombinatorics.hh"
#include "fem/functionspace.hh"
#include "fem/lagrangeshapefunctions.hh"
#include "fem/scalarspace.hh"

/**
 * \file
 * \brief Lagrange Finite Elements
 * \author Martin Weiser
 *
 */

namespace Kaskade
{
  /**
   * \brief A local to global mapper for scalar Lagrange bases.
   * 
   * \tparam ScalarType The underlying field type (usually double).
   * \tparam GV The grid view on which the FE space is defined.
   * \tparam restricted 
   */
  template <class ScalarType, class GV, bool restricted=false>
  class LagrangeMapperImplementation
  {
  public:
    typedef ScalarType Scalar;

    typedef GV     GridView;
    typedef typename GridView::Grid     Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename Grid::template Codim<0>::Entity Cell;

    typedef typename LagrangeShapeFunctionSetContainer<typename Grid::ctype,Grid::dimension,Scalar,restricted>::value_type ShapeFunctionSet;
    
    
    /**
     * \brief The index set obtained from gridView().
     */
    IndexSet const& indexSet() const { return is; }
    
    GridView const& gridView() const { return gv; }

    /**
     * \param indexSet_ must contain a subset (or all) of the cells in the grid.
     * \param order_ .
     */
    LagrangeMapperImplementation(GridView const& gridView_, int order_):
      ord(order_), is(gridView_.indexSet()), gv(gridView_)
    {}

    ShapeFunctionSet const& shapeFunctions(Cell const& cell) const
    {
      return lagrangeShapeFunctionSet<typename Grid::ctype,Grid::dimension,Scalar,restricted>(cell.type(),ord);
    }

    ShapeFunctionSet const& lowerShapeFunctions(Cell const& cell) const
    {
      return lagrangeShapeFunctionSet<typename Grid::ctype,Grid::dimension,Scalar,restricted>(cell.type(),ord-1);
    }
    
    ShapeFunctionSet const& emptyShapeFunctionSet() const
    {
      return lagrangeShapeFunctionSet<typename Grid::ctype,Grid::dimension,Scalar,restricted>(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension),-1);
    }

    int order() const { return ord; }

    typedef ScalarConverter<Cell> Converter;

    /**
     * @brief A class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For Lagrange elements, this realizes just the identity.
     */
    class Combiner {
    public:
      template <class Sequence>
      Combiner(Sequence const& s): n(s.size()) {}

      /**
       * @brief In-place computation of \f$ A \leftarrow A K \f$.
       * Since \f$ K \f$ is the identity, this is a no-op.
       */
      template <class Matrix>
      void rightTransform(Matrix& /* A */) const {}

      /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
      template <int n, int m>
      void rightTransform(std::vector<VariationalArg<Scalar,n,m> >& /* v */) const {}

      /// In-place computation of \f$ A \leftarrow K^+ A \f$.
      template <class Matrix>
      void leftPseudoInverse(Matrix& /* A */) const {}

      /// Implicit conversion to a sparse matrix.
      /// This is just the identity.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >() const
      {
        Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > K(n,n,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::random);
        for (int i=0; i<n; ++i)
          K.incrementrowsize(i);
        K.endrowsizes();
        for (int i=0; i<n; ++i)
          K.addindex(i,i);
        K.endindices();
        for (int i=0; i<n; ++i)
          *K[i].begin() = 1;
        return K;
      }

    private:
      int n;
    };

  private:
    int               ord;
    IndexSet const&   is;
    GridView          gv;
  };


  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  template <class ScalarType, class GV>
  class DiscontinuousLagrangeMapperImplementation: public LagrangeMapperImplementation<ScalarType,GV>
  {
    typedef LagrangeMapperImplementation<ScalarType,GV> Base;
    typedef typename Base::Cell Cell;

  public:
    DiscontinuousLagrangeMapperImplementation(GV const& gridView, int order):
      Base(gridView,order)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
      int const ord = this->order();
      int const dim = GV::Grid::dimension;


      if (gt.dim() == dim) // dofs live only on codim 0 entities (cells).
      {
        if ( gt.isSimplex() )
        {
          if (dim==1) return ord+1;
          if (dim==2) return (ord+1)*(ord+2)/2;
          if (dim==3) return (ord+1)*(ord+2)*(ord+3)/6;
          assert("Unknown dimension"==0);
        }
        if ( gt.isCube() ) return pow(ord+1,dim);
        if ( gt.isPyramid() || gt.isPrism() || gt.isNone() ) assert( "Not implemented"==0);
      }

      // subentities with codim > 0 do not carry local degrees of freedom
      return 0;
    }

    // Returns the geometry type, subentity number in cell and subentity
    // codimension for the subentity to which the dof is
    // associated. Here we have a discontinuous discretization, where
    // all dofs are associated to the cells.
    template <class ShapeFunction, class Dummy>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
                     Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Dummy&) const
    {
      gt = cell.type();
      subentity = 0;
      codim = 0;
      indexInSubentity = n;
    }
  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for FEFunctionSpace s of piecewise polynomials of order Order.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class DiscontinuousLagrangeMapper:
  public UniformScalarMapper<DiscontinuousLagrangeMapperImplementation<ScalarType,GV> >
  {
    typedef DiscontinuousLagrangeMapperImplementation<ScalarType,GV> Implementation;
    typedef UniformScalarMapper<Implementation> Base;
    typedef DiscontinuousLagrangeMapper<ScalarType,GV> Self;

  public:
    typedef GV GridView;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,ScalarType> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = -1;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     */
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };


    DiscontinuousLagrangeMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {}
  };

  //---------------------------------------------------------------------

  template <class ScalarType, class GV, class ShapeFunctionFilter=ScalarSpaceDetail::AllShapeFunctions>
  class ContinuousLagrangeMapperImplementation
  : public LagrangeMapperImplementation<ScalarType,GV,std::is_same<ShapeFunctionFilter,ScalarSpaceDetail::RestrictToBoundary>::value>, public ShapeFunctionFilter
  {
    typedef LagrangeMapperImplementation<ScalarType,GV,std::is_same<ShapeFunctionFilter,ScalarSpaceDetail::RestrictToBoundary>::value> Base;
    typedef typename Base::Cell Cell;

  public:
    typedef typename GV::IndexSet IndexSet;
    typedef ScalarType Scalar;
    typedef GV GridView;

    ContinuousLagrangeMapperImplementation(GV const& gridView, int order, ShapeFunctionFilter shapeFunctionFilter=ShapeFunctionFilter())
    : Base(gridView,order), ShapeFunctionFilter(shapeFunctionFilter)
    {}


    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
      int const ord = this->order();

      if(gt.isSimplex())
      {
        if (gt.dim()==0) return 1;
        if (gt.dim()==1) return ord-1;
        if (gt.dim()==2) return std::max(0,(ord-2)*(ord-1)/2);
        if (gt.dim()==3) return std::max(0,(ord-3)*(ord-2)*(ord-1)/6);
        assert("Unknown dimension"==0);
      }
      if(gt.isCube()) return std::max(0.0,pow(ord-1,gt.dim()));
      if(gt.isPyramid() || gt.isPrism()) assert("Not  implemented"==0);

      assert("Unknown geometry type"==0);

      return -1;
    }

    // Returns the geometry type, subentity number in cell and subentity
    // codimension for the subentity to which the dof is
    // associated. Here we have a continuous discretization, where all
    // dofs are associated to the entities on which the shape function
    // node is located.
    template <class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
                     Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      int const dim = Cell::dimension;
      int const ord = this->order();
      IndexSet const& is = this->indexSet();

      int dummy;
      std::tie(dummy,codim,subentity,indexInSubentity) = sf.location();
      gt = is.geomTypes(codim)[0]; // limited to just one geometry type....

      // local index is globally unique in the interior of the cell and
      // on the vertices. Otherwise we need to compute a globally unique
      // numbering based on the global numbering of incident vertices.
      if (codim>0 && codim<dim && ord>2) {
        // Obtain the local barycentric index of the shape function.
        std::array<int,dim+1> xi = barycentric(SimplexLagrangeDetail::tupleIndex<dim>(ord,n),ord);

        int nVertices = dim+1-codim;

        // Obtain a globally unique sorting of the barycentric coordinates, and extract
        // the corresponding permutation and selection of the barycentric coordinates
        // in our subentity.
        GlobalBarycentricPermutation<dim> gbp(is,cell);
        int pi[nVertices];
        gbp.barycentricSubsetPermutation(subentity,codim,pi);
        
        // Since here the node is not located on a lower dimensional
        // subentity, we know that all barycentric indices are at least 1. Therefore we
        // may restrict ourselves to the interior nodes by subtracting 1
        // from all barycentric indices. Note that the "order" associated
        // with these interior nodes shrinks by 1+dimension(our subentity),
        // i.e. dim+1-codim !
        int idx[nVertices];
        for (int i=0; i<nVertices; ++i)
          idx[i] = xi[pi[i]]-1;

        // Finally we compute the local index of this permuted node within the 
        // interior of our subentity.
        indexInSubentity = SimplexLagrangeDetail::local(idx,dim-codim,ord-nVertices);
      }

      // consider restrictions to the boundary here
      this->treatBoundary(data,this->gridView(),cell,codim,subentity);
    }

  };

  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for globally continuous FEFunctionSpace s of piecewise polynomials.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class ContinuousLagrangeMapper:
  public UniformScalarMapper<ContinuousLagrangeMapperImplementation<ScalarType,GV> >
  {
    typedef ContinuousLagrangeMapperImplementation<ScalarType,GV> Implementation;
    typedef UniformScalarMapper<Implementation> Base;

  public:
    typedef ScalarType Scalar;
    typedef GV GridView;
    typedef typename Base::ShapeFunctionSet ShapeFunctionSet;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = 0;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     **/
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<ContinuousLagrangeMapper>,m> type;
    };

    /**
     * \brief Constructor.
     * \param gridView the grid view on which to define the space, usually a leaf grid view
     * \param order polynomial ansatz order of shape functions (> 0)
     */
    ContinuousLagrangeMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {
      assert(order >= 1);
    }
  };
}

#endif
