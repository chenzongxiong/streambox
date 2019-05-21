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

#ifndef HIERARCHICSPACE_HH
#define HIERARCHICSPACE_HH

#include <numeric>
#include <tuple>
#include <type_traits>

#include "boost/multi_array.hpp"

#include "dune/grid/common/capabilities.hh"

#include "fem/converter.hh"
#include "fem/firstless.hh"
#include "fem/fixdune.hh"
#include "fem/functionspace.hh"
#include "fem/gridcombinatorics.hh"
#include "fem/hierarchicshapefunctions.hh"
#include "fem/scalarspace.hh"
#include "fem/sign.hh"
#include "utilities/power.hh"

/**
 * @file
 * @brief Hierarchic Finite Elements
 * @author Martin Weiser
 *
 */

namespace Kaskade
{
  /**
   * \brief A local to global mapper implementation for scalar hierarchic bases.
   *
   * \tparam ScalarType the scalar type
   * \tparam GV the grid view type
   */
  template <class Scalar_, class GV, bool restricted>
  class HierarchicMapperImplementationBase
  {
  public:
    typedef Scalar_                                   RT;
    typedef Scalar_                                   Scalar;
    typedef GV                                        GridView;
    typedef typename GridView::Grid                   Grid;
    typedef typename GridView::IndexSet               IndexSet;
    typedef typename Grid::template Codim<0>::Entity  Cell;

  public:
    static constexpr int dim = Grid::dimension;

    typedef typename HierarchicShapeFunctionSetContainer<typename Grid::ctype,dim,Scalar,restricted>::value_type ShapeFunctionSet;

    /**
     * \param indexSet must contain a subset (or all) of the cells in the grid.
     * \param order .
     */
    HierarchicMapperImplementationBase(GridView const& gridView_, int order):
      ord(order), gv(gridView_), is(gridView_.indexSet())
    {}

    ShapeFunctionSet const& emptyShapeFunctionSet() const
    {
      return hierarchicShapeFunctionSet<typename Grid::ctype,dim,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension),-1);
    }

    IndexSet const& indexSet() const { return is; }

    GridView const& gridView() const { return gv; }

    /**
     * \brief Returns the nominal order of the shape function set.
     */
    int order() const { return ord; }
    
    /**
     * A class mapping local shape function values and gradients to
     * global shape function values and gradients.  For scalar
     * hierarchic shape functions this is easy, only the gradients are
     * affected by the affine transformation from the reference element
     * to an actual cell.
     */
    typedef ScalarConverter<Cell> Converter;

    /**
     * A class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For hierarchic elements, this realizes just a sign change,
     * i.e. \f$ K \f$ is a diagonal matrix with entries +1 and -1.
     */
    class Combiner {
    public:
      Combiner(typename UniformScalarMapperCombinerArgument<boost::compressed_pair<bool,int> >::type const& s_): s(s_) {}

      /**
       * \brief In-place computation of \f$ A \leftarrow A K \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void rightTransform(Matrix& A) const {
        for (int j=0; j<A.M(); ++j) {
          if (sign(j)<0)
            for (int i=0; i<A.N(); ++i)
              A[i][j] = -A[i][j];
        }
      }

      /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
      template <int n, int m>
      void rightTransform(std::vector<VariationalArg<Scalar,n,m> >& v) const {
        for (int i=0; i<v.size(); ++i)
          if (sign(i)<0) {
            v[i].value = -v[i].value;
            v[i].derivative = -v[i].derivative;
          }
      }

      /**
       * \brief In-place computation of \f$ A \leftarrow K^+ A \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void leftPseudoInverse(Matrix& A) const {
        for (int i=0; i<A.N(); ++i)
          if (sign(i)<0)
            for (int j=0; j<A.M(); ++j)
              A[i][j] = - A[i][j];
      }

      /// Implicit conversion to a sparse matrix.
      /// This is just a signed permutation matrix.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >() const
        {
        int n = s.size();
        Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > K(n,n,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::random);
        for (int i=0; i<n; ++i)
          K.incrementrowsize(i);
        K.endrowsizes();
        for (int i=0; i<n; ++i)
          K.addindex(i,i);
        K.endindices();
        for (int i=0; i<n; ++i)
          *K[i].begin() = sign(i);
        return K;
        }

    private:
      int sign(int i) const { return s[i].second().second(); }

      typename UniformScalarMapperCombinerArgument<boost::compressed_pair<bool,int> >::type s;
    };

  private:
    int             ord;
    GridView        gv;
    IndexSet const& is;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \tparam Scalar scalar type
   * \tparam GV grid view
   */
  template <class Scalar, class GV, bool restricted=false>
  class HierarchicMapperImplementation: public HierarchicMapperImplementationBase<Scalar,GV,restricted>
  {
    typedef HierarchicMapperImplementationBase<Scalar,GV,restricted> Base;

  public:
    HierarchicMapperImplementation(GV const& gridView, int order):
      Base(gridView,order),ord(order)
    {}

    typename Base::ShapeFunctionSet const& shapeFunctions(typename Base::Cell const& cell, int ordr=0) const
    {
      return hierarchicShapeFunctionSet<typename GV::Grid::ctype,GV::Grid::dimension,Scalar>(cell.type(),ord);
    }

    typename Base::ShapeFunctionSet const& lowerShapeFunctions(typename Base::Cell const& cell) const
    {
      return hierarchicShapeFunctionSet<typename GV::Grid::ctype,GV::Grid::dimension,Scalar>(cell.type(),ord-1);
    }

  private:
    int ord;
  };

  //---------------------------------------------------------------------

  template <class ScalarType, class GV>
  class DiscontinuousHierarchicMapperImplementation: public HierarchicMapperImplementation<ScalarType,GV>
  {
    typedef HierarchicMapperImplementation<ScalarType,GV> Base;
    typedef typename Base::Cell Cell;

  public:
    DiscontinuousHierarchicMapperImplementation(GV const& gridView, int order):
      Base(gridView,order)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
      int n = 0;
      if (gt.dim() == GV::Grid::dimension) {
        // dofs live only on codim 0 entities (cells).
        for (int p=0; p<=this->order(); ++p)
          n += HierarchicDetail::size(GV::Grid::dimension,p);
      }

      return n;
    }

    template <class Cell, class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& /*sf*/, int n,
        Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      codim = 0;
      gt = cell.type();
      subentity = 0;
      indexInSubentity = n;
      //sign = 1;
    }

  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for FEFunctionSpace s of piecewise polynomials of order Order.
   *
   * \tparam Scalar scalar type
   * \tparam GV grid view
   */
  template <class Scalar, class GV>
  class DiscontinuousHierarchicMapper:
  public UniformScalarMapper<DiscontinuousHierarchicMapperImplementation<Scalar,GV>,boost::compressed_pair<bool,int> >
  {
    typedef DiscontinuousHierarchicMapperImplementation<Scalar,GV> Implementation;
    typedef UniformScalarMapper<Implementation,boost::compressed_pair<bool,int> > Base;
    typedef DiscontinuousHierarchicMapper<Scalar,GV> Self;
  public:
    typedef GV GridView;
    typedef int ConstructorArgument;
    typedef HierarchicSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    static int const continuity = -1;

    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };


    DiscontinuousHierarchicMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {}
  };

  //---------------------------------------------------------------------

  template <class Scalar, class GV, class ShapeFunctionFilter=ScalarSpaceDetail::AllShapeFunctions>
  class  ContinuousHierarchicMapperImplementation: public HierarchicMapperImplementation<Scalar,GV,ScalarSpaceDetail::isRestriction<ShapeFunctionFilter>()>, public ShapeFunctionFilter
  {
    typedef HierarchicMapperImplementation<Scalar,GV,ScalarSpaceDetail::isRestriction<ShapeFunctionFilter>()> Base;
    typedef typename Base::Cell Cell;
    static constexpr int dim = GV::dimension;

  public:
    ContinuousHierarchicMapperImplementation(GV const& gridView, int order, ShapeFunctionFilter shapeFunctionFilter = ShapeFunctionFilter()):
      Base(gridView,order), ShapeFunctionFilter(shapeFunctionFilter)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.s
    int dofOnEntity(Dune::GeometryType gt) const
    {
      using namespace HierarchicDetail;

      int codim = GV::Grid::dimension-gt.dim();

      int n = 0;
      for (int p=0; p<=this->order(); ++p)
        n += cdimSize(GV::Grid::dimension,p,codim) / nSubentities(GV::Grid::dimension,codim);

      return n;
    }

    template <class Cell, class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
        Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      using namespace HierarchicDetail;
      typename Base::IndexSet const& is = this->indexSet();

      // Compute global ids of cell vertices.
      int vIds[dim+1];
      vertexids(dim,0,0,vIds);
      for (int j=0; j<=dim; ++j)
        vIds[j] = this->indexSet().subIndex(cell,j,dim);

      int p, k;
      std::tie(p,codim,subentity,k) = sf.location();
      gt = is.geomTypes(codim)[0]; // limited to just one geometry type....

      // Compute tuple index corresponding to vertex permutation.
      std::tie(k,data.second()) = sfPermutation(dim,p,codim,subentity,k,vIds);

      // Compute index in subentity taking all preceding lower order
      // shapefunctions into account.
      indexInSubentity = 0;
      for (int o=0; o<p; ++o)
        indexInSubentity += cdimSize(dim,o,codim)/nSubentities(dim,codim);
      indexInSubentity += k;

      // consider restrictions to the boundary here
      this->treatBoundary(data,this->gridView(),cell,codim,subentity);
    }
  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for globally continuous FEFunctionSpace s of piecewise polynomials.
   *
   * \tparam Scalar scalar type
   * \tparam GV grid view
   */
  template <class Scalar, class GV>
  class ContinuousHierarchicMapper:
  public UniformScalarMapper<ContinuousHierarchicMapperImplementation<Scalar,GV>,boost::compressed_pair<bool,int> >
  {
    typedef ContinuousHierarchicMapperImplementation<Scalar,GV> Implementation;
    typedef UniformScalarMapper<Implementation,boost::compressed_pair<bool,int>> Base;
    typedef ContinuousHierarchicMapper<Scalar,GV> Self;

  public:
    typedef GV GridView;
    typedef int ConstructorArgument;
    typedef HierarchicSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    static int const continuity = 0;

    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };

    ContinuousHierarchicMapper(GridView const& gridView, int order)
    :  Base(Implementation(gridView,order)) {
      assert(order >= 0);
    }
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  template <class Scalar, class GV, bool restricted=false>
  class HierarchicExtensionMapperImplementation: public HierarchicMapperImplementationBase<Scalar,GV,restricted>
  {
    typedef HierarchicMapperImplementationBase<Scalar,GV,restricted> Base;

  public:
    HierarchicExtensionMapperImplementation(GV const& gridView, int order):
      Base(gridView,order)
    {}

    typename Base::ShapeFunctionSet const& shapeFunctions(typename Base::Cell const& cell,int ordr=0) const
    {
      return hierarchicExtensionShapeFunctionSet<typename GV::Grid::ctype,GV::Grid::dimension,Scalar>(cell.type(),this->order());
    }

    typename Base::ShapeFunctionSet const& lowerShapeFunctions(typename Base::Cell const& cell) const
    {
      return hierarchicExtensionShapeFunctionSet<typename GV::Grid::ctype,GV::Grid::dimension,Scalar>(cell.type(),-1);
    }
  };


  template <class Scalar, class GV>
  class DiscontinuousHierarchicExtensionMapperImplementation
  : public HierarchicExtensionMapperImplementation<Scalar,GV>
  {
    typedef HierarchicExtensionMapperImplementation<Scalar,GV> Base;
    typedef typename Base::Cell Cell;

  public:
    DiscontinuousHierarchicExtensionMapperImplementation(GV const& gridView, int order):
      Base(gridView,order)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
      // dofs live only on codim 0 entities (cells).
      return gt.dim()==GV::Grid::dimension? HierarchicDetail::size(GV::Grid::dimension,this->order()): 0;
    }

    template <class Cell, class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& /*sf*/, int n,
        Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      codim = 0;
      gt = cell.type();
      subentity = 0;
      indexInSubentity = n;
      //sign = 1;
    }

  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for FEFunctionSpace s of piecewise polynomials of given order.
   *
   * \tparam Scalar scalar type
   * \tparam GV grid view
   */
  template <class Scalar, class GV>
  class DiscontinuousHierarchicExtensionMapper:
  public UniformScalarMapper<DiscontinuousHierarchicExtensionMapperImplementation<Scalar,GV>,boost::compressed_pair<bool,int> >
  {
    typedef DiscontinuousHierarchicExtensionMapperImplementation<Scalar,GV> Implementation;
    typedef UniformScalarMapper<Implementation,boost::compressed_pair<bool,int> > Base;
    typedef DiscontinuousHierarchicExtensionMapper<Scalar,GV> Self;
  public:

    typedef GV GridView;
    typedef int ConstructorArgument;
    typedef HierarchicSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    static int const continuity = -1;

    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };


    DiscontinuousHierarchicExtensionMapper(GV const& gridView, int order):
      Base(Implementation(gridView,order))
    {}
  };

  //---------------------------------------------------------------------

  template <class Scalar, class GV, class ShapeFunctionFilter = ScalarSpaceDetail::AllShapeFunctions>
  class  ContinuousHierarchicExtensionMapperImplementation
  : public HierarchicExtensionMapperImplementation<Scalar,GV,ScalarSpaceDetail::isRestriction<ShapeFunctionFilter>()>, public ShapeFunctionFilter
  {
    typedef HierarchicExtensionMapperImplementation<Scalar,GV,ScalarSpaceDetail::isRestriction<ShapeFunctionFilter>()> Base;
    typedef typename Base::Cell Cell;

  public:
    ContinuousHierarchicExtensionMapperImplementation(GV const& gridView, int order, ShapeFunctionFilter shapeFunctionFilter = ShapeFunctionFilter()):
      Base(gridView,order), ShapeFunctionFilter(shapeFunctionFilter)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
// std::cerr << "   ContinuousHierarchicExtensionMapperImplementation::dofOnEntity\n";   
      using namespace HierarchicDetail;

      int codim = GV::Grid::dimension-gt.dim();
      return cdimSize(GV::Grid::dimension,this->order(),codim) / nSubentities(GV::Grid::dimension,codim);
    }

    template <class Cell, class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
        Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      using namespace HierarchicDetail;

      // Compute global ids of cell vertices.
      int const d = Cell::dimension;
      int vIds[d+1];
      vertexids(d,0,0,vIds);
      for (int j=0; j<=d; ++j)
        vIds[j] = this->indexSet().subIndex(cell,j,d);

      int nominalOrder, k;
      std::tie(nominalOrder,codim,subentity,k) = sf.location();
      gt = this->indexSet().geomTypes(codim)[0]; // limited to just one geometry type....

      // Compute tuple index corresponding to vertex permutation.
      std::tie(indexInSubentity,data.second()) = sfPermutation(d,nominalOrder,codim,subentity,k,vIds);

      // consider restrictions to the boundary here
      this->treatBoundary(data,this->gridView(),cell,codim,subentity);
    }
  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for globally continuous FEFunctionSpace s of piecewise polynomials.
   *
   * \tparam Scalar scalar type
   * \tparam GV grid view
   */
  template <class Scalar, class GV>
  class ContinuousHierarchicExtensionMapper:
  public UniformScalarMapper<ContinuousHierarchicExtensionMapperImplementation<Scalar,GV>,boost::compressed_pair<bool,int>>
  {
    typedef ContinuousHierarchicExtensionMapperImplementation<Scalar,GV> Implementation;
    typedef UniformScalarMapper<Implementation,boost::compressed_pair<bool,int>> Base;
    typedef ContinuousHierarchicExtensionMapper<Scalar,GV> Self;

  public:
    typedef GV GridView;
    typedef int ConstructorArgument;
    typedef HierarchicSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    static int const continuity = 0;

    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };

    ContinuousHierarchicExtensionMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {
      assert(order >= 0);
    }
  };

} // End of namespace Kaskade

//-------------------------------------------------------------------

#endif
