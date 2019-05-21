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

#ifndef CONSTANTSPACE_HH
#define CONSTANTSPACE_HH

#include <numeric>
#include <limits>

#include "fem/lagrangeshapefunctions.hh"
#include "fem/functionspace.hh"
#include "fem/views.hh"

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace ConstantSpaceDetail {

    template <class T, int digits>
    struct fillBits { static T const value = 2*fillBits<T,digits-1>::value+1; };

    template <class T> struct fillBits<T,0> { static T const value = 0; };

    template <class T>
    struct integral_limits
    {
      /**
       * This should be the value of std::numeric_limits<T>::max() as a
       * compile time constant.
       */
      static T const max = fillBits<T,std::numeric_limits<T>::digits>::value;
    };


  } // End of namespace ConstantSpaceDetail
  /**
   * \endcond
   */

  //---------------------------------------------------------------------

  /**
   * A FE mapper class for FE spaces of globally constant
   * functions. These can be interpreted as scalar variables as well,
   * but the interpretation as constant functions allows for easy
   * integration into the usual assembly framework. Since the FE space
   * has dimension 1, there is exactly one degree of freedom.
   */
  template <class ScalarType, class GV>
  class ConstantMapper
  {
    typedef ConstantMapper<ScalarType,GV> Self;

  public:
    typedef ScalarType RT __attribute__((deprecated));
    typedef ScalarType Scalar;

    typedef GV GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;

    static int const dim = Grid::dimension;
    static int const maxLocalSize = 1;

    typedef typename Grid::template Codim<0>::Entity                            Cell;
    typedef LagrangeSimplexShapeFunctionSet<typename Grid::ctype,dim,ScalarType> ShapeFunctionSetImplementation;  // constant functions are the same for any polytope
    typedef typename LagrangeShapeFunctionSetContainer<typename Grid::ctype,dim,ScalarType>::value_type ShapeFunctionSet;
    typedef int ConstructorArgument;
    typedef RangeView<std::vector<size_t>::const_iterator> GlobalIndexRange;

    typedef std::pair<size_t,int>                  IndexPair;
    typedef std::vector<IndexPair>::const_iterator SortedIndexIterator;
    typedef RangeView<SortedIndexIterator>         SortedIndexRange;
  
    /**
     * \brief Continuity index. Since constant functions are \f$ C^\infty \f$,
     * even analytic, this is infinity (or as close as int comes to
     * infinity).
     */
    // I'd rather use numeric_limits<int>::max(), but this is not a
    // constant expression...
    static int const continuity = ConstantSpaceDetail::integral_limits<int>::max;
    
    /**
     * \brief Whether the ansatz functions have global support (i.e. lead to dense matrices).
     */
    static bool const globalSupport = true;

    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };

    ConstantMapper(GridView const& gridView, int /* dummy */)
    : indexSet(gridView.indexSet()), globIdx(1), sortedIdx(1) {
      globIdx[0] = 0;
    sortedIdx[0] = std::make_pair(0,0);
    }

    virtual ~ConstantMapper() {}

    ShapeFunctionSet const& shapefunctions(Cell const& cell) const
    {
      return lagrangeShapeFunctionSet<typename Grid::ctype,dim,RT>(Dune::GeometryType(Dune::GeometryType::simplex,dim),0);  // constant functions are the same for any polytope
    }

    ShapeFunctionSet const& shapefunctions(size_t) const
    {
      Dune::GeometryType gt; 
      gt.makeSimplex(dim); // constant functions are the same for any polytope
      return lagrangeShapeFunctionSet<typename Grid::ctype,dim,RT>(gt,0);
    }

    /**
     * \brief Returns the maximal polynomial order of shape functions encountered in any cell.
     */
    int maxOrder() const
    {
      return 0;
    }

    /**
     * \brief Returns an empty range just for initialization purposes, since
     * RangeView is not default constructible.
     */
    GlobalIndexRange initGlobalIndexRange() const 
    {
      return GlobalIndexRange(globIdx.begin(),globIdx.begin());
    }

    /**
     * \brief Returns a const sequence containing the global indices of the
     * shape functions associated to this cell. 
     * 
     * Global indices start at 0 and are consecutive.
     */
    GlobalIndexRange globalIndices(typename Grid::template Codim<0>::Entity const& /* cell */) const
    {
      return GlobalIndexRange(globIdx.begin(),globIdx.end());
    }

    /**
     * \brief Returns a const sequence containing the global indices of the
     * shape functions associated to this cell. 
     * 
     * Global indices start at 0 and are consecutive.
     */
    GlobalIndexRange globalIndices(size_t) const 
    {
      return GlobalIndexRange(globIdx.begin(),globIdx.end());
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
    SortedIndexRange sortedIndices(Cell const& ) const 
    {
      return SortedIndexRange(sortedIdx.begin(),sortedIdx.end());
    }
    
    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(size_t) const 
    {
      return SortedIndexRange(sortedIdx.begin(),sortedIdx.end());
    }

    /**
     * \brief Returns the number of global degrees of freedom managed.
     */
    size_t size() const { return 1; }

    /**
     * \brief A class mapping local shape function values and derivatives to global shape function values and derivatives.  
     * 
     * For scalar constant shape functions this is absolutely trivial :).
     */
    class Converter
    {
    public:
      Converter() {}

      Converter(Cell const& /* cell */) {}

      void moveTo(Cell const& /* cell */) { }

      void setLocalPosition(Dune::FieldVector<typename Grid::ctype,Grid::dimension> const& /* x */) {}

      /// Applies the transformation \f$ \psi(x) \f$ to shape function value.
      Dune::FieldMatrix<RT,1,1> global(Dune::FieldMatrix<RT,1,1> const& sf) const { return sf; }


      /**
       * \brief Applies the transformation \f$ \psi \f$ to shape function value, derivative, and Hessian, returning a filled VariationalArg.
       */
      template <class Scalar>
      VariationalArg<Scalar,dim,1> global(VariationalArg<Scalar,dim,1>  const& sf, int deriv=1) const
      {
        VariationalArg<RT,dim,1> phi;
        phi.value = 1;
        phi.derivative = 0;
        phi.hessian = 0;
        return phi;
      }
        
      Dune::FieldMatrix<RT,1,1> local(Dune::FieldMatrix<RT,1,1> const& glob) const { return glob; }
    };

    /**
     * @brief A class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For constant elements, this realizes just the identity.
     */
    class Combiner {
    public:
      Combiner(Cell const& cell) {}

      /**
       * @brief In-place computation of \f$ A \leftarrow A K \f$.
       * Since \f$ K \f$ is the identity, this is a no-op.
       */
      template <class Matrix>
      void rightTransform(Matrix& /* A */) const {}

      /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
      template <int n, int m>
      void rightTransform(std::vector<VariationalArg<RT,n,m> >& /* v */) const {}

      /// In-place computation of \f$ A \leftarrow K^+ A \f$.
      template <class Matrix>
      void leftPseudoInverse(Matrix& /* A */) const {}

      /// Implicit conversion to a sparse matrix.
      /// This is just 1.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<RT,1,1> >() const
      {
        Dune::BCRSMatrix<Dune::FieldMatrix<RT,1,1> > K(1,1,Dune::BCRSMatrix<Dune::FieldMatrix<RT,1,1> >::random);
        K.incrementrowsize(0);
        K.endrowsizes();
        K.addindex(0,0);
        K.endindices();

        *K[0].begin() = 1;
        return K;
      }

    };

    /**
     * Returns a combiner for the given cell.
     */
    Combiner combiner(Cell const& cell, size_t index) const { return Combiner(cell); }



    /**
     * @brief Recomputes the internal information after grid changes.
     * Since there is no internal state of this trivial mapper, this
     * method does exactly nothing.
     */
    void update() {}

  private:
    IndexSet const& indexSet;
    std::vector<size_t> globIdx;
    std::vector<IndexPair> sortedIdx;
  };
} /* end of namespace Kaskade */

#endif
