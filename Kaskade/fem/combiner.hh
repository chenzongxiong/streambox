/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef COMBINER_HH
#define COMBINER_HH

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

namespace Kaskade
{
  /**
   * \ingroup fem
   * \brief A base class for implementing combiners with diagonal structure.
   * 
   * \tparam Scalar a scalar field type, usually double.
   * 
   * This implements the \ref Combiner concept: a matrix \f$ K \f$ mapping a subset
   * of global degrees of freedom (those given by globalIndices()) to
   * local degrees of freedom (shape functions).
   *
   * In some finite element spaces, the matrix \f$ K \f$ is diagonal. This base class
   * eases the implementation
   */
  template <class Scalar = double>
  class DiagonalCombiner {
  public:
    /**
     * Constructor.
     * 
     * \param n the number of shape functions on this cell
     * 
     * The matrix \f$ K \f$ is initialized to the identity. Derived classes shall
     * overwrite the orient member storing the diagonal entries.
     */
    DiagonalCombiner(int n): orient(n,1.0)
    {
    }
    
    /**
     * \brief In-place computation of \f$ A \leftarrow A K \f$.
     * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
     */
    template <class Matrix>
    void rightTransform(Matrix& A) const {
      // multiply from right -> modify columns
      assert(A.M()==orient.size());
      for (int i=0; i<A.N(); ++i)
        for (int j=0; j<A.M(); ++j)
          A[i][j] *= orient[j];
    }
    
    /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
    template <int n, int m>
    void rightTransform(std::vector<VariationalArg<Scalar,n,m> >& v) const {
      assert(v.size()==orient.size());
      for (int i=0; i<v.size(); ++i) {
        v[i].value *= orient[i];
        v[i].derivative *= orient[i];
        // TODO: hessians
      }
    }
    
    /**
     * \brief In-place computation of \f$ A \leftarrow K^+ A \f$.
     * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
     */
    template <class Matrix>
    void leftPseudoInverse(Matrix& A) const {
      // Since K^{-1} = K, this is simple...
      // multiply from left -> modify rows
      assert(A.N()==orient.size());
      for (int i=0; i<A.N(); ++i)
        for (int j=0; j<A.M(); ++j)
          A[i][j] *= orient[i];
    }
    
    /// Implicit conversion to a sparse matrix.
    /// This is just the diagonal.
    operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >() const
    {
      int n = orient.size();
      Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > K(n,n,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::random);
      for (int i=0; i<n; ++i)
        K.incrementrowsize(i);
      K.endrowsizes();
      for (int i=0; i<n; ++i)
        K.addindex(i,i);
      K.endindices();
      for (int i=0; i<n; ++i)
        *K[i].begin() = orient[i];
      return K;
    }
    
  protected:
    // contains the diagonal of K
    std::vector<Scalar> orient;
  };
  
}


#endif
