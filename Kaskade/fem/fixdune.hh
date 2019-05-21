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

#if !defined(FIXDUNE_HH)
#define FIXDUNE_HH

#include <type_traits>

#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/geometry/type.hh"
#include "dune/geometry/referenceelements.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/bcrsmatrix.hh"


/**
 * \file fixdune.hh
 * \brief This file contains various utility functions that augment the basic functionality of Dune.
 */


namespace Dune
{

  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (s,x) \mapsto sx \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator*(T s, Dune::FieldVector<T,n> x)
  {
    x *= s;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (x,s) \mapsto sx \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator*(Dune::FieldVector<T,n> x, T s)
  {
    x *= s;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (s,x) \mapsto sx \f$
   */
  template <class T, int n, class S,
            class enable = typename std::enable_if<std::is_arithmetic<S>::value>::type>
  Dune::FieldVector<T,n> operator*(Dune::FieldVector<T,n> x, S& s)
  {
    x *= static_cast<T>(s);
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (x,s) \mapsto sx \f$
   */
  template <class T, int n, class S,
            class enable = typename std::enable_if<std::is_arithmetic<S>::value>::type>
  Dune::FieldVector<T,n> operator*(S& s, Dune::FieldVector<T,n> x)
  {
    x *= static_cast<T>(s);
    return x;
  }
  
  /// Division of vector by scalar
  template <class T, int n>
  Dune::FieldVector<T,n> operator/(Dune::FieldVector<T,n> x, T s)
  {
    x /= s;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Vector addition \f$ (x,y) \mapsto x+y \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator+(Dune::FieldVector<T,n> x, Dune::FieldVector<T,n> const& y)
  {
    x += y;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Vector negation \f$ x \mapsto - x \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator-(Dune::FieldVector<T,n> const& x)
  {
    Dune::FieldVector<T,n> y;
    for (int i=0; i<n; ++i)
      y[i] = -x[i];
    return y;
  }
  
  
  /// Vector conversion
  template <class T, int n>
  Dune::FieldVector<T,n> asVector(Dune::FieldMatrix<T,n,1> const& x)
  {
    Dune::FieldVector<T,n> y;
    for (int i=0; i<n; ++i)
      y[i] = x[i][0];
    return y;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief outer vector product \f$ (x,y) \mapsto x y^T \f$.
   */
  template <class T, int n>
  Dune::FieldMatrix<T,n,n> outerProduct(Dune::FieldVector<T,n> const& x, Dune::FieldVector<T,n> const& y) {
    Dune::FieldMatrix<T,n,n> A;
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
        A[i][j] = x[i]*y[j];
    return A;
  }
  
  
#ifndef DUNE_ALBERTA_ALGEBRA_HH
  /**
   * \ingroup linalgbasic
   * \brief vector product \f$ (x,y) \mapsto x \times y \f$.
   */
  template <class T>
  Dune::FieldVector<T,3> vectorProduct(Dune::FieldVector<T,3> const& x, Dune::FieldVector<T,3> const& y)
  {
    Dune::FieldVector<T,3> t;
    for (int i=0; i<3; ++i)
      t[i] = x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3];
    return t;
  }
#endif
  
  /**
   * \ingroup linalgbasic
   * \brief Scalar multiplication \f$ (s,A) \mapsto sA \f$.
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator*(T s, Dune::FieldMatrix<T,n,m> x)
  {
    x *= s;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Scalar division \f$ (A,s) \mapsto (1/s)*A \f$.
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator/(Dune::FieldMatrix<T,n,m> A, T s)
  {
    A *= 1/s;
    return A;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Matrix-vector multiplication \f$ (A,x) \mapsto Ax \f$
   */
  template <class T, class S, int n, int m>
  Dune::FieldVector<S,n> operator*(Dune::FieldMatrix<T,n,m> const& A,
                                   Dune::FieldVector<S,m> const& x)
  {
    Dune::FieldVector<S,n> b(0);
    //A.umv(x,b);
    for (int i=0; i<n; ++i)
      for (int j=0; j<m; ++j)
        b[i] += A[i][j]*x[j];
    return b;
  }

  /**
   * \ingroup linalgbasic
   * \brief Matrix addition \f$ (A,B) \mapsto A+B \f$
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator+(Dune::FieldMatrix<T,n,m> const& A,
                                     Dune::FieldMatrix<T,n,m> const& B)
  {
    Dune::FieldMatrix<T,n,m> C;
    for (int i=0; i<n; ++i)
      C[i] = A[i]+B[i]; // vector addition of rows
    return C;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Matrix subtraction \f$ (A,B) \mapsto A-B \f$
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator-(Dune::FieldMatrix<T,n,m> const& A,
                                     Dune::FieldMatrix<T,n,m> const& B)
  {
    Dune::FieldMatrix<T,n,m> C;
    for (int i=0; i<n; ++i)
      C[i] = A[i]-B[i]; // vector addition of rows
    return C;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Matrix negation \f$ A \mapsto - A \f$
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator-(Dune::FieldMatrix<T,n,m> const& A)
  {
    Dune::FieldMatrix<T,n,m> C;
    for (int i=0; i<n; ++i)
      C[i] = -A[i]; // vector negation of rows
    return C;
  }
  
  /// \internal 
  namespace FixDuneDetail
  {
    template <class T, int n, int m, int k, int l>
    struct Mm;
    
    // standard matrix-matrix product
    template <class T, int n, int m, int l>
    struct Mm<T,n,m,m,l>
    {
      static auto product(Dune::FieldMatrix<T,n,m> const& A, Dune::FieldMatrix<T,m,l> const& B)
      {
        Dune::FieldMatrix<T,n,l> AB;
        for (int i=0; i<n; ++i)
          for (int j=0; j<l; ++j) 
          {
            T x = 0;
            for (int k=0; k<m; ++k)
              x += A[i][k]*B[k][j];
            AB[i][j] = x;
          }
          return AB;
      }
    };
    
    template <class T, int n, int m>
    struct Mm<T,1,1,n,m>
    {
      static auto product(Dune::FieldMatrix<T,1,1> const& A, Dune::FieldMatrix<T,n,m> B)
      {
        B *= A[0][0];
        return B;
      }
    };

    template <class T, int n, int m>
    struct Mm<T,n,m,1,1>
    {
      static auto product(Dune::FieldMatrix<T,n,m> A, Dune::FieldMatrix<T,1,1> const& B)
      {
        A *= B[0][0];
        return A;
      }
    };
    template <class T>
    struct Mm<T,1,1,1,1>
    {
      static auto product(Dune::FieldMatrix<T,1,1> A, Dune::FieldMatrix<T,1,1> const& B)
      {
        A *= B[0][0];
        return A;
      }
    };
  }
  /// \endinternal
  
  /**
   * \ingroup linalgbasic
   * \brief Matrix multiplication \f$ (A,B) \mapsto AB \f$
   */
  template <class T, int n, int m, int k, int l>
  auto operator*(Dune::FieldMatrix<T,n,m> const& A, Dune::FieldMatrix<T,k,l> const& B)
  {
    return FixDuneDetail::Mm<T,n,m,k,l>::product(A,B);
  }
}




template <class Scalar, int n>
Dune::FieldVector<Scalar,n*n> vectorize(Dune::FieldMatrix<Scalar,n,n> const& A)
{
  Dune::FieldVector<Scalar,n*n> result(0.);
  
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      result[i+j*n] = A[i][j];
    
    return result;
}

template <int n, class Scalar>
Dune::FieldMatrix<Scalar,n,n> unvectorize(Dune::FieldVector<Scalar,n*n> const& v)
{
  Dune::FieldMatrix<Scalar,n,n> result(0.);
  
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      result[i][j] = v[i+j*n];
    
    return result;
}




/**
 * \ingroup linalgbasic
 * \brief Matrix contraction (Frobenius product) \f$ (A,B) \mapsto A:B = \sum_{i,j} A_{ij} B_{ij} \f$
 */
template <class T, int n, int m>
T contraction(Dune::FieldMatrix<T,n,m> const& A, Dune::FieldMatrix<T,n,m> const& B) 
{
  T x = 0;
  for (int i=0; i<n; ++i)
    for (int j=0; j<m; ++j)
      x += A[i][j] * B[i][j];
  return x;
}

/**
 * \ingroup linalgbasic
 * \brief Matrix trace \f$ A \mapsto \mathrm{tr}\,A = \sum_{i} A_{ii} \f$
 *
 * Note: Instead of trace(A*B) use the equivalent, but faster contraction(A,B).
 */
template <class T, int n>
T trace(Dune::FieldMatrix<T,n,n> const& A) 
{
  T x = 0;
  for (int i=0; i<n; ++i)
    x += A[i][i];
  return x;
}

/**
 * \ingroup linalgbasic
 * \brief Matrix determinant \f$ A \mapsto \mathrm{det}\,A  \f$
 * \warning DEPRECATED, use A.determinant() directly, as provided by Dune.
 */
template <class T, int n>
T determinant(Dune::FieldMatrix<T,n,n> const& A) 
{
  return A.determinant();
}

/**
 * \ingroup linalgbasic
 * \brief Matrix diagonal as a vector.
 */
template <class T, int n>
Dune::FieldVector<T,n> diag(Dune::FieldMatrix<T,n,n> const& A)
{
  Dune::FieldVector<T,n> d;
  for (int i=0; i<n; ++i)
    d[i] = A[i][i];
  return d;
}

/**
 * \ingroup linalgbasic
 * \brief Returns the identity matrix of size n \f$ \mapsto I \f$
 */
template <class T, int n>
Dune::FieldMatrix<T,n,n> unitMatrix() 
{
  Dune::FieldMatrix<T,n,n> I(0);
  for (int i=0; i<n; ++i)
    I[i][i] = 1;
  return I;
}

/**
 * \ingroup linalgbasic
 * \brief Matrix transposition \f$ A \mapsto A^T \f$
 */
template <class T, int n, int m>
Dune::FieldMatrix<T,m,n> transpose(Dune::FieldMatrix<T,n,m> const& A) 
{
  Dune::FieldMatrix<T,m,n> At;
  for (int i=0; i<n; ++i)
    for (int j=0; j<m; ++j)
      At[j][i] = A[i][j];
  return At;
}

/**
 * \brief Returns the matrix type or its transposed type, depending on the given transpose flag.
 */
template <class Matrix, bool transposed>
struct Transpose {};

template <class Scalar, int n, int m, bool transposed>
struct Transpose<Dune::FieldMatrix<Scalar,n,m>,transposed>
{
  typedef Dune::FieldMatrix<Scalar,transposed? m:n,transposed? n:m> type;
};


/**
 * \ingroup linalgbasic
 * \brief Computes Z = X*Y.
 * X and Y need to be compatible, i.e. X.M()==Y.N()
 * has to hold. No aliasing may occur, i.e. &Z != &X and &Z != &Y must
 * hold. Z is reshaped as needed.
 */
template <class MatrixZ, class Matrix>
void MatMult(MatrixZ& z, Matrix const& x, Matrix const& y) 
{
  assert(((void*)&z != (void*)&x) && ((void*)&z != (void*)&y));
  assert(x.M() == y.N());
  
  z.setSize(x.N(),y.M());
  
  for (int i=0; i<z.N(); i++ ) {
      for (int j=0; j<z.M(); j++ ) {
          typename MatrixZ::block_type temp = 0.0;
          for (int k=0; k<x.M(); k++)
            temp += x[i][k]*y[k][j];
          z[i][j] = temp;
        }
    }
}

//---------------------------------------------------------------------

namespace Kaskade {
  
  /**
   * \ingroup linalgbasic
   * \brief A class for fixed size tensors of third order.
   *
   * This class is used to represent second derivatives of vector valued functions. In this
   * case, the convention is that \f[ h[i][j][k] = \frac{\partial f_i}{\partial x_j \partial x_k}. \f]
   *
   * \tparam T the entry type
   * \tparam n1 the first tensor dimension
   * \tparam n2 the second tensor dimension
   * \tparam n3 the third tensor dimension
   */
  // TODO: should this be implemented as Dune::FieldVector<Dune::FieldMatrix<T,n2,n3>,n1> ?
  template <class T, int n1, int n2, int n3>
  class Tensor3 {
  public:
    /**
     * \brief Constructor initializing all entries to the same value.
     */
    Tensor3(T const& init = T()) { for(auto& x: data) x = init; } // TODO: is this a performance problem (first default init, then assignment)?
    Dune::FieldMatrix<T,n2,n3> const& operator[] (int i) const { return data[i]; }
    Dune::FieldMatrix<T,n2,n3>      & operator[] (int i)       { return data[i]; }
    
    /**
     * \brief Addition update.
     */
    Tensor3<T,n1,n2,n3>& operator+=(Tensor3<T,n1,n2,n3> const& x)
    {
      for (int i=0; i<n1; ++i)
        data[i] += x[i];
      return *this;
    }
    
  private:
    std::array<Dune::FieldMatrix<T,n2,n3>,n1> data;
  };
  
  /**
   * \ingroup linalgbasic
   * \brief Multiplication of a tensor by a scalar.
   */
  template <class T, int n1, int n2, int n3>
  Tensor3<T,n1,n2,n3> operator*(T a, Tensor3<T,n1,n2,n3> const& x)
  {
    Tensor3<T,n1,n2,n3> ax;
    for (int i=0; i<n1; ++i)
      ax[i] = a*x[i];
    return ax;
  }

}

//---------------------------------------------------------------------

/**
 * \cond internals
 *
 * This namespace contains helper functions for implementing Hierarchic mappers.
 */
namespace FixDuneDetail {

  //TODO: in Dune 2.0 this is not needed anymore!
  template <int Codim>
  struct GetIndexOfSubEntity
  {
    template <class Cell, class IS>
    static int value(Cell const& cell, int codim, int entity, IS const& is)
    {
      if (codim==Codim)
        //       return is.template subIndex<Codim>(cell,entity);
        return is.subIndex(cell,entity,Codim);
      else
        return GetIndexOfSubEntity<Codim-1>::value(cell,codim,entity,is);
    }
  };

  template <>
  struct GetIndexOfSubEntity<0>
  {
    template <class Cell, class IS>
    static int value(Cell const& cell, int codim, int entity, IS const& is)
    {
      assert (codim==0);
      //     return is.template subIndex<0>(cell,entity);
      return is.subIndex(cell,entity,0);
    }
  };

  template <int Codim>
  struct GetNumberOfSubEntities
  {
    template <class Cell>
    static int value( Cell const& cell, int codim )
    {
      if (codim==Codim)
        return cell.template count<Codim>() ;
      else
        return GetNumberOfSubEntities<Codim-1>::value(cell, codim) ;
    }
  } ;

  template <>
  struct GetNumberOfSubEntities<0>
  {
    template <class Cell>
    static int value( Cell const& cell, int codim )
    {
      assert (codim==0);
      return cell.template count<0>();
    }
  };

} // End of namespace FixDuneDetail
/// \endcond ---------------------------------------------------------------------

/**
 * Computes the index of the k-th subentity of codimension c in the
 * index set is.  This is available in Dune only with c given
 * statically.
 *
 * TODO: use IndexSet::IndexType instead of size_t for Dune 1.1
 */
template <class IndexSet, class Cell>
size_t subIndex(IndexSet const& is, Cell const& cell, int codim, int subentity) 
{
  return FixDuneDetail::GetIndexOfSubEntity<Cell::dimension>::value(cell,codim,subentity,is);
}

/**
 *   \brief Computes the number of subentities of codimension c of a given entity.
 *
 *  This is available in Dune only with c given statically.
 *
 */
template <class Cell>
int count(Cell const& cell, int codim)
{
  return FixDuneDetail::GetNumberOfSubEntities<Cell::dimension>::value(cell, codim) ;
}


//  ---------------------------------------------------------------------

/**
 * \ingroup utilities
 * \brief Checks whether a point is inside or outside the reference element, and how much
 *
 * Computational geometry predicates are notoriously difficult to get right. In Dune, it may happen that
 * a point inside a cell is reported not to be contained in any of its children. Often this is due to
 * rounding errors. This function returns not just a binary decision subject to rounding errors, but
 * a floating point value that gives the magnitude of being inside or outside, such that the result can
 * be tested by arbitrary tolerances or sorted according to magnitude.
 *
 * \tparam LocalCoordinate the type of local coordinate vectors. Usually a Dune::FieldVector<ctype,dim>.
 *
 * \todo The current implementation is just for cubes and simplices. Pyramids and prisms are not covered.
 *
 * \return a value less or equal to zero if the point is inside, a value greater zero if outside
 */
template <class LocalCoordinate>
typename LocalCoordinate::field_type checkInside(Dune::GeometryType const& gt, LocalCoordinate const& x) {
  if (gt.isSimplex()) {
      // The reference simplex is the set defined by componentwise x_i>=0 and 1-sum x_i>=0. Thus we return
      // the negative of the minimum of these conditions.
      return -std::min(*std::min_element(x.begin(),x.end()),1-x.one_norm());
    } else if (gt.isCube()) {
      // The reference cube is the set defined by componentwise xi>=0 and xi_<=1.
      return -std::min(*std::min_element(x.begin(),x.end()),1-*std::max_element(x.begin(),x.end()));
    } else {
      // Other elements are not (yet) covered. Fallback to standard Dune implementation
      return Dune::ReferenceElements<typename LocalCoordinate::field_type,LocalCoordinate::dimension>::general(gt).checkInside(x)? 0: 1;
    }
}

/// copy between two block vectors of different element lengths
template <class Scalar, int n, int m, class Allocator1, class Allocator2>
void transferBlockVector(Dune::BlockVector<Dune::FieldVector<Scalar,n>,Allocator1> const& from, Dune::BlockVector<Dune::FieldVector<Scalar,m>,Allocator2>& to)
{
  to.resize(from.dim()*m);
  for(size_t i=0; i<from.N(); ++i)
    for(int j=0; j<n; ++j)
      to[(i*n+j)/m][(i*n+j)%m] = from[i][j];
}

///
template <class Scalar, int n, int m, class Allocator>
void bcrsPrint(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> const& A, std::ostream& os=std::cout)
{
  for(size_t i0=0; i0<A.N(); ++i0)
    for(size_t i1=0; i1<n; ++i1)
      {
        for(size_t j0=0; j0<A.M(); ++j0)
          for(size_t j1=0; j1<m; ++j1)
            if(A.exists(i0,j0)) os << A[i0][j0][i1][j1] << " ";
        os << std::endl;
      }
  os << std::endl;
}

template <class Scalar, int n, int m, class Allocator>
std::ostream& operator<<(std::ostream& os, Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> const& A)
{
  bcrsPrint(A,os);
  return os;
}

#endif
