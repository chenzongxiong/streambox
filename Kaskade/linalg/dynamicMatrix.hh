/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DYNAMICMATRIX_HH
#define DYNAMICMATRIX_HH

#include <tuple>
#include <vector>

#include "dune/common/dynvector.hh"
#include "dune/common/densematrix.hh"
#include "dune/common/fmatrix.hh"

/// \internal 
namespace Kaskade {
  template<class K> class DynamicMatrix;
  
  namespace DynamicMatrixDetail {
    
    // The row vector class accessing DynamicMatrix rows. It simply stores a pointer to the
    // first element, the number (of columns in the matrix) and the stride (number of rows, 
    // as the matrix stores elements column-major as in BLAS/LAPACK). It does not own or manage
    // the memory.
    // The interface adheres to the Dune::DenseVector requirements.
    template <class K>
    class StridedVector: public Dune::DenseVector<StridedVector<K>>
    {
      typedef Dune::DenseVector<StridedVector<K>> Base;
      
    public:
      typedef typename Base::size_type size_type;
      typedef typename Base::value_type value_type;
      
      StridedVector(StridedVector const& x): Dune::DenseVector<StridedVector<K>>(), p(x.p), n(x.n), stride(x.stride) {}
      
      /**
       * \brief Constructor.
       * \param p pointer to the first vector element
       * \param n number of elements in the vector
       * \param stride the offset between vector elements (a contiguously stored one has stride 1)
       */
      StridedVector(K* p_, size_type n_, size_type stride_): p(p_), n(n_), stride(stride_) {  }
      
      StridedVector& operator=(StridedVector const& v) = default;
      StridedVector& operator=(K const& x) { for (size_type i=0; i<n; ++i) vec_access(i) = x; return *this; }
      
      
      value_type&       vec_access(size_type i)       { return *(p+i*stride); }
      value_type const& vec_access(size_type i) const { return *(p+i*stride); }
      size_type         vec_size() const { return n; }
      
    private:
      K* p;
      size_type n;
      size_type stride;
    };
    
    // Helper class for extracting a pointer to the first entry (or the entry itself in case of elementary types
    template <class K>
    struct GetAddress
    {
      static K      * from(K      & k) { return &k; }
      static K const* from(K const& k) { return &k; }
    };
    
    template <class K, int n, int m>
    struct GetAddress<Dune::FieldMatrix<K,n,m>>
    {
      static K*       from(Dune::FieldMatrix<K,n,m>&       k) { return &k[0][0]; }
      static K const* from(Dune::FieldMatrix<K,n,m> const& k) { return &k[0][0]; }
    };
    
    template <class K, int n>
    struct GetAddress<Dune::FieldVector<K,n>>
    {
      static K*       from(Dune::FieldVector<K,n>&       k) { return &k[0]; }
      static K const* from(Dune::FieldVector<K,n> const& k) { return &k[0]; }
    };
    
    // Helper class for assigning Dune::FieldMatrix of different scalar type (which leads to segmentation faults as of Dune 2.4)
    template <class A, class B> 
    struct Copy
    {
      static void apply(A const& a, B& b)
      {
        b = a;
        assert(!std::isnan(b));
      }
    };
    
    template <class A, class B, int n, int m>
    struct Copy<Dune::FieldMatrix<A,n,m>, Dune::FieldMatrix<B,n,m>>
    {
      static void apply(Dune::FieldMatrix<A,n,m> const& a, Dune::FieldMatrix<B,n,m>& b)
      {
        for (int i=0; i<n; ++i)
          for (int j=0; j<n; ++j)
            Copy<A,B>::apply(a[i][j],b[i][j]);
      }
    };
    
    // Check wether the scalar entries are arranged directly in Lapack layout, such that
    // BLAS and LAPACK routines can work on them.
    template <class K>
    struct LapackLayout
    {
      static bool const value = std::is_floating_point<K>::value;
      static int  const rows  = 1;                                     // scalar rows
      static int  const cols  = 1;                                     // scalar columns
    };
    
    template <class K>
    struct LapackLayout<std::complex<K>>
    {
      static bool const value = LapackLayout<K>::value;
      static int  const rows  = LapackLayout<K>::rows; 
      static int  const cols  = LapackLayout<K>::cols; 
    };
    
    template <class K, int n, int m>
    struct LapackLayout<Dune::FieldMatrix<K,n,m>>
    {
      static bool const value = m==1 && LapackLayout<K>::value;
      static int  const rows  = n*LapackLayout<K>::rows;
      static int  const cols  = m*LapackLayout<K>::cols;  
   };
    
    template <class K, int n>
    struct LapackLayout<Dune::FieldVector<K,n>>
    {
      static bool const value = LapackLayout<K>::value;
      static int  const rows  = n*LapackLayout<K>::rows;
      static int  const cols  = LapackLayout<K>::cols;  
    };
    
    // matrix-vector multiplication calling BLAS
    void mv(int n, int m, double const* A, int lda, double const* x, double* y);
    void mv(int n, int m, float  const* A, int lda, float  const* x, float * y);
  }
}

namespace Dune
{

  template <class K>
  struct DenseMatVecTraits<Kaskade::DynamicMatrix<K>>
  {
    typedef Kaskade::DynamicMatrix<K> derived_type;
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K> row_type;
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K const> const_row_type;

    typedef row_type       row_reference;
    typedef const_row_type const_row_reference;

    typedef std::vector<K>                      container_type;
    typedef K                                   value_type;
    typedef typename container_type::size_type  size_type;
  };

  template <class K>
  struct DenseMatVecTraits<Kaskade::DynamicMatrixDetail::StridedVector<K>>
  {
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K> derived_type;
    typedef derived_type row_type;

    typedef K                                                                    value_type;
    typedef typename std::vector<typename std::remove_const<K>::type>::size_type size_type;
  };

  template< class K >
  struct FieldTraits<Kaskade::DynamicMatrix<K>>
  {
    typedef typename FieldTraits<K>::field_type field_type;
    typedef typename FieldTraits<K>::real_type  real_type;
  };
}
/// \endinternal

namespace Kaskade {
  /**
   * \ingroup linalgbasic
   * \brief A LAPACK-compatible dense matrix class with shape specified at runtime.
   * 
   * The memory layout is BLAS-like column major. For K = Dune::FieldMatrix<double,1,1>,
   * classical LAPACK and BLAS routines can be applied.
   * 
   * In contrast to Dune::DynamicMatrix this works with contiguous memory. 
   * In contrast to Dune::Matrix this does not reallocate on resize (unless the number of entries grows).
   * 
   * Note that the inherited "solve" method does not work for K=Dune::FieldMatrix<??,1,1>, probably due
   * to a Dune bug around common/densematrix.hh:808 (wrong argument type - pivots are matrix entries, not 
   * vector entries...).
   * 
   * \todo: Check use of DynamicMatrix<double> without FieldMatrix.
   */
  template<class K>
  class DynamicMatrix: public Dune::DenseMatrix<DynamicMatrix<K>>
  {
    typedef Dune::DenseMatrix<DynamicMatrix<K>> Base;
    using Self = DynamicMatrix<K>;
    typedef typename Base::row_reference row_reference;
    typedef typename Base::const_row_reference const_row_reference;
    typedef DynamicMatrixDetail::StridedVector<K> row_base;
    typedef DynamicMatrixDetail::StridedVector<K const> const_row_base;
    
  public:
    typedef typename Base::size_type  size_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::row_type   row_type;

    /**
     * \brief Creates a 0 x 0 matrix.
     */
    DynamicMatrix (): rs(0), cs(0), dat(1) {} // allocate at least one entry to prevent segmentation faults 

    /**
     * \brief Copy constructor.
     */
    DynamicMatrix(Self const& a) = default;
    
    template <class OtherK>
    friend class DynamicMatrix;
    
    /**
     * \brief Copies from a matrix with different (but compatible) entry type.
     */
    template <class OtherK>
    DynamicMatrix(DynamicMatrix<OtherK> const& A)
    : DynamicMatrix(A.N(),A.M())
    {
      copy(A.dat,dat);
    }

    /**
     * \brief Move constructor.
     */
    DynamicMatrix (Self&& a) = default;

    /**
     * \brief Creates a r x c matrix, initialized with given value.
     */
    DynamicMatrix (size_type r, size_type c, value_type v = value_type()) :
      rs(r), cs(c), dat(r*c,v)
    {}
    
    
    
    /**
     * \brief Copy assignment.
     */
    DynamicMatrix& operator=(DynamicMatrix const& a) = default;

    /**
     * \brief Move assignment.
     */
    DynamicMatrix& operator=(DynamicMatrix&& a) = default;
    
    /**
     * \brief Assignment from matrix with other scalar type.
     */
    template <class OtherK>
    DynamicMatrix& operator=(DynamicMatrix<OtherK> const& A)
    {
      resize(A.N(),A.M());
      copy(A.dat,dat);
      return *this;
    }

    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     * 
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::DynamicMatrix interface.
     */
    void resize (size_type r, size_type c)
    {
      rs = r;
      cs = c;
      dat.resize(r*c);
    }

    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     * 
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::Matrix interface.
     */
    void setSize (size_type r, size_type c)
    {
      resize(r,c);
    }

    void fill (K val)
    {
      for (int i=0;i<rs*cs;i++) dat[i]=val;
    }
    
    template<class X, class Y>
    void mv (const X& x, Y& y) const
    {
      using namespace DynamicMatrixDetail;
      using L = LapackLayout<K>;
      
      if (L::value)
        // If the scalar entries are arranged in LAPACK memory layout, we just call BLAS.
        DynamicMatrixDetail::mv(rs*L::rows,cs*L::cols,data(),lda(),
                                GetAddress<typename X::value_type>::from(x[0]),GetAddress<typename Y::value_type>::from(y[0]));
      else
        // As of Dune 2.4, the DenseMatrix::mv implementation is buggy and prevents the use of 
        // non-scalar entry types. We provide a fixed version here. (2016-02) Apparently corrected
        // in later Dune version.
        for (size_type i=0; i<rs; ++i)
        {
          y[i] = 0; // this was initialized by value_type(0), which can be a matrix...
          for (size_type j=0; j<cs; ++j)
            y[i] += (*this)[i][j] * x[j];
        }
    }
    
    /**
     * \todo docme
     */
    using Base::operator=;
    
    /**
     * \brief The leading dimension (for calling BLAS/LAPACK).
     */
    size_type lda() const { return rs*DynamicMatrixDetail::LapackLayout<K>::rows; }
    
    /**
     * \brief Raw access to data (for calling BLAS/LAPACK).
     * \warning Unless the entry type is Dune::FieldMatrix<Scalar,n,1>, the order of matrix entries is NOT column-major.
     */
    auto data() { return DynamicMatrixDetail::GetAddress<K>::from(dat[0]); }

    /**
     * \brief Raw access to data (for calling BLAS/LAPACK).
     * \warning Unless the entry type is Dune::FieldMatrix<Scalar,n,1>, the order of matrix entries is NOT column-major.
     */
    auto data() const { return DynamicMatrixDetail::GetAddress<K>::from(dat[0]); }

    // make this thing a matrix in the sense of the Dune interface
    size_type mat_rows() const { return rs; }
    size_type mat_cols() const { return cs; }
    row_reference mat_access(size_type i) { return row_reference(&dat[i],cs,rs); }
    const_row_reference mat_access(size_type i) const { return const_row_reference(&dat[i],cs,rs); }
    
  private:
    size_type rs, cs;
    std::vector<K> dat;
    
    // The assignment Dune::FieldMatrix<float,3,3> = Dune::FieldMatrix<double,3,3> leads to a segmentation fault 
    // possibly after interminate recursion as of Dune 2.4. Thus we provide a fix. (2016-02)
    template <class OtherK>
    void copy(std::vector<OtherK> const& a, std::vector<K>& b)
    {
      for (size_t i=0; i<a.size(); ++i)
        DynamicMatrixDetail::Copy<OtherK,K>::apply(a[i],b[i]);
    }
  };
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the singular value decomposition \f$ A = U \Sigma V^T \f$.
   *
   * If \f$ A \f$ is a \f$ n \times m \f$ matrix, then the SVD gives orthonormal \f$ U, V \f$ of size 
   * \f$ n\times n \f$ and \f$ m\times m\f$, respectively, and a \f$ n\times m\f$ matrix \f$ \Sigma \f$ with
   * only nonnegative diagonal entries.
   * 
   * \return the tuple (U,sigma,V), where sigma is a std::vector of the \f$ \min(n,m) \f$ diagonal entries of \f$ \Sigma \f$.
   */
  template <class Scalar>
  std::tuple<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>,std::vector<Scalar>,DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>>
  svd(DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>> const& A);
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the solution of \f$ Ax = b \f$.
   *
   * \return the tuple (U,sigma,V), where sigma is a std::vector of the \f$ \min(n,m) \f$ diagonal entries of \f$ \Sigma \f$.
   */
  template <class Scalar>
  Dune::DynamicVector<Dune::FieldVector<Scalar,1>> gesv(DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>> const& A,
                                                        Dune::DynamicVector<Dune::FieldVector<Scalar,1>> const& b);
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the inverse \f$ A^{-1} \f$.
   * \param A the input matrix to be inverted, is overwritten with its inverse
   */
  template <class Scalar>
  void invert(DynamicMatrix<Scalar>& A);
}

#endif
