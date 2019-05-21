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

#ifndef CRSUTIL_HH
#define CRSUTIL_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include "fem/fixdune.hh"

//---------------------------------------------------------------------

namespace Kaskade
{
  template <class OutIter>
  OutIter vectorToSequence(double x, OutIter iter)
  {
    *iter = x;
    return ++iter;
  }

  template <class OutIter>
  OutIter vectorToSequence(float x, OutIter iter)
  {
    *iter = x;
    return ++iter;
  }

  
  template <class K, int size, class OutIter>
  OutIter vectorToSequence(Dune::FieldVector<K,size> const& v, OutIter iter)
  {
    for (size_t i=0; i<size; ++i) 
    {
      #ifndef NDEBUG
      if (std::isinf(v[i]))
      {
        std::cerr << __FILE__ << ':' << __LINE__ << ":VectorToSequence::call  v[" << i << "] = inf\n" << std::flush;
        throw -111;
      }
      if (std::isnan(v[i]))
      {
        std::cerr << __FILE__ << ':' << __LINE__ << ":VectorToSequence::call  v[" << i << "] = nan\n" << std::flush;
        throw -111;
      }
      #endif
      iter = vectorToSequence(v[i],iter);
    }
    return iter;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief copies the entries of a vector sequentially to the output iterator.
   *
   * The Vector type can be Dune::BlockVector<B,A>, Dune::FieldVector<K,size>, 
   * or a finite element function. The type of the scalar entries of the
   * vector must be convertible to the value type of the output
   * iterator.
   */
  template <class B, class A, class OutIter>
  OutIter vectorToSequence(Dune::BlockVector<B,A> const& v, OutIter iter)
  {
    for (int i=0; i<v.N(); ++i)
      iter = vectorToSequence(v[i],iter);
    return iter;
  }
  
  template <class K, class OutIter>
  OutIter vectorToSequence(std::vector<K> const& v, OutIter iter)
  {
    for (auto const& x: v)
      iter = vectorToSequence(x,iter);
    return iter;
  }
  

  //---------------------------------------------------------------------



  template <class InIter>
  InIter vectorFromSequence(double& x, InIter iter)
  {
    x = *iter;
    return ++iter;
  }

  template <class InIter>
  InIter vectorFromSequence(float& x, InIter iter)
  {
    x = *iter;
    return ++iter;
  }

  template <class K, int size, class InIter>
  InIter vectorFromSequence(Dune::FieldVector<K,size>& v, InIter iter)
  {
    for (size_t i=0; i<v.N(); ++i)
      iter = vectorFromSequence(v[i],iter);
    return iter;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief copies the elements obtained from the input iterator sequentially to the entries of the vector.
   *
   * The value type of the input iterator must be convertible to the scalar
   * entries of the vector. The input iterator is advanced by the number
   * of scalar entries in the vector.
   */
  template <class B, class A, class InIter>
  InIter vectorFromSequence(Dune::BlockVector<B,A>& v, InIter in)
  {
    for (int i=0; i<v.N(); ++i)
      in = vectorFromSequence(v[i],in);
    return in;
  }
  
  template <class K, class InIter>
  InIter vectorFromSequence(std::vector<K>& v, InIter iter)
  {
    for (auto& x: v)
      iter = vectorFromSequence(x,iter);
    return iter;
  }
  
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * @brief A class template that supports converting certain Dune
   * matrices into the coordinate (triplet) format.
   */
  template <class Matrix>
  struct Matrix_to_Triplet {};

  template <class K, int N, int M, class Allocator>
  struct Matrix_to_Triplet<Dune::BCRSMatrix<Dune::FieldMatrix<K,N,M>,Allocator> >
  {
  private:
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<K,N,M>,Allocator> Matrix;
    typedef typename Matrix::ConstRowIterator RI;
    typedef typename Matrix::ConstColIterator CI;

  public:
    /**
     * @brief converts a Dune matrix into triplet format
     *
     * The coordinates are written to the output iterators i (rows) and
     * j (columns), while the matrix entries are written to the output
     * iterator z. The offsets startrow and startcol are added to the
     * indices. This is useful if this matrix is just a block in a
     * larger blocked matrix, or if the output shall be suitable for
     * Matlab (base 1 indices).
     *
     * The output iterators are required to provide enough space to take
     * all the matrix elements. The required number can be obtained
     * beforehand from nnz().
     *
     * If onlyLowerTriangle is true, no superdiagonal entries of the
     * matrix are written. In this case, a symmetric row/column blocking
     * is assumed.
     *
     * If symmetric is true, no superdiagonal entries of the matrix are
     * accessed, but all subdiagonal entries are mirrored. The symmetric
     * flag is not transitive, i.e. it is assumed that all matrix
     * elements (even if they are matrices themselves) are stored
     * completely.
     */
    template <class OutIteratorI, class OutIteratorD>
    static void call(Matrix const& a,
		     OutIteratorI& i, OutIteratorI& j,
		     OutIteratorD& z,
		     int startrow, int const startcol,
		     bool onlyLowerTriangle = false,
		     bool symmetric = false)
    {
      for (RI row=a.begin(); row!=a.end(); ++row) 
        for (CI col=row->begin() ; col!=row->end(); ++col) 
	{
	  if (!onlyLowerTriangle || (row.index()>=col.index()))
	    Matrix_to_Triplet<typename Matrix::block_type>::call(*col,i,j,z,startrow+row.index()*N,startcol+col.index()*M,
								 (row.index()==col.index()) && onlyLowerTriangle);
	  if (symmetric && !onlyLowerTriangle && row.index()>col.index())
	    Matrix_to_Triplet<typename Matrix::block_type>::call(*col,j,i,z,startcol+row.index()*N,startrow+col.index()*M);
        }
    }

    /**
     * @brief returns the number of structurally nonzero scalars in the matrix
     *
     * If onlyLowerTriangle is true, no superdiagonal entries of the
     * matrix are counted. In this case, a symmetric row/column blocking
     * is assumed.
     *
     * If symmetric is true, no superdiagonal entries are counted (since
     * they are probably not stored), but the subdiagonal ones are
     * counted twice. The symmetric flag is not transitive, i.e. all
     * counted matrix elements (even if they are matrices themselves)
     * are fully counted.
     */
    static size_t nnz(Matrix const& a,
                      bool onlyLowerTriangle = false,
                      bool symmetric = false)
    {
      size_t size = 0;
      RI rend = a.end();
      for (RI row=a.begin(); row!=rend; ++row)
      {
        CI cend = row->end();
        for (CI col=row->begin(); col!=cend; ++col)
        {
          if (!onlyLowerTriangle || (row.index()>=col.index()))
            size += Matrix_to_Triplet<typename Matrix::block_type>::nnz(*col,(row.index()==col.index()) && onlyLowerTriangle);

          if (symmetric && !onlyLowerTriangle && row.index()>col.index())
            size += Matrix_to_Triplet<typename Matrix::block_type>::nnz(*col);
        }
      }
      return size;
    }
  };

  template <class Entry, class Index> // forward declaration
  class NumaBCRSMatrix;
  
  template <class K, int N, int M, class Index>
  struct Matrix_to_Triplet<NumaBCRSMatrix<Dune::FieldMatrix<K,N,M>,Index> >
  {
  private:
    typedef NumaBCRSMatrix<Dune::FieldMatrix<K,N,M>,Index> Matrix;

  public:
    /**
     * @brief converts a Dune matrix into triplet format
     *
     * The coordinates are written to the output iterators i (rows) and
     * j (columns), while the matrix entries are written to the output
     * iterator z. The offsets startrow and startcol are added to the
     * indices. This is useful if this matrix is just a block in a
     * larger blocked matrix, or if the output shall be suitable for
     * Matlab (base 1 indices).
     *
     * The output iterators are required to provide enough space to take
     * all the matrix elements. The required number can be obtained
     * beforehand from nnz().
     *
     * If onlyLowerTriangle is true, no superdiagonal entries of the
     * matrix are written. In this case, a symmetric row/column blocking
     * is assumed.
     *
     * If symmetric is true, no superdiagonal entries of the matrix are
     * accessed, but all subdiagonal entries are mirrored. The symmetric
     * flag is not transitive, i.e. it is assumed that all matrix
     * elements (even if they are matrices themselves) are stored
     * completely.
     */
    template <class OutIteratorI, class OutIteratorD>
    static void call(Matrix const& a,
                     OutIteratorI& i, OutIteratorI& j,
                     OutIteratorD& z,
                     int startrow, int const startcol,
                     bool onlyLowerTriangle = false,
                     bool symmetric = false)
    {
      auto rend = a.end();
      for (auto row=a.begin(); row!=rend; ++row) 
      {
	auto cend = row->end();
        for (auto col=row->begin() ; col!=cend; ++col) 
	{
          if (!onlyLowerTriangle || (row.index()>=col.index()))
            Matrix_to_Triplet<typename Matrix::block_type>::call(*col,i,j,z,startrow+row.index()*N,startcol+col.index()*M,
								 (row.index()==col.index()) && onlyLowerTriangle);
          if (symmetric && !onlyLowerTriangle && row.index()>col.index())
            Matrix_to_Triplet<typename Matrix::block_type>::call(*col,j,i,z,startcol+row.index()*N,startrow+col.index()*M);
        }
      }
    }

    /**
     * @brief returns the number of structurally nonzero scalars in the matrix
     *
     * If onlyLowerTriangle is true, no superdiagonal entries of the
     * matrix are counted. In this case, a symmetric row/column blocking
     * is assumed.
     *
     * If symmetric is true, no superdiagonal entries are counted (since
     * they are probably not stored), but the subdiagonal ones are
     * counted twice. The symmetric flag is not transitive, i.e. all
     * counted matrix elements (even if they are matrices themselves)
     * are fully counted.
     */
    static size_t nnz(Matrix const& a,
                      bool onlyLowerTriangle = false,
                      bool symmetric = false)
    {
      size_t size = 0;
      auto rend = a.end();
      for (auto row=a.begin(); row!=rend; ++row)
      {
        auto cend = row->end();
        for (auto col=row->begin(); col!=cend; ++col)
        {
          if (!onlyLowerTriangle || (row.index()>=col.index()))
            size += Matrix_to_Triplet<typename Matrix::block_type>::nnz(*col,(row.index()==col.index()) && onlyLowerTriangle);

          if (symmetric && !onlyLowerTriangle && row.index()>col.index())
            size += Matrix_to_Triplet<typename Matrix::block_type>::nnz(*col);
        }
      }
      return size;
    }
  };

  template <class K, int n, int m>
  struct Matrix_to_Triplet<Dune::FieldMatrix<K,n,m> >
  {
    template <class OutIteratorI, class OutIteratorD>
    static void call(Dune::FieldMatrix<K,n,m> const& a,
		     OutIteratorI& i, OutIteratorI& j,
		     OutIteratorD& z,
		     int startrow, int const startcol,
		     bool onlyLowerTriangle = false)
    {
      if (onlyLowerTriangle) 
      {
        assert(n==m);
        for (int row=0; row<n; ++row)
          for (int col=0; col<=row; ++col) {
            *i++ = row + startrow;
            *j++ = col + startcol;
            *z++ = a[row][col];
          }
      } 
      else
        for (int row=0; row<n; ++row)
          for (int col=0; col<m; ++col) {
            *i++ = row + startrow;
            *j++ = col + startcol;
            *z++ = a[row][col];
          }
    }

    template <class OutIteratorI, class OutIteratorD>
    static void callTransposed(Dune::FieldMatrix<K,n,m> const& a,
        OutIteratorI& i, OutIteratorI& j,
        OutIteratorD& z,
        int startrow, int const startcol)
    {
      for (int row=0; row<n; ++row)
        for (int col=0; col<m; ++col) {
          *i++ = col + startrow;
          *j++ = row + startcol;
          *z++ = a[row][col];
        }
    }

    static size_t nnz(Dune::FieldMatrix<K,n,m> const& a,
        bool onlyLowerTriangle = false)
    {
      if (onlyLowerTriangle) {
        assert(n==m);
        return n*(n+1)/2;
      } else return n*m;
    }
  };


  /**
   * @brief converts a matrix to the coordinate (triplet) format.
   *
   * The Matrix type needs to satisfy the Dune matrix interface.
   */
  template <class Matrix, class OutIteratorI, class OutIteratorD>
  void matrix_to_triplet(Matrix const& a,
      OutIteratorI i, OutIteratorI j,
      OutIteratorD z)
  {
    Matrix_to_Triplet<Matrix>::call(a,i,j,z,0,0);
  }

//   /**
//    * \brief Sorts the sparse matrix nonzero positions into row buckets.
//    * \param nrows number of rows
//    * \param ridx vector of row indices
//    * \param cidx vector of column indices
//    * \param[out] rowColumnIndices stores the column indices for each row
//    */
//   template <class SparseInt>
//   void getBCRSIndicesFromTriplet(SparseInt nrows, std::vector<SparseInt> const& ridx, std::vector<SparseInt> const& cidx, 
//                                  std::vector<std::vector<SparseInt>>& rowColumnIndices)
//   {
//     rowColumnIndices.clear();
//     rowColumnIndices.resize(nrows);
//     for(SparseInt i=0; i<ridx.size(); ++i) 
//       rowColumnIndices[ridx[i]].push_back(cidx[i]);
//   }
  
  /**
   * \brief Sorts the sparse matrix nonzero positions into row buckets.
   * 
   * From the row and column indices of a sparse matrix with scalar entries, e.g., in triplet format,
   * obtain the row and column indices of a sparse matrix with nr-by-nc block-entries, as in
   * Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,nr,nc>>. The indices are stored by row.
   * 
   * \param nrows number of rows of BCRS target matrix (i.e. rows of triplet source / nr)
   * \param ridx vector of row indices
   * \param cidx vector of column indices
   * \param nr row size of sparse matrix entries
   * \param nc column size of sparse matrix entries
   */
  template <class SparseInt>
  std::vector<std::vector<SparseInt>> getBCRSIndicesFromTriplet(SparseInt nrows, std::vector<SparseInt> const& ridx, std::vector<SparseInt> const& cidx,
                                                                int nr=1, int nc=1)
  {
    std::vector<std::vector<SparseInt>> rowColumnIndices(nrows);
    for(SparseInt i=0; i<ridx.size(); ++i)                        // sort entries into rows
      rowColumnIndices[ridx[i]/nr].push_back(cidx[i]/nc);         // respect nr-by-nc blocking
    for (auto& r: rowColumnIndices)                               // sort column entries in ascending order
    {                                                             // and remove multiple entries
      std::sort(begin(r),end(r));                                 // (e.g. created by nr,nc > 1
      r.erase(std::unique(begin(r),end(r)),end(r));
    }
    return rowColumnIndices;
  }
  
  //---------------------------------------------------------------------


  template <class B, class A>
  std::ostream& operator << (std::ostream& out, Dune::BCRSMatrix<B,A> const& a)
  {
    typedef typename Dune::BCRSMatrix<B,A>::ConstRowIterator RI;
    typedef typename Dune::BCRSMatrix<B,A>::ConstColIterator CI;

    int const prec = 6;
    int const width = prec + 3;

    out << "[\n";
    for (RI row=a.begin(); row!=a.end(); ++row) {
      int oldidx = -1;
      for (CI col=row->begin() ; col!=row->end(); ++col) {
        out.width((col.index()-oldidx-1)*(width+1)); out << "";
        out.setf(std::ios_base::fixed,std::ios_base::floatfield);
        out.width(width); out.precision(prec); out << *col << ' ';
        oldidx = col.index();
      }
      out << '\n';
    }
    out << "\n]";

    return out;
  }

  //---------------------------------------------------------------------

} // namespace Kaskade
#endif
