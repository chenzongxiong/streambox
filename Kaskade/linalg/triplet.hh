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

#ifndef TRIPLET_HH
#define TRIPLET_HH

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <utility>

#include "dune/common/fmatrix.hh"
#include "dune/istl/bcrsmatrix.hh"

#include "linalg/crsutil.hh"
#include "linalg/umfpack_tools.hh"

namespace Kaskade
{
  /// \internal
  // forward declaration
  template <class Entry, class Index>
  class NumaBCRSMatrix;
  /// \endinternal

  /**
   * \ingroup linalgbasic
   * \brief Sparse Matrix in triplet format.
   *
   * A simple sparse matrix format consisting of scalar-valued (row,col,value)
   * entries. Use this as the least common denominator of sparse linear algebra
   * interfaces.
   * 
   * Multiple entries with same row and col are semantically
   * summed up.
   * \tparam Scalar (double or float, field_type)
   * \tparam SparseIndexInt (int,long)
   * 
   * long has to be used as SparseIndexInt if the number of nonzeros of the
   * sparse matrix exceeds INT_MAX 2,147,483,647 or if the umfpack  direct solver is used
   * and the internal workspace exceeds this number, what will happen much earlier.
   */
  template<class Scalar,class SparseIndexInt=int>
  class MatrixAsTriplet
  {
  public:
    /// STL-compliant typedefs
    typedef Scalar value_type;
    typedef typename std::vector<Scalar>::iterator iterator;
    typedef typename std::vector<Scalar>::const_iterator const_iterator;
    typedef typename std::vector<SparseIndexInt>::iterator index_iterator;
    typedef typename std::vector<SparseIndexInt>::const_iterator const_index_iterator;

    /// Constructor
    MatrixAsTriplet() {}

    /// Constructor that allocates memory directly
    explicit MatrixAsTriplet(size_t s) : ridx(s), cidx(s), data(s){}

    /// Constructor from a Dune::BCRSMatrix
    template <int n, int m, class Allocator>
    explicit MatrixAsTriplet(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> const& other)
    {
      fillFromBCRSInterface<n,m>(other);
    }

    /// Constructor from a NumaBCRSMatrix
    template <int n, int m, class Index>
    explicit MatrixAsTriplet(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Index> const& other)
    {
      fillFromBCRSInterface<n,m>(other);
    }

    /**
     * \brief Move constructor taking raw data.
     *
     * The size of row, col, and value
     * has to be identical. All entries in row and col need to be
     * nonnegative. The supplied containers are empty after
     * construction.
     */
    MatrixAsTriplet(std::vector<SparseIndexInt>&& row, std::vector<SparseIndexInt>&& col,
                    std::vector<Scalar>&& value)
    {
      assert(row.size()==col.size() && row.size()==value.size());
      assert(*std::min_element(row.begin(),row.end()) >=0);
      assert(*std::min_element(col.begin(),col.end()) >=0);

      ridx = std::move(row);
      cidx = std::move(col);
      data = std::move(value);
    }

    /// Constructs a matrix of size (nrows,ncols), filled with value. If value==0.0, then it is an empty matrix
    MatrixAsTriplet(SparseIndexInt nrows, SparseIndexInt ncols, Scalar const& value)
    {
      if(value != 0.0)
      {
        reserve(nrows*ncols);
        resize(0);
        for(size_t i=0; i<nrows; ++i)
          for(size_t k=0; k<ncols; ++k)
          {
            ridx.push_back(i);
            cidx.push_back(k);
            data.push_back(value);
          }
      }
    }

    /// Constructs a diagonal matrix of size (nrows,ncols). The diagonal is shifted down shift rows (if shift is negative, it is shifted up)
    MatrixAsTriplet(SparseIndexInt nrows, SparseIndexInt ncols, SparseIndexInt shift, Scalar const& value)
    {
      if(value != 0.0)
      {
        resize(0);
        for(size_t i=0; i<ncols; ++i)
        {
          size_t rowidx=i+shift;
          if(rowidx >= 0 && rowidx < nrows)
          {
            ridx.push_back(rowidx);
            cidx.push_back(i);
            data.push_back(value);
          }
        }
      }
    }

    /// Constructs a diagonal matrix of size (nrows,ncols). The diagonal is shifted down shift rows (if shift is negative, it is shifted up)
    MatrixAsTriplet(SparseIndexInt nrows, SparseIndexInt ncols, SparseIndexInt shift, std::vector<Scalar> const& value)
    {
      resize(0);
      for(int i=0; i<ncols; ++i)
      {
        int rowidx=i+shift;
        if(rowidx >= 0 && rowidx < nrows)
        {
          ridx.push_back(rowidx);
          cidx.push_back(i);
          data.push_back(value[i]);
        }
      }
    }

    /// Copy constructor
    MatrixAsTriplet(MatrixAsTriplet const& other) = default;

    /// Move constructor
    MatrixAsTriplet(MatrixAsTriplet&& other) = default;

    /// deprecated, use clear instead, will be deleted after 2014-08-31
    void flush() __attribute__((deprecated))
    {
      clear();
      std::cerr << "MatrixAsTriplet::flush() is deprecated, use clear() instead.\n";
    }

    void clear()
    {
      cidx.clear();
      ridx.clear();
      data.clear();
    }


    /**
     * \brief Shift the row and column indices such that \f$ A_{00} \f$ is (structurally) nonzero.
     */
    void setStartToZero()
    {
      SparseIndexInt r0 = *std::min_element(ridx.begin(),ridx.end());
      SparseIndexInt c0 = *std::min_element(cidx.begin(),cidx.end());
      if(r0!=0 || c0!=0)
        shiftIndices(-r0,-c0);
    }

    /// Returns number of rows (computes them, if not known)
    SparseIndexInt nrows() const
    {
      if(ridx.size()==0) return 0;
      return *std::max_element(ridx.begin(),ridx.end())+1;
    }

    /**
     * \brief Returns number of rows (computes them, if not known)
     * This is a Dune ISTL compatibility method.
     */
    SparseIndexInt N() const { return nrows(); }

    /// Returns number of cols (computes them, if not known)
    SparseIndexInt ncols() const
    {
      if(cidx.size()==0) return 0;
      return *std::max_element(cidx.begin(),cidx.end())+1;
    }

    /**
     * \brief Returns number of columns (computes them, if not known)
     * This is a Dune ISTL compatibility method.
     */
    SparseIndexInt M() const { return ncols(); }

    /**
     * \brief Resizes the memory.
     *
     * Existing values are preserved if s >= nnz(),
     * and the new entries are undefined. For s<nnz(), all the values
     * are undefined.
     */
    void resize(size_t s) 
    {
      ridx.resize(s); 
      cidx.resize(s); 
      data.resize(s); 
    }

    /// Reserve memory
    void reserve(size_t s) 
    {
      ridx.reserve(s); 
      cidx.reserve(s); 
      data.reserve(s); 
    }

    /**
     * Number of non-zeros. This is the number of entries in the triplet
     * format and can be larger than the number of nonzero entries in
     * the represented sparse matrix.
     */
    size_t nnz() const { return ridx.size(); }

    /// Matrix-vector multiplication: out += alpha * (*this) * in
    /**
     * ?????: Was ist denn nr? number of blocks/rows??
     */
    template <class X, class Y>
    void axpy(Y& out, X const& in, Scalar alpha=1.0,int nr=1) const
    {
      assert(ridx.size()==cidx.size() && ridx.size()==data.size());
      assert(in.size() >= ncols()*nr);
      assert(out.size() >= nrows()*nr);
      SparseIndexInt cols=ncols();
      SparseIndexInt rows=nrows();
      for(int j=0; j<nr; ++j)
        for(size_t i=0; i<ridx.size(); ++i)
        {
          assert(ridx[i]+j*rows < out.size());
          assert(cidx[i]+j*cols < in.size());
          out[ridx[i]+j*rows] += alpha*in[cidx[i]+j*cols]*data[i];
        }
    }

    /// Matrix-vector multiplication
    /**
     * \f$ x = x + \alpha A^{T}y \f$
     */
    template <class X, class Y>
    void atxpy(Y& y, X const& x, Scalar alpha=1.0) const
    {
      assert(ridx.size()==cidx.size() && ridx.size()==data.size());
      assert(y.size() == ncols());
      assert(x.size() == nrows());
      for(size_t i=0; i<data.size(); ++i) y[cidx[i]] += alpha*data[i]*x[ridx[i]];
    }

    /// Matrix-vector multiplication
    /**
     * \f$ x = x + \alpha A^{T}y \f$
     */
    template <class X, class Y>
    void usmtv(Scalar const alpha, X const& x, Y& y) const
    {
      atxpy(y,x,alpha);
    }


    /// scaled matrix-vector multiplication
    template <class X, class Y>
    void usmv(Scalar const alpha, X const& x, Y& y) const
    {
      axpy(y, x, alpha);
    }


    /// Matrix-vector multiplication: out = (*this) * in
    template <class X, class Y>
    void ax(Y& out, X const& in) const
    {
      assert(ridx.size()==cidx.size() && ridx.size()==data.size());
      assert(in.size() >= ncols());
      out.resize(0);
      out.resize(nrows(),0.0);
      for(size_t i=0; i<ridx.size(); ++i)
      {
        assert(ridx[i] < out.size());
        assert(cidx[i] < in.size());
        out[ridx[i]] += in[cidx[i]]*data[i];
      }
    }

    /// matrix-vector multiplication
    template <class X, class Y>
    void mv(X const& x, Y& y) const
    {
      axpy(y,x);
    }

    /// Shifts matrix indices. Can be used together with += to concatenate sparese matrices
    void shiftIndices(SparseIndexInt row, SparseIndexInt col)
    {
      assert(ridx.size()==cidx.size() && ridx.size()==data.size());
      for(size_t i=0; i< ridx.size(); ++i)
      {
        ridx[i] += row;
        cidx[i] += col;
        assert(ridx[i] >= 0 && cidx[i] >= 0);
      }
    }

    /**
     * \brief Adds the given value to the specified matrix entry.
     * The entry is created in case it has been structurally zero before.
     */
    void addEntry(SparseIndexInt row, SparseIndexInt col, Scalar const& value)
    {
      ridx.push_back(row);
      cidx.push_back(col);
      data.push_back(value);
    }

    /// Matrix addition: (*this) += m  (works also for matrices with non-matching sparsity pattern)
    MatrixAsTriplet& operator+=(MatrixAsTriplet const& m)
    {
      if(ridx.size()==0)
      {
        ridx=m.ridx;
        cidx=m.cidx;
        data=m.data;
        return *this;
      }
      ridx.insert(ridx.end(),m.ridx.begin(),m.ridx.end());
      cidx.insert(cidx.end(),m.cidx.begin(),m.cidx.end());
      data.insert(data.end(),m.data.begin(),m.data.end());
      
      // conversion to compresses row and back to triplet to remove multiple entries
      
      std::vector<SparseIndexInt> Ap, Ai;
      std::vector<Scalar> Az;
      umfpack_triplet_to_col(nrows(), ncols(), ridx, cidx, data, Ap, Ai, Az);
      {
        std::vector<SparseIndexInt> nullInt1, nullInt2;
        std::vector<Scalar> nullDouble;
        cidx.swap(nullInt1);
        ridx.swap(nullInt2);
        data.swap(nullDouble);
      }
      umfpack_col_to_triplet(Ap, cidx);
      data.insert(data.begin(),&Az[0],&Az[Ap.back()]);
      ridx.insert(ridx.begin(),&Ai[0],&Ai[Ap.back()]);
      
      return *this;
    }

    /// Scaling
    MatrixAsTriplet& operator*=(Scalar const& scalar)
    {
      for(auto& v: data)
        v *= scalar;

      return *this;
    }

    /// Scaling
    MatrixAsTriplet& operator/=(Scalar const& scalar)
    {
      return (*this) *= 1/scalar;
    }

    /// Assignment
    MatrixAsTriplet& operator=(MatrixAsTriplet const& m) = default;

    /// Move assignment
    MatrixAsTriplet& operator=(MatrixAsTriplet&& other) = default;

    void scaleRows(std::vector<Scalar>const& scaling)
    {
      for(size_t i=0; i<ridx.size();++i)
        data[i] *= scaling[ridx[i]];
    }

    /// Transposition
    MatrixAsTriplet& transpose()
    {
      // simple: just exchange row and column indices.
      std::swap(ridx,cidx);
      return *this;
    }

    /// Deletes all sub-diagonal entries.
    void deleteLowerTriangle()
    {
      size_t count(0);
      for(size_t i=0; i<ridx.size(); ++i)
      {
        if(ridx[i] <= cidx[i]) count++;
      }
      std::vector<SparseIndexInt> r2,c2;
      std::vector<Scalar> d2;
      r2.reserve(count);
      c2.reserve(count);
      d2.reserve(count);
      for(size_t i=0; i<ridx.size(); ++i)
      {
        if(ridx[i] <= cidx[i])
        {
          r2.push_back(ridx[i]);
          c2.push_back(cidx[i]);
          d2.push_back(data[i]);
        }
      }
      std::swap(r2,ridx);
      std::swap(c2,cidx);
      std::swap(d2,data);
    }

    /** 
     * \brief Output as Matlab source code
     */
    std::ostream& print(std::ostream& out = std::cout, double threshold=0.0) const
    {
      out << "[" << std::endl;
      for(size_t i=0; i<ridx.size(); ++i) if(std::fabs(data[i]) > threshold)
        out << ridx[i] << "," << cidx[i] << "," << data[i] << ";" << std::endl;
      out << "]" << std::endl;
      return out;
    }

    /// transform into a set of column vectors
    void toColumns(std::vector<std::vector<Scalar> >& colsvector) const
    {
      SparseIndexInt rows=nrows();
      colsvector.resize(ncols());
      for(size_t i=0; i<colsvector.size();++i)
      {
        colsvector[i].resize(0);
        colsvector[i].resize(rows,0.0);
      }
      for(size_t i=0; i<data.size(); ++i)
        colsvector[cidx[i]][ridx[i]]=data[i];
    }

    void toRows(std::vector<std::vector<Scalar> >& rows) const
    {
      SparseIndexInt cols = ncols();
      rows.resize(nrows());
      std::for_each(rows.begin(),rows.end(),
          [&cols](std::vector<Scalar>& row)
          {
            row.resize(0);
            row.resize(cols,0.0);
          });
      for(size_t i=0; i<data.size(); ++i) rows[ridx[i]][cidx[i]] = data[i];
    }

    void toVector(std::vector<Scalar>& colsvector)
    {
      SparseIndexInt rows=nrows();
      colsvector.resize(0);
      colsvector.resize(nrows()*ncols(),0.0);
      for(size_t i=0; i<data.size(); ++i)
        colsvector[cidx[i]*rows+ridx[i]]=data[i];
    }
    /// add a column

    void addColumn(std::vector<Scalar>& colsvector, size_t position)
    {
      for(size_t i=0; i<colsvector.size();++i)
      {
        ridx.push_back(i);
        cidx.push_back(position);
        data.push_back(colsvector[i]);
      }
    }


    /// add to a dense matrix
    template<class Mat>
    void addToMatrix(Mat& mat)
    {
      for(size_t i=0; i<data.size(); ++i)
        mat[ridx[i]][cidx[i]] += data[i];
    }

    void erase(SparseIndexInt i)
    {
      ridx.erase(std::begin(ridx) + i);
      cidx.erase(std::begin(cidx) + i);
      data.erase(std::begin(data) + i);
    }

    MatrixAsTriplet lumped()
    {
      SparseIndexInt noRows = nrows(), nnz = data.size();
      std::vector<Scalar> diagonalEntries(noRows,0);

      // sum up row entries
      for(SparseIndexInt row=0; row<noRows; ++row)
      {
        SparseIndexInt i=0;
        while(i<data.size())
        {
          if(ridx[i] == row) diagonalEntries[row] += data[i];
          else ++i;
        }
      }

      assert(nrows()==0 && ncols()==0); // matrix should be empty now

      MatrixAsTriplet<Scalar,SparseIndexInt> result;
      for(SparseIndexInt i=0; i<noRows; ++i) result.addEntry(i,i,diagonalEntries[i]);
      return result;
    }

    /// get iterator on data field
    iterator begin(){ return data.begin(); }

    /// get const iterator on data field
    const_iterator begin() const { return data.begin(); }

    /// get iterator on end of data field
    iterator end() { return data.end(); }

    /// get const iterator on end of data field
    const_iterator end() const { return data.end(); }

    /// get number of entries in sparse matrix
    size_t size() const { return data.size(); }

    /// write to text file
    void toFile(const char* filename, unsigned int precision=10) const
    {
      std::ofstream myfile;
      myfile.precision(precision);
      myfile.open(filename);
      for(size_t i=0; i<data.size(); ++i)
        myfile << ridx[i] << " " << cidx[i] << " " << data[i] << std::endl;
      myfile.close();
    }  
    
    /// copy to Dune::BCRSMatrix
    template <int bcrsN=1>
    std::unique_ptr<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,bcrsN,bcrsN> > > toBCRS() const
    {
//      std::cout << "triplet: " << std::endl;
//      print(std::cout);
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,bcrsN,bcrsN> > BCRSMat;
      
      // get column indices for each row
      SparseIndexInt noRows = nrows()/bcrsN,
                     noCols = ncols()/bcrsN;
      std::vector<std::vector<SparseIndexInt> > bcrsids(noRows);
      getBCRSIndices<bcrsN>(noRows, ridx, cidx, bcrsids);
//      std::cout << "brcsids: " << std::endl;
//      for(size_t i=0; i<bcrsids.size(); ++i)
//        for(size_t j=0; j<bcrsids[i].size(); ++j)
//          std::cout << i << "," << j << ": " << bcrsids[i][j] << std::endl;

      SparseIndexInt nonz = 0;
      for (auto const& s : bcrsids) nonz += s.size();
//      std::cout << "nRows: " << noRows << ", nCols: " << noCols << ", nnz: " << nonz << "/" << nnz() << std::endl;
      std::unique_ptr<BCRSMat> result(new BCRSMat(noRows,noCols,nonz,BCRSMat::row_wise));
//      std::cout << "initialized bcrs matrix" << std::endl;
      // create sparsity pattern
      for (typename BCRSMat::CreateIterator row=result->createbegin(); row!=result->createend(); ++row)
        for(SparseIndexInt col = 0; col < bcrsids[row.index()].size(); ++col) row.insert(bcrsids[row.index()][col]);

      // fill matrix
      for (SparseIndexInt i=0; i<data.size(); ++i)
      {
        size_t blockRow = ridx[i]/bcrsN, blockCol = cidx[i]/bcrsN;
//        std::cout << "block: " << blockRow << ", " << blockCol << std::endl;
//        std::cout << "in block: " << ridx[i] - bcrsN*blockRow << ", " << cidx[i] - bcrsN*blockCol << std::endl;
        (*result)[blockRow][blockCol][ridx[i]-bcrsN*blockRow][cidx[i]-bcrsN*blockCol] = data[i];
      }
//      std::cout << "bcrs" << std::endl;
//      for (int k=0; k<result->N(); ++k)
//        for (typename BCRSMat::ConstColIterator ca=(*result)[k].begin(); ca!=(*result)[k].end(); ++ca)
//          std::cout << k << ", " << ca.index() << ": " << *ca << std::endl;


      return result;
    }
    
    bool isSymmetric() const
    {
      for(size_t i=0; i<data.size(); ++i)
      {
        // check if symmetric entry exists
        auto tmp = findEntry(cidx[i],ridx[i]);
        if(tmp.second == false) return false;
        else if(std::fabs(data[i] - data[tmp.first]) > 1e-6) return false;

#ifdef TESTOUTPUT
        if(std::fabs(data[i]-data[tmp.first])/std::max(std::fabs(data[i]),std::fabs(data[tmp.second])) > 1e-9)
        {
          std::cout << "id: " << i << " <-> " << tmp.first << std::endl;
          std::cout << "relative error: " << std::fabs(data[i]-data[tmp.first])/std::max(std::fabs(data[i]),std::fabs(data[tmp.second])) << std::endl;
          std::cout << "absolute error: " << std::fabs(data[i]-data[tmp.first]) << std::endl;
          return false;
        }
#endif

      }

      return true;
    }

    /**
     * \return std::pair<size_t,bool> retVal. retVal.first: index. retVal.second: true, if entry was found
     */
    std::pair<size_t,bool> findEntry(SparseIndexInt row, SparseIndexInt col) const
    {
      for(size_t i=0; i<data.size(); ++i) if(ridx[i]==row && cidx[i]==col) return std::make_pair(i,true);
      return std::make_pair(0,false);
    }

    
    /// row indices
    std::vector<SparseIndexInt>  ridx;
    /// column indices
    std::vector<SparseIndexInt> cidx;
    /// data
    std::vector<Scalar>  data;

  private:
    template <int bcrsN, class SparseInt>
    void getBCRSIndices(SparseInt nrows, std::vector<SparseInt> const& ridx, std::vector<SparseInt> const& cidx, std::vector<std::vector<SparseInt> >& rowColumnIndices) const
    {
      rowColumnIndices.clear();
      rowColumnIndices.resize(nrows);
      for(SparseInt i=0; i<ridx.size(); ++i) rowColumnIndices[ridx[i]/bcrsN].push_back(cidx[i]/bcrsN);

      for(std::vector<SparseInt>& s : rowColumnIndices)
      {
        std::sort(s.begin(),s.end());
        auto last = std::unique(s.begin(),s.end());
        s.erase(last,s.end());
      }
    }

    template <int n, int m, class Matrix>
    void fillFromBCRSInterface(Matrix const& other)
    {
      reserve(other.nonzeroes()*n*m);

      auto rend = other.end();
      for (auto r=other.begin(); r!=rend; ++r)  // iterate over rows
      {
        SparseIndexInt ri = r.index() * n;      // block row start index in matrix with scalar entries
        auto cend = r->end();
        for (auto c=r->begin(); c!=cend; ++c)   // iterate over columns
        {
          SparseIndexInt ci = c.index() * m;    // block col start index in matrix with scalar entries
          auto const& entry = *c;
           
          for(size_t k=0; k<n; ++k)             // iterate over block entries
            for(size_t l=0; l<m; ++l)
            {
              ridx.push_back(ri+k);
              cidx.push_back(ci+l);
              data.push_back(entry[k][l]);
            }
        }
      }
    }
    
    template <class Iterator>
    struct Read
    {
      Read(Iterator& in_): in(in_) {}

      template <class Element>
      void operator()(Element& e) const { in = vectorFromSequence(e,in); }

    private:
      Iterator& in;
    };

    template <class Iterator>
    struct Write
    {
      Write(Iterator& out_): out(out_) {}

      template <class Element>
      void operator()(Element const& e) const { out = vectorToSequence(e,out); }

    private:
      Iterator& out;
    };
  };

  /**
   * \ingroup linalgbasic
   * \brief removes subdiagonal entries from the matrix entries stored in elementary triplet format
   */
  template<class Scalar>
  void deleteLowerTriangle(std::vector<int>& ridx, std::vector<int>& cidx, std::vector<Scalar>& data)
  {
    size_t count(0);
    for(size_t i=0; i<ridx.size(); ++i)
      if(ridx[i] <= cidx[i])
        count++;

    std::vector<int> r2,c2;
    std::vector<Scalar> d2;
    r2.reserve(count);
    c2.reserve(count);
    d2.reserve(count);

    for(size_t i=0; i<ridx.size(); ++i)
      if(ridx[i] <= cidx[i]) {
        r2.push_back(ridx[i]);
        c2.push_back(cidx[i]);
        d2.push_back(data[i]);
      }

    std::swap(r2,ridx);
    std::swap(c2,cidx);
    std::swap(d2,data);
  }

  template <typename Scalar, typename SparseIndexInt>
  std::ostream& operator<<(std::ostream &s, MatrixAsTriplet<Scalar,SparseIndexInt> const& mat)
  {
    return mat.print(s,0);
  }

  namespace AssemblyDetail
  {
    /**
     * \cond internals
     */
    template <class IdxOutIter, class DataOutIter>
    struct BlockToTriplet
    {
      BlockToTriplet(size_t firstRowBlock_, size_t firstColumnBlock_, std::vector<size_t> const& rowOff_, std::vector<size_t> const& colOff_,
          IdxOutIter& ri_, IdxOutIter& ci_, DataOutIter& xi_, bool onlyLowerTriangle_):
            firstRowBlock(firstRowBlock_), firstColumnBlock(firstColumnBlock_), rowOff(rowOff_), colOff(colOff_),
            ri(ri_), ci(ci_), xi(xi_), onlyLowerTriangle(onlyLowerTriangle_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        // Check if block is in requested range
        if (inRange(mb.rowId,mb.colId))
	  Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ri,ci,xi,
								rowOff[MatrixBlock::rowId-firstRowBlock],colOff[MatrixBlock::colId-firstColumnBlock],
								mb.rowId==mb.colId && onlyLowerTriangle,
								mb.symmetric);
        // For mirror blocks, the transposed block needs to be written
        // if the transposed block is in the requested
        // range. Transposition is implicitly achieved by swapping
        // column and row index output iterators.
        if (MatrixBlock::mirror && inRange(mb.colId,mb.rowId) && onlyLowerTriangle==false)
          Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ci,ri,xi,
              colOff[MatrixBlock::rowId-firstColumnBlock],rowOff[MatrixBlock::colId-firstRowBlock],
              false,
              mb.symmetric);
      }

      bool inRange(size_t r, size_t c) const
      {
        return r>=firstRowBlock && r<firstRowBlock+rowOff.size() && c>=firstColumnBlock && c<firstColumnBlock+colOff.size();
      }


      size_t firstRowBlock, firstColumnBlock;
      std::vector<size_t> const& rowOff;
      std::vector<size_t> const& colOff;
      IdxOutIter& ri;
      IdxOutIter& ci;
      DataOutIter& xi;
      bool onlyLowerTriangle;
    };

    /**
     * \brief Specialize this template for your matrix type in order to use it with
     * VariationalFunctionalAssembler::template get<Matrix>() or MatrixRepresentedOperator::template get<Matrix>()
     *
     * Its apply function takes the following parameters:
     * \param block MatrixBlockArray (heterogeneous collection of NumaBCRSMatrices) holding the assembled data. See Matrix_to_Triplet for an example how to access the data.
     * \param firstRowBlock first row in MatrixBlockArray block
     * \param firstColumnBlock first column in MatrixBlockArray block,
     * \param rowOff std::vector of offsets for following blocks
     * \param colOff std::vector of offsets for following blocks
     * \param extractOnlyLowerTriangle if true, only the lower triangular part of the source data is copied
     * \param nnz number of scalar nonzero entries
     * \param nrows number of scalar rows 
     * \param ncols number of scalar columns
     */
    template <class Matrix> struct Fill;
    template <class Matrix> struct FillPointer;

    /**
     *  \brief Allows to use Kaskade::MatrixAsTriplet with VariationalFunctionalAssembler and
     *  classes satisfying the Kaskade::MatrixRepresentedOperator concept.
     */
    template <class Scalar, class SparseIndexInt>
    struct Fill<MatrixAsTriplet<Scalar,SparseIndexInt>>
    {
      typedef MatrixAsTriplet<Scalar,SparseIndexInt> Matrix;

      template <class MatrixBlockArray>
      static Matrix apply(MatrixBlockArray const& block, size_t firstRowBlock, size_t firstColumnBlock, std::vector<size_t> const& rowOff,
                          std::vector<size_t> const& colOff, bool extractOnlyLowerTriangle, size_t nnz, size_t, size_t)
      {
        Matrix result(nnz);

        typename Matrix::index_iterator rowIterator = result.ridx.begin(), 
                                        columnIterator = result.cidx.begin();
        typename Matrix::iterator dataIterator = result.data.begin();
        for_each(block,BlockToTriplet<typename Matrix::index_iterator,typename Matrix::iterator>(firstRowBlock,firstColumnBlock,rowOff,colOff,
                                                                                                 rowIterator, columnIterator, dataIterator, extractOnlyLowerTriangle));

        return result;
      }
    };

    /**
     *  \brief Allows to use Kaskade::MatrixAsTriplet with VariationalFunctionalAssembler and
     *  classes satisfying the Kaskade::MatrixRepresentedOperator concept.
     */
    template <class Scalar, class SparseIndexInt>
    struct FillPointer<MatrixAsTriplet<Scalar,SparseIndexInt> >
    {
      typedef MatrixAsTriplet<Scalar,SparseIndexInt> Matrix;

      template <class MatrixBlockArray>
      static std::unique_ptr<Matrix> apply(MatrixBlockArray const& block, size_t firstRowBlock, size_t firstColumnBlock, std::vector<size_t> const& rowOff,
                                           std::vector<size_t> const& colOff, bool onlyLowerTriangle, size_t nnz, size_t, size_t)
      {
        std::unique_ptr<Matrix> result(new Matrix(nnz));
        
        typename Matrix::index_iterator rowIterator = result->ridx.begin(),
        columnIterator = result->cidx.begin();
        typename Matrix::iterator dataIterator = result->data.begin();
        for_each(block,BlockToTriplet<typename Matrix::index_iterator,typename Matrix::iterator>(firstRowBlock,firstColumnBlock,rowOff,colOff,
                                                                                                 rowIterator, columnIterator, dataIterator, onlyLowerTriangle));
        return result;
      }
    };

    /**
     *  \brief Allows to use Dune::BCRSMatrix with VariationalFunctionalAssembler and
     *  classes satisfying the Kaskade::MatrixRepresentedOperator concept.
     */
    template <class Scalar, class Allocator>
    struct Fill<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Allocator> >
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Allocator> Matrix;

      template <class MatrixBlockArray>
      static Matrix apply(MatrixBlockArray const& block, size_t firstRowBlock, size_t firstColumnBlock, std::vector<size_t> const& rowOff,
                          std::vector<size_t> const& colOff, bool onlyLowerTriangle, size_t nnz, size_t nrows, size_t ncols)
      {
        // read as triplet
        typedef MatrixAsTriplet<Scalar,size_t> Triplet;
        Triplet triplet( Fill<Triplet>::apply(block, firstRowBlock, firstColumnBlock, rowOff, colOff, onlyLowerTriangle, nnz, nrows, ncols) );

        // get column indices for each row
       auto bcrsids = getBCRSIndicesFromTriplet(nrows, triplet.ridx, triplet.cidx);

        Matrix result(nrows,ncols,nnz,Matrix::row_wise);
        // create sparsity pattern
        for (typename Matrix::CreateIterator row=result.createbegin(); row!=result.createend(); ++row)
          for(size_t col = 0; col < bcrsids[row.index()].size(); ++col) 
            row.insert(bcrsids[row.index()][col]);

        // fill matrix
        for (size_t i=0; i<nnz; ++i)
          result[triplet.ridx[i]][triplet.cidx[i]] = triplet.data[i];
        return result;
      }
    };

    /**
     *  \brief Allows to use Dune::BCRSMatrix with VariationalFunctionalAssembler and
     *  classes satisfying the Kaskade::MatrixRepresentedOperator concept.
     */
    template <class Scalar, class Allocator, int nr, int nc>
    struct FillPointer<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,nr,nc>,Allocator> >
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,nr,nc>,Allocator> Matrix;

      template <class MatrixBlockArray>
      static std::unique_ptr<Matrix> apply(MatrixBlockArray const& block, size_t firstRowBlock, size_t firstColumnBlock, std::vector<size_t> const& rowOff,
                                           std::vector<size_t> const& colOff, bool onlyLowerTriangle, size_t nnz, size_t nrows, size_t ncols)
      {
        // make sure that the number of scalar rows and cols is a multiple of the block entry sizes
        assert(nrows%nr == 0 && ncols%nc == 0);
        
        // read as triplet
        typedef MatrixAsTriplet<Scalar,size_t> Triplet;
        Triplet triplet( Fill<Triplet>::apply(block, firstRowBlock, firstColumnBlock, rowOff, colOff, onlyLowerTriangle, nnz, nrows, ncols) );
      
        // Check for NaN in the entries.
#ifndef NDEBUG
        for (size_t i=0; i<triplet.data.size(); ++i) 
          if(std::isnan(triplet.data[i]))
          {
            std::cout << "nan found: " << i << " -> " << triplet.ridx[i] << ", " << triplet.cidx[i] << std::endl;
            exit(-1);
          }
#endif          

        // get column indices for each row
        auto bcrsids = getBCRSIndicesFromTriplet(nrows/nr, triplet.ridx, triplet.cidx, nr, nc);
        size_t bcrsNnz = 0;           // compute number of occupied nr x nc blocks in target matrix
        for (auto const& r: bcrsids)
          bcrsNnz += r.size();

        std::unique_ptr<Matrix> result(new Matrix(nrows/nr,ncols/nc,bcrsNnz,Matrix::row_wise));
        // create sparsity pattern
        for (auto row=result->createbegin(); row!=result->createend(); ++row)
          for(size_t col = 0; col < bcrsids[row.index()].size(); ++col) 
            row.insert(bcrsids[row.index()][col]);

        // fill matrix, entering scalars into the corresponding entry block
        for (size_t i=0; i<nnz; ++i)
          (*result)[triplet.ridx[i]/nr][triplet.cidx[i]/nc][triplet.ridx[i]%nr][triplet.cidx[i]%nc] = triplet.data[i];
        
        return result;
      }
    };
    
    /**
     *  \brief Allows to use NumaBCRSMatrix with VariationalFunctionalAssembler and
     *  classes satisfying the Kaskade::MatrixRepresentedOperator concept.
     */
    template <class Scalar, class Index>
    struct Fill<NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index>>
    {
      typedef NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index> Matrix;

      template <class MatrixBlockArray>
      static Matrix apply(MatrixBlockArray const& block, size_t firstRowBlock, size_t firstColumnBlock, std::vector<size_t> const& rowOff,
                          std::vector<size_t> const& colOff, bool onlyLowerTriangle, size_t nnz, size_t nrows, size_t ncols)
      {
        // read as triplet
        typedef MatrixAsTriplet<Scalar,size_t> Triplet;
        Triplet triplet( Fill<Triplet>::apply(block, firstRowBlock, firstColumnBlock, rowOff, colOff, onlyLowerTriangle, nnz, nrows, ncols) );
        return Matrix(triplet,false);
      }
    };

    
    /**
     * \endcond
     */
  } // end of namespace AssemblyDetail
} // end of namespace Kaskade

#endif
