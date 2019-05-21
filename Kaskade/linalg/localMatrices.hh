/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LOCALMATRICES_HH
#define LOCALMATRICES_HH

#include <vector>

#include "linalg/threadedMatrix.hh"

namespace Kaskade
{
  
  /**
   * \ingroup linalgbasic
   * \brief Providing a matrix or array interface to LAPACK-ordered entries.
   * 
   * The class represents dense or diagonal matrices, depending on the template
   * argument diagonal.
   * 
   * \tparam Entry the type of matrix entries
   * \tparam diagonal if true, only diagonal
   * \tparam SortedRowIdx a range of row indices (global row, local row)
   * \tparam SortedColIdx a range of column indices (global col, local col)
   */
  template <class Entry, bool diagonal, class SortedRowIdx, class SortedColIdx>
  class LocalMatrix 
  {
  public:
    static bool const lumped = diagonal; // backwards compatibility, remove once NumaBCRSMatrix::scatter has been updated
    
    LocalMatrix(SortedRowIdx const& ridx, SortedColIdx const& cidx, Entry* data_)
    : ridx_(ridx), cidx_(cidx), data(data_)
    {
    }
    
    /**
     * \brief A sequence of (global row, local row) indices.
     */
    SortedRowIdx ridx() const { return ridx_; }
    
    /**
     * \brief A sequence of (global col, local col) indices.
     */
    SortedColIdx cidx() const { return cidx_; }
    
    /**
     * \brief The number of entries stored.
     */
    size_t size() const
    {
      return diagonal? ridx().size(): ridx().size()*cidx().size();
    }
    
    /**
     * \brief Resets the data pointer.
     */
    void relocate(Entry* newData)
    {
      data = newData;
    }
    
    /**
     * \brief Access the matrix entries
     */
    Entry& operator()(int row, int col)
    {
      // For diagonal matrix blocks, only the diagonal of local matrices are computed
      // and accessed. Store these not as a diagonal of a full matrix, but in a 
      // contiguous vector.
      if (diagonal)
      {
        assert(row==col);
        return data[row];
      }
      else
        // LAPACK dense storage scheme.
        return data[row*cidx().size()+col];
    }
    
    Entry const& operator()(int row, int col) const
    {
      // For diagonal matrix blocks, only the diagonal of local matrices are computed
      // and accessed. Store these not as a diagonal of a full matrix, but in a 
      // contiguous vector.
      if (diagonal)
      {
        assert(row==col);
        return data[row];
      }
      else
        return data[row*cidx().size()+col];
    }
    
  private:
    SortedRowIdx ridx_;
    SortedColIdx cidx_;
    
    // a pointer to the raw storage for the matrix entries
    Entry* data;
  };
  
  /**
   * \ingroup linalgbasic
   * \brief A structure for holding a sequence of several local matrices to be filled sequentially
   *        and to be scattered jointly.
   * 
   * This realizes a container of dense or diagonal matrices (depending on the \arg diagonal template parameter),
   * with the critical feature that the entries for all the matrices are stored in a contiguous memory block. This
   * should improve locality of memory access during assembly and scattering of local matrices.
   * 
   * \tparam Entry the type of matrix entries
   * \tparam diagonal if true, only diagonal
   * \tparam SortedRowIdx a range of row indices (global row, local row)
   * \tparam SortedColIdx a range of column indices (global col, local col)
   */
  template <class Entry, bool diagonal, class SortedRowIdx, class SortedColIdx, class IndexType=std::size_t>
  class LocalMatrices
  {
  public:
    using value_type = LocalMatrix<Entry,diagonal,SortedRowIdx,SortedColIdx>;
    
    /**
     * \param maxStorage_ the desired maximum size of the occupied memory. This is no hard limit, merely a hint. Defaults to 256kB.
     */
    LocalMatrices(NumaBCRSMatrix<Entry,IndexType>& globalMatrix_, size_t maxStorage_=256*1024)
    : globalMatrix(globalMatrix_), maxStorage(maxStorage_)
    {
      // Allocate memory as desired. The vector will only grow if there is a local matrix with more entries than defined here.
      localData.reserve(maxStorage/sizeof(Entry));
    }
    
    /**
     * \brief Destructor.
     * 
     * The destructor scatters all remaining local matrices into the global matrix.
     */
    ~LocalMatrices()
    {
      scatter();
    }
    
    // Storage for local stiffness matrices.
    std::vector<value_type> localMatrices;
        
    /**
     * \brief Appends another (zero-initialized) local matrix.
     * 
     * In case that appending the new local matrix exceeds the maximum storage size, the existing local matrices
     * are scattered into the global matrix before the new matrix is created.
     * 
     * \param ridx a range of (global idx, local idx) pairs sorted by global index defining the local matrix rows
     * \param cidx the same for the columns
     */
    void push_back(SortedRowIdx const& ridx, SortedColIdx const& cidx)
    {
      // for diagonal matrices, only the diagonal is stored in a contiguous vector
      size_t size = diagonal? ridx.size(): ridx.size()*cidx.size();
      assert(ridx.size()==ridx.size() || !diagonal);
      
      // Appending the new matrix exceeds the storage capacity. Scatter the existing matrices to their target and 
      // remove them.
      if (localData.capacity() < localData.size()+size)
        scatter();
      
      // Append the new local matrix
      Entry* pos = &*localData.insert(localData.end(),size,Entry(0));
      localMatrices.emplace_back(ridx,cidx,pos);
    }
    
    value_type      & back()       { return localMatrices.back(); }
    value_type const& back() const { return localMatrices.back(); }
      
    /**
     * \brief Scatters the local matrices into the global one and removes them from the local container.
     */
    void scatter()
    {
      globalMatrix.scatter(begin(localMatrices),end(localMatrices));
      clear();
    }
    
    /**
     * \brief clears all the data, leaving an empty state
     */
    void clear() 
    {
      localData.clear();
      localMatrices.clear();
    }
    
    /**
     * \brief reports the size of the local matrices storage in bytes
     * 
     * This can be used to limit the number of local matrices such that their memory fits into the CPU cache.
     */
    size_t storageSize() const
    {
      return localData.size() * sizeof(Entry);
    }
    
    /**
     * \brief reports the maximal desired local matrices storage size 
     */
    size_t storageSizeLimit() const 
    {
      return maxStorage;
    }
    
  private:
    // This allows to hold a single array into which the local matrices are scattered, and
    // thus improves memory locality (for better cache hit rates and less false sharing).
    // Contains the local matrices one after the other, for diagonal matrices, only storage
    // for the diagonal is provided.
    std::vector<Entry>               localData;      // entries of elemental matrices
    NumaBCRSMatrix<Entry,IndexType>& globalMatrix;
    size_t                           maxStorage;     // maximal desired storage space in bytes
  };
  
}

#endif