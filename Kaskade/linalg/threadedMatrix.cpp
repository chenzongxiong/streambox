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

#include <algorithm>

#include "dune/grid/config.h"

#include "linalg/threadedMatrix.hh"

namespace Kaskade
{
  namespace ThreadedMatrixDetail
  {
    template <class Index>
    size_t CRSChunkPatternCreator<Index>::balanceForward(size_t const covered, size_t const nnz, int chunks, 
                                                         std::vector<typename CRSChunkPatternCreator<Index>::IndexArray>& moveRows)
    {
      // include the given rows into our chunk by prepending them
      this->firstRow -= moveRows.size();                         // we own more rows
      cols.insert(cols.begin(),moveRows.size(),IndexArray());    // create space for the new rows in front
      std::move(moveRows.begin(),moveRows.end(),cols.begin());   // move rows to our chunk
      moveRows.clear();                                          
      
      // calculate how many rows to cut off at the end 
      size_t target = (this->node()+1)*nnz/chunks - covered;     // number of nonzeroes we should have
      auto pos = cols.begin();
      size_t own = 0;
      for ( ; pos != cols.end() && own <= target; ++pos)         // step forward until we've seen at least as many entries as we should own
        own += pos->size();                                      // add up the number of entries so far
      
      std::move(pos,cols.end(),std::back_inserter(moveRows));   // move the remaining rows to the next chunk.
      
      cols.erase(pos,cols.end());
      this->lastRow -= moveRows.size();
      
      // return the number of entries retained so far
      return covered + own;
    }
    
    template <class Index>
    size_t CRSChunkPatternCreator<Index>::balanceBackward(size_t covered, size_t const nnz, int chunks, 
                                                          std::vector<typename CRSChunkPatternCreator<Index>::IndexArray>& moveRows)
    {
      // include the given rows into our chunk by appending them
      this->lastRow += moveRows.size();
      std::move(moveRows.begin(),moveRows.end(),std::back_inserter(cols));
      moveRows.clear();
      
      size_t target = covered - this->node()*nnz/chunks;
      auto pos = cols.begin();
      size_t own = nonzeroes();
      for ( ; pos!=cols.end() && own>target; ++pos)
        own -= pos->size();
      
      std::move(cols.begin(),pos,std::back_inserter(moveRows));   // move the remaining rows to the next chunk.
      
      cols.erase(cols.begin(),pos);
      this->firstRow += moveRows.size();
      
      return covered - own;
    }

    // explicit instantiaion of method
    template size_t CRSChunkPatternCreator<size_t>::balanceForward(size_t const covered, size_t const nnz, int chunks, 
                                                                   std::vector<CRSChunkPatternCreator<size_t>::IndexArray>& moveRows);
    template size_t CRSChunkPatternCreator<int>::balanceForward(size_t const covered, size_t const nnz, int chunks, 
                                                                std::vector<CRSChunkPatternCreator<int>::IndexArray>& moveRows);
    template size_t CRSChunkPatternCreator<size_t>::balanceBackward(size_t const covered, size_t const nnz, int chunks, 
                                                                    std::vector<CRSChunkPatternCreator<size_t>::IndexArray>& moveRows);
    template size_t CRSChunkPatternCreator<int>::balanceBackward(size_t const covered, size_t const nnz, int chunks, 
                                                                 std::vector<CRSChunkPatternCreator<int>::IndexArray>& moveRows);
    
    
    template <class Index>
    size_t CRSChunkPattern<Index>::nonzeroes() const
    {
      if (this->symmetric())
      {
        size_t nnz = 0;
        for (Index r=0; r<this->last()-this->first(); ++r)
        {
          nnz += 2*(colStarts[r+1]-colStarts[r]); // count entries twice
          if (colStarts[r+1]>colStarts[r] && cols[colStarts[r+1]-1]==r+this->first()) 
            --nnz; // existing diagonal entry has been counted twice - subtract
        }
        return nnz;
      }
      else
        return storage();
    }
    // explicit instantiaion of method
    template size_t CRSChunkPattern<size_t>::nonzeroes() const;
    template size_t CRSChunkPattern<int>::nonzeroes() const;
    
    
    // -----------------------------------------------------------------------------------
    
    template <class Index>
    CRSChunkPattern<Index>::CRSChunkPattern(CRSChunkPatternCreator<Index> const& creator)
    : CRSChunkPatternInfo<Index>(creator)
    , colStarts(NumaAllocator<size_t>(creator.node())), cols(NumaAllocator<Index>(creator.node()))
    {
      // Compute start of row indices.
      colStarts.resize(this->last()-this->first()+1,0);
      for (Index ri=this->first(); ri!=this->last(); ++ri)
        colStarts[ri-this->first()+1] = creator.row(ri).size();
      std::partial_sum(colStarts.begin(),colStarts.end(),colStarts.begin());
      
      // Get space for all nonzeros
      cols.reserve(colStarts.back());
      
      // Copy the column indices
      for (Index ri=this->first(); ri!=this->last(); ++ri)
      {
        auto const& row = creator.row(ri);
        cols.insert(cols.end(),row.begin(),row.end());
      }
      
      // investigate sparsity. Remember we may have an empty chunk with zero rows.
      size_t nRows = this->last()-this->first();
      nnzPerRow = nRows>0? nonzeroes() / nRows
                         : 0;
    }
    // explicit instantiaion of method
    template CRSChunkPattern<size_t>::CRSChunkPattern(CRSChunkPatternCreator<size_t> const& creator);
    template CRSChunkPattern<int>::CRSChunkPattern(CRSChunkPatternCreator<int> const& creator);
  }

  // -----------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------

  template <class Index>
  void NumaCRSPatternCreator<Index>::balance()
  {
    std::vector<typename ChunkCreator::IndexArray> moveRows;
    size_t nnz = nonzeroes(); // total number of nonzeroes
    size_t covered = 0;       // nonzeroes covered in previous chunks
    
    // Move surplus elements towards the end.
    for (int i=0; i<nodes(); ++i)
      covered = creators[i].balanceForward(covered,nnz,nodes(),moveRows);
    assert(moveRows.empty() && covered==nnz);
    
    // Now every chunk range [0,n[ has no more entries than it should, but 
    // maybe too few. Moving surplus elements from later chunks to previous
    // chunks corrects this.
    for (int i=nodes()-1; i>=0; --i)
      covered = creators[i].balanceBackward(covered,nnz,nodes(),moveRows);
    assert(moveRows.empty() && covered==0);
  }
  // explicit instantiaion of method
  template void NumaCRSPatternCreator<size_t>::balance();
  template void NumaCRSPatternCreator<int>::balance();


  
}