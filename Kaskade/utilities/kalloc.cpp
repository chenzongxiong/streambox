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


#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef HAVE_NUMA
#include "numa.h"
#endif


#include "fem/firstless.hh"
#include "detailed_exception.hh"
#include "kalloc.hh"

#ifdef Cygwin
#include <sstream>
std::string to_string(int i)
{
  std::stringstream s;  s << i; return s.str();
}
#else
using std::to_string;
#endif

namespace Kaskade
{
  namespace
  {
    // rounds up n to the next multiple of a
#if 0  //hym: this not gonna be used
    size_t roundUp(size_t n, size_t a)
    {
      return n/a + (n%a ? 1 : 0);
    }
#endif
    void advance(void*& p, ptrdiff_t n)
    {
      char* q = static_cast<char*>(p);
      q += n;
      p = q;
    }

    // obtain system memory
    void* sysAlloc(size_t n, int node)
    {
      #ifdef HAVE_NUMA
      void* p = numa_alloc_onnode(n,node);
      if (p == nullptr)
        throw DetailedException("Memory allocation of " + to_string(n) + " bytes on NUMA node " + to_string(node) + "failed.\n",
                                __FILE__,__LINE__);
      return p;
      #else
      return std::malloc(n);
      #endif
    }

    // release system memory
    void sysFree(void* p, size_t n)
    {
      #ifdef HAVE_NUMA
      numa_free(p,n);
      #else
      std::free(p);
      #endif
    }

    std::string memSize(size_t n)
    {
      static char const* unit[] = { " B ", " KB", " MB", " GB" };
      int k = 0;
      while (n>9999 && k<3)
      {
        ++k;
        n /= 1024;
      }
      return to_string(n) + unit[k];
    }

  }

  Kalloc::Chunk::Chunk(void* mem_)
  : mem(mem_)
  {}

  bool Kalloc::Chunk::operator==(Chunk const& c) const
  {
    return mem == c.mem;
  }

  // small object arrays not yet implemented --> at least align bytes per allocation
  Kalloc::Kalloc(int node_, int align_, size_t blocksize_, bool checking_)
  : node(node_), align(align_), blocksize(blocksize_), large(blocksize), mid(0), nBlocks(0), nRangeBlocks(0), checkLent(checking_)
  {
    assert(blocksize % align == 0);
    assert(blocksize >= 4096);

    // mid-size pages are maintained in buckets of multiples of the alignment size
    // the largest bucket has index bucket(large)
    freeChunks.resize(bucket(large)+1);
    nLentRanges.resize(freeChunks.size(),0);

    totalAllocations = 0;
    globalAllocations = 0;
    freedLocal = 0;
    freedGlobal = 0;
    bucketMiss = 0;
  }

  Kalloc::Kalloc(Kalloc&& kalloc)
  : node(kalloc.node), align(kalloc.align), blocksize(kalloc.blocksize), large(kalloc.large), mid(kalloc.mid),
    checkLent(kalloc.checkLent)
  {
#ifndef BOOST_DISABLE_THREADS
     boost::lock_guard<boost::mutex> lock(kalloc.mutex);
#endif

     freeChunks = std::move(kalloc.freeChunks);
     nLentRanges = std::move(kalloc.nLentRanges);

     nBlocks = kalloc.nBlocks;
     systemBlocks = std::move(kalloc.systemBlocks);
     nRangeBlocks = kalloc.nRangeBlocks;

     lent = std::move(kalloc.lent);

    totalAllocations = 0;
    globalAllocations = 0;
    freedLocal = 0;
    freedGlobal = 0;
    bucketMiss = 0;
  }

  Kalloc::~Kalloc()
  {
    // perform sanity checks
    if (nBlocks > 0)
      for (auto const& entry: systemBlocks)
        sysFree(entry.first,entry.second);
  }

  void* Kalloc::alloc(size_t n)
  {
    // first lock the data structures
#ifndef BOOST_DISABLE_THREADS
    boost::lock_guard<boost::mutex> lock(mutex);
#endif

    // then perform allocation
    return allocUnlocked(n);
  }

  void* Kalloc::allocUnlocked(size_t n)
  {
    size_t request = n;

    // exclude trivial case
    if (n==0)
      return nullptr;

   ++totalAllocations;

    // check in which range domain the request falls
    if (n > large)
    {
      ++nBlocks;
      void* p = sysAlloc(n,node);
      systemBlocks[p] = n;
      if (checkLent)
        lent[p] = n;
      return p;
    }
    else if (n>=mid)
    {
      // look for the range bucket: round up the size to align
      int const b = bucket(n);
      assert(static_cast<size_t>(b) < freeChunks.size()); //hym
      n = bucketsize(b);

      // check for a free range being available
      if (freeChunks[b].empty())
      {
        ++bucketMiss;

        // nope. look for larger ranges to split
        void* p;
        size_t m;
        std::tie(p,m) = findFreePage(b+1);

        // carve out the requested range and add it to the free range list
        insert(p,b);
        advance(p,n);
        m -= n;

        // distribute the remaining memory to the free range buckets
        distribute(p,m);
      }

      // extract the range from the free blocks
      ++nLentRanges[b];
      void* range = extract(b);
      if (checkLent)
        lent[range] = request;
      return range;
    }
    else
      abort(); // not implemented yet - mid==0

    return nullptr; // never get here
  }

  std::pair<void*,size_t> Kalloc::findFreePage(int b)
  {
    // look for available ranges
    for (int c=b; static_cast<size_t>(c)<freeChunks.size(); ++c) //hym
      if (!freeChunks[c].empty())
        return std::make_pair(extract(c),bucketsize(c));

    // none found - request new block from system
    ++nBlocks;
    ++nRangeBlocks;
    void* p = sysAlloc(blocksize,node);
    systemBlocks[p] = blocksize;
    return std::make_pair(p,blocksize);
  }

  void Kalloc::reserve(size_t n, size_t k)
  {
    // only prepare if this falls in the range size
    if (n>large || n<mid || n==0)
      return;

    int const b = bucket(n);
    n = bucketsize(b);

    // check whether enough memory chunks of required size are already available
    if (freeChunks[b].size() >= k)
      return;

    // reserve that much space for chunks in the list
    freeChunks[b].reserve(k);

    k = k - freeChunks[b].size(); // we need that much

    void* p;
    size_t m;

    // obtain new ranges until we've satisfied the demand
    while (k>0)
    {
      std::tie(p,m) = findFreePage(bucket(k*n));

      // carve out at most k chunks
      while (k>0 && m>=n)
      {
        insert(p,b);
        advance(p,n);
        m -= n;
        --k;
      }

      // distribute the remaining memory to the free range buckets
      distribute(p,m);
    }
  }

  void Kalloc::free(void* chunk, size_t n)
  {
    // first lock the data structures
#ifndef BOOST_DISABLE_THREADS
    boost::lock_guard<boost::mutex> lock(mutex);
#endif

    // then free memory
    freeUnlocked(chunk,n);
  }

  void Kalloc::freeUnlocked(void* chunk, size_t n)
  {
    // exclude trivial case
    if (n==0 || chunk==nullptr)
      return;

    if (checkLent)
    {
      auto it =  lent.find(chunk);
      if (it==lent.end()) {
        std::cerr << "free(" << (size_t)(chunk) << "," << n << ") does not correspond to any lent chunk!\n";
        std::cerr << "chunks are: ";
      }
      assert(it!=lent.end());
      if (it->second != n)
        std::cerr << "************ at " << chunk << " given back " << n << " but should be " << it->second << " ***************\n";
    }

    // check in which chunk domain the request falls
    if (n > large)
    {
      auto it = systemBlocks.find(chunk);
      assert(it != systemBlocks.end());
      systemBlocks.erase(it);
      sysFree(chunk,n);
      ++freedGlobal;
      --nBlocks;
    }
    else if (n>=mid)
    {
      ++freedLocal;

      // look for the bucket where this chunk belongs
      int b = bucket(n);

      // put the chunk back into the pool
      insert(chunk,b);
      --nLentRanges[b];
    }
    else
    {
    }
  }

  std::ostream& Kalloc::print(std::ostream& out) const
  {
    // now lock the data structures
#ifndef BOOST_DISABLE_THREADS
    boost::lock_guard<boost::mutex> lock(mutex);
#endif

    out << "----------------------------------------\n";
    out << "Kalloc memory manager on node " << node << "\n";
    out << "total blocks: " << nBlocks << " (" << memSize(nBlocks*blocksize) << ")\n";
    out << "range blocks: " << nRangeBlocks << " (" << memSize(nRangeBlocks*blocksize) << ") \n";
    out << "----------------------------------------\n";
    out << "bucket  range size   free   lent    memory\n";
    out << "----------------------------------------\n";

    // right-aligned output
    out << std::right;
    size_t nTotalFreeRanges = 0, nTotalLentRanges = 0, totalRangeMemory = 0;
    for (int i=0; static_cast<size_t>(i)<freeChunks.size(); ++i) //hym
    {
      out << std::setw(4) << i << std::setw(5) << i << std::setw(8) << memSize(bucketsize(i)) << std::setw(8) << freeChunks[i].size()
          << std::setw(7) << nLentRanges[i] << std::setw(10) << memSize(bucketsize(i)*freeChunks[i].size()) << "\n";
      nTotalFreeRanges += freeChunks[i].size();
      nTotalLentRanges += nLentRanges[i];
      totalRangeMemory += bucketsize(i)*(freeChunks[i].size()+nLentRanges[i]);
    }
    out << "----------------------------------------\n";
    out << "total            " << std::setw(8)<< nTotalFreeRanges << std::setw(7) << nTotalLentRanges << std::setw(10) << memSize(totalRangeMemory) << "\n";
    out << "----------------------------------------\n";
    std::cerr << "total allocations: " << totalAllocations << "  global: " << nBlocks << "  freedLocal: " << freedLocal << "  freedGlobal: " << freedGlobal << "  free managed range memory: " << memSize(freeRangedMemory()) << "  bucket miss rate: " << (float)bucketMiss/totalAllocations << "\n";
    out << "----------------------------------------\n";

    return out;
  }


  // removes the last chunk from the given bucket
  void* Kalloc::extract(int b)
  {
    assert(0<=b && static_cast<size_t>(b)<freeChunks.size() && !freeChunks[b].empty()); //hym
    void* chunk = freeChunks[b].back().mem;
    assert(chunk != nullptr);
    freeChunks[b].pop_back();

    return chunk;
  }

  // inserts the chunk at the end of the given bucket
  void Kalloc::insert(void* chunk, int b)
  {
    assert(chunk != nullptr);
    assert(0<=b && static_cast<size_t>(b)<freeChunks.size()); //hym
    freeChunks[b].push_back(Chunk(chunk));
  }

  // compute the size of ranges stored in the given bucket
  size_t Kalloc::bucketsize(int b) const
  {
    if (b==0) return align;
    if (b%2)  return (1<<(b+1)/2)*align;
    return 3*(1<<(b/2-1))*align;
  }

  // compute the smallest bucket the ranges of which can store n bytes (best fit)
  int Kalloc::bucket(size_t n) const
  {
    // This iterative approach is probably a little bit naive, but appears to be
    // fast enough. It is definitely not the bottleneck as of 2014-03-11.
    int b = 0;
    while (bucketsize(b) < n)
      ++b;

    return b;
  }

  size_t Kalloc::freeRangedMemory() const
  {
    size_t mem = 0;
    for (int i=0; static_cast<size_t>(i)<freeChunks.size(); ++i) //hym
      mem += freeChunks[i].size()*bucketsize(i);
    return mem;
  }

  void Kalloc::coalesce()
  {
    // This implements a kind of brute force approach - take all ranges, sort them by address,
    // and merge adjacent ones. The performance is proably suboptimal - make sure this is called
    // rarely.

    // This stores for all ranges the address and size
    std::vector<std::pair<char*,size_t>> ranges;
    for (int i=0; static_cast<size_t>(i)<freeChunks.size(); ++i) //hym
    {
      size_t bs = bucketsize(i);
      std::transform(freeChunks[i].begin(),freeChunks[i].end(),std::back_inserter(ranges),
                     [=](Chunk const& chunk) { return std::make_pair(static_cast<char*>(chunk.mem),bs); });
    }

    // sort by address
    std::sort(ranges.begin(),ranges.end(),FirstLess());

    // look for runs of consecutive ranges
    for (int i=0; static_cast<size_t>(i)<ranges.size(); ++i) //hym
    {
      int j = i;
      while (static_cast<size_t>(j+1)<ranges.size() && ranges[j].first+ranges[j].second==ranges[j+1].first) //hym
        ++j;
      // now j points to the last range in the consecutive run
      if (j>i)
        merge(ranges.begin()+i,ranges.begin()+j+1);
      i = j; // start over with the next run
    }
  }

  template <class Iter>
  void Kalloc::merge(Iter first, Iter last)
  {
    size_t n = 0; // total size of chunks to be merged
    std::cerr << "merging: ";

    // find all mentioned chunks in the buckets and extract them
    for (Iter i=first; i!=last; ++i)
    {
      int const bs = i->second; // bucket size
      int const b = bucket(bs);  // bucket number
      std::cerr << "[" << std::hex << (size_t)(i->first) << "," <<  bs << "," << std::dec <<  b << "] ";
      assert(i->first);
      assert(0<=b && static_cast<size_t>(b)<freeChunks.size()); //hym
      auto it = std::find(freeChunks[b].begin(),freeChunks[b].end(),Chunk(i->first));
      assert(it != freeChunks[b].end());
      freeChunks[b].erase(std::find(freeChunks[b].begin(),freeChunks[b].end(),Chunk(i->first))); // don't need to check for not found - we know it must be there
      n += bs;
    }
    std::cerr << "\n";

    // scatter the joined memory to our chunk buckets
    distribute(first->first,n);
  }

  void Kalloc::distribute(void* p, size_t m)
  {
    while (m>0)
    {
      // find a suitable bucket: first look for one *into* which the range fits
      int c = std::min(bucket(m),static_cast<int>(freeChunks.size()-1));
      int cs = bucketsize(c);

      // if the fit is not exact, the range is smaller - put it into the
      // next smaller bucket (with smaller ranges) and distribute the
      // remaining  range
      if (static_cast<size_t>(cs)>m) //hym
      {
        --c;
        cs = bucketsize(c);
      }
      if (static_cast<size_t>(cs)>m) //hym
        assert(static_cast<size_t>(cs)<=m); //hym

      // insert the fraction of the range and proceed
      insert(p,c);
      advance(p,cs);
      m -= cs;
    }
  }
}
