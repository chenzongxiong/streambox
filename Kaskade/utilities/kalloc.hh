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

#ifndef KALLOC_HH
#define KALLOC_HH

#include <map>
#include <vector>

#ifndef BOOST_DISABLE_THREADS
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#endif



namespace Kaskade 
{
  
  /**
   * \ingroup threading
   * \brief A simple memory manager for NUMA systems.
   * 
   * The memory manager is not intended to be a drop-in replacement for the standard malloc (and relies itself on standard malloc).
   * It is intended to implement standard allocators for containers that need to store their data on a particular NUMA node.
   * 
   * The memory is managed in two groups.
   * - Very large memory ranges (larger than half of the blocksize) are allocated and released directly via libnuma. These ranges
   *   are aligned at page boundaries (usually 4k).
   * - Mid and small size memory ranges are carved out from libnuma memory blocks and maintained in buckets of blocks of same size. 
   *   These ranges are aligned at the given alignment value, usually at cache lines of 64 bytes. 
   */
  class Kalloc
  {
    // Meta-data for chunks of memory managed.
    // Currently this just contains a void pointer. Do we need this for future flexibility,
    // i.e. do we want to store some more metadata of chunks? Linked list of free chunks instead
    // of currently relying on system new for managing chunk pointers?
    // If not, we could remove this and simply work with void* internally.
    class Chunk
    {
    public:
      Chunk(void* mem_);
      
      bool operator==(Chunk const& c) const;
      
      void* mem; // pointer to memory chunk 
    };
    
  public:
    /**
     * \brief Constructor.
     * \param node The NUMA node on which to allocate memory.
     * \param align The alignment of the memory to return. Has to be a power of 2, and significantly below blocksize. Defaults to common cache line size of 64 bytes. 
     * \param blocksize The size of memory blocks to be requested from libnuma. Has to be a multiple of \f$ 4096 \f$. Defaults to 2^21, which is two MB.
     * \param checking If true, performs additional usage sanity checks for detecting memory management bugs.
     */
    Kalloc(int node_, int align_=64, size_t blocksize_=2*1024*1024, bool checking=
#if defined(NDEBUG)
 false
#else
 true
#endif
    );
    
    /**
     * \brief Move Constructor.
     */
    Kalloc(Kalloc&& kalloc);  
    
    /**
     * \brief Destructor.
     */
    ~Kalloc();
    
    /**
     * \brief Allocates n bytes of memory.
     * 
     * If n exceeds the alignment size, the returned memory range is guaranteed to be aligned 
     * as to the value specified on construction.
     * 
     * The method is thread-safe.
     * 
     * \param n requested memory size in bytes.
     */
    void* alloc(size_t n);
    
    /**
     * \brief Releases memory range.
     * 
     * The method is thread-safe.
     * 
     * \param p pointer returned previously by alloc
     * \param n memory size as requrested in the corresponding call to alloc
     */
    void free(void* p, size_t n);
    
    /**
     * \brief Tells the allocator that subsequently several blocks of the same size will be requested.
     * 
     * This is a hint to the allocator that can improve its performance.
     * 
     * \param n the size of the memory blocks to be requested
     * \param k the number of the memory blocks to be requested
     */
    void reserve(size_t n, size_t k);
    
    /**
     * \brief Reports the alignment size.
     */
    size_t alignment() const { return align; }
    
    /**
     * \brief Prints memory management statistics to the given stream.
     */
    std::ostream& print(std::ostream& out) const;

  protected:
    /**
     * \brief Allocates n bytes of memory.
     * 
     * This method is NOT thread-safe. Use this only if you know exactly there's only one thread accessing the 
     * allocator. 
     */
    void* allocUnlocked(size_t n);
    
    /**
     * \brief Releases memory range.
     * 
     * This method is NOT thread-safe. Use this only if you know exactly there's only one thread accessing the 
     * allocator. 
     */
    void freeUnlocked(void* p, size_t n);
    
  private:
    // Conceptually constant parameters
    int node;                   // NUMA node on which to allocate 
    int align;                  // memory alignment
    size_t blocksize;           // NUMA memory block size
#ifndef BOOST_DISABLE_THREADS    
    mutable boost::mutex mutex; // mutex for protection against concurrent access
#endif
    
    size_t large;               // maximum size of ranges stored in lists
    size_t mid;                 // minimum size of mid ranges
    
    
    // The free memory ranges we manage
    std::vector<std::vector<Chunk>> freeChunks; // buckets of free memory chunks 
    
    // Statistic data for reporting and on-line performance optimization, e.g., range coalescence
    size_t                 nBlocks;        // number of memory blocks obtained from the system
    size_t                 nRangeBlocks;   // number of memory blocks used for smaller ranges
    std::vector<size_t>    nLentRanges;    // the number of lent memory ranges from that bucket
    std::map<void*,size_t> systemBlocks;   // keeps track of all memory blocks obtained from the system
    
    void* extract(int bucket);
    void insert(void* range, int bucket);
    size_t bucketsize(int bucket) const;
    int bucket(size_t n) const;
    
    size_t freeRangedMemory() const;       // computes the amount of free memory available in range lists
    void coalesce();                       // finds adjacent free ranges and  merges them
    template <class Iter>
    void merge(Iter first, Iter last);     // merges consecutive ranges
    void distribute(void* mem, size_t n);  // subidivides the given memory into chunks of suitable size and sort into buckets
    std::pair<void*,size_t> findFreePage(int b);  // find a bucket (at least b) with a free page, allocating one from system if needed
    
    size_t totalAllocations, globalAllocations, freedLocal, freedGlobal, bucketMiss; // statistics
    
    bool checkLent;                        // if true, enables keeping track of all memory ranges handed out
    std::map<void*,size_t> lent;           // for debugging only: can keep track of all memory ranges handed out
  };
  
  /**
   * \brief A class for local memory management.
   * 
   * Use this only locally, as it is not thread-safe. If unsure, use Kalloc.
   */
  class KallocUnlocked: public Kalloc
  {
  public:
    /**
     * \brief Constructor.
     * \param node The NUMA node on which to allocate memory.
     * \param align The alignment of the memory to return. Has to be a power of 2, and significantly below blocksize. Defaults to common cache line size of 64 bytes. 
     * \param blocksize The size of memory blocks to be requested from libnuma. Has to be a multiple of \f$ 4096 \f$. Defaults to 2^21, which is two MB.
     * \param checking If true, performs additional usage sanity checks for detecting memory management bugs.
     */
    KallocUnlocked(int node_, int align_=64, size_t blocksize_=2*1024*1024, bool checking=true): Kalloc(node_,align_,blocksize_,checking) {}
    
    /**
     * \brief Allocates n bytes of memory.
     * 
     * If n exceeds the alignment size, the returned memory range is guaranteed to be aligned 
     * as to the value specified on construction.
     * 
     * The method is NOT thread-safe.
     * 
     * \param n requested memory size in bytes.
     */
    void* alloc(size_t n) { return allocUnlocked(n); }
    
    /**
     * \brief Releases memory range.
     * 
     * The method is NOT thread-safe.
     * 
     * \param p pointer returned previously by alloc
     * \param n memory size as requrested in the corresponding call to alloc
     */
    void free(void* p, size_t n) { return freeUnlocked(p,n); }
  };
}

#endif

