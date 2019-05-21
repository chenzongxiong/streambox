/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef THREADING_HH
#define THREADING_HH

#include <functional>
#include <future>
#include <map>
#include <queue>
#include <utility>
#include <vector>

#define BOOST_THREAD_PROVIDES_INTERRUPTIONS

#ifndef BOOST_DISABLE_THREADS
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#endif

//#include <boost/timer/timer.hpp>  // xzl: this seems to conflict with boost/progress.hpp

#include "utilities/kalloc.hh"

namespace Kaskade
{
#ifndef BOOST_DISABLE_THREADS
  /**
   * \ingroup threading
   * \brief A global lock for the Dune::QuadratureRules factory, which is not thread-safe as of 2015-01-01.
   */
  extern boost::mutex DuneQuadratureRulesMutex;

  /**
   * \ingroup threading
   * \brief A global lock for the Dune::GenericReferenceElement singletons, which are not thread-safe.
   */
  extern boost::mutex refElementMutex;
#endif

  //---------------------------------------------------------------------------

  /**
   * \ingroup threading
   * \brief Computes partitioning points such that the sum of weights in each partition is roughly the same.
   *
   * Let \f$ w_i, 0\le i < N \f$ denote the weights x[i]. On exit, x has size \f$ k=n+1 \f$ with values \f$ x_i \f$
   * such that \f$ x_0 = 0 \f$, \f$ x_{k-1} = N\f$, and
   * \f[ \sum_{j=x_i}^{x_{i+1}-1} w_j \approx \frac{1}{k} \sum_{j=0}^{N-1} w_j. \f]
   *
   * \param[in,out] x the array of (nonnegative) weights
   * \param[in]     n the desired number of partitions (positive).
   */
  void equalWeightRanges(std::vector<size_t>& x, size_t n);

  /**
   * \ingroup threading
   * \brief Computes partitioning points of ranges for uniform weight distributions.
   *
   * The functionality is similar to \ref equalWeightRanges, but for uniform weight distribution, the
   * partitioning points can be computed directly.
   *
   * \tparam Index an integral type
   * \param i the number of the range for which the starting point is to be computed. Has to be in [0,n]
   * \param n the number of ranges
   * \param m the total number of entries
   */
  template <class BlockIndex, class Index>
  Index uniformWeightRangeStart(BlockIndex i, BlockIndex n, Index m)
  {
    assert(i>=0 && i<=n && n>0 && m>0);
    return (i*m)/n;
  }

  /**
   * \ingroup threading
   * \brief Computes the range in which an index is to be found when partitioned for uniform weights.
   *
   * \tparam Index an integral type
   * \param j the index for which the range number is to be computed. Has to be in [0,m[
   * \param n the number of ranges
   * \param m the total number of entries
   */
  template <class Index>
  Index uniformWeightRange(Index j, Index n, Index m)
  {
    assert(j>=0 && j<m && n>0 && m>0);

    // We're looking for i such that floor(i*m/n) <= j < floor((i+1)*m/n), i.e. index j has
    // to be in the half-open range given by the partitioning points of the range. Now this
    // implies i*m/n-1 <= j < (i+1)*m/n and also j*n/m-1 < i <= (j+1)*n/m. As i is integer,
    // this implies floor(j*n/m) <= i <= floor((j+1)*n/m). Typically, n/m is small, in which
    // case this closed interval is likely to contain just one natural number - the result.
    Index low = (j*n)/m;       // floor is implied by integer arithmetics
    Index high = ((j+1)*n)/m;  // floor is implied by integer arithmetics

    if (low==high)
    {
      assert(uniformWeightRangeStart(low,n,m)<=j && j<uniformWeightRangeStart(low+1,n,m));
      return low;
    }

    // The implied inequality chains above are not sharp estimates, therefore the interval
    // [low,high] can contain several natural numbers. This is always the case if n>m, i.e., we have
    // more ranges than entries and several ranges are empty. But it can also happen if j*n/m is
    // just below the next integral number. Now we have to walk through all natural numbers in
    // the interval to find the correct range. TODO: is there a direct way?
    Index i = low;
    while (i<high && !(uniformWeightRangeStart(i,n,m)<=j && j<uniformWeightRangeStart(i+1,n,m)))
      ++i;

    assert(uniformWeightRangeStart(i,n,m)<=j && j<uniformWeightRangeStart(i+1,n,m));
    return i;
  }

  //---------------------------------------------------------------------------

  #ifndef BOOST_DISABLE_THREADS
  /**
   * \ingroup threading
   * \brief A concurrent fifo queue.
   *
   * The queue is not copyable as it is intended to hold non-copyable std:packaged_task objects.
   */
  template <class T>
  class ConcurrentQueue {
  public:
    /**
     * \brief Constructs an empty queue.
     */
    ConcurrentQueue()
    {
    }

    /**
     * \brief Moves a queue.
     */
    ConcurrentQueue(ConcurrentQueue&& q)
    {
      boost::lock_guard<boost::mutex> lock(q.mutex);
      queue = std::move(q.queue);
    }

    /**
     * \brief Assignment.
     */
    ConcurrentQueue& operator=(ConcurrentQueue const& q)
    {
      boost::lock_guard<boost::mutex> lock(q.mutex);
      queue = q.queue;
      return *this;
    }

    bool empty() const
    {
      boost::lock_guard<boost::mutex> lock(mutex);
      return queue.empty();
    }

    size_t size() const
    {
      boost::lock_guard<boost::mutex> lock(mutex);
      return queue.size();
    }

    /**
     * \brief Stores an element at the end of the queue.
     */
    void push_back(T&& t)
    {
      {
      boost::lock_guard<boost::mutex> lock(mutex);
      queue.push(std::move(t));
      } // release lock before waking up consumers
      filled.notify_one();
    }

    /**
     * \brief Retrieves the foremost element.
     *
     * This method blocks if the queue is empty and waits for data to become available.
     */
    T pop_front()
    {
      boost::unique_lock<boost::mutex> lock(mutex);

      // Wait for data to become available.
      while (queue.empty())
        filled.wait(lock);

      // extract and remove data. When assignment throws an exception, the queue
      // remains unmodified
      T t = std::move(queue.front());
      queue.pop();
      return t;
    }

  private:
    std::queue<T> queue;
    mutable boost::mutex mutex; // has to be mutable because of empty/size
    boost::condition_variable filled;
  };
#endif

  //----------------------------------------------------------------------------
#ifndef BOOST_DISABLE_THREADS
  /**
   * \ingroup threading
   * \brief Abstract interface for tasks to be scheduled for concurrent execution.
   */
  typedef std::packaged_task<void()> Task;

  /**
   * \ingroup threading
   * \brief Abstract waitable job ticket for submitted tasks.
   */
  typedef std::future<void> Ticket;
#else
  /**
   * \ingroup threading
   * \brief Abstract interface for tasks to be scheduled for concurrent execution. This is just a workaround for missing std::packaged_task on GCC/Windows.
   */
  typedef std::function<void()> Task;

  /**
   * \ingroup threading
   * \brief Abstract waitable job ticket for submitted tasks. This is just a workaround for missing std::future on GCC/Windows.
   */
  struct Ticket
  {
    void get() const {}
    void wait() const {}
  };
#endif



  //----------------------------------------------------------------------------

  /**
   * \ingroup threading
   * \brief Implementation of thread pools suitable for parallelization of (more or less) memory-bound algorithms (not only) on NUMA machines.
   *
   * This class maintains two thread pools satisfying different needs.
   *
   * The threads in the first (global) pool can be moved by the operating system freely between nodes and CPUs, and should be used (by submitting tasks via \ref run)
   * whenever there is no particular need to execute the task on a particular NUMA node, i.e. if the task is not memory-bandwidth-bound or works on
   * data that is not located on a particular NUMA node. The number of global threads defaults to twice the number of available CPUs (such that all CPUs can be
   * busy even if some threads wait on a mutex), but is at least 4 unless limited to a smaller number on construction of the thread pool, see \ref instance.
   *
   * The threads in the second (NUMA) pool are pinned to nodes (the OS is allowed to move them between CPUS on the same node), and should be used
   * (by submitting tasks via \ref runOnNode) whenever the tasks are memory-bandwidth-bound and the data resides on a particular NUMA node.
   * Note that as the threads are locked to the nodes, the operating system cannot move the threads. This is intended to keep
   * the thread close to its data in order to have local memory access. On the other hand, on multi-user machines it can lead to
   * several threads competing for the same CPU while other CPUs are idle. Use the node-locked threads only if locality of memory
   * access is of top priority.
   *
   * Relying on the first-touch policy for controlling memory locality is not guaranteed to keep threads and data close to each other.
   * First, the allocator may re-use memory blocks previously touched and released by a different thread, leading to remote data
   * access. Second, the operating system may decide to move a thread to a different node without knowledge of which thread
   * is memory-bound or compute-bound.
   *
   * Caveat: When submitting tasks recursively to the task pool, it is easy to create a deadlock. While the top level tasks wait for the
   *         completion of the lower level tasks, the lower level task is not processed as all threads in the pool are occupied by the
   *         top level tasks. Take care not to submit more recursive tasks than there are worker threads available.
   */
  class NumaThreadPool {
  public:
    /**
     * \brief Returns a globally unique thread pool instance.
     *
     * On the very first call of this method, the singleton thread pool is created with a limitation of the number of global (unpinned)
     * threads limited by the given number of maxThreads. Later calls with different value of maxThreads do not change the number of threads.
     * In case the number of global threads shall be limited throughout, call this method at the very start of the program.
     *
     * \param maxThreads an upper bound for the number of global threads to create.
     *
     * As it makes little sense to have multiple thread pools fight for physical ressources, a single instance should be employed.
     */
    static NumaThreadPool& instance(int maxThreads = std::numeric_limits<int>::max());

    /**
     * \name System information
     * @{
     */

    /**
     * \brief Reports the number of NUMA nodes (i.e., memory interfaces/CPU sockets)
     */
    int nodes() const
    {
      return nNode;
    }

    /**
     * \brief Reports the total number of CPUs (usually a multiple of nodes)
     */
    int cpus() const
    {
      return nCpu;
    }

    /**
     * \brief Reports the number of CPUs on the given node (usually the same for all nodes).
     */
    int cpus(int node) const
    {
      return cpuByNode[node].size();
    }

    /**
     * \brief Reports the maximal number of CPUs on one node.
     */
    int maxCpusOnNode() const
    {
      return maxCpusPerNode;
    }

    /**
     * @}
     */

    /**
     * \name Task submission
     * @{
     */

    /**
     * \brief Schedules a task to be executed on an arbitrary CPU.
     *
     * Waiting for tasks to be completed can be done by obtaining a waitable Ticket object from the task
     * before submitting it.
     *
     * \param task the task object to be executed.
     */
    Ticket run(Task&& task);

    /**
     * \brief Schedules a task to be executed on a CPU belonging to the given NUMA node.
     *
     * Waiting for tasks to be completed can be done by obtaining a waitable Ticket object from the task
     * before submitting it.
     *
     * \param node the number of the node on which to execute the task. 0 <= node < nodes().
     * \param task the task object to be executed. The object has to live at least until the call to its call operator returns.
     */
    Ticket runOnNode(int node, Task&& task);

    /**
     * @}
     */

    /**
     * \name NUMA-aware memory management
     * @{
     */

    /**
     * \brief Returns the allocator used for the given node.
     */
    Kalloc& allocator(int node);

    /**
     * \brief Allocates memory on a specific node.
     *
     * Note: this is comparatively slow and should only be used for allocating large chunks of memory
     * to be managed locally. The memory has to be released by a subsequent call to deallocate.
     */
    void* allocate(size_t n, int node);

    /**
     * \brief frees a chunk of memory previously allocated
     */
    void deallocate(void* p, size_t n, int node);

    /**
     * \brief Reports the alignment size of allocator at given NUMA node.
     */
    size_t alignment(int node) const;

    /**
     * \brief Tells the allocator to prepare for subsequent allocation of several memory blocks of same size.
     * \param n the requested size of the memory blocks
     * \param k the number of memory blocks that will be requested
     * \param node on which NUMA node
     *
     * Use this as a hint for the allocator that allows to improve its performance.
     */
    void reserve(size_t n, size_t k, int node);

    /**
     * @}
     */

  private:

    // private constructor to be called by instance()
    NumaThreadPool(int maxThreads);
    ~NumaThreadPool();

#ifndef BOOST_DISABLE_THREADS
    Ticket runOnQueue(ConcurrentQueue<Task>& queue, Task&& task);
#endif


    int nCpu, nNode, maxCpusPerNode;
#ifndef BOOST_DISABLE_THREADS
public: // xzl debugging
    std::vector<ConcurrentQueue<Task>>    nodeQueue;  // task queues for passing to worker threads on nodes
private:
    ConcurrentQueue<Task>                 globalQueue;
    boost::thread_group                   threads;    // worker threads
#endif
    std::vector<int>                      nodeByCpu;  // NUMA topology info
    std::vector<std::vector<int>>         cpuByNode;  // NUMA topology info
    std::map<void*,std::pair<size_t,int>> memBlocks;  // NUMA memory allocations
    std::vector<Kalloc>                   nodeMemory; // NUMA memory management
  };

  /**
   * \ingroup threading
   * \brief A parallel for loop that executes the given functor in parallel on different CPUs
   *
   * \tparam Func a functor with an operator ()(int i, int n)
   *
   * \param f the functor to call
   * \param maxTasks the maximal number of tasks to create
   *
   * The given function object is called in parallel for values (i,n) with i ranging from 0 to n-1. The number n is
   * the same for all calls and guaranteed not to exceed maxTasks (but may be smaller).
   *
   * The function returns after the last task has been completed, therefore the computational
   * effort for tasks 0,...,n-1 should be roughly equal, no matter which value n has. For an optimal performance,
   * maxTasks should be either a small multiple of the number of CPUs in the system (such that ideally all CPUs are
   * busy the whole time) or much larger than that (such that an imbalance has only a small impact).
   */
  template <class Func>
  void parallelFor(Func const& f, int maxTasks = std::numeric_limits<int>::max())
  {
    NumaThreadPool& pool = NumaThreadPool::instance();
    int nTasks = std::min(4*pool.cpus(),maxTasks);

    std::vector<Ticket> tickets(nTasks);
    for (int i=0; i<nTasks; ++i)
      tickets[i] = pool.run(Task([i,&f,nTasks] { f(i,nTasks); }));

    for (auto& ticket: tickets)
    {
      ticket.wait();
    }
  }

  /**
   * \ingroup threading
   * \brief A parallel for loop that executes the given functor in parallel on different NUMA nodes
   *
   * \tparam Func a functor with an operator ()(int i, int n)
   *
   * \param f the functor to call
   * \param maxTasks the maximal number of tasks to create
   *
   * The given function object is called in parallel for values (i,n) with i ranging from 0 to n-1. The number n is
   * min(nNodes,maxTasks) for all calls to f. Every task is executed on the node given by arguemnt i.
   *
   * The function returns after the last task has been completed, therefore the computational
   * effort for tasks 0,...,n-1 should be roughly equal, no matter which value n has.
   */
  template <class Func>
  void parallelForNodes(Func const& f, int maxTasks = std::numeric_limits<int>::max())
  {
    NumaThreadPool& pool = NumaThreadPool::instance();
    int nTasks = std::min(pool.nodes(),maxTasks);
    std::vector<Ticket> tickets(nTasks);
    for (int i=0; i<nTasks; ++i)
      tickets[i] = pool.runOnNode(i,Task([i,&f,&nTasks] { f(i,nTasks); }));
    for (auto& ticket: tickets)
      ticket.wait();
  }

  //----------------------------------------------------------------------------

  /**
   * \cond internals
   */
  namespace ThreadingDetail
  {
    class NumaAllocatorBase
    {
    public:
      NumaAllocatorBase(int node);

      size_t max_size() const;

      void* allocate(size_t n);

      void deallocate(void* p, size_t n);

      int node() const
      {
        return nod;
      }

    private:
      int nod;
    public:  // xzl: for debugging
      Kalloc* allocator; // xzl: ctor will point this to existing threadpool's NUMA allocator
    };
  }
  /**
   * \endcond
   */

  /**
   * \ingroup threading
   * \brief An STL allocator that uses memory of a specific NUMA node only.
   */
  template <class T>
  class NumaAllocator
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef T& reference;
    typedef T const& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    // make sure that on copying, the target copy has its data in the same
    // NUMA memory region as the source by default -- this improves data locality.
    typedef std::true_type propagate_on_container_copy_assignment;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type propagate_on_container_swap;
    
    template <class U>
    struct rebind
    {
      typedef NumaAllocator<U> other;
    };

    /**
     * \brief Construct an allocator for allocating on the given NUMA node.
     *
     * \param node the NUMA node on which to allocate the memory. This has to be less than NumaThreadPool::instance().nodes().
     *             For negative values, the memory is allocated in interleaved mode.
     */
    NumaAllocator(int node): alloc(node)
    {
    }

    /* xzl: the following three appears necessary, otherwise rebinding will run into errors 
     * cf: alloc.cpp
     */
    NumaAllocator(const NumaAllocator & other) : alloc(other.node()) { }
    
    template <class U>
    NumaAllocator(const NumaAllocator<U> & other) : alloc(other.node()) { }
    
    ~NumaAllocator() {}
    
    /**
     * \brief Reports the node on which we allocate.
     */
    int node() const
    {
      return alloc.node();
    }

    pointer address( reference x ) const
    {
      return &x;
    }

    const_pointer address( const_reference x ) const
    {
      return &x;
    }

    size_type max_size() const
    {
      return alloc.max_size() / sizeof(T);
    }

    /**
     * \brief Allocates the requested amount of memory.
     *
     * If \arg n == 0, a null pointer is returned.
     *
     * \param n number of objects of type T
     */
    pointer allocate( size_type n, std::allocator<void>::const_pointer /* hint */ = 0 )
    {
      if (n>0)
        return static_cast<pointer>(alloc.allocate(n*sizeof(T)));
      else
        return nullptr;
    }

    void deallocate(pointer p, size_type n)
    {
      if (p)
        alloc.deallocate(static_cast<void*>(p),n*sizeof(T));
    }

    template< class U, class... Args >
    void construct( U* p, Args&&... args )
    {
      ::new((void*)p) U(std::forward<Args>(args)...);
    }

    template <class U>
    void destroy( U* p )
    {
      p->~U();
    }

    /**
     * \brief comparison for equality
     *
     * Allocators compare equal, if they allocate on the same NUMA node.
     */
    template <class U>
    bool operator==(NumaAllocator<U> const& other) const
    {
      return node()==other.node();
    }

    template <class U>
    bool operator!=(NumaAllocator<U> const& other) const
    {
      return !(node() == other.node());
    }

//  private:
  public:			// xzl: for debugging
    ThreadingDetail::NumaAllocatorBase alloc; // xzl: storing a ptr to threadpool's existing allocator
  };

  //----------------------------------------------------------------------------

  /**
   * \brief A utility class implementing appropriate copy semantics for boost mutexes.
   *
   * \todo: design appropriate semantics for move construction/assignment
   */
  class Mutex
  {
  public:
    /**
     * \brief Default constructor.
     *
     * The constructed object's mutex is unlocked.
     */
    Mutex() = default;

    /**
     * \brief Copy constructor.
     *
     * The constructed object's mutex is unlocked, independently of the state of m.
     */
    Mutex(Mutex const& m) {}

    /**
     * \brief Assignment.
     *
     * Our mutex doesn't change its state.
     */
    Mutex& operator=(Mutex const& m) { return *this; }

    #ifndef KASKADE_SEQUENTIAL
    /**
     * \brief provides access to the mutex to perform the locking.
     */
    boost::mutex& get() { return mutex; }
    #endif

  private:
    #ifndef KASKADE_SEQUENTIAL
    boost::mutex mutex;
    #endif
  };

}

#endif
