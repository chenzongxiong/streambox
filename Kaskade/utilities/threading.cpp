/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <new>
#include <numeric>  // needed for clang++

#ifdef HAVE_NUMA
#include "numa.h"
#endif

// #undef NDEBUG

#include "utilities/threading.hh"

namespace Kaskade
{
#ifndef BOOST_DISABLE_THREADS
  boost::mutex DuneQuadratureRulesMutex;
  boost::mutex refElementMutex;
#endif

  //----------------------------------------------------------------------------

  void equalWeightRanges(std::vector<size_t>& w, size_t nRequest)
  {
    size_t totalWeight = std::accumulate(w.begin(),w.end(),0);
    size_t n = std::min(nRequest,w.size());

    std::vector<size_t> x(n+1);
    // first entry is always 0
    x[0] = 0;

    size_t weight = 0; // count the weight in partition i-1
    int k = 1;
    while (weight<totalWeight)
    {
      assert(k<=(int)n);

      x[k] = x[k-1]+1;           // at least one entry per partition
      weight += w[x[k-1]];       // and its weight is counted

      // add more entries to the current partition until we're at the desired weight
      while (weight < (k*totalWeight)/n || (x[k] < w.size() && w[x[k]]==0))
      {
        weight += w[x[k]];
        ++x[k];
        assert(x[k]<=w.size());
      }

      ++k;
    }
    assert(x[k-1] == w.size());

    // maybe we ended up with fewer partitions than intended. This can happen
    // if the weight distribution is highly nonuniform.
    x.erase(x.begin()+k,x.end());

    // maybe we ended up with fewer partitions than requested. This happens, e.g.,
    // if few items are given in the input. Specify the missing partitions as empty.
    x.insert(x.end(),nRequest-x.size()+1,x.back());

    // return the partitions in x
    std::swap(w,x);
  }


  //----------------------------------------------------------------------------

  namespace {


    //----------------------------------------------------------------------------
#ifndef BOOST_DISABLE_THREADS
    // Code to be executed by worker threads
    class Worker
    {
    public:
      // Create thread pinned on given node (unless node is negative, then no pinning is done)
      Worker(int node_, ConcurrentQueue<Task>& tasks_)
      : node(node_), tasks(tasks_)
      {
      }

      void operator()()
      {

//        assert(boost::this_thread::interruption_enabled());

#ifdef HAVE_NUMA
        if (numa_available()>=0 && node>=0)
        {
          numa_run_on_node(node);   // bind thread to the NUMA node, leaving choice of CPU to the operating system
          numa_set_preferred(node); // bind memory allocation from this thread to the memory of this node
        }
#endif

//        assert(boost::this_thread::interruption_enabled());

        Task t = tasks.pop_front();  //  get commission

        while (true)                 // run forever, i.e. until explicit interruption
        {
          t();                       // do the work
          t = tasks.pop_front();     // get new commission
        }
      }

    private:
      int node;                      // the node this worker thread runs on
      ConcurrentQueue<Task>& tasks;  // the queue from which tasks are extracted
    };
#endif

  } // end of anonymous namespace


  //----------------------------------------------------------------------------

  namespace ThreadingDetail
  {

    NumaAllocatorBase::NumaAllocatorBase(int node_): nod(node_)
    {
#ifdef HAVE_NUMA
      if (nod >= 0)
        allocator = &NumaThreadPool::instance().allocator(nod);
      else
        allocator = nullptr;
#endif
    }

    size_t NumaAllocatorBase::max_size() const
    {
#ifdef HAVE_NUMA
      long free;
      if (nod>=0)
        return numa_node_size(nod,&free);
      else
      {
        std::allocator<char> a;
        return a.max_size();
      }
#else
      std::allocator<char> a;
      return a.max_size();
#endif
    }

    void* NumaAllocatorBase::allocate(size_t n)
    {
#ifdef HAVE_NUMA
      void* mem = nod>=0? allocator->alloc(n): numa_alloc_interleaved(n);
      if (mem==nullptr) {
	std::cerr << "cannot serve memory request (size=" << n << " B) on node " << nod << "\n";
	throw std::bad_alloc();
      }
#else
      void* mem = std::malloc(n);
#endif

      return mem;
    }

    void NumaAllocatorBase::deallocate(void* p, size_t n)
    {

#ifdef HAVE_NUMA
      if (nod>=0)
        allocator->free(p,n);
      else
        numa_free(p,n);
#else
      std::free(p);
#endif
    }

  }

  //----------------------------------------------------------------------------


  NumaThreadPool& NumaThreadPool::instance(int maxThreads)
  {
    static NumaThreadPool pool(maxThreads);
    return pool;
  }


  Ticket NumaThreadPool::run(Task&& task)
  {
#ifdef BOOST_DISABLE_THREADS
  task();
  return Ticket();
#else
    return runOnQueue(globalQueue,std::move(task));
#endif

  }

  Ticket NumaThreadPool::runOnNode(int node, Task&& task)
  {
#ifdef BOOST_DISABLE_THREADS
  task();
  return Ticket();
#else
    assert(0<=node && node<nodes());
    return runOnQueue(nodeQueue[node],std::move(task));
#endif
  }

#ifndef BOOST_DISABLE_THREADS
  Ticket NumaThreadPool::runOnQueue(ConcurrentQueue<Task>& queue, Task&& task)
  {
    Ticket tick = task.get_future();
  #if defined(KASKADE_SEQUENTIAL)
    // execute the task immediately
    task();
  #else
    // submit it to the worker thread
    queue.push_back(std::move(task));
  #endif
    return tick;
  }
#endif


  NumaThreadPool::NumaThreadPool(int maxThreads)
  {
#ifndef BOOST_DISABLE_THREADS
  nCpu = boost::thread::hardware_concurrency();
#else
  nCpu = 1;
#endif

#if defined(HAVE_NUMA)
    // initialize libnuma
    if (numa_available()<0) // no numa available
    {
      nNode = 1;

      // associate cpus and nodes trivially
      nodeByCpu.resize(nCpu,0);

      cpuByNode.resize(nNode);
      for (int c=0; c<nCpu; ++c)
        cpuByNode[0].push_back(c);
    }
    else // we have numa - look out for number of nodes and cpus association to nodes
    {
      nNode = numa_max_node()+1;
      nCpu = numa_num_configured_cpus();

      nodeByCpu.resize(nCpu);
      cpuByNode.resize(nNode);

      // for each node, obtain the cpus in this node
      for (int c=0; c<nCpu; ++c)
      {
        int n = numa_node_of_cpu(c);
        nodeByCpu[c] = n;
        cpuByNode[n].push_back(c);
      }

    }
#else // we don't know nothing about NUMA - probably we're on a single-socket-multi-core system or pureley sequential
      nNode = 1;

      // associate cpus and nodes trivially: only one node, containing all the cpus
      nodeByCpu.resize(nCpu,0);

      cpuByNode.resize(nNode);
      for (int c=0; c<nCpu; ++c)
        cpuByNode[0].push_back(c);
#endif

    nodeMemory.reserve(nNode);
    for (int i=0; i<nNode; ++i)
      nodeMemory.push_back(Kalloc(i,64,4*1024*1024,false));

#ifndef BOOST_DISABLE_THREADS
    nodeQueue.resize(nNode);

    // On each node create one thread per cpu listening for the node task queue.
    for (int n=0; n<nNode; ++n)
      for (unsigned int i=0; i<cpuByNode[n].size(); ++i)
        threads.create_thread(Worker(n,nodeQueue[n]));

    // Create as many unpinned global threads as allowed.
    for (int i=0; i<std::min(maxThreads,std::max(4,2*nCpu)); ++i)
      threads.create_thread(Worker(-1,globalQueue));
#endif

    // compute maximal number of CPUs on any node
    maxCpusPerNode = 0;
    for (auto const& cpus: cpuByNode)
      maxCpusPerNode = std::max(maxCpusPerNode,static_cast<int>(cpus.size()));
  }

  Kalloc& NumaThreadPool::allocator(int node)
  {
    assert(0<=node && node<nNode);
    return nodeMemory[node];
  }

  void* NumaThreadPool::allocate(size_t len, int node)
  {
    assert(0<=node && node<nNode);
    void* mem = nodeMemory[node].alloc(len);

    if (mem==nullptr)
      throw std::bad_alloc();

    return mem;
  }

  void NumaThreadPool::deallocate(void* mem, size_t n, int node)
  {
    // look for memory block to get its length (required by numa)
    nodeMemory[node].free(mem,n);
  }

  size_t NumaThreadPool::alignment(int node) const
  {
    assert(0<=node && node<(int64_t)nodeMemory.size());
    return nodeMemory[node].alignment();
  }

  void NumaThreadPool::reserve(size_t n, size_t k, int node)
  {
    assert(0<=node && node<(int64_t)nodeMemory.size());
    nodeMemory[node].reserve(n,k);
  }

  NumaThreadPool::~NumaThreadPool()
  {
#ifndef BOOST_DISABLE_THREADS
    // tell all worker threads to stop
    threads.interrupt_all();
    // wait for worker threads to terminate
    threads.join_all();
#endif

    // clean up the memory
    for (auto mem: memBlocks)
      nodeMemory[mem.second.second].free(mem.first,mem.second.first);
  }

  //----------------------------------------------------------------------------


}
