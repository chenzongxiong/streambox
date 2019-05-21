#include <numaif.h>
#include <set>
#include <vector>
#include "utilities/threading.hh"
#include "log.h"

#define SIZE_K (1024)
#define SIZE_M (1024 * 1024)
#define SIZE_G (1024 * 1024 * 1024)

using namespace std;
using namespace Kaskade;

int get_numa_node(void *ptr) {
	int numa_node = -1;
	get_mempolicy(&numa_node, NULL, 0, (void*)ptr, MPOL_F_NODE | MPOL_F_ADDR);
	return numa_node;
}

void threadpool() {
	NumaThreadPool& pool = NumaThreadPool::instance();
	I("total cpu %d node %d\n", pool.cpus(), pool.nodes());
}

void parfor() {

	parallelFor([&](size_t block, size_t nBlocks)
			{
			I("hello world\n");
			}, 4);

	parallelFor([&](size_t block, size_t nBlocks)
			{
			I("hello world1\n");
			}, 4);

	parallelFor([&](size_t block, size_t nBlocks)
			{
			I("hello world2\n");
			}, 4);

}

void kalloc_numa()
{
	// create a dedicated allocator (different from the default ones that are
	// already created)
	Kalloc alloc(0);
	//  alloc.print(cout);
	//  alloc.alloc(128);
	alloc.alloc(1 * SIZE_M);
	//  alloc.alloc(96 * SIZE_M);
	alloc.print(cout);
	NumaThreadPool& pool = NumaThreadPool::instance();
	// the default numa memory allocator
	NumaAllocator<int> nalloc(pool.nodes() - 1);
	auto p = nalloc.allocate(1024 * SIZE_M);
	assert(p);
	vector<int, NumaAllocator<int>> as(1024, 0, nalloc);
	as.reserve(512 * SIZE_M);
	cout << "numa node is " << get_numa_node(&(as[0]));
	cout << " ----- nalloc ----- ";
	nalloc.alloc.allocator->print(cout);
}

int main(int argc, char **argv) {
	// global memory pool as part of the thread pool (nodeMemory)
	// NumaAllocator per se has low memory overhead -- only saving a pointer and
	// a node id
	NumaAllocator<int> alloc(0);
	NumaAllocator<long> allocb(0); /* expected by vector<bool> */
	//  std::vector<int,> a();
	std::vector<int, NumaAllocator<int>> a(1024, 0, NumaAllocator<int>(0));
	std::vector<int, NumaAllocator<int>> assss(1024, 0, NumaAllocator<int>(0));
	vector<int, NumaAllocator<int>> a1(1024, 0, alloc);
	vector<int, NumaAllocator<int>> as(alloc);
	using bool_allocator = NumaAllocator<unsigned long>;
	//using bool_allocator = NumaAllocator<bool>;  
	//vector<bool, bool_allocator> boolvec();   /* working */  
	vector<bool, bool_allocator>(bool_allocator(0));
	//  vector<bool, NumaAllocator<bool>> boolvec1(100, true);   

	/* why isn't this working? */
	//  vector<bool, NumaAllocator<bool>> boolvec2(100, true, NumaAllocator<bool>(0));   

	//  vector<int, NumaAllocator<int>> myvec;

	//using rbtreenode = std::_Rb_tree_node<int>; 
	//using set_allocator = NumaAllocator<rbtreenode>;
	using set_allocator = NumaAllocator<int>;
	set_allocator talloc(0);
	//NumaAllocator<int> int_allocator(talloc);
	set<int, std::less<int>, set_allocator> a11(talloc);
	//  set<int, std::less<int>, set_allocator> a3(talloc);
	//  set<int, std::less<int>, NumaAllocator<rbtreenode>> a4();
	threadpool();
	kalloc_numa();
	return 0;
}



