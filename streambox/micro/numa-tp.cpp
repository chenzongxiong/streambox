/*
 * numa-tp.cpp
 *
 *  Created on: Sep 21, 2016
 *      Author: xzl
 *
 *
 */

#include <numaif.h>
#include <vector>
#include "utilities/threading.hh"
#include "CTPL/ctpl_stl.h"

using namespace std;
using namespace Kaskade;

const int max_tasks = 1024;


/* global, when out of scope (upon program ends), will wait until
 * all outstanding tasks exits.
 */

ctpl::thread_pool p(2);


void burn_cpu() {
	while (1) {
		for (unsigned int i = 0; i < 0xffffffff; i++)
			;
		sched_yield();
	}
}

void burn_cpu_awhile() {
	for (unsigned int i = 0; i < 0xffffffff; i++)
		;
}


void parfor() {

	cout << "par for....\n";

  parallelFor([&](size_t block, size_t nBlocks)
      {
        burn_cpu();
      }, max_tasks);
}

void _sleep(int id) {
	cout << "---- sleep starts\n";
  std::this_thread::sleep_for(std::chrono::seconds(3));
  cout << "---- sleep is done\n";
}

//! sleeps for one second and returns 1
auto sleep_func = [](){ _sleep(1); };

/* if future from async is not moved, upon destruction it will
 * block unitl async is done.
 */
void test_async() {
	cout << "async w/o getting future\n";
	std::async(std::launch::async, sleep_func);
	cout << "after async\n";

	cout << "async w getting future\n";
	auto f = std::async(std::launch::async, sleep_func);
	cout << "after async\n";
	f.get();
	cout << "got future\n";
}

/* whether numa threadpool supports "fire and forget" */
void test_numa_async() {

	NumaThreadPool& pool = NumaThreadPool::instance();

	// w/o taking the returned future. the spanwed task will
	// continue to "free run" (underneath, it's a packaged_task)
	cout << "before run on node\n";
	pool.runOnNode(0, Task(sleep_func));
	cout << "end run on node\n";
	sleep(5);

	// taking the returned future
	cout << "before run on node\n";
	auto f = pool.runOnNode(0, Task(sleep_func));
	cout << "end run on node\n";
	f.wait();
	cout << "got future\n";
}

void test_ctpl_async() {

	auto f = p.push(_sleep);
	cout << "task pushed \n";
	f.get();
	cout << "got future\n";

	cout << "test fire and forget\n";
	p.push(_sleep);
  cout << "task pushed \n";
}

void test_numa_parallel() {

	cout << "run on numa nodes....\n";

	NumaThreadPool& pool = NumaThreadPool::instance();
	const int bundle_count = pool.cpus() * 2;
	const int bundles_per_node = bundle_count / pool.nodes();

	vector<Ticket> tickets(bundle_count);

  for (int i = 0; i < bundles_per_node; i++) {
 	  for (int nodeid = 0; nodeid < pool.nodes(); nodeid++) {
 	  	tickets[nodeid * bundles_per_node + i] =
				pool.runOnNode(nodeid, Task([=] {
					/* ----- lambda starts ----- */
 	  				 cout << "task starts\n";
 	  				 burn_cpu_awhile();
						 cout << "task done\n";
						 return 1;
					 }));
				 /* ----- lambda ends ----- */
 	  	}
 	 }

  for (auto && t : tickets)
  	t.wait();
  cout << "all done\n";

}

int main(int argc, char **argv) {
//	parfor();
//	test_numa_parallel();
//	test_async();
	test_ctpl_async();
	cout << "returned from test func\n";
//	test_numa_async();
}
