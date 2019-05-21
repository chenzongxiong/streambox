/*
 * mt-numa.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: xzl
 *
  use cmake:
  make mt-numa.bin -j20

  obsoleted.
  g++-5 -std=c++11 -O2 -o mt-numa.bin mt-numa.cpp -lpthread \
 	 -lboost_date_time -lboost_system \
 	 -lpcre -lnuma -Wall -I../ -I../Kaskade

 */

#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <string.h>

#include <fstream>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <string>
//#include <regex>
#include <pcre.h>
#include <numa.h>
//#include <thread>
#include <atomic>
#include <set>

#include <cstdlib>
#include <ctime>
#include <random>

#include <unordered_set>

#include "tbb/concurrent_unordered_set.h"

#include "utilities/threading.hh" /* for numa allocator */

#include "log.h"

using namespace std;
using namespace boost::posix_time;

/* choose type of tests */
//#define TEST_STR	/* strstr or regex on string buffers */
#define TEST_SET		/* write to hashtable */

/* will evenly distr these threads to numa nodes.
 *
 * thread i --> node (i % node_num)
 *
 * warning if #threads > #cores */
//#define NUM_THREADS     	12
#define NUM_THREADS     	24
//#define NUM_THREADS     	56

#define BUFFER_SIZE 			(500 * 1024 * 1024)
//#define BUFFER_SIZE 			(500 * 1024)   // dbg
//#define BUFFER_SIZE 			(1024 * 1024)   // dbg
#define MAX_NUMA_NODES						8

vector<vector<char *>> buffers(MAX_NUMA_NODES); /* buffers[nid][buf_id] */

/* unordered set on diff NUMA nodes */
using T = uint64_t;
using AllocatorT = Kaskade::NumaAllocator<T>;
//using AllocatorT = tbb::tbb_allocator<T>;
using SetT = set<T, less<T>, AllocatorT>;

//using UnoderedSetT = unordered_set<T, hash<T>, equal_to<T>, AllocatorT>;
//using UnoderedSetT = tbb::concurrent_unordered_set<T, hash<T>, equal_to<T>, AllocatorT>;  // on specific NUMA node
using UnoderedSetT = tbb::concurrent_unordered_set<T>;		// overall NUMA nodes
vector<vector<UnoderedSetT *>> usets(MAX_NUMA_NODES);

struct thr_param {
	int 				tid;
	int					nid_exec; /* the node the thr should be running on */
	int					nid_mem;	/* the node of the buffer consumed by the thr */
	const char *buf;
	UnoderedSetT *set;
	uint64_t 		sz;
};

int num_nodes;
int max_node;
thr_param params[NUM_THREADS];

atomic<uint64_t> total_bytes (0);

uint64_t pcre_regex(const char *buf, uint64_t size) {

	pcre *reCompiled;
	pcre_extra *pcreExtra;
	int pcreExecRet;
	int subStrVec[30];
	const char *pcreErrorStr;
	int pcreErrorOffset;
	char *aStrRegex;
	const char *psubStrMatchStr;
	int j;
	uint64_t count = 0;
	const char * testStrings = buf;

	aStrRegex = (char *)"\\d\\d\\d\\d";
	//std::cout << "Regex to use: " << aStrRegex << std::endl;

	// First, the regex string must be compiled.
	reCompiled = pcre_compile(aStrRegex, 0, &pcreErrorStr, &pcreErrorOffset, NULL);


	// pcre_compile returns NULL on error, and sets pcreErrorOffset & pcreErrorStr
	if(reCompiled == NULL) {
		std::cout << "ERROR: Could not compile " << aStrRegex << " " << pcreErrorStr << std::endl;
		exit(1);
	}

	// Optimize the regex
	pcreExtra = pcre_study(reCompiled, 0, &pcreErrorStr);
	/* pcre_study() returns NULL for both errors and when it can not optimize the regex.  The last argument is how one checks for
	   errors (it is NULL if everything works, and points to an error string otherwise. */
	if(pcreErrorStr != NULL) {
		std::cout << "ERROR: Could not study " << aStrRegex << " " << pcreErrorStr << std::endl;
		exit(1);
	}

	int subject_length = size;
	int start_offset = 0;

	/* find all matched substring*/
	while (1) {
		//std::cout << "string " << testStrings << std::endl;

		pcreExecRet = pcre_exec(reCompiled,
				pcreExtra,
				testStrings,
				subject_length,
				start_offset,           // Start looking at this point
				0,                      // OPTIONS
				subStrVec,              //subStrVec[0] is the position of the first matched substring
				30);                    // Length of subStrVec
		//std::cout << "pcreExecRet is:" << pcreExecRet << std::endl;

		/* no matched string*/
		if(pcreExecRet < 0) {
			goto out;
		} else {
			//printf("Result: We have a match!\n");
			count += pcreExecRet;
//			printf("submatches %d. len %d \n", pcreExecRet, subStrVec[1] -subStrVec[0]);

#if 0
			/* extract the matched substring */
			for(int i=0; i<pcreExecRet; i++) {
/*
				char *substring_start = testStrings + subStrVec[2*i];
				int substring_length = subStrVec[2*i+1] - subStrVec[2*i];
*/


//				printf("%2d: %.*s\n", i, substring_length, substring_start);
				//pcre_get_substring(testStrings, subStrVec, pcreExecRet, i, &(psubStrMatchStr));
				//std::cout << "matched string: " << psubStrMatchStr << std::endl;
			}
			// Free up the substring
//			pcre_free_substring(psubStrMatchStr);
#endif

			//subStrVec[2 * i] is the position of the matched substring i
//			std::cout << "subStrVec: " << (testStrings + subStrVec[0]) << std::endl;
			start_offset = subStrVec[2 * (pcreExecRet - 1) + 1];
		}
		//std::cout << std::endl;
	}
out:
	pcre_free(reCompiled);

	// Free up the EXTRA PCRE value (may be NULL at this point)
	if(pcreExtra != NULL) {
#ifdef PCRE_CONFIG_JIT
		pcre_free_study(pcreExtra);
#else
		pcre_free(pcreExtra);
#endif
	}

//		std::cout << "one regex thread finished, count: " << count << std::endl;
	return count;
}

/* return: count */
uint64_t str_str(const char *buf, uint64_t size) {

//	std::cout << "buffer size" <<  strlen(buf) << std::endl;

	uint64_t count = 0;
	const char * needle = "China";
	const int len = strlen(needle);

	char *p = strstr((char *)buf, needle);
	while(p != NULL) {
		p = p + len;
		count ++;
		p = strstr(p, needle);
	}
//	std::cout << "one str_str thread finished, count: " << count << std::endl;
	return count;
}

void* thr_func(void* t){
	xzl_assert(t);
	thr_param* p = (thr_param *)t;

	/* dbg */
#ifdef TEST_STR
	I("exec thr %d on node %d. buf %.2f MB (0x%lx) on node %d",
			p->tid, p->nid_exec, (double)p->sz/1024/1024, p->buf, p->nid_mem);
#else
	I("exec thr %d on node %d. buf %.2f MB (0x%lx) set (0x%lx) on node %d",
			p->tid, p->nid_exec, (double)p->sz/1024/1024, (unsigned long)(p->buf),
			(unsigned long)(p->set), p->nid_mem);
#endif

	numa_run_on_node(p->nid_exec); /* should cause migration */

	while (1) {

#ifdef TEST_STR
		str_str(p->buf, p->sz);
//		pcre_regex(p->buf, p->sz);
		total_bytes.fetch_add(p->sz);
#else
		{
			/* interpret the char buffer as an array of T. insert each T */
			long count = p->sz / sizeof(T);
			for (long i = 0; i < count; i++) {
//				p->set->insert(((T *)(p->buf))[i]);
				T const & val = ((T *)(p->buf))[i];
				p->set->insert(val);

//				p->set->insert(std::rand());	// debugging

				if (i % 1000 == 999) {
					total_bytes.fetch_add(sizeof(T) * 1000);
				}
			}
//			I("count = %ld, final size = %lu", count, p->set->size());
//			p->set->clear(); /* drop */
		}
#endif

		/* check numa node */
		int cpu = sched_getcpu();
		int node = numa_node_of_cpu(cpu);
		assert(node == p->nid_exec);
	}

	pthread_exit(NULL);
}

inline int thr_to_exec_node(int i) {
	return (i % num_nodes);
}

/* what is the node of the buffer accessed by thr i? */
inline int thr_to_mem_node(int i) {
//	return 0;							// traffic jam
	return (i % num_nodes); // right mapping
//	return ((i + 1) % num_nodes); // bad mapping
//	return ((i + 2) % num_nodes); // also bad mapping
}

/* fill thr params.
 * this assigns buffers to thrs */
void prep_params() {
	for (int i = 0; i < NUM_THREADS; i++) {
		/* assign numa buffer to thr. can assign "wrong values" here */
		int nid_mem = thr_to_mem_node(i);

		params[i].tid = i;
		params[i].nid_exec = thr_to_exec_node(i);
		params[i].nid_mem = thr_to_mem_node(i);
		params[i].sz = BUFFER_SIZE;

		xzl_assert(!buffers[nid_mem].empty());
		params[i].buf = buffers[nid_mem].back();
		buffers[nid_mem].pop_back();

#if 0 // each thr has its own set
		xzl_assert(!usets[nid_mem].empty());
		params[i].set = usets[nid_mem].back();
		usets[nid_mem].pop_back();
#endif

#if 0 // each numa node shares one set
		xzl_assert(!usets[nid_mem].empty());
//		params[i].set = usets[nid_mem].back();
//		usets[nid_mem].pop_back();
		params[i].set = usets[nid_mem][0];
#endif

#if 1
		// all threads on all numa nodes share one set
		xzl_assert(!usets[nid_mem].empty());
		params[i].set = usets[0][0];
#endif
	}
}

int main (int argc, char **argv)
{
	int rc;
	int i;

	pthread_t threads[100];
	pthread_attr_t attr;

	void *status;
	char input_file[100];
	const char *default_input_file = "/ssd/1g.txt";

	max_node = numa_max_node();
	num_nodes = max_node + 1;

	if (argc < 2) {
		strncpy(input_file, default_input_file, 100);
		W("warning: input file path default to %s", input_file);
	} else {
		strncpy(input_file, argv[1], 100);
		I("input file is %s", input_file);
	}

	/* create per thread buffer, evenly on all numa nodes */
	for (int i = 0; i < NUM_THREADS; i++) {
		int nid = thr_to_mem_node(i);
		char *p = (char *)numa_alloc_onnode(BUFFER_SIZE, nid);
		xzl_assert(p);
		buffers[nid].push_back(p);
	}

#ifdef TEST_SET
	/* create per thread unordered set, on all numa nodes */
	for (int i = 0; i < NUM_THREADS; i++) {
		int nid = thr_to_mem_node(i);
//		usets[nid].push_back(new UnoderedSetT(AllocatorT(nid)));	// on given NUMA node
		usets[nid].push_back(new UnoderedSetT());		// NUMA unaware
	}
#endif

	/* debug */
	printf("buffers for total %d numa nodes: ", num_nodes);
	for (int i = 0; i < num_nodes; i++) {
//	for (auto & nbuf : buffers) {
		printf("%lu ", buffers[i].size());
	}
	printf("each %.2f MB", (double)BUFFER_SIZE/1024/1024);
	printf("\n");

#ifdef TEST_STR
	/* fill all buffers with file contents */
	int fd = open(input_file, O_RDONLY);
	xzl_assert(fd);

	struct stat finfo;
	fstat(fd, &finfo);
	xzl_assert(finfo.st_size >= BUFFER_SIZE);

	for (auto & nbuf : buffers) {
		for (auto & buf : nbuf) {
			auto sz = pread(fd, buf, BUFFER_SIZE, 0);
			xzl_assert(sz == BUFFER_SIZE);
			/* zero-terminate each buffer */
			buf[sz - 1] = '\0';
			/* dbg */
			assert((int64_t)strlen(buf) == sz - 1);
		}
	}
#else
	/* fill all buffers with random integers */
	srand(time(NULL));
	const long count = BUFFER_SIZE / sizeof(T);
	for (auto & nbuf : buffers) {
		for (auto & buf : nbuf) {
			for (long i = 0; i < count; i++) {
				((T *)(buf))[i] = std::rand();
			}
		}
	}
#endif

//	/* random number for testing hashtable */
//	std::random_device rd;
//	std::mt19937_64 gen(rd());

	prep_params();

	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for(int i = 0; i < NUM_THREADS; i++) {
		rc = pthread_create(&threads[i], &attr,
				thr_func, /* func to exec */
				(void *)(&params[i])
				);
		xzl_assert(!rc);
	}

	ptime last_ts = boost::posix_time::microsec_clock::local_time();
	uint64_t last_bytes = total_bytes;

	/* periodic check */
	while (1) {
		std::this_thread::sleep_for(std::chrono::milliseconds(2000));

		/* ... check ... */
		auto b = total_bytes.load();
		ptime end = boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration diff = end - last_ts;
		double interval_microsec = (double) diff.total_microseconds();
		double diff_bytes = b - last_bytes;

		if (diff_bytes == 0) { /* no thread report */
			I(".");
		} else {
			I("tput %.4f GB/s",
				diff_bytes / (1024*1024*1024) / (interval_microsec / 1e6 ));
		}

		last_ts = end;
		last_bytes = b;
	}

	// free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for( i=0; i < NUM_THREADS; i++ ){
		rc = pthread_join(threads[i], &status);
		xzl_assert(rc);
	}

/*
		for(i = 0; i < j; i++){
			if(buff[i]){
				numa_free(buff[i], size);
			}
		}

  XXX free sets? XXX
*/

	return 0;
}



