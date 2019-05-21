/*
 *  compile: 
 
 g++-5 -std=c++11 -O2 -o multi-thread-regex-numa.bin multi-thread-regex-numa-6t.cpp -lpthread \
 -lre2 -lboost_date_time -lpcre -lnuma
 
 */

/*
 * This is used to test NUMA-aware threads
 * 1.k2 configuration:  
 *	node 0 cpus: 0 2 4 6 8 10 12 14 16 18 20 22
 *      node 1 cpus: 1 3 5 7 9 11 13 15 17 19 21 23
 * In this code, we only allocate memory from Node 0
 *
 * 2.Have to define one of there types:
 * 	-- LOCAL_ONLY:  use core 0, 2, 4...
 *      -- REMOTE_ONLY: use core 1, 3, 5...
 *	-- LOCAL_REMOTE: use core 0, 1, 2, 3...
 *
 * 3.Here the NUM_THREADS has to be : 2, 4, 6, 8, 10, or 12. 
 */


/*
 * this is used to test the scalability of regular expression libraries
 * Before running:
 * 1. modify input file name
 * 2. set NUM_THREADS, which is the max # of threads you want to use.  
 */

/*
 * ref: https://www.tutorialspoint.com/cplusplus/cpp_multithreading.htm
 */

/*
 * 1. For std::regex, using -O2 can achive 10X faster performance than using -O0
 * 2. Both std::regex and re2:regex can scall well in multi-thread, but re2::regex is much faster than std::regex
 * 3. If we using 96 cores, the re2 throughput can reach around 20GB/s.
 * 4. Performance numbers:
 *	Uing 1 core
 *		-- re2: 0.59 GB/s
 *		-- std::regex: 0.003 GB/s (using -O0)
 *		-- std::regex: 0.03 GB/s (using -O2)
 *	Using 5 cores
 *		-- re2: 2.2 GB/s
 *		-- std::regex: 0.014 GB/s (using -O0)
 *		-- std::regex: 0.14 GB/s (using -O2)
 *	Using 10 cores
 *		-- re2: 3.5 GB/s
 *		-- std::regex: 0.028 GB/s (using -O0)
 *		-- std::regex: 0.28 GB/s (using -O2)
 */

#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <fstream>
//#include <boost/regex.hpp>
#include <re2/re2.h>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <string>
#include <regex>
//#include <boost/regex.hpp>
#include <pcre.h>
#include <numa.h>
#include <thread>
#include <atomic>

using namespace std;
using namespace boost::posix_time;

#define NUM_THREADS     12

//#define REMOTE 1 //using core: 0,1,2,3,4,5,6,7,8,9,10,11
		 // else useing core: 0,2,4,6,8,10.12,14,16,18,20,22
/*
#ifdef REMOTE
static atomic<int> thread_num (-1);
#else
static atomic<int> thread_num (-2);
#endif
*/

//XXX Have to chose one of them 
//#define LOCAL_ONLY
//#define REMOTE_ONLY
#define LOCAL_REMOTE

#ifdef LOCAL_ONLY // 0, 2, 4 ...
	static atomic<int> thread_num (-2);
#endif

#ifdef REMOTE_ONLY //1,3,5...
	static atomic<int> thread_num (-1);
#endif

#ifdef LOCAL_REMOTE // 0,1,2,3..
	static atomic<int> thread_num (-1);
#endif

int size = 0;

void *wait_on_core(void *t)
{
	int i, s;
	long tid;

	cpu_set_t cpuset;
	pthread_t thread;

	thread = pthread_self();
	CPU_ZERO(&cpuset);
	//core 0, 2, 4, 8... 22 are on node 0
	//thread_num ++;
#ifdef REMOTE
	thread_num.fetch_add(1);
#else
	thread_num.fetch_add(2);
#endif

	CPU_SET(thread_num, &cpuset);
	s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if(s != 0){
		std::cout << "setaffinity faile!!!" << std::endl;
		exit(0);
	}

	sleep(3);
	std::thread::id this_id = std::this_thread::get_id();
	std::cout << "On core: " << sched_getcpu() << std::endl;

	pthread_exit(NULL);
}

void bind_to_core(){
	cpu_set_t cpuset;
	pthread_t thread;

	thread = pthread_self();
	CPU_ZERO(&cpuset);

#ifdef LOCAL_ONLY
	thread_num.fetch_add(2);
#endif

#ifdef REMOTE_ONLY
	thread_num.fetch_add(2);
#endif

#ifdef LOCAL_REMOTE
	thread_num.fetch_add(1);
#endif

	CPU_SET(thread_num, &cpuset);
	int s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if(s != 0){
		std::cout << "setaffinity faile!!!" << std::endl;
		exit(0);
	}
//	std::cout << "On core: " << sched_getcpu() << std::endl;
//	std::cout << sched_getcpu() << std::endl;
}


void *wait(void *t)
{
	int i;
	long tid;

	tid = (long)t;

	sleep(1);
	cout << "Sleeping in thread " << endl;
	cout << "Thread with id : " << tid << "  ...exiting " << endl;
	pthread_exit(NULL);
}

void * pcre_regex(void * t){

	bind_to_core();

	pcre *reCompiled;
	pcre_extra *pcreExtra;
	int pcreExecRet;
	int subStrVec[30];
	const char *pcreErrorStr;
	int pcreErrorOffset;
	char *aStrRegex;
	const char *psubStrMatchStr;
	int j;	
	int count = 0;
	char * testStrings = (char *)t;

	//aStrRegex = "\\d\\d\\d\\d";
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

	int subject_length = strlen(testStrings);
	int start_offset = 0;
	
	/* find all matched substring*/ 
	while (1) { 
		//std::cout << "string " << testStrings << std::endl;
	
		pcreExecRet = pcre_exec(reCompiled,
				pcreExtra,
				testStrings,
				subject_length,
				start_offset,                      // Start looking at this point
				0,                      // OPTIONS
				subStrVec,              //subStrVec[0] is the position of the first matched substring
				30);                    // Length of subStrVec
		//std::cout << "pcreExecRet is:" << pcreExecRet << std::endl;

		/* no matched string*/
		if(pcreExecRet < 0) {
			//std::cout << "NO match" << std::endl;
			//return 0;
			goto out;
		} else {
			//printf("Result: We have a match!\n");
			//std::cout << "Result: We have a match!" << std::endl; 
			count += pcreExecRet;
//			printf("submatches %d. len %d \n", pcreExecRet, subStrVec[1] -subStrVec[0]);

#if 1			
			// PCRE contains a handy function to do the above for you:
			for(int i=0; i<pcreExecRet; i++) {
/*
				char *substring_start = testStrings + subStrVec[2*i];
				int substring_length = subStrVec[2*i+1] - subStrVec[2*i];
*/
				

//				printf("%2d: %.*s\n", i, substring_length, substring_start);

				//pcre_get_substring(testStrings, subStrVec, pcreExecRet, i, &(psubStrMatchStr));
				//std::cout << "matched string: " << psubStrMatchStr << std::endl;
			}
#endif

			//subStrVec[0] is the position of matched substring
//			std::cout << "subStrVec: " << (testStrings + subStrVec[0]) << std::endl;

//			testStrings += subStrVec[0];
			start_offset = subStrVec[1];
//			testStrings += strlen(psubStrMatchStr);

			// Free up the substring
//			pcre_free_substring(psubStrMatchStr);
		} 
		//std::cout << std::endl;
	}
out:
	//printf("done. count =%d\n", count);
	
	// Free up the regular expression.
	pcre_free(reCompiled);

	// Free up the EXTRA PCRE value (may be NULL at this point)
	if(pcreExtra != NULL) {
#ifdef PCRE_CONFIG_JIT
		pcre_free_study(pcreExtra);
#else
		pcre_free(pcreExtra);
#endif
	}
        //std::cout << "cout is " << count << std::endl;	
	pthread_exit(NULL);
}

void * re2_regex(void * t){
	char *buff = (char *)t;
	std::string matched;
	RE2 re("(\\d\\d\\d\\d)");
	//RE2 re("(((http[s]{0,1}|ftp)://)[a-zA-Z0-9\\.\\-]+\\.([a-z]|[A-Z]|[0-9]|[/.]|[~]|[-])*)");
	int count = 0;
	re2::StringPiece input(buff);
	while(RE2::FindAndConsume(&input, re, &matched)) {
		count ++;
	}
	//std::cout << "one re2_regex thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}

void * std_regex(void * t){
	char * buff = (char *)t;
	std::string str(buff);
	std::regex ex("(\\d\\d\\d\\d)");
	std::smatch res;
	int count = 0;
	std::string::const_iterator searchStart(str.cbegin());
	while( std::regex_search(searchStart, str.cend(), res, ex)){
		count ++;
		searchStart += res.position() + res.length();
	}
	//std::cout << "one std_regex thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}

void * std_regex_2(void * t){
	char * buff = (char *)t;
	std::string str(buff);
	std::regex ex("(\\d\\d\\d\\d)");

	std::sregex_iterator it(str.begin(), str.end(), ex);
	std::sregex_iterator reg_end;
	int count = 0;
	for (; it != reg_end; ++it) {
		count ++;
	}
	std::cout << "one std_regex thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}

/*
void * str_str(void *t){
	bind_to_core();
	char *buff = (char *)t;
	char p; 
	int cnt = 0;
	for(int i = 0; i < size; i++){
		buff[i] = 'a';
		cnt ++;
	}
	
	//std::cout << "count is " << cnt << std::endl;
	pthread_exit(NULL);
}
*/

/*
void * str_str(void * t){
	const char *buff = (char *)t;
	const char *target = "very";
	const char *result = buff;
	int count = 0;
	while((result = std::strstr(result, target)) != NULL){
		++result;
		count ++;
	}
	//std::cout << "count is " << count << std::endl;
	pthread_exit(NULL);
}
*/

void * str_str(void * t){
	bind_to_core();
	char * buff = (char *)t;
	//std::cout << "buffer size" <<  strlen(buff) << std::endl;
	char * p = strstr(buff, "very");
	volatile int count = 0;
	
	while(p != NULL){
		//p = p + 4;
		p++;
		count ++;
		p = strstr(p, "very");
	}

	//std::cout << "one str_str thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}

int main ()
{
	int rc;
	int i,j,k;
	//pthread_t threads[NUM_THREADS];
	pthread_t threads[100];
	//char *buff[100];
	//char *buff;
	//int size;
	pthread_attr_t attr;
	void *status;

//	int fd = open("/tmp/large1.txt", O_RDONLY); //100MB
//	int fd = open("/shared/stream_data_set/movies.txt", O_RDONLY);
//	int fd = open("/tmp/1g.txt.part00", O_RDONLY);
//	int fd = open("/tmp/movies.txt.part00.part00.part00.part00.part00.part00", O_RDONLY);

//	int fd = open("/tmp/1gb-2.txt", O_RDONLY);
//	int fd = open("/tmp/10mb.txt", O_RDONLY);
//	int fd = open("/ssd/word_10MB.txt", O_RDONLY);
		int fd = open("/ssd/1g.txt", O_RDONLY);
//	int fd = open("/home/xzl/tmp/word_10MB.txt", O_RDONLY);
	
	struct stat finfo; 
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);
	//char *buff = (char *) malloc(finfo.st_size);

	char *buff = (char *) numa_alloc_onnode(finfo.st_size, 0); //only allocate memory from node 0
	size = pread(fd, buff, finfo.st_size, 0);
	std::cout << "Input size: "<< size/1024/1024 << "MB" << std::endl;

	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for(int j = 2; j <= NUM_THREADS; j = j + 2) {
#ifdef LOCAL_ONLY
		std::cout << "using LOCAL_ONLY cores: ";
		for(int m = 0; m < j; m++){
			std::cout << 2 * m << " "; 
		}
		std::cout << std::endl;
#endif

#ifdef REMOTE_ONLY
		std::cout << "using REMOTE_ONLY cores: ";
		for(int m = 0; m < j; m++){
			std::cout << 2*m +1 << " ";
		}
		std::cout << std::endl;
#endif

#ifdef LOCAL_REMOTE
		std::cout << "using LOCAL_REMOTE cores: ";
		for(int m = 0; m < j; m++){
			std::cout << m << " ";
		}
		std::cout << std::endl;
#endif
/*
		for(i = 0; i < j; i++){
			buff[i] = (char *) numa_alloc_onnode(finfo.st_size, 0); //only allocate memory from node 0
			size = pread(fd, buff[i], finfo.st_size, 0);
		}
		std::cout << "Input size: "<< size/1024/1024 << "MB" << std::endl;
*/
		j = 6;
		//for(k = 0; k < 20; k++){
		while(1){
#ifdef LOCAL_ONLY //0,2,4,6,
			thread_num.store(-2);
#endif

#ifdef REMOTE_ONLY //1,3,5
			thread_num.store(-1);
#endif

#ifdef LOCAL_REMOTE //0,1,2,3
			thread_num.store(-1);
#endif
			ptime begin = boost::posix_time::microsec_clock::local_time();
			//std::cout << "begin is: " << to_simple_string(begin) << std::endl;	
			//for( i=0; i < NUM_THREADS; i++ ){
			for(i = 0; i < j; i++){
				//cout << "main() : creating thread, " << i << endl;
				//rc = pthread_create(&threads[i], &attr, wait, (void *)buff );
				//rc = pthread_create(&threads[i], &attr, re2_regex, (void *)buff );
				//rc = pthread_create(&threads[i], &attr, std_regex, (void *)buff );
				rc = pthread_create(&threads[i], &attr, str_str, (void *)buff);
				//rc = pthread_create(&threads[i], &attr, std_regex_2, (void *)buff );
				//rc = pthread_create(&threads[i], &attr, pcre_regex, (void *)buff[i] );
				if (rc){
					cout << "Error:unable to create thread," << rc << endl;
					exit(-1);
				}
			}

			// free attribute and wait for the other threads
			pthread_attr_destroy(&attr);
			for(i = 0; i < j; i++){
				//for( i=0; i < NUM_THREADS; i++ ){
				rc = pthread_join(threads[i], &status);
				if (rc){
					cout << "Error:unable to join," << rc << endl;
					exit(-1);
				}
				//cout << "Main: completed thread id :" << i ;
				//cout << "  exiting with status :" << status << endl;
			}

			//end time
			ptime end  = boost::posix_time::microsec_clock::local_time();
			//std::cout << "end is: " << to_simple_string(end) << std::endl;
			boost::posix_time::time_duration diff = end - begin;
		/*	double interval_sec = (double) diff.total_milliseconds() / 1000;
			std::cout << "interval_sec = " << interval_sec << std::endl;
			//std::cout << "Throughput is: " << NUM_THREADS * (size / interval_sec / (1024*1024*1024))<< "GB/s" << std::endl;
			std::cout << "Bandwidth is: " << j * (size / interval_sec / (1024*1024*1024))<< "GB/s" << std::endl;
		*/
			//double interval_msec = (double) diff.total_milliseconds();
			double interval_microsec = (double) diff.total_microseconds();
			//std::cout << "interval_microsec: " << interval_microsec << std::endl;
			//std::cout << "Throughput is: " << (j * (size / interval_microsec / (1024*1024*1024))) / (1000 * 1000) << "GB/s" << std::endl;			
			//std::cout << "size is: " << size << "; interval_microsec is: " << interval_microsec << std::endl;
			std::cout << "Throughput is: " << j * (size/interval_microsec) * (1000*1000) / (1024*1024*1024) << "GB/s" << std::endl;

		}//end of for(k)
		std::cout << "--------" << std::endl;

/*
		for(i = 0; i < j; i++){
			if(buff[i]){
				numa_free(buff[i], size);				
			}
		}
*/
	}//for(j)

	cout << "Main: program exiting." << endl;
		
	if(buff){
		numa_free(buff, finfo.st_size);
	}
	 	
	pthread_exit(NULL);
}


