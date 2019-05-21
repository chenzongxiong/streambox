/*
   compile:

   g++-5 -std=c++11 -O2 -o multi-thread-regex.bin \
    multi-thread-regex.cpp \
   -lpthread -lboost_regex -lre2 -lboost_date_time

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
#include <boost/regex.hpp>
#include <re2/re2.h>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <string>
#include <regex>
#include <boost/regex.hpp>

using namespace std;
using namespace boost::posix_time;

//#define NUM_THREADS     96
#define NUM_THREADS     40

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

void * re2_regex(void * t){
	char *buff = (char *)t;
	std::string matched;
	//RE2 re("(\\d\\d\\d\\d)");
	RE2 re("(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)");
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
	std::cout << "one std_regex thread finished, count is" << count << std::endl;
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


void * str_str(void * t){
	char * buff = (char *)t;
	char * p = strstr(buff, "very");
	int count = 0;

	while(p != NULL){
		p = p + 4;
		count ++;
		p = strstr(p, "very");
	}

	std::cout << "one str_str thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}


int main ()
{
	int rc;
	int i;
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
	void *status;

	const char * fname =
			"/ssd/1g.txt";
//			"/ssd/word_100MB.txt"

	int fd = open(fname, O_RDONLY);
	struct stat finfo; 
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	else
		printf("open %s okay\n", fname);

	fstat(fd, &finfo);
	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);
	std::cout << "size(bytes):" << size << std::endl;

	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	while(1){
		//start time
		ptime begin = boost::posix_time::microsec_clock::local_time();

		for( i=0; i < NUM_THREADS; i++ ){
			//cout << "main() : creating thread, " << i << endl;
			//rc = pthread_create(&threads[i], &attr, wait, (void *)buff );
			rc = pthread_create(&threads[i], &attr, re2_regex, (void *)buff );
//			rc = pthread_create(&threads[i], &attr, std_regex, (void *)buff );
//			rc = pthread_create(&threads[i], &attr, str_str, (void *)buff );
			//rc = pthread_create(&threads[i], &attr, std_regex_2, (void *)buff );
			if (rc){
				cout << "Error:unable to create thread," << rc << endl;
				exit(-1);
			}
		}

		// free attribute and wait for the other threads
		pthread_attr_destroy(&attr);
		for( i=0; i < NUM_THREADS; i++ ){
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
		boost::posix_time::time_duration diff = end - begin;
		double interval_sec = (double) diff.total_milliseconds() / 1000;
		//std::cout << "interval_sec is " << interval_sec << std::endl;
		//std::cout << "size * NUM_THREAD is " << (unsigned long long ) (size *  NUM_THREADS) << std::endl;
		//std::cout << "Bandwidth is: " << (size * NUM_THREADS) / interval_sec / (1024*1024*1024)<< "GB/s" << std::endl;
		std::cout << "Bandwidth is: " << NUM_THREADS * (size / interval_sec / (1024*1024*1024))<< "GB/s" << std::endl;
		//std::cout << "Bandiwdth is " << (((double)((size * NUM_THREADS)/(1024 * 1024 * 1024)) * 1000) / (double) diff.total_milliseconds()) << " GB/s" << std::endl;

	} //end while(1)
	
	cout << "Main: program exiting." << endl;
	pthread_exit(NULL);
}


