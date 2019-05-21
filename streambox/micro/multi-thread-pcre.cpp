/*
 *  compile: 
 
 g++-5 -std=c++11 -O2 -o multi-thread-pcre.bin multi-thread-pcre.cpp -lpthread \
 -lre2 -lboost_date_time -lpcre
 
 
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

using namespace std;
using namespace boost::posix_time;

#define NUM_THREADS     24

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

				pcre_get_substring(testStrings, subStrVec, pcreExecRet, i, &(psubStrMatchStr));
				std::cout << "matched string: " << psubStrMatchStr << std::endl;
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
        std::cout << "cout is " << count << std::endl;	
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


void * str_str(void * t){
	char * buff = (char *)t;
	char * p = strstr(buff, "very");
	int count = 0;

	while(p != NULL){
		p = p + 4;
		count ++;
		p = strstr(p, "very");
	}

	//std::cout << "one str_str thread finished, count is" << count << std::endl;
	pthread_exit(NULL);
}


int main ()
{
	int rc;
	int i;
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
	void *status;

//	int fd = open("/tmp/large1.txt", O_RDONLY); //100MB
//	int fd = open("/shared/stream_data_set/movies.txt", O_RDONLY);
//	int fd = open("/tmp/1g.txt.part00", O_RDONLY);
//	int fd = open("/tmp/movies.txt.part00.part00.part00.part00.part00.part00", O_RDONLY);

	//int fd = open("/tmp/1gb.txt", O_RDONLY);
	int fd = open("/tmp/10mb.txt", O_RDONLY);
//	int fd = open("/ssd/word_10MB.txt", O_RDONLY);
//	int fd = open("/home/xzl/tmp/word_10MB.txt", O_RDONLY);
	
	struct stat finfo; 
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);
	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);
	std::cout << "Input size: "<< size/1024/1024 << "MB" << std::endl;

	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	while(1){
		//start time
		ptime begin = boost::posix_time::microsec_clock::local_time();

		for( i=0; i < NUM_THREADS; i++ ){
			//cout << "main() : creating thread, " << i << endl;
			//rc = pthread_create(&threads[i], &attr, wait, (void *)buff );
			//rc = pthread_create(&threads[i], &attr, re2_regex, (void *)buff );
			//rc = pthread_create(&threads[i], &attr, std_regex, (void *)buff );
			//rc = pthread_create(&threads[i], &attr, str_str, (void *)buff );
			//rc = pthread_create(&threads[i], &attr, std_regex_2, (void *)buff );
			rc = pthread_create(&threads[i], &attr, pcre_regex, (void *)buff );
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


