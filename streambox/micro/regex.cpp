/*
	g++-5 -std=c++11 regex.cpp -g -O2 -o regex.bin -lboost_regex -lre2 -lboost_date_time
*/
#include <iostream>
#include <regex>
#include <string>
#include <fstream>
#include <boost/regex.hpp>
#include <re2/re2.h>
#include "boost/date_time/posix_time/posix_time.hpp"

using namespace std;
using namespace boost::posix_time;

/* http://stackoverflow.com/questions/5620235/cpp-regular-expression-to-validate-url */
int main1(const char * str)
{
  std::string url (str);
  unsigned counter = 0;

  std::regex url_regex (
    //R"((?:[a-z][\w-]+:(?:\/{1,3}|[a-z0-9%])|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}\/)(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|[^\s`!()\[\]{};:'".,<>?������������]))",
    R"(^(([^:\/?#]+):)?(//([^\/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?)",
    //R"(ftp|http|https)://(w+.)*(w*)/([wd]+/{0,1})+)", 
    //R"xxx(^(([^:/?#]+):)?(//([^/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?)xxx",
    std::regex::extended
  );
  std::smatch url_match_result;

  std::cout << "Checking: " << url << std::endl;

  //if (std::regex_match(url, url_match_result, url_regex)) {
   if (std::regex_search(url, url_match_result, url_regex)) {
   for (const auto& res : url_match_result) {
      std::cout << counter++ << ": " << res << std::endl;
    }
  } else {
    std::cerr << "Malformed url." << std::endl;
  }

  return EXIT_SUCCESS;
}

int main2(const char * st){
    std::regex exp("(\\b\\S*\\b)");
    std::smatch res;
    //std::string str = "first second third forth";
    std::string str = "http://davinci.marc.gatech.edu/~tad/dennis/no-cense.htm http://gutenberg.net.au";

    std::string::const_iterator searchStart( str.cbegin() );
    while ( regex_search( searchStart, str.cend(), res, exp ) )
    {
        std::cout << ( searchStart == str.cbegin() ? "" : " " ) << res[0] << std::endl;
        searchStart += res.position() + res.length();
    }
    std::cout << std::endl;
    return 0;
}

int main3(void){
    std::regex exp("(\\b\\S*\\b)");
    std::smatch res;
    std::string str = " hello world http://davinci.marc.gatech.edu/~tad/dennis/no-cense.htm http://gutenberg.net.au";

    unsigned counter = 0;
    std::regex url_regex (
        R"(^(([^:\/?#]+):)?(//([^\/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?)",
        std::regex::extended
    );
    std::smatch url_match_result;

    std::string::const_iterator searchStart( str.cbegin() );
    while ( regex_search( searchStart, str.cend(), res, exp ) ){
    	//std::cout << ( searchStart == str.cbegin() ? "" : " " ) << res[0] << std::endl;
	std::cout << res[0] << std::endl;
 	//res[0] one of two URL
	counter = 0;
	std::string url (res[0]);
	if (std::regex_search(url, url_match_result, url_regex)) {
		for (const auto& rs : url_match_result) {
			std::cout << counter++ << ": " << rs << std::endl;
		}
	}

	searchStart += res.position() + res.length();
    }
    std::cout << std::endl;
    return 0;
}

int main4(void){
    std::ifstream input("/tmp/small.txt");
    //for(std::string line; getline(input, line);){
    //	std::cout << line << std::endl;
    //}
    
    std::regex exp("(\\b\\S*\\b)");
    std::smatch res;

    unsigned counter = 0;
    std::regex url_regex (
		    R"(^(([^:\/?#]+):)?(//([^\/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?)",
		    std::regex::extended
		    );
    std::smatch url_match_result;

    std::string line;
    while(std::getline(input, line)){
	    //std::cout << line << std::endl;
	    std::string::const_iterator searchStart( line.cbegin() );
	    while ( std::regex_search( searchStart, line.cend(), res, exp ) ){
		    //std::cout << ( searchStart == str.cbegin() ? "" : " " ) << res[0] << std::endl;
		    //std::cout << res[0] << std::endl;
		    //res[0] one of two URL
		    counter = 0;
		    std::string url (res[0]);
		    if (std::regex_match(url, url_match_result, url_regex)) {
			    for (const auto& rs : url_match_result) {
			    	    std::cout << counter++ << ": " << rs << std::endl;
			    }
			    //std::cout << std::string url_string (url_match_result) << std::endl;
		    }

		    searchStart += res.position() + res.length();
	    }
	    std::cout << std::endl;
    
    }

    return 0;
}

int main5(){
	std::string url = "hello world http://davinci.marc.gatech.edu/~tad/dennis/no-cense.htm hello worl";
	std::string url2 = "http://www.davinci.marc.gatech.edu/~tad/dennis/no-cense.htm";

	//boost::regex re("(ftp|http|https)://(w+.)*(w*)/([wd]+/{0,1})+");
	boost::regex re( "(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)" );
	if (!boost::regex_match(url2, re))
			std::cout << "Your URL is not formatted correctly!" << std::endl;
	        else	
			std::cout << "The regexp " << re << " is invalid!" << std::endl;
}

int main6(){
    //std::ifstream input("/home/miaohy/temp/large.txt");
    //std::ifstream input("/tmp/small.txt");
    std::ifstream input("/tmp/large.txt");
    std::regex exp("(\\b\\S*\\b)");
    std::smatch res;
    
    boost::regex re( "(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)" );
    unsigned counter = 0;

    std::string line;
    while(std::getline(input, line)){
	    std::string::const_iterator searchStart( line.cbegin() );
	    while ( std::regex_search( searchStart, line.cend(), res, exp ) ){
		    //counter = 0;
		    std::string url (res[0]);
		    
		    if (boost::regex_match(url, re)){
		    	counter ++;
		    	//std::cout << "valid url: " << url << std::endl;
		    }
		    searchStart += res.position() + res.length();
	    }
    }
    std::cout << "main 6 done, count is :" << counter << std::endl;
    return 0;
}

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
// find all RUL
int main7(){

	int fd = open("/tmp/large.txt", O_RDONLY);
	struct stat finfo;
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);

	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);

	buff[size - 1] = '\0';

	cout << "entire buff is " << size << std::endl;

//	std::regex ex("(\\d\\d\\d\\d)");
	std::regex ex(R"(\d\d\d\d)");
	std::regex exp("(\\b\\S*\\b)");
	boost::regex re( "(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)" );
	//std::cmatch m;
	std::smatch res;
	int count = 0;
	std::string str (buff);
	std::string::const_iterator searchStart( str.cbegin() );
	while ( std::regex_search( searchStart, str.cend(), res, exp ) ){
		std::string url (res[0]);
		if (boost::regex_match(url, re)){
			count ++;
		}
		searchStart += res.position() + res.length();
	}
	std::cout << "count is " << count << std::endl;


    return 0;
}

//hym: strstr single thread
int main8 ()
{
 //int fd = open("/tmp/small.txt", O_RDONLY);
 int fd = open("/tmp/large.txt", O_RDONLY);
 struct stat finfo;
 if(fd < 0)
         std::cout << "open fail!" << std::endl;
 fstat(fd, &finfo);
 
 char *buff = (char *) malloc(finfo.st_size);
 int size = pread(fd, buff, finfo.st_size, 0); 
 std::cout << size << std::endl;

 int count;

 while(1){
 	char *p = strstr(buff, "very");
 	count = 0;
	ptime begin = boost::posix_time::microsec_clock::local_time();
 	
	while(p != NULL){
 		p++;
		p = strstr(p, "very");
		count ++;
 	}
	
	ptime end  = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration diff = end - begin;
	double interval_sec = (double) diff.total_milliseconds() / 1000;
	std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;

 	
 }
 return 0;
}

//hym: re2 single thread
int main9(){
	cout << " ----------- RE2 --------- \n";

	int fd = open("/tmp/large.txt", O_RDONLY);
	//int fd = open("/tmp/small.txt", O_RDONLY);
	struct stat finfo;
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);

	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);

	std::cout << size << std::endl;
	std::string matched, match1;
	//RE2 re("(\\d\\d\\d\\d)");

	// works
//	RE2 re("(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)");

	//'matched' is the entire URL string
	RE2 re("(((http[s]{0,1}|ftp)://)[a-zA-Z0-9\\.\\-]+\\.([a-z]|[A-Z]|[0-9]|[/.]|[~]|[\\-])*)");

	int count;
	//while(1){
		count = 0;
		ptime begin = boost::posix_time::microsec_clock::local_time();

		re2::StringPiece input(buff);
		while(RE2::FindAndConsume(&input, re, &matched, &match1)) {
			count ++;
			std::cout << matched << match1 << std::endl;
		}
		ptime end  = boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration diff = end - begin;
		double interval_sec = (double) diff.total_milliseconds() / 1000;		
		std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;
		std::cout << count << std::endl;
	//}
	return 0;
}

// std::regex : only find one match -> wrong
int main10(){

	cout << " ----------- std::regex --------- \n";

	int fd = open("/tmp/large.txt", O_RDONLY);
	struct stat finfo;
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);

	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);

	buff[size - 1] = '\0';

	cout << "entire buff is " << size << std::endl;

//	std::regex ex("(\\d\\d\\d\\d)");
	std::regex ex(R"(\d\d\d\d)");
	std::cmatch m;
	
	int count; 
	while(1){
		count = 0;
		ptime begin = boost::posix_time::microsec_clock::local_time();

		std::regex_search(buff, m, ex); //only can find one match
		for(unsigned int i = 0; i < m.size(); i++){
			count ++;
		}

		ptime end  = boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration diff = end - begin;
		double interval_sec = (double) diff.total_milliseconds() / 1000;
		std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;
		std::cout << count << std::endl;
	}

}

//hym: regex can find all matches -> using const_iterator
int main11(){
        int fd = open("/tmp/large.txt", O_RDONLY);
        struct stat finfo;
        if(fd < 0)
                std::cout << "open fail!" << std::endl;
        fstat(fd, &finfo);

        char *buff = (char *) malloc(finfo.st_size);
        int size = pread(fd, buff, finfo.st_size, 0); 

        std::cout << size << std::endl;

        //std::regex ex("(\\d\\d\\d\\d)");
        std::regex ex(R"(\d\d\d\d)");
	std::cmatch m;

        std::string str(buff);
        //std::smatch match_result;
        std::smatch res;

        int count; 
        while(1){
                count = 0;
                ptime begin = boost::posix_time::microsec_clock::local_time();
                std::string::const_iterator searchStart(str.cbegin());
                while( std::regex_search(searchStart, str.cend(), res, ex)){
                        count ++;                   
                        searchStart += res.position() + res.length();
                        //std::cout << "count is " << count << std::endl;
                }   
                    
                ptime end  = boost::posix_time::microsec_clock::local_time();
                boost::posix_time::time_duration diff = end - begin;
                double interval_sec = (double) diff.total_milliseconds() / 1000;
                std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;
                std::cout << "count in std::regex is: " << count << std::endl;
        }   

}


//hym: regex can find all matches -> using sregex_iterator
int main12(){
        int fd = open("/tmp/large.txt", O_RDONLY);
        struct stat finfo;
        if(fd < 0)
                std::cout << "open fail!" << std::endl;
        fstat(fd, &finfo);

        char *buff = (char *) malloc(finfo.st_size);
        int size = pread(fd, buff, finfo.st_size, 0); 

        std::cout << size << std::endl;

        //std::regex ex("(\\d\\d\\d\\d)");
        std::regex ex(R"(\d\d\d\d)");
	std::cmatch m;

        std::string str(buff);
        //std::smatch match_result;
        std::smatch res;

        int count; 
        while(1){
                count = 0;
                ptime begin = boost::posix_time::microsec_clock::local_time();
                std::string::const_iterator searchStart(str.cbegin());
  /*              while( std::regex_search(searchStart, str.cend(), res, ex)){
                        count ++;                   
                        searchStart += res.position() + res.length();
                        //std::cout << "count is " << count << std::endl;
                }   
  */                
  		std::sregex_iterator it(str.begin(), str.end(), ex);
		std::sregex_iterator reg_end;
		for (; it != reg_end; ++it) {
			count ++;
		}

                ptime end  = boost::posix_time::microsec_clock::local_time();
                boost::posix_time::time_duration diff = end - begin;
                double interval_sec = (double) diff.total_milliseconds() / 1000;
                std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;
                std::cout << "count in std::regex is: " << count << std::endl;
        }   

}

/*
 *hym: test re2 latency, and compare the result with
 *ref:http://sljit.sourceforge.net/regex_perf.html
 *result: our implementation seems worse than there mplementation
 *         https://docs.google.com/spreadsheets/d/1l6rgwdKpaum4wdMlAU7KWxPhPtyg9z38iDOmGVqb4kk/edit#gid=0
 *reason: we c++, they c
 */
int main13(){
	cout << " ----------- RE2 Latency --------- \n";

	int fd = open("/tmp/bench/mtent12.txt", O_RDONLY);
	//int fd = open("/tmp/small.txt", O_RDONLY);
	struct stat finfo;
	if(fd < 0)
		std::cout << "open fail!" << std::endl;
	fstat(fd, &finfo);

	char *buff = (char *) malloc(finfo.st_size);
	int size = pread(fd, buff, finfo.st_size, 0);

	std::cout << size << std::endl;
	std::string matched;
	//RE2 re("(\\d\\d\\d\\d)");	
	//RE2 re("(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)");
	//RE2 re("(Twain)");
	//RE2 re("(([A-Za-z]awyer|[A-Za-z]inn)\s)");
	//RE2 re(("(?i)Twain)");
	//RE2 re("([a-z]shing)");
	//RE2 re("([\"'][^\"']{0,30}[?!\\.][\"'])");
	std::string re[] = {
		"(Twain)",
		"((?i)Twain)",
		"([a-z]shing)",
		"(Huck[a-zA-Z]+|Saw[a-zA-Z]+)",
		"(\\b\\w+nn\\b)",
		"([a-q][^u-z]{13}x)",
		"(Tom|Sawyer|Huckleberry|Finn)",
		"((?i)Tom|Sawyer|Huckleberry|Finn)",
		"(.{0,2}(Tom|Sawyer|Huckleberry|Finn))",
		"(.{2,4}(Tom|Sawyer|Huckleberry|Finn))",
		"(Tom.{10,25}river|river.{10,25}Tom)",
		"([a-zA-Z]+ing)",
		"(\\s[a-zA-Z]{0,12}ing\\s)",
		"(([A-Za-z]awyer|[A-Za-z]inn)\\s)",
		"([\"'][^\"']{0,30}[?!\\.][\"'])"
	};
	int count;
	int i;
	while(1){
	   for(i = 0; i < 15; i++ ){
		count = 0;
		ptime begin = boost::posix_time::microsec_clock::local_time();

		re2::StringPiece input(buff);
		while(RE2::FindAndConsume(&input, (RE2)re[i], &matched)) {
			count ++;
			//std::cout << matched << std::endl;
		}
		ptime end  = boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration diff = end - begin;
		double interval_sec = (double) diff.total_milliseconds() / 1000;		
		std::cout << re[i] <<":\n  time: " << diff.total_milliseconds() << "ms" << std::endl;
		//std::cout << "Bandwidth is: " << size / interval_sec / (1024*1024*1024) << "GB/s" << std::endl;
		std::cout << "  matched: " <<count << std::endl;
	
	   } //for
	   std::cout << endl;
	}// while(1)
	return 0;
	
}

int main()
{
    std::string target("baaaby");
    std::smatch sm;

    std::regex re1("a(a)*b");
    std::regex_search(target, sm, re1);
    std::cout << "entire match: " << sm[0] << '\n'
              << "submatch #1: " << sm[1] << '\n';

    std::regex re2("a(a*)b");
    std::regex_search(target, sm, re2);
    std::cout << "entire match: " << sm[0] << '\n'
              << "submatch #1: " << sm[1] << '\n';

    {
    	/* std::regex only returns 1st match */
			std::regex ex3(R"(\d\d\d\d)");
			std::string target("1245 9870");
			std::regex_search(target, sm, ex3);
			std::cout << "entire match: " << sm[0] << '\n'
								<< "submatch #1: " << sm[1] << '\n';
    }

    /* xzl: return all matches
     * ref: http://www.informit.com/articles/article.aspx?p=2064649&seqNum=5
     * */
    {
    	printf("------------ xzl: test sregex_iterator\n");

    	/* std::regex only returns 1st match */
			std::regex ex3(R"(\d\d\d\d)");
			std::string target("1245 9870");

			sregex_iterator it(target.begin(), target.end(), ex3);
			sregex_iterator reg_end;

			for (; it != reg_end; ++it) {
			     std::cout << "Substring found: ";
			     std::cout << it->str() << ", Position: ";
			     std::cout << it->position() << std::endl;
			}
    }

    /* xzl: same as above, but using char * interface */
    {
    	printf("------------ xzl: test cregex_iterator\n");

			std::regex ex3(R"(\d\d\d\d)");
			const char *target = "1245 9870";

			cregex_iterator it(target, target + strlen(target), ex3);
			cregex_iterator reg_end;

			for (; it != reg_end; ++it) {
			     std::cout << "Substring found: ";
			     std::cout << it->str() << ", Position: ";
			     std::cout << it->position() << std::endl;
			}
    }


    //const char *str = "http://pge.rastko.net";
    //const char *str = "http://localhost.com/path\\?hue\\=br\\#cool";
    //const char *str = "http://gutenberg.net.au";
    const char *str = "hello world http://davinci.marc.gatech.edu/~tad/dennis/no-cense.htm hello worl";

    //main1(str);
    //main2(str);
    //main3();
    //main4();
    //main5();
    //main6(); //find all valid url
    //main7();
    //main8(); //strstr()
    main9(); //re2()
    //main10(); //std::regex -> only find one match
    //main11(); //std::regex -> find all matches: using const_iterator
    //main12(); //std::regex -> find all matches: using sregex_iterator 
//    main13();
}


