/*
 * json-parser.cpp
 *
 *  Created on: Jan 23, 2017
 *      Author: xzl
 *
 *      g++-5 json-parser.cpp -o json-parser.bin -std=c++11 -g
 *
 *			g++-5  ../measure.c json-parser.cpp -o json-parser.bin -std=c++11 -g -O2
 *
 *      https://github.com/nlohmann/json/tree/master
 */

#include <sstream>
#include <iostream>
#include <fstream>      // std::ifstream
#include </ssd/local/include/json.hpp>		/* from the modern c++ json parser, /ssd/local/include */
#include <boost/progress.hpp> /* progress bar */
//#include "json.hpp"
#include "../measure.h"

// for convenience
using json = nlohmann::json;
using namespace std;

void single_json() {
	// read a JSON file
//	std::ifstream i("/ssd/twit_data/stream_122716.json");
	std::ifstream i("/tmp/test.json");
	json j;
	i >> j;

	// write prettified JSON to another file
	std::ofstream o("/tmp/pretty.json");
	o << std::setw(4) << j << std::endl;
}

void multi_json_from_file() {

	std::string line;

//	std::ifstream infile("/tmp/test.json");
	std::ifstream infile("/ssd/twit_data/stream_122716.json");
	std::ofstream o("/tmp/pretty.json");

//	auto j3 = json::parse("{ \"happy\": true, \"pi\": 3.141 }");
//	cout << std::setw(4) << j3 << std::endl;

	int cnt = 0;
	while (std::getline(infile, line, '\n')) {
//		cout << "read a newline" << endl;
//		cout << line << endl;
//		cout << "line size:" << line.length() << endl;

//		int a = line.length();
//		printf("line size %d\n", a);

//		cout << "last char:" << line[a - 1];
//		cout << line;

//		o << line;
//		stringstream ss(line);
//		ss >> j;

//		auto j3 = json::parse("{ \"happy\": true, \"pi\": 3.141 }");
		auto j3 = json::parse(line);
//		cout << j.dump(4);

		/* write to file */
		o << std::setw(4) << j3 << std::endl;

		/* write to stdout */
//		cout << std::setw(4) << j3 << std::endl;
		cnt ++;
	}
	printf("%d json objects written\n", cnt);
}


void multi_json_extract(const char *in_fname, const char *out_fname) {

	std::string line;

//	std::ifstream infile("/ssd/twit_data/stream_122716.json");
	std::ifstream infile(in_fname);
//	std::ofstream o("/tmp/pretty.json");
	std::ofstream o(out_fname);

//	auto j3 = json::parse("{ \"happy\": true, \"pi\": 3.141 }");
//	cout << std::setw(4) << j3 << std::endl;

	int cnt = 0, max_cnt = 100*1000;
	boost::progress_display show_progress(max_cnt);
	k2_measure("start");

	while (std::getline(infile, line, '\n')) {
		auto j3 = json::parse(line);

//		cout << j3["timestamp_ms"] << endl;
//		cout << j3["text"] << endl;

//		o << j3["timestamp_ms"] << endl;
		o << j3["text"] << endl;


		cnt ++;
		show_progress += 1;
		if (cnt == max_cnt)
			break;
	}
	k2_measure("end");
	k2_measure_flush();
	printf("%d json objects written\n", cnt);
}



int main(int argc, char **argv)
{
//	single_json();
//	multi_json_from_file();
	multi_json_extract(argv[1], argv[2]);
	return 0;
}
