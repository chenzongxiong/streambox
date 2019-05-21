/*
 * test-common.h
 *
 *  Created on: Oct 7, 2016
 *      Author: xzl
 */

#ifndef TEST_TEST_COMMON_H_
#define TEST_TEST_COMMON_H_

//#define xstr(s) str(s)
//#define str(s) #s
#define xstr(s) #s //hym: avoid conflict with .str()

/* expect: Pipeline* p = Pipeline::create(NULL) */
#define source_transform(T) \
	  auto T##_output = dynamic_cast<PCollection *>(p->apply1(&T)); \
	  T##_output->_name = xstr(T##_out)
	  //T##_output->_name = str(T##_out)

/* only works for apply1() */
#define connect_transform(T1, T2) \
    auto T2##_output = dynamic_cast<PCollection *>(T1##_output->apply1(&T2)); \
	T2##_output -> _name = xstr(T2##_out)
	//T2##_output -> _name = str(T2##_out)

#define source_transform_1to2(T) \
	  vector<PCollection *> T##_output = p->apply2(&T); \
	  T##_output[0]->_name = xstr(T##_out0);  \
		T##_output[1]->_name = xstr(T##_out1)
	  //T##_output->_name = str(T##_out)

#define connect_transform_1to2(S, O0, O1) \
	auto O0##_output = dynamic_cast<PCollection *>(S##_output[0]->apply1(&O0)); \
	O0##_output -> _name = xstr(O0##_out); \
	auto O1##_output = dynamic_cast<PCollection *>(S##_output[1]->apply1(&O1)); \
	O1##_output -> _name = xstr(O1##_out) \

#define connect_transform_2to1(T1, T2, S) \
    auto S##_output = dynamic_cast<PCollection *>(T1##_output->apply1(&S)); \
    S##_output -> _name = xstr(S##_out); \
    auto S##_output1 = dynamic_cast<PCollection *>(T2##_output->apply1(&S)); \
    xzl_assert(S##_output == S##_output1);  /* same output */

void print_config(void);

using namespace std;
struct pipeline_config {
	unsigned long records_per_interval;
	/// unsigned long target_tput; /* krec/s */
    uint64_t target_tput;
	unsigned long record_size; /* e.g. string range etc. */
    unsigned long window_size;
    unsigned long window_count;
    unsigned long campaigns;
	std::string input_file;
	unsigned long cores;
};

void parse_options(int ac, char *av[], pipeline_config* config);

#include "linux-sizes.h"

#endif /* TEST_TEST_COMMON_H_ */
