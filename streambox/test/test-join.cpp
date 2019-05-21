/*
#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<typeinfo>

#include "Create.h"
#include "Read.h"
*/
#include "core/Pipeline.h"
//#include "TransformEvaluator.h"
//#include "TransformEvaluate.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Mapper/SimpleMapperEvaluator.h"
#include "Join/JoinEvaluator1.h"
#include "Sink/RecordBitmapBundleSinkEvaluator.h"
#include <tuple>
#include "test-common.h" /*
 *                    +--- SimpleMapper -- (KVPair)--\
 * Unbound -(long)->  +                             Join --> (KVPai
 *                    +--- SimpleMapper -- (KVPair)--/
 *
 */
 /*
 template<class T>
 auto operator<<(std::ostream& os, const T& t) -> decltype(t.print(os), os)
 {
      t.print(os);
      return os;
 }
*/
template<class T>
using BundleT = RecordBitmapBundle<T>;

#ifdef DEBUG
pipeline_config config = {
		.records_per_interval = 100,
		.target_tput = (20 * 100),
		.record_size = sizeof(long),
		.input_file = "/ssd/test-digit.txt",
};
#else
pipeline_config config = {
		.records_per_interval = 1000 * 1000,
		//.target_tput = (10 * 1000 * 1000),
		.target_tput =  2000 * 1000,   //k2: 2000 should be okay
		.record_size = sizeof(long),
		//.input_file = "/ssd/1g.txt",
		.input_file = "/ssd/test-digit.txt",
};
#endif

int main(int ac, char *av[])
{
	parse_options(ac, av, &config);
	print_config();
	// create a new Bounded transform
	//UnboundedInMem<long, BundleT> unbound ("[unbounded]", 50, 0);
	UnboundedInMem<long, BundleT> unbound ("[unbounded]",
  				config.input_file.c_str(),
  				config.records_per_interval , /* records per wm interval */
  				config.target_tput,	      /* target tput (rec/sec) */
  				config.record_size,
				0 //XXX session_gap_ms = ????
  			);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);
/*
	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	SimpleMapper<long, pair<long, long>> mapper0 ("[mapper0]");
	SimpleMapper<long, pair<long, long>> mapper1 ("[mapper1]");

	connect_transform(unbound, mapper0);
  	connect_transform(unbound, mapper1);
*/
	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	SimpleMapper<long, pair<long, long>> mapper0 ("[mapper0]");
	auto c11 = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";
	mapper0.set_side_info(1); //left side 

	SimpleMapper<long, pair<long, long>> mapper1 ("[mapper1]");
	auto c12 = dynamic_cast<PCollection *>(o2[1]->apply1(&mapper1));
	c12->_name = "mapper1_out";
	mapper1.set_side_info(2); //right side


	Join<pair<long,long>> join("[join]", seconds(1));
	auto c2 = dynamic_cast<PCollection *>(c11->apply1(&join));
	auto c3 = dynamic_cast<PCollection *>(c12->apply1(&join));
	c2 = c2;
	assert(c2 == c3); // same output.
	join.set_side_info(3);

	RecordBitmapBundleSink<pair<long, vector<long>>> sink("[sink]");
	auto c4 = dynamic_cast<PCollection *>(c3->apply1(&sink));
	c4 = c4;
	assert(c4); //just for fixing warning
	sink.set_side_info(4);

	RecordBitmapBundleSink<pair<long, vector<long>>> sink1("[sink1]");
	auto c5 = dynamic_cast<PCollection *>(c4->apply1(&sink1));
	c5 = c5;
	assert(c5); //just for fixing warning
	//assert(false && "set RecordBitmapBundleSink side_info to 4!!");
	sink1.set_side_info(5);
	
	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);

	return 0;
}

