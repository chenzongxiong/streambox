#include <re2/re2.h>
#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Mapper/GrepMapperEvaluator.h"
#include "Sink/RecordBitmapBundleSinkEvaluator.h"
#include "test-common.h"
/*
 *      UnboundedInMem
 *           |  RecordBundle/RecordBitmapBundle<string_range>
 *           V
 *        GrepMapper
 *           |  RecordBundle/RecordBitmapBundle<string>
 *           V
 *       RecordBitmapBundleSink<string>
 *
 *  NB: this needs more work. GrepMapper() can only output RecordBitmapBundle.
 *  Consider using testWnidowedGrep()
 */

/* -- define the bundle type we use -- */
template<class T>
using BundleT = RecordBundle<T>;
//using BundleT = RecordBitmapBundle<T>;
#ifdef DEBUG
pipeline_config config = {
		.records_per_interval = 1000,
		.target_tput = (20 * 1000),
		.record_size = 1000,
		.input_file = "/ssd/1g.txt",
};
#else
pipeline_config config = {
		.records_per_interval = 4096,
		.target_tput = (7 * 1000 * 1000),
		.record_size = 1024,
		.input_file = "/ssd/1g.txt",
};
#endif

void testGrep()
{	
	UnboundedInMem<string_range, BundleT>
	unbound("unbounded-inmem",
		config.input_file.c_str(),
		config.records_per_interval, 	/* records per wm interval */
		config.target_tput, 	 	/* target tput (rec/sec) */
		config.record_size,
		0  				/* session_gap_ms */
	);

  // create a new pipeline
  Pipeline* p = Pipeline::create(NULL);

  PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
  unbound_output->_name = "unbound_output";

  GrepMapper<string_range, creek::string, RecordBundle>
  			mapper (R"(\d\d\d\d)", "[grep-mapper]"); // for std::regex

  RecordBundleSink<creek::string> sink("sink");

  connect_transform(unbound, mapper);
  connect_transform(mapper, sink);

  EvaluationBundleContext eval(1, config.cores);
  eval.runSimple(p);
}

int main(int ac, char *av[])
{
	parse_options(ac, av, &config);
	print_config();

#if 0
	RE2::Options opt;
	opt.set_max_mem(256 << 20);  /* does not help */
#endif

	testGrep();
}
