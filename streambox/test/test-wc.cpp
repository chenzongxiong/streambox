#include "core/Pipeline.h"

/* No need to include the transform headers.
 * Since the -Evaluator.h file already included the corresponding
 * transform header.
 */
//#include "Unbounded.h"
#include "Source/UnboundedInMemEvaluator.h"
//#include "WordCountMapper.h"
#include "Mapper/WordCountMapperEvaluator.h"
//#include "WinGBK.h"
#include "Win/WinGBKEvaluator.h"
//#include "WindowKeyedReducer.h"
//#include "WindowKeyedReducerEvaluator.h"
#include "WinKeyReducer/WinKeyReducerEval.h"
#include "Sink/WindowsBundleSinkEvaluator.h"

#include "test-common.h"

/*
 *      UnboundedInMem
 *            | BundleT<string_range>
 *            V
 *       WordCountMapper
 *            | BundleT<KVPair<string,long>>
 *            V
 *          WinGBK
 *            | WindowsKeyedBundle<KVPair<string,long>>
 *            V
 *      WindowKeyedReducer (stateful)
 *            | WindowsBundle<KVPair<string,long>>
 *            V
 *				WindowsBundle
 *
 *	 where @BundleT can be specified below.
 */


/* -- define the bundle type we use -- */
template<class T>
using BundleT = RecordBundle<T>;
//using BundleT = RecordBitmapBundle<T>;

/* default config. can override on cmdline options */
#ifdef DEBUG
pipeline_config config = {
		.records_per_interval = 1000,
		.target_tput = (20 * 1000),
		.record_size = 1000,
		.input_file = "./inputdata.dat",
//		.cores = std::thread::hardware_concurrency() - 1,
};
#else
pipeline_config config = {
		.records_per_interval = (1000 * 1000),
		.target_tput = (3700 * 1000),
		.record_size = 100,
		.input_file = "./inputdata.dat",
//		.cores = std::thread::hardware_concurrency() - 1,
};
#endif

void testWordCount()
{

	UnboundedInMem<string_range, BundleT>
	unbound("unbounded-inmem",
		config.input_file.c_str(),
		config.records_per_interval , 	/* records per wm interval */
		config.target_tput, 		/* target tput (rec/sec) */
		config.record_size,
		0  				/* session_gap_ms */
	);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	WordCountMapper<string_range, pair<creek::string,long>, BundleT> mapper ("[wc-mapper]");

	WinGBK<pair<creek::string, long>, BundleT,
		WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(1));
	WinKeyReducer<pair<creek::string,long>,  /* pair in */
		WinKeyFragLocal_Std, WinKeyFrag_Std, /* kv d/s */
		pair<creek::string,long>,  /* pair out */
		WindowsBundle /* bundle out */
			> reducer("[reducer]");
	WindowsBundleSink<pair<creek::string,long>> sink("sink");

	connect_transform(unbound, mapper);
	connect_transform(mapper, wgbk);
	connect_transform(wgbk, reducer);
	connect_transform(reducer, sink);
	//  connect_transform(wgbk, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

//#include <jemalloc/jemalloc.h>

int main(int ac, char *av[])
{
	parse_options(ac, av, &config);

	//  malloc_stats_print(NULL, NULL, NULL);  // can test jemalloc
	print_config();
	testWordCount();
}

