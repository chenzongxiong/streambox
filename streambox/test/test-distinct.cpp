#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Mapper/WordCountMapperEvaluator.h"
#include "Win/FixedWindowIntoEvaluator.h"
//#include "WindowedSumEvaluator.h"
#include "WinSum/WinSum_mergeset.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
#include "test-common.h"
/*
 *      UnboundedInMem
 *           |  BundleT<string_range>
 *           V
 *     	 WordCountMapper (line mode)
 *           |  BundleT<KVstring>
 *           V
 *      FixedWindowInto
 *           |  WindowsBundle<creek::string>
 *           V
 *       WindowedSum (stateful)
 *           |  RecordBundle<unordered_set<creek::string>>
 *           V
 *      RecordBundleSink
 *
 * where @BundleT can be specified below.
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
	.record_size = 100,  /* string range, will be split into urlid */
	.input_file = "/ssd/urlid2.txt"
};
#else
pipeline_config config = {
	.records_per_interval = (1000 * 1000),
	.target_tput = (2000 * 1000),
	.record_size = 100,  /* string range, will be split into urlid */
	.input_file = "/ssd/train-out.txt"
};
#endif


void testDist()
{
	UnboundedInMem<string_range, BundleT>
		unbound("unbounded-inmem",
				config.input_file.c_str(),
				config.records_per_interval , 	/* records per wm interval */
				config.target_tput, 		/* target tput (rec/sec) */
				config.record_size,  		/* used to estimate file size only */
				0  				/* session_gap_ms */
		       );

#ifdef USE_FOLLY_HASHMAP
	using Set = folly::AtomicHashMap<uint64_t, uint64_t>;
#elif defined(USE_CUCKOO_HASHMAP)
	//  using Set = cuckoohash_map<creek::string, int>;
	using Set = cuckoohash_map<creek::string, int, CityHasher<creek::string>>;
#elif defined(USE_TBB_DS)
	//	  using Set = tbb::concurrent_unordered_set<creek::string, creek::FBStringHash>;
	//  using Set = cuckoohash_map<creek::string, ValueContainerT>;
	//  using Set = tbb::concurrent_unordered_set<creek::string>;
	using Set = creek_set_array::SetArray;
#else
	using Set = unodered_set<creek::string>;
#endif

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	WordCountMapper<string_range, creek::string, BundleT, wc_mapper::LINE>
		mapper ("[split-mapper]");
	FixedWindowInto<creek::string, BundleT> fwi ("window", seconds(1));
	WinSum_mergeset<creek::string, shared_ptr<Set>> agg ("agg",
			1 /* >1 for sliding window, slow right now */);
	RecordBundleSink<shared_ptr<Set>> sink("sink");

	connect_transform(unbound, mapper);
	connect_transform(mapper, fwi);
	connect_transform(fwi, agg);
	connect_transform(agg, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}


int main(int ac, char *av[])
{
	parse_options(ac, av, &config);

	print_config();
	testDist();
}
