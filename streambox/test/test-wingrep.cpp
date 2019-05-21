#include "config.h"
#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Win/FixedWindowIntoEvaluator.h"
#include "Mapper/WindowedGrepMapperEvaluator.h"
//#include "WindowedSumEvaluator.h"
#include "WinSum/WinSum_mergevector.h"
#include "WinSum/WinSumEval.h"
#include "Sink/RecordBundleSinkEvaluator.h"
#include "test-common.h"
/*
 *      UnboundedInMem
 *           |  BundleT<string_range>
 *           V
 *      FixedWindowInto
 *           |  WindowsBundle<string_range>
 *           V
 *     WindowedGrepMapper
 *           |  WindowsBundle<string>
 *           V
 *       WindowedSum (stateful)
 *           |  RecordBundle<vector<string>>
 *           V
 *     RecordBundleSink
 *
 * @BundleT can be specified below.
 */


/* -- define the bundle type we use -- */
template<class T>
using BundleT = RecordBundle<T>;
//using BundleT = RecordBitmapBundle<T>;

#ifdef DEBUG
pipeline_config config = {
		.records_per_interval = SZ_1M,
		.target_tput = SZ_500K,
		.record_size = 100,
		.input_file = "/ssd/9g.txt",
};
#else
pipeline_config config = {
		.records_per_interval = SZ_1M,
		.target_tput = SZ_16M,
		.record_size = 1000,
		.input_file = "/ssd/9g.txt",
};
#endif

void testWindowedGrep()
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
	unbound_output->_name = "unbound_output";

	FixedWindowInto<string_range, BundleT> fwi ("window", seconds(1));
	//auto c2 = dynamic_cast<PCollection *>(o2->apply1(&fwi));
	//c2->_name = "win_out";
	connect_transform(unbound, fwi);

	/* four digits as year */
	//WindowedGrepMapper<> grep (R"(\d\d\d\d)", "mapper"); // std::regex
	//WindowedGrepMapper<> grep (R"(http://)", "mapper");		// re2

	/* url */
	WindowedGrepMapper<> grep (
			"(((http[s]{0,1}|ftp)://)[a-zA-Z0-9\\.\\-]+\\.([a-z]|[A-Z]|[0-9]|[/.]|[~]|[\\-])*)",
			"[grep-mapper]"
			);

#if 0
	WindowedGrepMapper<> grep (
			//R"([-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*))",
			//"(http|https)://([^/ :]+):?([^/ ]*)(/?[^ #?]*)\\x3f?([^ #]*)#?([^ ]*)",
			"(((http[s]{0,1}|ftp)://)[a-zA-Z0-9\\.\\-]+\\.([a-z]|[A-Z]|[0-9]|[/.]|[~]|[\\-])*)",
			"mapper"); /* will match empty string? */

	// need "extended"
	//	      R"(^(([^:\/?#]+):)?(//([^\/?#]*))?([^?#]*)(\?([^#]*))?(#(.*))?)",
#endif

	connect_transform(fwi, grep);

	//#if USE_TBB_DS
	//	  using Vector = tbb::concurrent_vector<creek::string>;
	//#else
	//	  using Vector = vector<creek::string>;
	//#endif

	using Vector = creek::concurrent_vector<creek::string>;

	/* can't afford copying aggregation results (as vector) back & forth */
	WinSum_mergevector<creek::string, shared_ptr<Vector>> agg ("agg",
			1 /* >1 for sliding window, slow right now */);
	RecordBundleSink<shared_ptr<Vector>> sink("sink");

	connect_transform(grep, agg);
	connect_transform(agg, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

int main(int ac, char *av[])
{
	parse_options(ac, av, &config);
	print_config();

	testWindowedGrep();
}
