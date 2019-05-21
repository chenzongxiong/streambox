#include "core/Pipeline.h"
#include "Source/NetworkLatencySource.h"
#include "Values.h"
#include "test-common.h" 
#include "Source/NetworkLatencySourceEvaluator.h"
#include "Win/WinGBKEvaluator.h"
//#include "WindowKeyedReducer.h"
//#include "WindowKeyedReducerEvaluator.h"
#include "WinKeyReducer/WinKeyReducerEval.h"
#include "Sink/WindowsBundleSinkEvaluator.h"

template<class T>
using BundleT = RecordBundle<T>;

#ifdef DEBUG
pipeline_config config = {
		.records_per_interval = 10 * 1000,
		.target_tput = (200 * 1000),
		.record_size = sizeof(struct srcdst_rtt),
		.input_file = "/ssd/region_Raw_PingmeshData.result",
};
#else
pipeline_config config = {
		.records_per_interval = 500 * 1000,
		.target_tput = 800 * 1000,
		.record_size = sizeof(struct srcdst_rtt),
		.input_file = "/ssd/region_Raw_PingmeshData.result",
};
#endif

int main(int ac, char *av[]){

	parse_options(ac, av, &config);

	print_config();

	NetworkLatencySource netsource(
				"[netsource]",
  				config.input_file.c_str(),
  				config.records_per_interval , /* records per wm interval */
  				config.target_tput,	      /* target tput (rec/sec) */
				config.record_size,
				0 //XXX session_gap_ms = ????
			);

	Pipeline* p = Pipeline::create(NULL);

	PCollection *netsource_output = dynamic_cast<PCollection *>(p->apply1(&netsource));
	netsource_output->_name = "src_out";

	WinGBK<pair<creek::ippair, long>, BundleT,
		WinKeyFragLocal_Simple> wgbk ("[wingbk]", seconds(1));
	WinKeyReducer<pair<creek::ippair, long>,
		WinKeyFragLocal_Simple,
		WinKeyFrag_SimpleMT, pair<creek::ippair, long>, WindowsBundle> reducer("[reducer]");

	/* can compile but shoudln't making WindowKeyedFragmentSimpleMT input: can't easily
	 * iterate over it
	 */
	//	WinGBK<pair<creek::ippair, long>, BundleT,
	//		WindowKeyedFragmentUnsafeSimple> wgbk ("[wingbk]", seconds(1));
	//	WindowKeyedReducer<pair<creek::ippair, long>,
	//		WindowKeyedFragmentSimpleMT,
	//		WindowKeyedFragmentSimpleMT> reducer("[reducer]");

	WindowsBundleSink<pair<creek::ippair, long>> sink("sink");
	
	connect_transform(netsource, wgbk);
	connect_transform(wgbk, reducer);
	connect_transform(reducer, sink); 
	//  connect_transform(wgbk, sink); 

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
	return 0;
}
