#include "core/Pipeline.h"
#include "Source/Unbounded.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Win/WinGBKEvaluator.h"
#include "WinKeyReducer/WinKeyReducerEval.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
#include <fstream>
#include "LinearRoad/LinearRoadParser.h"
#include "LinearRoad/LinearRoadFilter.h"
#include "LinearRoad/LinearRoadTollNotification.h"

#include "test-common.h"

typedef uint64_t Timestamp;
using NanoSeconds = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;

template<class T>
using BundleT = RecordBundle<T>;

pipeline_config config = {
  .records_per_interval = 2500000,
  .target_tput = std::numeric_limits<uint64_t>::max(),
  // .target_tput = 10000,
  // set record_size to be 128, bug about mutex can be reproducible.
  .record_size = sizeof(LinearRoadRecord),
  .input_file = "/home/zxchen/cardatapoints.out0",
  // .input_file = "./crashdata.dat",
  .cores = 7,
};

Timestamp getTimestamp()
{
  return std::chrono::duration_cast<NanoSeconds>(
    Clock::now().time_since_epoch())
    .count();
}

void testInput() {
  UnboundedInMem<string_range, BundleT>
    unbound("[unbounded-inmem]",
            config.input_file.c_str(),
            config.records_per_interval , 	/* records per wm interval */
            config.target_tput, 		/* target tput (rec/sec) */
            config.record_size,
            0,  				/* session_gap_ms */
            TXT_LINES
      );




  // create a new pipeline
	Pipeline* p = Pipeline::create(NULL);
	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	// parse dataset line by line into LinearRoadParser
	LinearRoadParser<string_range, LinearRoadRecord, RecordBundle> parser("[linear_road_parser]");
    LinearRoadFilter<LinearRoadRecord, LinearRoadRecord, RecordBundle> filter("[linear_road_filter]");
    // LinearRoadTollNotification<LinearRoadRecord, LinearRoadTollRecord, RecordBundle> toll("[linear_road_toll]");
    LinearRoadTollNotification<LinearRoadRecord, pair<uint64_t, long>, RecordBundle> mapper("[linear_road_toll]");
    // create a 1 seconds window, group by key
	WinGBK<pair<uint64_t, long>, BundleT,
           WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(1));

    // reduce aggregation
    WinKeyReducer<pair<uint64_t, long>,  /* pair in */
            WinKeyFragLocal_Std, WinKeyFrag_Std, /* kv d/s */
            pair<uint64_t, long>,  /* pair out */
            WindowsBundle /* bundle out */
            > reducer("[reducer]");

    WindowsBundleSink<pair<uint64_t, long>> sink("[sink]");

	connect_transform(unbound, parser);
    connect_transform(parser, filter);
    connect_transform(filter, mapper);
    connect_transform(mapper, wgbk);
    connect_transform(wgbk, reducer);
    connect_transform(reducer, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);

}

int main(int ac, char *av[]) {
  parse_options(ac, av, &config);
  print_config();

  Timestamp begin = getTimestamp();
  testInput();
  Timestamp end = getTimestamp();


  double elapsed_time = double(end - begin) / (1024 * 1024 * 1024);
  std::cout << "time elapses: " << elapsed_time << std::endl;
}
