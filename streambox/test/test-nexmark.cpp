#include <fstream>
#include <map>
#include <atomic>

#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Nexmark/NexmarkRecord.hpp"
#include "Nexmark/NexmarkParser.hpp"
#include "Nexmark/NexmarkFilter.hpp"
#include "Nexmark/NexmarkAggregation.hpp"
#include "Nexmark/NexmarkSink.hpp"
#include "Win/WinGBKEvaluator.h"
// #include "Win/SlidingWindowInto.hpp"
#include "Win/FixedWindowInto.h"
#include "WinKeyReducer/WinKeyReducerEval.h"
#include "Sink/WindowsBundleSinkEvaluator.h"

#include "Values.h"
#include "test-common.h"

typedef uint64_t Timestamp;
using NanoSeconds = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;
using namespace std;

template<class T>
using BundleT = RecordBundle<T>;

pipeline_config config = {
    .records_per_interval = 2500000,
    .target_tput = std::numeric_limits<uint64_t>::max(),
    .record_size = 32,
    .window_size = 2,
    .window_count = 10000,
    .campaigns = 1,
    .input_file = "/home/zxchen/nexmark_test_data.bin",
    .cores = 8,
};


Timestamp getTimestamp() {
    return std::chrono::duration_cast<NanoSeconds>(
        Clock::now().time_since_epoch())
        .count();
}

using KVPair = pair<uint64_t, uint64_t>;

void testInput() {

    UnboundedInMem<string_range, BundleT>
        unbound("[unbounded-inmem]",
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

    NexmarkParser<string_range, NexmarkRecord, RecordBundle> parser("[nexmark_parser]");
    NexmarkFilter<NexmarkRecord, NexmarkRecord, RecordBundle> filter("[nexmark_filter]");
    NexmarkAggregation<NexmarkRecord, KVPair, BundleT> mapper("[nexmark_aggregation]");

	WinGBK<KVPair, BundleT, WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(config.window_size));

    // reduce aggregation
    WinKeyReducer<KVPair,  /* pair in */
                  WinKeyFragLocal_Std,
                  WinKeyFrag_Std, /* kv d/s */
                  KVPair,  /* pair out */
                  RecordBundle
                  > reducer("[reducer]");


    NexmarkAggregation<KVPair, KVPair, BundleT> mapper2("[nexmark_aggregation2]");
	WinGBK<KVPair, BundleT, WinKeyFragLocal_Std> wgbk2("[wingbk2]", seconds(config.window_size));

    // reduce aggregation
    WinKeyReducer<KVPair,  /* pair in */
                  WinKeyFragLocal_Std,
                  WinKeyFrag_Std, /* kv d/s */
                  KVPair,  /* pair out */
                  WindowsBundle /* bundle out */
                  // RecordBundle
                  > reducer2("[reducer2]");
    WindowsBundleSink<pair<uint64_t, uint64_t>> sink("[sink]");

	connect_transform(unbound, parser);
    connect_transform(parser, filter);
    connect_transform(filter, mapper);
	connect_transform(mapper, wgbk);
    connect_transform(wgbk, reducer);
    connect_transform(reducer, mapper2);
    connect_transform(mapper2, wgbk2);
    connect_transform(wgbk2, reducer2);
    connect_transform(reducer2, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);

}

// stream
// .map("record.price = (record.price * 89) / 100;")
// .select({"auction", "price"})
// .execute();
void testQ1() {
    UnboundedInMem<string_range, BundleT>
        unbound("[unbounded-inmem]",
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

    NexmarkParser<string_range, NexmarkRecord, RecordBundle> parser("[nexmark_parser]");
    NexmarkAggregation<NexmarkRecord, NexmarkRecord, BundleT> mapper("[nexmark_mapper]");
    RecordBundleSink<NexmarkRecord> sink("[sink]");

	connect_transform(unbound, parser);
    connect_transform(parser, mapper);
    connect_transform(mapper, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

// stream
// .filter("auction" == 1007L || "auction" == 2001L || "auction" == 2019L || "auction" ==  )
// .select({"auction", "price"})
// .execute();
void testQ2() {
    UnboundedInMem<string_range, BundleT>
        unbound("[unbounded-inmem]",
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

    NexmarkParser<string_range, NexmarkRecord, RecordBundle> parser("[nexmark_parser]");
    NexmarkFilter<NexmarkRecord, NexmarkRecord, RecordBundle> filter("[nexmark_filter]");
    NexmarkAggregation<NexmarkRecord, NexmarkRecord, BundleT> mapper("[nexmark_mapper]");
    RecordBundleSink<NexmarkRecord> sink("[sink]");

	connect_transform(unbound, parser);
    connect_transform(parser, filter);
    connect_transform(filter, mapper);
    connect_transform(mapper, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}


// Stream
// .groupBy("auction")
// .window(SlidingProcessingTimeWindow(windowSizeSec, windowPeriodSec))
// .aggregate(Count())
// .window(SlidingProcessingTimeWindow(windowSizeSec, windowPeriodSec))
// .aggregate(Max("count"))
// .execute();
void testQ3() {
    UnboundedInMem<string_range, BundleT>
        unbound("[unbounded-inmem]",
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

    NexmarkParser<string_range, NexmarkRecord, RecordBundle> parser("[nexmark_parser]");

    NexmarkAggregation<NexmarkRecord, KVPair, BundleT> mapper("[nexmark_mapper]");

	WinGBK<KVPair, BundleT, WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(config.window_size));

    // reduce aggregation
    WinKeyReducer<KVPair,  /* pair in */
                  WinKeyFragLocal_Std,
                  WinKeyFrag_Std, /* kv d/s */
                  KVPair,  /* pair out */
                  RecordBundle
                  > reducer("[reducer]");

    NexmarkAggregation<KVPair, KVPair, BundleT> mapper2("[nexmark_aggregation2]");
    WinGBK<KVPair, BundleT, WinKeyFragLocal_Std> wgbk2("[wingbk2]", seconds(config.window_size));

    // reduce aggregation
    WinKeyReducer<KVPair,  /* pair in */
                  WinKeyFragLocal_Std,
                  WinKeyFrag_Std, /* kv d/s */
                  KVPair,  /* pair out */
                  WindowsBundle /* bundle out */
                  > reducer2("[reducer2]");
    WindowsBundleSink<KVPair> sink("[sink]");

	connect_transform(unbound, parser);
    connect_transform(parser, mapper);
	connect_transform(mapper, wgbk);
    connect_transform(wgbk, reducer);
    connect_transform(reducer, mapper2);
    connect_transform(mapper2, wgbk2);
    connect_transform(wgbk2, reducer2);
    connect_transform(reducer2, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

// stream
// .window(SlidingProcessingTimeWindow(windowSizeSec, windowPeriodSec))
// .aggregate(Max("price"))
// .execute();
void testQ4() {
    UnboundedInMem<string_range, BundleT>
        unbound("[unbounded-inmem]",
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

    NexmarkParser<string_range, NexmarkRecord, RecordBundle> parser("[nexmark_parser]");
    // SlidingWindowInto<NexmarkRecord, BundleT> sliding_win("sliding_win",
    //                                                       seconds(config.window_size),
    //                                                       seconds(config.window_size));

    FixedWindowInto<NexmarkRecord, BundleT> fixed_win("fixed_win",
                                                      seconds(config.window_size));
    // TODO: implemente MAX aggregation
    // WinGBK<NexmarkRecord, BundleT, WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(config.window_size));

    // NexmarkAggregation<NexmarkRecord, NexmarkRecord, BundleT> mapper("[nexmark_mapper]");

    // WinKeyReducer<NexmarkRecord,  /* pair in */
    //               WinKeyFragLocal_Std,
    //               WinKeyFrag_Std, /* kv d/s */
    //               NexmarkRecord,  /* pair out */
    //               // WindowsBundle /* bundle out */
    //               RecordBundle
    //               > reducer("[nexmark_reducer]");

    // RecordBundleSink<NexmarkRecord> sink("[sink]");

	connect_transform(unbound, parser);
    // connect_transform(parser, sliding_win);
    // connect_transform(filter, mapper);
    // connect_transform(mapper, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}


ofstream statistics_file;
ofstream latency_file;
// boost::posix_time::ptime zxchen_start; //global variable
Timestamp latencyStart;
std::atomic<bool> first (true);

// int zxchen_first = 1;
// int first_end = 1;

std::map<Window, ptime, Window> window_keeper;
// ptime czx_start = min_date_time;

// #ifdef MEASURE_40M_TUPLES
std::string papi_conf_file = "../papi/papi_conf_global.cfg";
// on gpu-1: Cycles, Instructions, Branches, Branches Mispred, L1D, L1I, dTLB, iTLB
// TLB_PRESET
// BRANCH_PRESET
// DCACHE_MISS_PRESET
// ICACHE_MISS_PRESET
// TOTAL_INSTR_CYCLE_PRESET
// std::string papi_seq = "L1_CACHE_MISS_PRESET";
std::string papi_seq = "L2_CACHE_MISS_PRESET";
// std::string papi_seq = "CACHE_MISS_PRESET";
// std::string papi_seq = "TLB_PRESET";
// std::string papi_seq = "BRANCH_PRESET";
// std::string papi_seq = "TOTAL_INSTR_CYCLE_PRESET";
// #endif



int main(int ac, char *av[]) {
    parse_options(ac, av, &config);
    print_config();

    std::string filename = "./results/test_yahoo_records_" + std::to_string(config.records_per_interval) + "_cores_" + std::to_string(config.cores) + "_window_size_" + std::to_string(config.window_size) + "_campaigns_" + std::to_string(config.campaigns) + ".csv";
    std::string latency_fname = "./results/latency/test_yahoo_records_" + std::to_string(config.records_per_interval) + "_cores_" + std::to_string(config.cores) + "_window_size_" + std::to_string(config.window_size) + "_campaigns_" + std::to_string(config.campaigns) + ".txt";
    statistics_file.open(filename.c_str());

    latency_file.open(latency_fname.c_str());

    Timestamp begin = getTimestamp();
    // testInput();
    // testQ1();
    // testQ2();
    testQ3();
    // testQ4();
    Timestamp end = getTimestamp();


    double elapsed_time = double(end - begin) / (1024 * 1024 * 1024);
    std::cout << "time elapses: " << elapsed_time << std::endl;
}
