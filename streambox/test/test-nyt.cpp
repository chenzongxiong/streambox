#include <map>
#include <fstream>

#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Win/WinGBKEvaluator.h"
#include "WinKeyReducer/WinKeyReducerEval.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
#include "Values.h"
#include "NewYorkTaxi/NYTParser.h"
#include "NewYorkTaxi/NYTParserEvaluator.h"
#include "NewYorkTaxi/NYTQuery.h"
#include "Win/FixedWindowInto.h"
#include "test-common.h"


typedef uint64_t Timestamp;
using NanoSeconds = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;
using namespace std;

template<class T>
using BundleT = RecordBundle<T>;

pipeline_config config = {
    .records_per_interval = 1000000,
    .target_tput = std::numeric_limits<uint64_t>::max(),
    .record_size = sizeof(NYTRecord),
    .window_size = 1,
    .window_count = 10000,
    .campaigns = 1,
    // .input_file = "/home/zchen/nydata/data.bin",
    .input_file = "/home/zongxiong/nyt/data.bin",
    // .input_file = "/tmp/data.bin",
    .cores = 8,
};

Timestamp getTimestamp() {
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
                0  				/* session_gap_ms */
            );

    // create a new pipeline
	Pipeline* p = Pipeline::create(NULL);
	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	NYTParser<string_range, NYTRecord, RecordBundle> parser("[nyt_parser]");
    NYTQuery<NYTRecord, std::pair<uint64_t, uint64_t>, RecordBundle> query("[nyt_query]");
    // create a window, group by key
	WinGBK<std::pair<uint64_t, uint64_t>, RecordBundle,
           WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(config.window_size));
    // reduce
    WinKeyReducer<std::pair<uint64_t, uint64_t>,  /* pair in */
                  WinKeyFragLocal_Std, WinKeyFrag_Std, /* kv d/s */
                  std::pair<uint64_t, uint64_t>,  /* pair out */
                  WindowsBundle /* bundle out */
                  > reducer("[reducer]");

    WindowsBundleSink<std::pair<uint64_t, uint64_t>> sink("[sink]");

    connect_transform(unbound, parser);
    connect_transform(parser, query);
    connect_transform(query, wgbk);
    connect_transform(wgbk, reducer);
    connect_transform(reducer, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}


ofstream statistics_file;
ofstream latency_file;

std::map<Window, ptime, Window> window_keeper;

// #ifdef MEASURE_40M_TUPLES
std::string papi_conf_file = "../papi/papi_conf_global.cfg";
// on gpu-1: Cycles, Instructions, Branches, Branches Mispred, L1D, L1I, dTLB, iTLB
// TLB_PRESET
// BRANCH_PRESET
// DCACHE_MISS_PRESET
// ICACHE_MISS_PRESET
// TOTAL_INSTR_CYCLE_PRESET
// std::string papi_seq = "L1_CACHE_MISS_PRESET";
// std::string papi_seq = "TLB_PRESET";
// std::string papi_seq = "BRANCH_PRESET";
std::string papi_seq = "TOTAL_INSTR_CYCLE_PRESET";
// #endif


int main(int ac, char *av[]) {
    parse_options(ac, av, &config);
    print_config();

    std::string filename = "./test1.csv";
    std::string latency_fname = "./test2.csv";
    statistics_file.open(filename.c_str());
    latency_file.open(latency_fname.c_str());

    Timestamp begin = getTimestamp();
    testInput();
    Timestamp end = getTimestamp();


    double elapsed_time = double(end - begin) / (1024 * 1024 * 1024);
    std::cout << "time elapses: " << elapsed_time << std::endl;

    return 0;
}
