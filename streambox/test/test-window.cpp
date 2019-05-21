#include <fstream>
#include <map>

#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
#include "Win/CountBasedTumblingWindow.h"
#include "WinKeyReducer/WinKeyReducer.h"
#include "NewYorkTaxi/NYTWindowedSum.h"
#include "NewYorkTaxi/NYTRecord.h"
#include "NewYorkTaxi/NYTParser.h"
#include "NewYorkTaxi/NYTParserEvaluator.h"
#include "NewYorkTaxi/NYTQuery.h"
#include "Values.h"
#include "test-common.h"

typedef uint64_t Timestamp;
using NanoSeconds = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;
using namespace std;

template<class T>
using BundleT = RecordBundle<T>;

template<class T>
using WindowsBundleT = WindowsBundle<T>;

pipeline_config config = {
    .records_per_interval = 1000000,
    .target_tput = std::numeric_limits<uint64_t>::max(),
    .record_size = sizeof(NYTRecord),
    .window_size = 2,
    .campaigns = 1,
    .input_file = "/home/zongxiong/nyt/data.bin",
    .cores = 8,
};


// Timestamp getTimestamp() {
//     return std::chrono::duration_cast<NanoSeconds>(
//         Clock::now().time_since_epoch())
//         .count();
// }

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

    NYTParser<string_range, NYTRecord, BundleT> parser("[nyt_parser]");
    NYTQuery<NYTRecord, std::pair<uint64_t, uint64_t>, RecordBundle> query("[nyt_query]");
    CountBasedTumblingWindow<std::pair<uint64_t, uint64_t>, BundleT> cbtw("[count_window]", 1000);

    // reduce aggregation
    // NYTWindowedSum<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> reducer("[reducer]");
    WinKeyReducer<pair<uint64_t, uint64_t>,  /* pair in */
                  WinKeyFragLocal_Std,
                  WinKeyFrag_Std, /* kv d/s */
                  pair<uint64_t, uint64_t>,  /* pair out */
                  WindowsBundle /* bundle out */
                  > reducer("[reducer]");

    WindowsBundleSink<std::pair<uint64_t, uint64_t>> sink("[sink]");

    connect_transform(unbound, parser);
    connect_transform(parser, query);
    connect_transform(query, cbtw);
    connect_transform(cbtw, reducer);
    connect_transform(reducer, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);

}


ofstream statistics_file;
ofstream latency_file;
// boost::posix_time::ptime zxchen_start; //global variable
// int zxchen_first = 1;
// int first_end = 1;

std::map<Window, ptime, Window> window_keeper;
// ptime czx_start = min_date_time;

int main(int ac, char *av[]) {
    parse_options(ac, av, &config);
    print_config();

    std::string filename = "./results/test_yahoo_records_" + std::to_string(config.records_per_interval) + "_cores_" + std::to_string(config.cores) + "_window_size_" + std::to_string(config.window_size) + "_campaigns_" + std::to_string(config.campaigns) + ".csv";
    std::string latency_fname = "./results/latency/test_yahoo_records_" + std::to_string(config.records_per_interval) + "_cores_" + std::to_string(config.cores) + "_window_size_" + std::to_string(config.window_size) + "_campaigns_" + std::to_string(config.campaigns) + ".txt";
    statistics_file.open(filename.c_str());

    latency_file.open(latency_fname.c_str());

    // Timestamp begin = getTimestamp();
    testInput();
    // Timestamp end = getTimestamp();


    // double elapsed_time = double(end - begin) / (1024 * 1024 * 1024);
    // std::cout << "time elapses: " << elapsed_time << std::endl;
}
