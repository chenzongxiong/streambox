#include <Mapper/YahooMapper.h>
#include <Win/FixedWindowInto.h>
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
#include "Select/SimpleSelect.h"
#include "Source/YahooBenchmarkSource.h"
#include "Win/FixedWindowIntoEvaluator.h"
//#include "WindowedSumEvaluator.h"
#include "WinSum/WinSum_mergeset.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
#include "test-common.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <WinSum/WinSum_mergevector.h>

#include "Yahoo/YahooParserEvaluator.h"

/*
 *      UnboundedInMem
 *            | BundleT<string_range>
 *            V
 *      YahooRecordParser
 *            | BUndleT<YahooRecord>
 *            V
 *         YahootMapper
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

pipeline_config config = {
    .records_per_interval = 1250000,
    .target_tput = std::numeric_limits<uint64_t>::max(),
//    .record_size = 136,
    .record_size = 78,
    .input_file =
    "/home/zxchen/data_test/Data.txt",
    .cores = 7,
};

void testYahooBenchmark()
{
    UnboundedInMem<string_range, BundleT> unbound(
        "[unbound]",
        config.input_file.c_str(),
        config.records_per_interval , /* records per wm interval */
        config.target_tput,	      /* target tput (rec/sec) */
        config.record_size,
        0,
        TXT_YAHOO
        );

    // create a new pipeline
    Pipeline* p = Pipeline::create(NULL);

    // Source
    // PCollection *yahooSource_output = dynamic_cast<PCollection *>(p->apply1(&yahooSource));
    // yahooSource_output->_name = "src_out";
    PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
    unbound_output->_name = "src_out";

    YahooParser<string_range, YahooRecord, RecordBundle> parser("[yahoo_parser]");

    // // With group by
    // YahooMapper<Event, pair<creek::string,long>, BundleT> mapper ("[yahoo-mapper]", campaigns);
    // WinGBK<pair<creek::string, long>, BundleT, WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(1));
    // WinKeyReducer<pair<creek::string,long>,  // pair in
    //               WinKeyFragLocal_Std, WinKeyFrag_Std, // kv d/s
    //               pair<creek::string,long>,  // pair out
    //               WindowsBundle // bundle out
    //               > reducer("[reducer]");
    // WindowsBundleSink<pair<creek::string,long>> sink("sink");

    connect_transform(unbound, parser);
    // connect_transform(unbound, mapper);
    // connect_transform(mapper, wgbk);
    // connect_transform(wgbk, reducer);
    // connect_transform(reducer, sink);

    // Eval the pipeline
    EvaluationBundleContext eval(1, config.cores);
    eval.runSimple(p);
}

//#include <jemalloc/jemalloc.h>

int main(int ac, char *av[])
{
    parse_options(ac, av, &config);

    //  malloc_stats_print(NULL, NULL, NULL);  // can test jemalloc
    print_config();
    testYahooBenchmark();

    return 0;
}
