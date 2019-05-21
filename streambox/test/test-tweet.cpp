#include "core/Pipeline.h"
#include "Source/UnboundedInMemEvaluator.h"
#include "Mapper/WinMapperEvaluator.h"
//#include "WinSum/WindowedSumEvaluator.h"
#include "WinSum/WinSum_addlong.h"
#include "WinSum/WinSum_addlong_tvpair.h"
#include "Sink/WindowsBundleSinkEvaluator.h"
//#include "WinKeyReducer/WinKeyReducerEvalRecordBundle.h"
#include "WinKeyReducer/WinKeyReducerEval.h"

#include "test-common.h"
/*
 *      UnboundedInMem
 *            | BundleT<string_range>
 *            V
 *       WinMapper
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
	.record_size = 200, /* estimate buffer size */
	//		.input_file = "/ssd/twit_data/stream_122716.json",
	.input_file = "/ssd/tweets-textonly.txt",
};
#else
pipeline_config config = {
	.records_per_interval = (1000 * 1000),
	.target_tput = (2000 * 1000),  /* 2000 ~ 1sec */
	.record_size = 200,
	//.input_file = "/ssd/tweets-textonly-x8.txt",
	.input_file = "/ssd/twitter_download/filtered_tweets.txt",
};
#endif

void tweet_sentiment()
{
	UnboundedInMem<string_range, BundleT>
	unbound("unbounded-inmem",
		config.input_file.c_str(),
		config.records_per_interval , 	/* records per wm interval */
		config.target_tput, 		/* target tput (rec/sec) */
		config.record_size,
		0,  				/* session_gap_ms */
		//TXT_TWEETS
		TXT_LINES
	);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	WinMapper<string_range, long, BundleT> mapper("winmapper", seconds(1));
	WinSum_addlong<long, long> agg ("agg",
			2 /* >1 for sliding window, slow right now */);
	RecordBundleSink<long> sink("sink");

	connect_transform(unbound, mapper);
	connect_transform(mapper, agg);
	connect_transform(agg, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

/* test two streams (same pipelines) WITHOUT join */

void tweet_2streams()
{
	UnboundedInMem<string_range, BundleT>
	unbound("unbounded-inmem",
		config.input_file.c_str(),
		config.records_per_interval , 	/* records per wm interval */
		config.target_tput, 		/* target tput (rec/sec) */
		config.record_size,
		0,  				/* session_gap_ms */
		TXT_LINES,
		2 				/* num of output streams */
	);

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	source_transform_1to2(unbound);

	WinMapper<string_range, long, BundleT> mapper0("winmapper", seconds(1));
	//  WindowedSum<long, long> agg0 ("agg", 2 /* >1 for sliding window, slow now */);
	WinSum_addlong_tvpair<long, creek::tvpair> agg0 ("agg", 2 /* >1 for sliding window, slow now */);
	//	RecordBundleSink<long> sink0("sink");
	RecordBundleSink<creek::tvpair> sink0("sink");

	WinMapper<string_range, long, BundleT> mapper1("winmapper", seconds(1));
	//	WindowedSum<long, long> agg1 ("agg", 2 /* >1 for sliding window, slow now */);
	WinSum_addlong_tvpair<long, creek::tvpair> agg1 ("agg", 2 /* >1 for sliding window, slow now */);
	//	RecordBundleSink<long> sink1("sink");
	RecordBundleSink<creek::tvpair> sink1("sink");

	connect_transform_1to2(unbound, mapper0, mapper1);
	connect_transform(mapper0, agg0);
	connect_transform(agg0, sink0);

	connect_transform(mapper1, agg1);
	connect_transform(agg1, sink1);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);

#if 0
	vector<PCollection*> o2 = p->apply2(&unbound);
	o2[0]->_name = "src_out0";
	o2[1]->_name = "src_out1";

	/* left */
	WinMapper<string_range, long, BundleT> mapper0("winmapper", seconds(1));
	auto mapper0_out = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";

	WinSum<long, long> agg0 ("agg",
			2 /* >1 for sliding window, slow right now */);
	RecordBundleSink<long> sink0("sink");
	connect_transform(mapper0, agg0);
	connect_transform(agg0, sink0);

	/* right */
	WinMapper<string_range, long, BundleT> mapper0("winmapper", seconds(1));
	auto mapper0_out = dynamic_cast<PCollection *>(o2[0]->apply1(&mapper0));
	c11->_name = "mapper0_out";

	WinSum<long, long> agg0 ("agg",
			2 /* >1 for sliding window, slow right now */);
	RecordBundleSink<long> sink0("sink");
	connect_transform(mapper0, agg0);
	connect_transform(agg0, sink0);
#endif

}

/* test two same pipelines with join */
/* Somehow there's a failed assertion when connecting J->JD->JDD, regardless of whether they
 * use RecordBundle or RecordBitmapBundle (see br:test/windowedsum-bundletype).
 * However, J->JD->JDD works for test-join.cpp. No idea why.
 *
 * According to hym, his J/JD expects JDD.
 * This workaround allows a pipeline to use no JDD and make JD the sink (dropping
 * all bundles/puncs there).
 *
 */
#ifndef WORKAROUND_JOIN_JDD
#error "must have this"
#endif

static void test_2streams_join()
{
	UnboundedInMem<string_range, BundleT>
		unbound("unbounded-inmem",
				config.input_file.c_str(),
				config.records_per_interval , 	/* records per wm interval */
				config.target_tput, 		/* target tput (rec/sec) */
				config.record_size,
				0,  				/* session_gap_ms */
				TXT_LINES,
				2 				/* num of output streams */
		       );

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	source_transform_1to2(unbound);

	WinMapper<string_range, long, BundleT> mapper0("winmapper", seconds(1));
	WinSum_addlong_tvpair<long, creek::tvpair> agg0 ("agg", 2 /* >1 for sliding window, slow now */);
	agg0.set_side_info(SIDE_INFO_L);

	WinMapper<string_range, long, BundleT> mapper1("winmapper", seconds(1));
	WinSum_addlong_tvpair<long, creek::tvpair> agg1 ("agg", 2 /* >1 for sliding window, slow now */);
	agg1.set_side_info(SIDE_INFO_R);

	Join<pair<long,long>, RecordBundle, RecordBundle> join("[join]", seconds(1));
	join.set_side_info(SIDE_INFO_J);

	RecordBundleSink<pair<long, vector<long>>> sink("[sink]");
	sink.set_side_info(SIDE_INFO_JD);

	//RecordBundleSink<pair<long, vector<long>>> sink1("[sink1]");
	//sink1.set_side_info(SIDE_INFO_JDD);

	connect_transform_1to2(unbound, mapper0, mapper1);

	connect_transform(mapper0, agg0);
	connect_transform(mapper1, agg1);

	connect_transform_2to1(agg0, agg1, join);

	connect_transform(join, sink);
	//connect_transform(sink, sink1);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

#include "Source/UnboundedInMemEvaluator.h"
#include "Mapper/WordCountMapperEvaluator.h"
#include "Win/WinGBKEvaluator.h"

/* this forces WinKeyReducer to invoke an alternate eval that produces
 * RecordBundle instead of WindowedBundle...
 */
//#ifndef WORKAROUND_WINKEYREDUCER_RECORDBUNDLE
//#error "must have this"
//#endif

static void test_wc_tvpair()
{
	UnboundedInMem<string_range, BundleT>
		unbound("unbounded-inmem",
				config.input_file.c_str(),
				config.records_per_interval , /* records per wm interval */
				config.target_tput, 						/* target tput (rec/sec) */
				config.record_size,
				0  													/* session_gap_ms */
		       );

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	WordCountMapper<string_range, pair<creek::string,long>, BundleT> mapper ("[wc-mapper]");

	WinGBK<pair<creek::string, long>, BundleT,
		WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(1));
	WinKeyReducer_winbundle<pair<creek::string,long>,
		WinKeyFragLocal_Std,
		WinKeyFrag_Std>
			reducer("[tvreducer]");
	RecordBundleSink<creek::tvpair> sink("sink");

	connect_transform(unbound, mapper);
	connect_transform(mapper, wgbk);
	connect_transform(wgbk, reducer);
	connect_transform(reducer, sink);
	//  connect_transform(wgbk, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

static void test_score_and_wc()
{
	UnboundedInMem<string_range, BundleT>
		unbound("unbounded-inmem",
				config.input_file.c_str(),
				config.records_per_interval , 	/* records per wm interval */
				config.target_tput, 		/* target tput (rec/sec) */
				config.record_size,
				0,  				/* session_gap_ms */
				TXT_LINES,
				2 				/* num of output streams */
		       );

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	source_transform_1to2(unbound);

	/* left -- scoring */
	WinMapper<string_range, long, BundleT> mapper0("winmapper", seconds(1));
	WinSum_addlong_tvpair<long, creek::tvpair> agg0 ("agg", 1 /* >1 for sliding window, slow now */);
	agg0.set_side_info(SIDE_INFO_L);

	/* right -- wc */
	WordCountMapper<string_range, pair<creek::string,long>, BundleT> mapper1 ("[wc-mapper]");
	WinGBK<pair<creek::string, long>, BundleT,
		WinKeyFragLocal_Std> wgbk ("[wingbk]", seconds(1));
	WinKeyReducer<pair<creek::string,long>,
		WinKeyFragLocal_Std, WinKeyFrag_Std,
		creek::tvpair, RecordBundle> reducer("[tvreducer]");
	reducer.set_side_info(SIDE_INFO_R);

	/* join & sink */
	Join<pair<long,long>, RecordBundle, RecordBundle> join("[join]", seconds(1));
	join.set_side_info(SIDE_INFO_J);

	RecordBundleSink<pair<long, vector<long>>> sink("[sink]");
	sink.set_side_info(SIDE_INFO_JD);


	connect_transform_1to2(unbound, mapper0, mapper1);

	/* left */
	connect_transform(mapper0, agg0);
	/* right */
	connect_transform(mapper1, wgbk);
	connect_transform(wgbk, reducer);

	connect_transform_2to1(agg0, reducer, join);

	connect_transform(join, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

/* no aggregation. for testing
 * source -> winmapper -> sink
 */
void tweet_noagg()
{

	UnboundedInMem<string_range, BundleT>
		unbound("unbounded-inmem",
				config.input_file.c_str(),
				config.records_per_interval , 	/* records per wm interval */
				config.target_tput, 		/* target tput (rec/sec) */
				config.record_size,
				0,  				/* session_gap_ms */
				TXT_TWEETS
		       );

	// create a new pipeline
	Pipeline* p = Pipeline::create(NULL);

	PCollection *unbound_output = dynamic_cast<PCollection *>(p->apply1(&unbound));
	unbound_output->_name = "src_out";

	WinMapper<string_range, long, BundleT> mapper("winmapper", seconds(1));
	WindowsBundleSink<long> sink("sink");

	connect_transform(unbound, mapper);
	connect_transform(mapper, sink);

	EvaluationBundleContext eval(1, config.cores);
	eval.runSimple(p);
}

int main(int ac, char *av[])
{
	parse_options(ac, av, &config);

	print_config();
	//	tweet();
	//	tweet_2streams();
	//	test_2streams_join();
	//	test_wc_tvpair();
	test_score_and_wc();
}
