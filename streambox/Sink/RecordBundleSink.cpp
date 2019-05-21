#define K2_NO_DEBUG 1

#include "Yahoo/YahooRecord.h"
#include "Sink.h"
#include "RecordBundleSinkEvaluator.h"
#include "RecordBundleSinkEvaluatorJD.h"

#ifdef USE_TBB_DS
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#endif

#include "Source/UnboundedInMemEvaluator.h"

/* ----  specialization: different ways to dump bundle contents  --- */

/* an integer only. for tweet (sentiment score) */
template<>
void RecordBundleSink<long>::printBundle
	(const RecordBundle<long> & input_bundle) {
    W("got one bundle: ", );
#ifndef NDEBUG
    for (auto && rec : input_bundle.content) {
    	I("==== rec (long) ===== ");
    	// XXX dump the ts as well ??
    	EE("long: %ld", rec.data); // dbg, do "grep printBundle"
//    	abort();
    }
    I("----------------");
#endif
}

/* for tweet (tvpair sink) */
template<>
void RecordBundleSink<creek::tvpair>::printBundle
	(const RecordBundle<creek::tvpair> & input_bundle) {
    W("got one bundle: ", );
#ifndef NDEBUG
    for (auto && rec : input_bundle.content) {
    	I("==== rec (long) ===== ");
    	// XXX dump the ts as well ??
    	EE("tvpair: %ld %ld (window end %s",
    			rec.data.first, rec.data.second,
    			to_simple_string(rec.ts).c_str()); // dbg, do "grep printBundle"
//    	abort();
    }
    I("----------------");
#endif
}


/* for wc? */
template<>
void RecordBundleSink<vector<creek::string>>::printBundle
	(const RecordBundle<vector<creek::string>> & input_bundle) {
    W("got one bundle: ");
#ifndef NDEBUG
    for (auto && rec : input_bundle.content) {
    	I("==== rec (vector<string>) ===== ");
    	// XXX dump the ts as well ??
    	for (auto && s : rec.data)
    		cout << to_simplest_string1(rec.ts).str() << " " << s << endl;
    }
    I("----------------");
#endif
}

/* for distinct count */
#if defined(USE_CUCKOO_HASHMAP)
using Set = cuckoohash_map<creek::string, int, CityHasher<creek::string>>;
#elif defined(USE_TBB_HASHMAP)
//using Set = tbb::concurrent_unordered_set<creek::string>;
using Set = creek_set_array::SetArray;
#elif defined(USE_FOLLY_HASHMAP)
using Set = folly::AtomicHashMap<uint64_t, uint64_t>;
#else
#err "undefined"
#endif
using SetPtr = shared_ptr<Set>;

template<>
void RecordBundleSink<SetPtr>::printBundle
	(const RecordBundle<SetPtr> & input_bundle) {
#if 0
	for (auto && rec : input_bundle.content) {
		auto & ptr = rec.data;
		for (auto & s : (*ptr)) {
			cout << s << endl;
		}
	}
#endif
}

//#ifdef USE_FOLLY_HASHMAP
//using HashSet = folly::AtomicHashMap<uint64_t, uint64_t>;
//using HashSetPtr = shared_ptr<HashSet>;
//template<>
//void RecordBundleSink<HashSetPtr>::printBundle
//	(const RecordBundle<HashSetPtr> & input_bundle) {
//
//}
//#endif


//template<>
//void RecordBitmapBundleSink<pair<long, vector<long>>>::ExecEvaluator(
//		int nodeid, EvaluationBundleContext *c)
//{
//	RecordBitmapBundleSinkEvaluator<pair<long, vector<long>>> eval(nodeid);
//	eval.evaluate(this, c);
//}



/* -------------------- ExecEvaluator --------------------
 * out of line to avoid circular dependency
 */

template<class T>
void RecordBundleSink<T>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* modeled after RecordBitmapBundleSink<T>::ExecEvaluator. dispatch
	 * based on side info...
	 */
	if(this->get_side_info() == SIDE_INFO_JD) {
#ifdef WORKAROUND_JOIN_JDD
		RecordBundleSinkEvaluator<T> eval(nodeid); /* construct a normal eval */
#else
		RecordBundleSinkEvaluatorJD<T> eval(nodeid); /* side info JD -- the right way */
#endif
		eval.evaluate(this, c, bundle_ptr);
	} else {
		xzl_assert(this->get_side_info() == SIDE_INFO_JDD
							|| this->get_side_info() == SIDE_INFO_NONE);
		RecordBundleSinkEvaluator<T> eval(nodeid); /* default side info */
		eval.evaluate(this, c, bundle_ptr);
	}
}

#if 0
template<class T>
void RecordBundleSink<T>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c)
{
	RecordBundleSinkEvaluator<T> eval(nodeid);
	eval.evaluate(this, c);
}
#endif

/* ---- instantiation for concrete types --- */

template
void RecordBundleSink<vector<creek::string>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void RecordBundleSink<string_range>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* for grep */
template
void RecordBundleSink<creek::string>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* for Yahoo */
template
void RecordBundleSink<YahooRecord>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/*
template
void RecordBundleSink<string>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
*/
/* for wingrep*/
/*
template
void RecordBundleSink<shared_ptr<tbb::concurrent_vector<std::string>>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
*/
#ifdef USE_TBB_DS
/* temp workaround. XXX we should output vector<> */
template
void RecordBundleSink<tbb::concurrent_vector<creek::string>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void RecordBundleSink<shared_ptr<tbb::concurrent_vector<creek::string>>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void RecordBundleSink<tbb::concurrent_vector<YahooRecord>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
#endif

/* ---- for distinct count ---- */
template
//void RecordBundleSink<shared_ptr<tbb::concurrent_unordered_set<creek::string>>>::ExecEvaluator(
//		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
void RecordBundleSink<SetPtr>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

//#ifdef USE_FOLLY_HASHMAP
//template
//void RecordBundleSink<HashSetPtr>::ExecEvaluator(
//		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
//#endif

/* for tweet (sentiment score) */
template
void RecordBundleSink<long>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* for tweet (test tvpair) */
template
void RecordBundleSink<creek::tvpair>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* for tweet (join) */
template
void RecordBundleSink<pair<long, vector<long>>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

