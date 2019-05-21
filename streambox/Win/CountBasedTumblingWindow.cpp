#include "Win/CountBasedTumblingWindow.h"
#include "Win/CountBasedTumblingWindowEvaluator.h"
#include "NewYorkTaxi/NYTRecord.h"

template <typename T, template<class> class InputBundle>
void CountBasedTumblingWindow<T, InputBundle>::ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                                                             shared_ptr<BundleBase> bundle_ptr) {
	CountBasedTumblingWindowEvaluator<T, InputBundle> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}

/* -------instantiation concrete classes------- */
// template
// void CountBasedTumblingWindow<string_range, RecordBundle>::ExecEvaluator(int nodeid,
//                                                                          EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

// template
// void CountBasedTumblingWindow<creek::string, RecordBundle>::ExecEvaluator(int nodeid,
//                                                                           EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

// template
// void CountBasedTumblingWindow<YahooRecord, RecordBundle>::ExecEvaluator(int nodeid,
//                                                                         EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* todo: more types */
// nyt
// template
// void CountBasedTumblingWindow<NYTRecord, RecordBundle>::ExecEvaluator(int nodeid,
//                                                                       EvaluationBundleContext *c,
//                                                                       shared_ptr<BundleBase> bundle_ptr);
template
void CountBasedTumblingWindow<std::pair<uint64_t, uint64_t>, RecordBundle>::ExecEvaluator(int nodeid,
                                                                                          EvaluationBundleContext *c,
                                                                                          shared_ptr<BundleBase> bundle_ptr);
