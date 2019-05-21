#define K2_NO_DEBUG 1

#include "Values.h"
#include "core/EvaluationBundleContext.h"
#include "WinGBK.h"
#include "WinGBKEvaluator.h"

template <class KVPair,
					template<class> class InputBundleT,
					template<class> class WindowKeyedFragmentT
					>
void WinGBK<KVPair, InputBundleT, WindowKeyedFragmentT>::ExecEvaluator(int nodeid,
			EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr) {

#ifndef NDEBUG /* if evaluators get stuck ...*/
	static atomic<int> outstanding (0);
#endif
	/* instantiate an evaluator */
	WinGBKEvaluator<KVPair, InputBundleT, WindowKeyedFragmentT> eval(nodeid);
	//eval.evaluate(this, c);

#ifndef NDEBUG		// for debug
	outstanding ++;
#endif
	eval.evaluate(this, c, bundle_ptr);

#ifndef NDEBUG   // for debug
	outstanding --; I("end eval... outstanding = %d", outstanding);
#endif
}

/* template instantiation with concrete types. NB: MT-safe WinKeyFrag often overkill. */

/* ---- for netmon --- */
template
void WinGBK<pair<creek::ippair, long>, RecordBundle, WinKeyFragLocal_Simple>
	::ExecEvaluator(int nodeid,
									EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* ---- for wordcount --- */
template
//void WinGBK<pair<creek::string, long>, RecordBundle<pair<std::string, long>>>
void WinGBK<pair<creek::string, long>, RecordBundle, WinKeyFragLocal_Std>
	::ExecEvaluator(int nodeid,
									EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* ---- for yahoo --- */
template
void WinGBK<pair<uint64_t, long>, RecordBundle, WinKeyFragLocal_Std>
::ExecEvaluator(int nodeid,
				EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

/* todo: instantiate more types */
// nyt
template
void WinGBK<std::pair<uint64_t, uint64_t>, RecordBundle, WinKeyFragLocal_Std>
::ExecEvaluator(int nodeid,
				EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
