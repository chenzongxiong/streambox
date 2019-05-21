/*
 * WindowKeyedReducer-wc.cpp
 *
 */

#include "config.h"

#include "WinKeyReducer.h"
#ifdef WORKAROUND_WINKEYREDUCER_RECORDBUNDLE
#include "WinKeyReducerEvalRecordBundle.h"
#else
#include "WinKeyReducerEval.h"
#endif

#ifdef USE_TBB_DS  /* needed for output */
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#endif

#if 1 /* there's a generic on in WinKeyReducer.cpp */
/* necessary since instantiation later needs these defs */
template <typename KVPair,
					template<class> class InputWinKeyFragT,
					template<class> class InternalWinKeyFragT,
					typename KVPairOut,
					template<class> class OutputBundleT_
>
void WinKeyReducer<KVPair,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_>::ExecEvaluator
	(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* instantiate an evaluator */
#ifdef WORKAROUND_WINKEYREDUCER_RECORDBUNDLE
	WinKeyReducerEvalRecordBundle<KVPair,InputWindowKeyedFragmentT,InternalWindowKeyedFragmentT> eval(nodeid);
#else
	WinKeyReducerEval<KVPair,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_> eval(nodeid);
#endif
	eval.evaluate(this, c, bundle_ptr);
}
#endif

/* WindowResult is based on shared ptr.
 * ValueContainerT is decl in header.
 *
 * We only have to specialize combine(), since iterating partitioned HT is different
 */

/* NB: this is for a particular window state (not all windows) */

using MyWinKeyReducer = WinKeyReducer<
		pair<creek::string, long>,
		WinKeyFragLocal_Std,
		WinKeyFrag_Std,
		pair<creek::string, long>,
		WindowsBundle>;
using KVPair = pair<creek::string, long>;
using InputT = WinKeyFragLocal_Std<KVPair>;
using WindowResultPtr = shared_ptr<WinKeyFrag_Std<KVPair>>;

#include "WinKeyReducer-wc-common.h"

#if 0
/* for a particular window */
/* not in use: since eval do local reduce and then combine.
 * And it requires to iterate HT, not supported by creek::map */
template <>
WindowResultPtr const & WinKeyReducer<pair<creek::string, long>,
			WinKeyFragLocal_Std, WinKeyFrag_Std>::aggregate
		(WindowResultPtr * acc, InputT const & in)
{
#if 0
	assert(acc);
	assert(WindowEqual()((*acc)->w, in.w));

	for (auto && kvs : in.vals) {
		auto & key = kvs.first;
		auto & v_container = kvs.second;
		(*acc)->add_vcontainer_safe(key, v_container);
	}
#endif
	xzl_bug("not impl");
	return *acc;
}

/* XXX if @mine is empty, should we just swap pointers? */
template <>
WindowResultPtr const & WinKeyReducer<pair<creek::string, long>,
					WinKeyFragLocal_Std, WinKeyFrag_Std>::combine
	(WindowResultPtr & mine, WindowResultPtr const & others)
{
	xzl_assert(mine && others);
	xzl_assert(WindowEqual()(mine->w, others->w));

	for (auto && kvs : others->vals) {
		auto & key = kvs.first;
		auto & v_container = kvs.second;
		mine->add_vcontainer_safe(key, v_container);
	}

	return mine;
}

template<>
std::pair<creek::string, long> WinKeyReducer<std::pair<creek::string, long>,
WinKeyFragLocal_Std, WinKeyFrag_Std>::do_reduce
  (creek::string const & key, InternalValueContainerT const & vcontainer) {

  auto end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}

template<>
std::pair<creek::string, long> WinKeyReducer<std::pair<creek::string, long>,
WinKeyFragLocal_Std, WinKeyFrag_Std>::do_reduce_unsafe
  (creek::string const & key, InputValueContainerT const & vcontainer) {

  auto end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}
#endif

/* template instantiation with concrete types */

template
void MyWinKeyReducer::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

//template
//void WindowKeyedReducer<pair<creek::string, long>,WindowKeyedFragmentUnsafe,WindowKeyedFragmentStd>::OnNewUpstreamWatermark(int nodeid,
//		EvaluationBundleContext *c, ptime wm);


