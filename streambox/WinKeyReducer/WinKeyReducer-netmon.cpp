/*
 * WindowKeyedReducer-netmon.cpp
 *  Split from WindowKeyedReducer.cpp. specialized version.
 */

#include "config.h"

#include "WinKeyReducer.h"
#include "WinKeyReducerEval.h"

#ifdef USE_TBB_DS  /* needed for output */
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#endif

#if 1
/* necessary since instantiation later needs these defs.  */
template <typename KVPair,
					template<class> class InputWinKeyFragT,
					template<class> class InternalWinKeyFragT,
					typename KVPairOut,
					template<class> class OutputBundleT_
					>
void WinKeyReducer<KVPair,InputWinKeyFragT,
			InternalWinKeyFragT, KVPairOut, OutputBundleT_>::ExecEvaluator
	(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* instantiate an evaluator */
	WinKeyReducerEval<KVPair,InputWinKeyFragT,InternalWinKeyFragT,
			KVPairOut, OutputBundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}
#endif

/* ----------------------- ops specialization for netmon --------------------------------- */

using KVPair1 = pair<creek::ippair, long>;
using MyWinKeyReducer = WinKeyReducer<
		KVPair1,
		WinKeyFragLocal_Simple, WinKeyFrag_SimpleMT,
		KVPair1,
		WindowsBundle>;

/* for a particular window.
 * state: partitioned HT, simple vcontainer, MT safe */
using WindowResultPtr1 = shared_ptr<WinKeyFrag_SimpleMT<KVPair1>>;

/* XXX if @mine is empty, should we just swap pointers? */
template <>
WindowResultPtr1 const & MyWinKeyReducer
::combine(WindowResultT & mine, LocalWindowResultT const & others)
{
	xzl_assert(mine && others);
	xzl_assert(WindowEqual()(mine->w, others->w));

	/* go through all internal maps in @others, i.e. partitioned HT
	 * XXX this in fact incur unnecessary hash compute */
#if 0
	int cnt = 0;
	for (int i = 0; i < others->vals.num_maps; i++) {
		for (auto && kvs : others->vals.get_map(i)) {
			auto & key = kvs.first;
			auto & v_container = kvs.second;
			mine->add_vcontainer_safe(key, v_container);
		}
		cnt += others->vals.get_map(i).size();
	}
	EE("---->total keys %d", cnt);
#endif

#if 1
	for (auto const & map : others->vals.maps) {
		for (auto && kvs : map) {
			auto & key = kvs.first;
			auto & v_container = kvs.second;
			mine->add_vcontainer_safe(key, v_container);
		}
	}
#endif

#if 0
	for (auto && kvs : others->vals) {
		auto & key = kvs.first;
		auto & v_container = kvs.second;
		mine->add_vcontainer_safe(key, v_container);
	}
//	EE("---->total keys %lu", others->vals.size());
#endif

	return mine;
}

/* ------------------- for netmon. also sum all values -------------------  */

template<>
KVPair1 MyWinKeyReducer::do_reduce(creek::ippair const & key,
	InternalValueContainerT const & vcontainer)
{

  auto end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}

/* same impl but use ValueContainerUnsafe */
template<>
KVPair1 MyWinKeyReducer::do_reduce_unsafe
  (creek::ippair const & key, InputValueContainerT const & vcontainer) {

//  VV("vcontainer has %lu elements...", vcontainer.vals.size());
//  auto bcontainer_ptr = vcontainer.vals[0];
//  VV("bcontainer0 %ld %ld ...", (*bcontainer_ptr)[0], (*bcontainer_ptr)[1]);

  auto end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}

/* template instantiation with concrete types */
template
void MyWinKeyReducer::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

#if 0
template
void WindowKeyedReducer<pair<creek::ippair, long>,
			WindowKeyedFragmentSimpleMT,
			WindowKeyedFragmentSimpleMT>
		::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

//template
//void WindowKeyedReducer<pair<creek::ippair, long>,
//		WindowKeyedFragmentSimpleMT,
//		WindowKeyedFragmentSimpleMT>
//	::OnNewUpstreamWatermark(int nodeid,
//		EvaluationBundleContext *c, ptime wm);

#endif
