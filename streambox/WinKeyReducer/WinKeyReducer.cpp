#include "config.h"

#include "WinKeyReducer.h"
#include "WinKeyReducerEval.h"

#ifdef USE_TBB_DS  /* needed for output */
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#endif

template <typename KVPair,
					template<class> class InputWindowKeyedFragmentT,
					template<class> class InternalWindowKeyedFragmentT,
					typename KVPairOut,
					template<class> class OutputBundleT_
>
void WinKeyReducer<KVPair,InputWindowKeyedFragmentT,
		InternalWindowKeyedFragmentT, KVPairOut, OutputBundleT_>
::ExecEvaluator
	(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* instantiate an evaluator */
	WinKeyReducerEval<KVPair,InputWindowKeyedFragmentT,
			InternalWindowKeyedFragmentT, KVPairOut, OutputBundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}


//template <typename KVPair,
//					template<class> class InputWindowKeyedFragmentT,
//					template<class> class InternalWindowKeyedFragmentT>
//void WindowKeyedReducer<KVPair,InputWindowKeyedFragmentT,InternalWindowKeyedFragmentT>::OnNewUpstreamWatermark(
//		int nodeid, EvaluationBundleContext *c, ptime wm)
//{
//	/* instantiate an evaluator */
//	WinKeyReducerEval<KVPair,InputWindowKeyedFragmentT,InternalWindowKeyedFragmentT> eval(nodeid);
//	eval.OnNewUpstreamWatermark(wm, this, c);
//}

/* ---------------------------------------------------------- */

// #if 0
// // a default implementation -- sum all values
// // (belonging to the same key)
// // XXX combine this with the following?
// template<>
// std::pair<long, long> WindowKeyedReducer<std::pair<long, long>,WindowKeyedFragmentT>::do_reduce
//   (long const & key, ValueContainerT const & vcontainer) {

//   auto end = vcontainer.cend(); /* avoid calling it in each iteration */

//   long sum = 0;
//   for (auto it = vcontainer.cbegin(); it != end; ++it)
//     sum += *it;

//   return make_pair(key, sum);
// }
// #endif




/* ---------------------------------------------------------- */
