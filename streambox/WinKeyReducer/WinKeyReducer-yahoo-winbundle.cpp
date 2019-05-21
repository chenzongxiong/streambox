//
// Created by manuelrenz on 12.04.18.
//


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

// using MyWinKeyReducer = WinKeyReducer<
//         pair<uint64_t, long>,
//         WinKeyFragLocal_Std,
//         WinKeyFrag_Std,
//         pair<uint64_t, long>,
//         WindowsBundle>;
// using KVPair = pair<uint64_t, long>;
// using InputT = WinKeyFragLocal_Std<KVPair>;
// using WindowResultPtr = shared_ptr<WinKeyFrag_Std<KVPair>>;

using MyWinKeyReducer = WinKeyReducer<
        pair<uint64_t, uint64_t>,
        WinKeyFragLocal_Std,
        WinKeyFrag_Std,
        pair<uint64_t, uint64_t>,
        WindowsBundle>;
using KVPair = pair<uint64_t, uint64_t>;
using InputT = WinKeyFragLocal_Std<KVPair>;
using WindowResultPtr = shared_ptr<WinKeyFrag_Std<KVPair>>;


#include "WinKeyReducer-yahoo-common.h"

/* template instantiation with concrete types */

template
void MyWinKeyReducer::ExecEvaluator(int nodeid,
                                    EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
bool MyWinKeyReducer::ReportStatistics(PTransform::Statstics* stat);
