#include "config.h"

#include "WinKeyReducer/WinKeyReducer.h"
#ifdef WORKAROUND_WINKEYREDUCER_RECORDBUNDLE
#include "WinKeyReducer/WinKeyReducerEvalRecordBundle.h"
#else
#include "WinKeyReducer/WinKeyReducerEval.h"
#endif

#ifdef USE_TBB_DS  /* needed for output */
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#endif

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


/* WindowResult is based on shared ptr.
 * ValueContainerT is decl in header.
 *
 * We only have to specialize combine(), since iterating partitioned HT is different
 */

/* NB: this is for a particular window state (not all windows) */

using NYTWinKeyReducer = WinKeyReducer<
    std::pair<uint64_t, uint64_t>,
    WinKeyFragLocal_Std, WinKeyFrag_Std,
    std::pair<uint64_t, uint64_t>,
    WindowsBundle>;

using WindowResultPtr = shared_ptr<WinKeyFrag_Std<std::pair<uint64_t, uint64_t>>>;

#include "NYTWinKeyReducerCommon.h"
/* template instantiation with concrete types */
template
void NYTWinKeyReducer::ExecEvaluator(int nodeid,
                                     EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
bool NYTWinKeyReducer::ReportStatistics(PTransform::Statstics* stat);
