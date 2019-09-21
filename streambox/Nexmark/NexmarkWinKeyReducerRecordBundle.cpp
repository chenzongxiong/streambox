#include "config.h"

#include "WinKeyReducer/WinKeyReducer.h"
#include "WinKeyReducer/WinKeyReducerEval.h"

#if 1
/* necessary since instantiation later needs these defs. XXX factor out to a common header. but not in WinKeyReducer.h */
template <typename KVPair,
          template<class> class InputWinKeyFragT,
          template<class> class InternalWinKeyFragT,
          typename KVPairOut,
          template<class> class OutputBundleT_
          >
void WinKeyReducer<KVPair,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_>::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	WinKeyReducerEval<KVPair,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}
#endif


/* WindowResult is based on shared ptr.
 * ValueContainerT is decl in header.
 *
 * We only have to specialize combine(), since iterating partitioned HT is different
 */

/* NB: this is for a particular window state (not all windows) */

using NexmarkWinKeyReducer = WinKeyReducer<
    std::pair<uint64_t, uint64_t>,
    WinKeyFragLocal_Std, WinKeyFrag_Std,
    std::pair<uint64_t, uint64_t>,
    RecordBundle>;
using WindowResultPtr = shared_ptr<WinKeyFrag_Std<std::pair<uint64_t, uint64_t>>>;

#include "NexmarkWinKeyReducerCommon.hpp"

template
void NexmarkWinKeyReducer::ExecEvaluator(int nodeid,
                                         EvaluationBundleContext *c,
                                         shared_ptr<BundleBase> bundle_ptr);
