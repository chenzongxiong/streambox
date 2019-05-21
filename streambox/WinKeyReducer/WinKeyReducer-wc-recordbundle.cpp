
#include "config.h"

#include "WinKeyReducer.h"
#include "WinKeyReducerEval.h"

#if 1
/* necessary since instantiation later needs these defs. XXX factor out to a common header. but not in WinKeyReducer.h */
template <typename KVPair,
		template<class> class InputWinKeyFragT,
		template<class> class InternalWinKeyFragT,
		typename KVPairOut,
		template<class> class OutputBundleT_
>
void WinKeyReducer<KVPair,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_>::ExecEvaluator
		(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
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

using MyKVPairIn = std::pair<creek::string, long>;
using MyKVPairOut = creek::tvpair;
using MyWinKeyReducer = WinKeyReducer<
		MyKVPairIn, /* kv in */
		WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
		MyKVPairOut, /* kvout, note that this is different than in pair */
		RecordBundle /* output bundle */
>;
using KVPair = pair<creek::string, long>;
using InputT = WinKeyFragLocal_Std<KVPair>;
using WindowResultPtr = shared_ptr<WinKeyFrag_Std<KVPair>>;


#include "WinKeyReducer-wc-common.h"

template
void MyWinKeyReducer::ExecEvaluator(int nodeid,
																		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
