
#include "WinSum_addlong_tvpair.h"

// ----------------------------------------------------- //
// operators
// ----------------------------------------------------- //

/* ---- for tweet, the output is <ts, long> which reuses hash join logic ---
 * NB: we don't manipulate ts here; the evaluator will write to the emitted
 * records based on window ts.
 */
template <>
creek::tvpair const & WinSum_addlong_tvpair<long, creek::tvpair>::aggregate(creek::tvpair * acc,
                                                                  long const & in) {
    acc->second += in;
    return *acc;
}

template <>
creek::tvpair const & WinSum_addlong_tvpair<long, creek::tvpair>::aggregate_init(creek::tvpair * acc) {
    acc->second = 0;
    return *acc;
}

template <>
creek::tvpair const & WinSum_addlong_tvpair<long, creek::tvpair>::combine
        (creek::tvpair & mine, creek::tvpair const & others) {
    mine.second += others.second;
    return mine;
}

#include "WinSumTvEval.h"

#define MyWinSum WinSum_addlong_tvpair
#define MyWinSumEval WinSumTvEval

#include "WinSum-dispatch-eval.h"

/* instantiate concreate classes, otherwise won't link */

template
void
WinSum_addlong_tvpair<long, creek::tvpair>
::ExecEvaluator(int nodeid,
                EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);
