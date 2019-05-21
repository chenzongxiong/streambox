//
// same as addlong, but output tvpair. see WinSumTvEval.h for more.
//

#ifndef CREEK_WINSUM_TVPAIRCOUNT_H
#define CREEK_WINSUM_TVPAIRCOUNT_H

#include "WinSumBase.h"

template <typename InputT, typename OutputT>
class WinSum_addlong_tvpair : public WinSumBase<WinSum_addlong_tvpair, InputT, OutputT> {
public:
    WinSum_addlong_tvpair(string name, int multi = 1)
            : WinSumBase<WinSum_addlong_tvpair, InputT, OutputT>(name, multi) { }

    /* for aggregating one single input. to be specialized */
    static OutputT const & aggregate_init(OutputT * acc);
    static OutputT const & aggregate(OutputT * acc, InputT const & in);

    /* combine the evaluator's (partial) aggregation results (for a particular window)
     * to the tran's internal state.
     */
    static OutputT const & combine(OutputT & mine, OutputT const & others);

    void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                       std::shared_ptr<BundleBase> bundle) override;
};

#endif //CREEK_WINSUM_TVPAIRCOUNT_H
