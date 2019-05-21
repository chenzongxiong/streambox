#ifndef NYT_WINDOWED_SUM_H
#define NYT_WINDOWED_SUM_H

#include "config.h"

#include "Values.h"
#include "core/StatefulTransform.h"

template <
    typename InputT,
    typename OutputT>
class NYTWindowedSum
    : public StatefulTransform<NYTWindowedSum<InputT, OutputT>, InputT, OutputT> {

public:

    using TransformT = NYTWindowedSum<InputT, OutputT>;

    NYTWindowedSum(string name)
        : StatefulTransform<TransformT, InputT, OutputT>(name), multi(multi) {}

    /* multiplier for supporting sliding window.
     * = 1 for tumbling window.	 * >1 for sliding window (x factor for the sliding window size)
	 */
    const int multi;

    /* for aggregating one single input. to be specialized */
    static OutputT const & aggregate_init(OutputT * acc);
    static OutputT const & aggregate(OutputT * acc, InputT const & in);
    /* combine the evaluator's (partial) aggregation results (for a particular window)
     * to the tran's internal state.
     */
    static OutputT const & combine(OutputT & mine, OutputT const & others);

    void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                       shared_ptr<BundleBase> bundle) override;

};

#endif // NYT_WINDOWED_SUM_H
