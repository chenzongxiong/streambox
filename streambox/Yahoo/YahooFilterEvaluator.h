#ifndef YAHOOFILTEREVALUATOR_H
#define YAHOOFILTEREVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Yahoo/YahooRecord.h"
#include "Yahoo/YahooFilter.h"

/* convert a stream of records to a stream of <YahooRecord> records */
template <typename InputT, typename OutputT, template<class> class BundleT_>
class YahooFilterEvaluator
        : public SingleInputTransformEvaluator<YahooFilter<InputT,OutputT>,
                BundleT_<InputT>, BundleT_<OutputT>> {

    using TransformT = YahooFilter<InputT,OutputT>;
    using InputBundleT = BundleT_<InputT>;
    using OutputBundleT = BundleT_<OutputT>;

public:
    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            // trans->record_counter_.fetch_add(1, std::memory_order_relaxed);
            if(trans->do_map(*it)) {
                output_bundle->add_record(*it);
            }
        }
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        return true;
        // return false;
    }

    YahooFilterEvaluator(int node)
            : SingleInputTransformEvaluator<TransformT,
            InputBundleT, OutputBundleT>(node) { }

};

#endif //YAHOOFILTEREVALUATOR_H
