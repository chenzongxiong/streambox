#ifndef NEXMARK_FILTER_EVALUATOR_HPP
#define NEXMARK_FILTER_EVALUATOR_HPP

#include "core/SingleInputTransformEvaluator.h"
#include "Nexmark/NexmarkFilter.hpp"

/* convert a stream of records to a stream of <NexmarkRecord> records */
template <typename InputT,
          typename OutputT,
          template<class> class BundleT_>
class NexmarkFilterEvaluator
    : public SingleInputTransformEvaluator<NexmarkFilter<InputT, OutputT>,
                                           BundleT_<InputT>, BundleT_<OutputT>
                                           > {

    using TransformT = NexmarkFilter<InputT, OutputT>;
    using InputBundleT = BundleT_<InputT>;
    using OutputBundleT = BundleT_<OutputT>;

public:
    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            if(trans->do_map(*it)) {
                output_bundle->add_record(*it);
            }
        }
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        return true;
        // return false;
    }

    NexmarkFilterEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT,
                                        InputBundleT, OutputBundleT>(node) {}

};

#endif //NEXMARK_FILTER_EVALUATOR_HPP
