#ifndef NEXMARKAGGREGATIONEVALUATOR_H
#define NEXMARKAGGREGATIONEVALUATOR_H

#include "Values.h"
#include "core/SingleInputTransformEvaluator.h"
#include "Nexmark/NexmarkAggregation.hpp"


template <typename InputT,
          typename OutputT,
          template<class> class BundleT
          >
class NexmarkAggregationEvaluator : public SingleInputTransformEvaluator<
    NexmarkAggregation<InputT, OutputT, BundleT>,
    BundleT<InputT>, BundleT<OutputT>>
{

    using InputBundleT = BundleT<InputT>;
    using OutputBundleT = BundleT<OutputT>;
    using TransformT = NexmarkAggregation<InputT, OutputT, BundleT>;

public:
    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            // TransformT::do_map(*it, output_bundle); /* static rocks! */
            trans->do_map(*it, output_bundle);
        }
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        return false;
    }

    NexmarkAggregationEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) {}
};

#endif //NEXMARKAGGREGATIONEVALUATOR_H
