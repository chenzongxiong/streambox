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
            trans->do_map(*it, output_bundle);
        }


        // uint64_t maxPrice = 0;
        // NexmarkRecord *maxRecord = nullptr;
        // ptime ts;
        // for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        //     trans->do_map(*it, output_bundle);
        //     // if (maxPrice < it->data.price) {
        //     //     maxPrice = it->data.price;
        //     //     // maxRecord = &(it->data);
        //     //     ts = it->ts;
        //     // }
        //     Record<NexmarkRecord> const &rec = *it;
        //     if (maxPrice < rec.data.price) {
        //         maxPrice = rec.data.price;
        //     }
        // }
        // output_bundle->emplace_record(std::pair<uint64_t, uint64_t>(1, maxPrice), ts);
        // std::cout << "input_bundle->size: " << input_bundle->size() << std::endl;

        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        // return false;
        return true;
    }

    NexmarkAggregationEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) {}
};

#endif //NEXMARKAGGREGATIONEVALUATOR_H
