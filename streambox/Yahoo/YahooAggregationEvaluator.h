//
// Created by manuelrenz on 09.04.18.
//

#ifndef YAHOOAGGREGATIONEVALUATOR_H
#define YAHOOAGGREGATIONEVALUATOR_H

#include "Values.h"
#include "Yahoo/YahooAggregation.h"
#include "core/SingleInputTransformEvaluator.h"


template <typename InputT, typename OutputT,
        //					typename InputBundleT, typename OutputBundleT,
    template<class> class BundleT
    >
// wc_mapper::mode mode>
class YahooAggregationEvaluator
        : public SingleInputTransformEvaluator<
    // YahooAggregation<InputT, OutputT, BundleT, mode>,
    YahooAggregation<InputT, OutputT, BundleT>,
                //RecordBitmapBundle<InputT>, RecordBitmapBundle<OutputT>> {
                //      			InputBundleT, OutputBundleT> {
          BundleT<InputT>, BundleT<OutputT>>
{

using InputBundleT = BundleT<InputT>;
using OutputBundleT = BundleT<OutputT>;
// using TransformT = YahooAggregation<InputT, OutputT, BundleT, mode>;
using TransformT = YahooAggregation<InputT, OutputT, BundleT>;

//  using TransformT = YahooAggregation<InputT, OutputT, InputBundleT, OutputBundleT, mode>;
//  using InputBundleT = RecordBitmapBundle<InputT>;
//  using OutputBundleT = RecordBitmapBundle<OutputT>;

public:
bool evaluateSingleInput (TransformT* trans,
                          shared_ptr<InputBundleT> input_bundle,
                          shared_ptr<OutputBundleT> output_bundle) override {

// go through Records in input bundle (the iterator automatically
// skips "masked" Records.
#if 0
    std::cout << "thread id in yahoo aggregte: " << std::this_thread::get_id() << std::endl;
#endif
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        TransformT::do_map(*it, output_bundle); /* static rocks! */
    }
    TransformT::record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
    return true;
}

YahooAggregationEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }
};

#endif //YAHOOAGGREGATIONEVALUATOR_H
