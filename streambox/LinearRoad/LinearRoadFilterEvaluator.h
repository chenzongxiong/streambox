#ifndef LRB_FILTER_EVALUATOR_H
#define LRB_FILTER_EVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "LinearRoad/LinearRoadRecord.h"
#include "LinearRoad/LinearRoadFilter.h"

template <typename InputT,
  typename OutputT,
  template<class> class BundleT_>

  class LinearRoadFilterEvaluator
  : public SingleInputTransformEvaluator<LinearRoadFilter<InputT,OutputT>,
  BundleT_<InputT>, BundleT_<OutputT>> {

  using TransformT = LinearRoadFilter<InputT,OutputT>;
  using InputBundleT = BundleT_<InputT>;
  using OutputBundleT = BundleT_<OutputT>;

public:
  bool evaluateSingleInput (TransformT* trans,
                            shared_ptr<InputBundleT> input_bundle,
                            shared_ptr<OutputBundleT> output_bundle) override {

    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {

      if(trans->do_map(*it)) {
        trans->record_counter_.fetch_add(1, std::memory_order_relaxed);
        output_bundle->add_record(*it);
      }
    }
    // std::cout << "record_counter: " << trans->record_counter_.load(std::memory_order_relaxed) << std::endl;

    return true;
    // return false;
  }

LinearRoadFilterEvaluator(int node)
  : SingleInputTransformEvaluator<TransformT,
    InputBundleT, OutputBundleT>(node) { }

};

#endif //LRB_FILTER_EVALUATOR_H
