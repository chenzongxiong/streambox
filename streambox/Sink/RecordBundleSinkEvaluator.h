#ifndef RECORDBUNDLESINKEVALUATOR_H
#define RECORDBUNDLESINKEVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Sink/Sink.h"

/* InputT: the element type of the record bundle */
template <typename InputT>
class RecordBundleSinkEvaluator
    : public SingleInputTransformEvaluator<RecordBundleSink<InputT>,
                                           RecordBundle<InputT>, RecordBundle<InputT>> {

	using TransformT = RecordBundleSink<InputT>;
	using InputBundleT = RecordBundle<InputT>;
	using OutputBundleT = RecordBundle<InputT>;

public:

	RecordBundleSinkEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }

    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        TransformT::printBundle(*input_bundle);
        // for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        //     if(trans->do_map(*it)) {
        //         output_bundle->add_record(*it);
        //     }
        // }
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        return false; /* no output bundle */
    }
};

#endif /* RECORDBUNDLESINKEVALUATOR_H */
