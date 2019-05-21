#ifndef RECORDBITMAPBUNDLESINKEVALUATOR_H
#define RECORDBITMAPBUNDLESINKEVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Sink/Sink.h"

/* InputT: the element type of the record bundle */
template <typename InputT>
class RecordBitmapBundleSinkEvaluator
    : public SingleInputTransformEvaluator<RecordBitmapBundleSink<InputT>,
      RecordBitmapBundle<InputT>, RecordBitmapBundle<InputT>> {

	using TransformT = RecordBitmapBundleSink<InputT>;
	using InputBundleT = RecordBitmapBundle<InputT>;
	using OutputBundleT = RecordBitmapBundle<InputT>;

public:

	RecordBitmapBundleSinkEvaluator(int node)
	: SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }

  bool evaluateSingleInput (TransformT* trans,
        shared_ptr<InputBundleT> input_bundle,
        shared_ptr<OutputBundleT> output_bundle) override {

    //TransformT::printBundle(*input_bundle);
    TransformT::report_progress(* input_bundle);
    return false; /* no output bundle */
  }

};

#endif /* RECORDBITMAPBUNDLESINKEVALUATOR_H */
