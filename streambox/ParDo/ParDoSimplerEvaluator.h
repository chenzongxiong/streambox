#ifndef PARDOSIMPLEREVALUATOR_H
#define PARDOSIMPLEREVALUATOR_H

/* InputT, OutputT: the input/output elements */
template <typename InputT, typename OutputT, typename DoFnSimpler>
class ParDoSimplerEvaluator
		: public SingleInputTransformEvaluator<ParDo<InputT, OutputT>,
		  		RecordBundle<InputT>,RecordBundle<OutputT>> {

	using TransformT = ParDoSimpler<InputT, OutputT, DoFnSimpler>;
  using InputBundleT = RecordBundle<InputT>;
  using OutputBundleT = RecordBundle<OutputT>;


  ParDoSimplerEvaluator(int node)
  : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }

  bool evaluateSingleInput (TransformT* trans,
        shared_ptr<InputBundleT> input_bundle,
        shared_ptr<OutputBundleT> output_bundle) override {

    // Go over each element. Note that we pass in reference so that
    // they can be updated in place.
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        TransformT::fn.dofn(*it, *output_bundle);
    }

  	return true;
  }
};

#endif /* PARDOSIMPLEREVALUATOR_H */
