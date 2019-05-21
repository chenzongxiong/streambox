#ifndef WC_MAPPER_EVAL_H
#define WC_MAPPER_EVAL_H

#include "Values.h"
#include "Mapper/WordCountMapper.h"
#include "core/SingleInputTransformEvaluator.h"

template <typename InputT, typename OutputT,
	 //					typename InputBundleT, typename OutputBundleT,
	 template<class> class BundleT,
	 wc_mapper::mode mode>
	 class WordCountMapperEvaluator
	 : public SingleInputTransformEvaluator<
	   //      			WordCountMapper<InputT, OutputT, InputBundleT, OutputBundleT, mode>,
	   WordCountMapper<InputT, OutputT, BundleT, mode>,
	   //RecordBitmapBundle<InputT>, RecordBitmapBundle<OutputT>> {
	   //      			InputBundleT, OutputBundleT> {
	   BundleT<InputT>, BundleT<OutputT>>
{

	using InputBundleT = BundleT<InputT>;
	using OutputBundleT = BundleT<OutputT>;
	using TransformT = WordCountMapper<InputT, OutputT, BundleT, mode>;

	//  using TransformT = WordCountMapper<InputT, OutputT, InputBundleT, OutputBundleT, mode>;
	//  using InputBundleT = RecordBitmapBundle<InputT>;
	//  using OutputBundleT = RecordBitmapBundle<OutputT>;

	public:
	bool evaluateSingleInput (TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle) override {

		// go through Records in input bundle (the iterator automatically
		// skips "masked" Records.
		for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
			//        trans->do_map(*it, output_bundle);
			TransformT::do_map(*it, output_bundle); /* static rocks! */
		}
		return true;
	}

	WordCountMapperEvaluator(int node)
		: SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }
};


#endif // WC_MAPPER_EVAL_H
