#ifndef LR_PARSER_EVAL_H
#define LR_PARSER_EVAL_H

#include "Values.h"
#include "LinearRoad/LinearRoadParser.h"
#include "core/SingleInputTransformEvaluator.h"

template <typename InputT, typename OutputT,
            template<class> class BundleT>
class LinearRoadParserEvaluator
            : public SingleInputTransformEvaluator<
	   LinearRoadParser<InputT, OutputT, BundleT>, BundleT<InputT>,
       BundleT<OutputT>>
{

	using InputBundleT = BundleT<InputT>;
	using OutputBundleT = BundleT<OutputT>;
	using TransformT = LinearRoadParser<InputT, OutputT, BundleT>;

	public:
	bool evaluateSingleInput (TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle) override {

      for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        auto out_record = trans->do_map(*it);
        // std::cout << out_record.data.toString() << std::endl;
        output_bundle->add_record(out_record);
      }
      // std::cout << "size of input bundle is: " << input_bundle->size() << std::endl;
      trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
      trans->byte_counter_.fetch_add(input_bundle->size()*60, std::memory_order_relaxed);
      return true;
	}

  LinearRoadParserEvaluator(int node)
		: SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }
};


#endif // LR_PARSER_EVAL_H
