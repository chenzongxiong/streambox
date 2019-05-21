#ifndef GREPMAPPEREVALUATOR_H
#define GREPMAPPEREVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "GrepMapper.h"

template <typename InputT, typename OutputT, template<class> class BundleT_>
class GrepMapperEvaluator
    : public SingleInputTransformEvaluator<GrepMapper<InputT,OutputT,BundleT_>,
      BundleT_<InputT>, BundleT_<OutputT>> {

	using TransformT = GrepMapper<InputT,OutputT>;
	using InputBundleT = BundleT_<InputT>;
	using OutputBundleT = BundleT_<OutputT>;

public:
	bool evaluateSingleInput (TransformT* trans,
      shared_ptr<InputBundleT> input_bundle,
      shared_ptr<OutputBundleT> output_bundle) override {

//#ifndef INPUT_ALWAYS_ON_NODE0
//  	assert(this->_node == input_bundle->node);
//#endif

    // go through Records in input bundle (the iterator automatically
    // skips "masked" Records.

//	std::cout << "2-1.GrepMapperEvaluator.h  evaluateSingleInput" << std::endl;  
	uint64_t cnt = 0;
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
//    	TransformT::do_map(*it, output_bundle); /* do_map stateful */
    	cnt += trans->do_map(*it, output_bundle);
    }
    if (cnt)
    	return true;
    else
    	return false;
	}

	GrepMapperEvaluator(int node)
	  : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }
};

#endif /* GREPMAPPEREVALUATOR_H */
