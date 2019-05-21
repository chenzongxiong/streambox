#ifndef WINMAPPEREVALUATOR_H_
#define WINMAPPEREVALUATOR_H_

#include "core/SingleInputTransformEvaluator.h"
#include "Mapper/WinMapper.h"

template <typename InputT, typename OutputT, template<class> class BundleT_>
class WinMapperEvaluator
	: public SingleInputTransformEvaluator<WinMapper<InputT,OutputT,BundleT_>,
	  BundleT_<InputT>, WindowsBundle<OutputT>> {

	using TransformT = WinMapper<InputT,OutputT,BundleT_>;
	using InputBundleT = BundleT_<InputT>;
	using OutputBundleT = WindowsBundle<OutputT>; /* always produce WindowsBundle */

public:

	WinMapperEvaluator(int node) :
			SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(
					node) {
	}

	bool evaluateSingleInput(TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle) override {

    /* go through Records w/ iterator, which deals with the "masked" records.
     * let @WindowsBundle takes care of data layout.
     */
#if 0
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        // the time offset within a window
        long offset = ((*it).ts - trans->start).total_microseconds() \
            % (trans->window_size).total_microseconds();

        output_bundle->add_value(
              Window((*it).ts - microseconds(offset), trans->window_size),
              text_to_score((*it).data));
//              *it);
    }
#endif

  	uint64_t cnt = 0;

  	for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
    	cnt += trans->do_map(*it, output_bundle);
    }

    if (cnt)
    	return true;
    else
    	return false;
	}

private:
//	/* XXX should be in WinMapper.h. For code ease, leave it here */
//	long text_to_score(string_range const & str) {
//		/* XXX todo */
//		return 0;
//	}
};


#endif /* WINMAPPEREVALUATOR_H_ */
