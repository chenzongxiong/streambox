#ifndef WINDOWEDGREPMAPPEREVALUATOR_H
#define WINDOWEDGREPMAPPEREVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Mapper/WindowedGrepMapper.h"

template <typename InputT, typename OutputT>
class WindowedGrepMapperEvaluator
	: public SingleInputTransformEvaluator<WindowedGrepMapper<InputT,OutputT>,
    WindowsBundle<InputT>, WindowsBundle<OutputT>> {

	using TransformT = WindowedGrepMapper<InputT,OutputT>;
	using InputBundleT = WindowsBundle<InputT>;
	using OutputBundleT = WindowsBundle<OutputT>;

public:

	WindowedGrepMapperEvaluator(int node) :
			SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(
					node) {
	}

	bool evaluateSingleInput(TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle) override {

		uint64_t cnt = 0;
		for (const auto & w : input_bundle->vals) {
			const auto & win = w.first;
			const auto & pfrag = w.second;
			for (const auto & rec : pfrag->vals) {
				//cout << string(rec.data.data, rec.data.len);
				//I("record has %lu bytes", rec.data.len);
				cnt += trans->do_map(win, rec, output_bundle);
			}
		}

		if (cnt) {
			return true;
		} else {
			return false;
		}
	}
};

#endif /* WINDOWEDGREPMAPPEREVALUATOR_H */
