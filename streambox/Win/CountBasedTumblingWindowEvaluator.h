#ifndef COUNT_BASED_TUMBLING_WINDOW_EVALUATOR_H
#define COUNT_BASED_TUMBLING_WINDOW_EVALUATOR_H

#include "core/EvaluationBundleContext.h"
#include "core/TransformEvaluator.h"
#include "core/SingleInputTransformEvaluator.h"
#include "Win/CountBasedTumblingWindow.h"

#include "NewYorkTaxi/NYTRecord.h"

/* output bundle format is fixed */
template <class InputT, template<class> class InputBundleT_>
class CountBasedTumblingWindowEvaluator
    : public SingleInputTransformEvaluator<CountBasedTumblingWindow<InputT, InputBundleT_>,
                                           // InputBundleT_<InputT>, WindowsBundle<InputT>> {
                                           InputBundleT_<InputT>,
                                           WindowsKeyedBundle<std::pair<uint64_t, uint64_t>, WinKeyFragLocal_Std>> {

	// using KVPair = KVIn;
	// using KVPairOut = KVOut;
    using TransformT = CountBasedTumblingWindow<InputT, InputBundleT_>;
	// // using TransformT = WinKeyReducer<KVPair, InputWinKeyFragT, InternalWinKeyFragT, KVOut, OutputBundleT_>;
	// using InputBundleT = WindowsKeyedBundle<KVPair, InputWinKeyFragT>;


    using InputBundleT = InputBundleT_<InputT>;
    // using OutputBundleT = WindowsBundle<InputT>;
    using OutputBundleT = WindowsKeyedBundle<std::pair<uint64_t, uint64_t>, WinKeyFragLocal_Std>;


public:
    CountBasedTumblingWindowEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) {
    }

    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        // applied a trick here, misleading StatefulTransform to process it as timestamp
        Window win;
        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            trans->counter.fetch_add(1, std::memory_order_relaxed);
            // uint64_t maxCount = trans->MAX_COUNT;
            // bool exchanged = trans->counter.compare_exchange_strong(maxCount, 0);
            // if (exchanged) {
            //     std::cout << "exchanged at win id:" << trans->win_id
            //               << ", count: " << trans->counter
            //               << std::endl;
            //     trans->win_id.fetch_add(1, std::memory_order_relaxed);
            // }
            // win = Window(from_time_t(trans->win_id));
            size_t win_id = trans->counter / trans->MAX_COUNT; // transfer to seconds
            win = Window(from_time_t(win_id));

            output_bundle->add_record(win, *it);;
            // output_bundle->add_record(win, it->data);
        }

#if 0
        for (auto it = output_bundle->begin(); it != output_bundle->end(); ++ it) {
            auto win = it->first;
            auto vals = it->second.vals;
            for (auto val : vals) {
                auto k = val.first;
                auto v = val.second.vals[0]->size();
                EE("input.size: %lu, win %s, k: %lu, v: %lu",
                   input_bundle->size(),
                   to_simplest_string(win.window_start()).c_str(),
                   k, v);
            }
        }
#endif
        // return false;
        return true;
    }
};

#endif // COUNT_BASED_TUMBLING_WINDOW_EVALUATOR_H
