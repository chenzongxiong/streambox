#ifndef WINDOWSBUNDLESINKEVALUATOR_H_
#define WINDOWSBUNDLESINKEVALUATOR_H_

#include "core/SingleInputTransformEvaluator.h"
#include "Sink/Sink.h"
#include "Values.h"
#include <map>
#include <atomic>
extern std::map<Window, ptime, Window> window_keeper;
extern uint64_t latencyStart;
extern atomic<bool> first;

using NanoSeconds = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;


/* InputT: the element type of the record bundle */
template <typename InputT>
class WindowsBundleSinkEvaluator
    : public SingleInputTransformEvaluator<WindowsBundleSink<InputT>,
//      WindowsKeyedBundle<InputT>, WindowsBundle<InputT>> {
                                           WindowsBundle<InputT>, WindowsBundle<InputT>> {   // sufficient for wc & topK?

	using TransformT = WindowsBundleSink<InputT>;
	using InputBundleT = WindowsBundle<InputT>;
//	using InputBundleT = WindowsKeyedBundle<InputT>;
	//using InputBundleT = WindowsKeyedBundle<KVPair>;
	using OutputBundleT = WindowsBundle<InputT>;

    std::atomic<bool> firstLatencyEvaluation;
public:

	WindowsBundleSinkEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }

    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        //XXX TransformT::printBundle(*input_bundle);
        // TransformT::report_progress(* input_bundle);
        TransformT::printBundle(*input_bundle);
        // std::cout << "EvaluateSingleInput" << std::endl;
        if (first) {
            first = false;
            uint64_t latencyEnd = std::chrono::duration_cast<NanoSeconds>(
                Clock::now().time_since_epoch())
                .count();
            double diff = (latencyEnd - latencyStart) / 1e6;
            std::cout << "Latency: " << diff << std::endl;
        }

        int size = 0;
        for (auto && win_frag: input_bundle->vals) {
            auto && win = win_frag.first;
            auto && pfrag = win_frag.second;

            for (auto &kv : pfrag->vals) {
                auto k = kv.data.first;
                auto v = kv.data.second;
                EE("win.start: %s, k: %lu, v: %lu", to_simplest_string(win.window_start()).c_str(), k, v);
            }
            size += pfrag->vals.size();
        }

        if (size > 0) {
            trans->record_counter_.fetch_add(size, std::memory_order_relaxed);
        }

        return false; /* no output bundle */
    }
};

#endif /* WINDOWSBUNDLESINKEVALUATOR_H_ */
