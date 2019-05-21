#ifndef NYT_WINDOWEDSUM_EVALULATOR_H
#define NYT_WINDOWEDSUM_EVALULATOR_H

#include "core/EvaluationBundleContext.h"
#include "core/SingleInputTransformEvaluator.h"
#include "NewYorkTaxi/NYTWindowedSum.h"

template <class InputT, class OutputT, class TransformT>
class NYTWindowedSumEvaluator
    :public SingleInputTransformEvaluator<TransformT,
                                          WindowsBundle<InputT>, RecordBundle<OutputT>>
{

    using InputBundleT = WindowsBundle<InputT>;
    using OutputBundleT = RecordBundle<OutputT>;

public:
    NYTWindowedSumEvaluator(int node)
		: SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) { }

	/* We cannot close windows based on the transform's local watermark.
	 *
	 * Instead, we should close windows based on "external" wm -- the upstream
	 * and inflight.
	 *
	 * This is because the watermark is computed by taking the internal state
	 * into account. This gurantees that even when the internal state is purged,
	 * the downstream transforms will not see data older than the watermark.
	 *
	 * Internal state will held back the transform's watermark.
	 * See RefreshWatermark() for more.
	 *
	 * @output_bundles: a container that receives new output bundles (can be
	 * multiple)
	 * NB:
	 * 1. caller holds mtx_watermark
	 * 2. only save output bundles to a local vec;
	 *    does not deposit to output or schedule consumers
	 * 	  (prefer to do so w/o conlock)
	 *
	 */

private:
	/* default way to gen a record for a window's results */
	virtual Record<OutputT> const makeRecord(OutputT & value,
                                             Window const & win)
	{
		return Record<OutputT>(value, win.window_end());
	}

public:
    /* Flush internal state (partial aggregation results) given the external
     * wm.
     *
     * @up_wm: the external wm.
     *
     * return: the local wm to be updated after the flushing.
     */
    ptime flushStateFixed(TransformT* trans, const ptime up_wm,
                          vector<shared_ptr<OutputBundleT>>* output_bundles,
                          EvaluationBundleContext* c, bool purge = true)
    {
        return min_date_time;
    }


private:
    ptime flushStateSliding(TransformT* trans, const ptime up_wm,
                            vector<shared_ptr<OutputBundleT>>* output_bundles,
                            EvaluationBundleContext* c)
    {
        return min_date_time;
    }

    ptime flushState(TransformT* trans, const ptime up_wm,
                     vector<shared_ptr<OutputBundleT>>* output_bundles,
                     EvaluationBundleContext* c, bool purge = true) override
    {
        if (trans->multi == 1)
            return flushStateFixed(trans, up_wm, output_bundles, c, purge);
        else
            return flushStateSliding(trans, up_wm, output_bundles, c);
    }

    /* We cannot close windows based on the transform's local watermark.
     * Instead, we should close windows based on inputs and upstream watermarks.
     *
     * This is because the watermark is computed by taking the internal state
     * into account. This gurantees that even when the internal state is purged,
     * the downstream transforms will not see data older than the watermark.
     *
     * Internal state will held back the transform's watermark.
     * See RefreshWatermark() for more.
     *
     * NB:
     * 1. caller holds mtx_watermark
     * 2. only deposit bundles; does not spawn consumers.
     * 	  (prefer to do so w/o conlock)
     *
     * @return: # of consumers to be spawned.
     *
     */
    /* does not do any summation, directly purge the requested window
       NOT TESTED. if in use, have to modify @watermark to be @up_wm. see func
       above.  */
    void fireAsWindowsBundle(ptime watermark, TransformT* trans,
                             EvaluationBundleContext *c,
                             bool purge = false) {
        PValue* v = trans->getFirstInput();
        assert(v);
        PTransform *upstream = v->producer;
        assert(upstream);

        ptime min_ts_flight = max_date_time;
        {
            unique_lock<mutex> lock(trans->mtx_watermark);
            for (auto && b : trans->inflight_bundles) {
                if (b->min_ts < min_ts_flight)
                    min_ts_flight = b->min_ts;
            }
        }

        ptime partial_min_ts = \
            min(v->min_ts, min(upstream->GetWatermarkSafe(), min_ts_flight));

        typename TransformT::windows_map \
            wm (trans->RetrieveWindows(purge, partial_min_ts));

        if (wm.size() == 0)   // this happens when no window is being closed.
            return;

        /* in the internal state @windows_map, each window may correspond to
         * multiple fragments;
         * in an output @WindowsBundle, each window only has one fragment.
         * We need to flatten @windows_map
         */
        auto output_bundle = make_shared<WindowsBundle<InputT>>();
        for (const auto & w : wm) {
            const auto & win = w.first;
            const auto & vec_pfrag = w.second;

            output_bundle->vals[win] = make_shared<WindowFragment<InputT>>(win);

            /* flatten the internal @windows_map */
            for (auto & pfrag : vec_pfrag) {
                output_bundle->vals[win]->merge(*pfrag);
            }
        }

        auto out = trans->getFirstOutput();
        assert(out);
        out->depositOneBundle(output_bundle, this->_node);
        c->SpawnConsumer(out, this->_node);
    }


    /* do local aggregation and then combine to the transform.
     * XXX the use of AggResultT (MT safe) may be an overkill, since we don't need
     * MT safe here. */
	bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

// //		using AggResultT = map<Window, OutputT, Window>;
//         using AggResultT = typename TransformT::AggResultT;

// 		AggResultT aggres;

// 		/* -- heavy lifting is done w/o locking the transform -- */
// 		for (auto it = input_bundle->vals.begin(); it != input_bundle->vals.end(); it ++) {
// 			auto && win = it->first;
// 			auto && pfrag = it->second;

// 			OutputT sum;
// 			TransformT::aggregate_init(&sum);

// 			for (auto && rec : pfrag->vals) {
// 				TransformT::aggregate(&sum, rec.data);
// 			}

// 			assert (aggres.count(win) == 0 && "bug? duplicate windows in same bundle");

// 			aggres[win] = sum;
// 		}

// 		/* -- combine -- */
//         trans->AddAggregatedResult(aggres);
//	  TransformT::AddAggregatedResult(aggres);
        std::cout << "window sum evaluator is called" << std::endl;
		return false; // no output bundle
	}

};

#endif /* NYT_WINDOWEDSUM_EVALULATOR_H */
