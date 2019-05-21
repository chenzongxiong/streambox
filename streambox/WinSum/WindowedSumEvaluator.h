/*
 (mostly) obsoleted by WindowedSumEvaluator1.h
 
 Simply update the internal "windows" state of the transform.
 does not actually execute the summation (which is deferred).
*/

#ifndef WINDOWED_SUM_EVAL_H
#define WINDOWED_SUM_EVAL_H

#include "core/EvaluationBundleContext.h"
#include "core/TransformEvaluator.h"
#include "WinSum/WinSumBase.h"

template <class InputT, class OutputT = InputT>
class WindowedSumEvaluator
		: public TransformEvaulator<WinSumBase<InputT, OutputT>> {

    using TransformT = WinSumBase<InputT, OutputT>;
    using InputBundleT = WindowsBundle<InputT>;
    using OutputBundleT = RecordBundle<OutputT>;

public:
  WindowedSumEvaluator(int node) {
  	ASSERT_VALID_NUMA_NODE(node);
  	this->_node = node;
  }

//  void OnNewUpstreamWatermark (ptime up_wm,
//      TransformT* trans, EvaluationBundleContext* c) override
//  {
//    fireAsRecordBundle(up_wm, trans, c, true); // aggregation
////    fireAsWindowsBundle(watermark, trans, c, true);
//  }

#ifndef NDEBUG
  ptime last_punc_ts;
  bool once= true;
#endif

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
   * 	  (prefer to do so w/o wmlock)
   *
   * @return: # of consumers to be spawned.
   *
   */
  int fireAsRecordBundle(ptime up_wm, TransformT* trans,
      EvaluationBundleContext *c,
      bool purge = false, bool wmlocked = false)
  {

    // Compute a "partial watermark" for closing windows.
    // We only consider upstream wm (snapshot), input, and pending
    // (excluding internal state)
    PValue* v = trans->getFirstInput(); // XXX deal with multiple ins
    assert(v);
    auto out = trans->getFirstOutput();
    assert(out);
    assert(out->consumer);

#if 0
    PTransform *upstream = v->producer;
    assert(upstream);
#endif

    int consumers = 0;

    ptime min_ts_flight = max_date_time;
    {
    	if (!wmlocked) {
				unique_lock<mutex> lock(trans->mtx_watermark);
				for (auto && b : trans->inflight_bundles) {
						if (b->min_ts < min_ts_flight)
							min_ts_flight = b->min_ts;
				}
    	} else {
    		for (auto && b : trans->inflight_bundles) {
						if (b->min_ts < min_ts_flight)
							min_ts_flight = b->min_ts;
				}
    	}
    }

    ptime partial_min_ts =
#if 0
        min(v->min_ts, min(upstream->watermark, min_ts_flight));
#endif
    	min(v->min_ts, min(up_wm, min_ts_flight));

    typename TransformT::windows_map \
      winmap (trans->RetrieveWindows(purge, partial_min_ts));

    if (winmap.size() == 0) {  // this happens when no window is being closed.
    	goto update_wm;
    }

    {
			// each output record also must carry a ts.
			auto output_bundle = make_shared<OutputBundleT>(
					winmap.size() * 4,  /* just a guess of the # of records in the output bundle */
					this->_node);

			// Iterate windows, each of which is associated with a vector
			// of InputT vectors.
			// We only do cilk_for *within* a window. If we want to do so *across*
			// windows, the outmost container has to support random access (e.g.
			// std::map won't work.)
			OutputT sum;
			for (auto it = winmap.begin(); it != winmap.end(); it++) {
				TransformT::aggregate_init(sum);
				for (auto&& vec : it->second) {
						for (auto&& v : vec->vals) {
							TransformT::aggregate(sum, v.data);
						}
				}
	//      E("sum is %lu (over %lu x vec)", sum, it->second.size());
				I("add one aggregated outpout (over %lu x vec)", it->second.size());
				output_bundle->add_record(Record<OutputT>(sum,
						it->first.window_end()));
			}

			out->depositOneBundle(output_bundle, this->_node);
			consumers ++;
	//    c->SpawnConsumer(out, this->_node);
    }

update_wm:
    /* also update local wm and emit punc */
    assert(wmlocked && "XXX must lock here");
    trans->watermark = min(v->min_ts, min(up_wm, trans->min_ts));

    auto new_punc = make_shared<Punc>(trans->watermark, this->_node);
		out->depositOneBundle(new_punc, this->_node);
		consumers ++;

		return consumers;
  }

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

  /* Upon receiving a punc from upstream.
   *
   * caller must hold mtx_watermark.
   * return: true if the punc needs to be passed to downstream immediately
   *
   * Since there's no inflight bundles, no actual work is done
   */
  bool process_punc(TransformT* trans, EvaluationBundleContext* c,
  		shared_ptr<Punc> punc, PValue *out)
  {

#ifndef NDEBUG /* extra check: punc must be monotonic */
  	if (once) {
  		once = false;
  	} else {
  		assert(last_punc_ts < punc->min_ts && "bug: punc regression");
  	}
  	last_punc_ts = punc->min_ts;
#endif

		W("%s (node %d) get punc: %s", trans->name.c_str(), this->_node,
				to_simple_string(punc->min_ts).c_str());

		assert(punc->min_ts > trans->watermark);
		return false;
  }

  /* lightweight: only "saving" incoming bundles into internal state. Leave heavy
   * aggregation work to the punc processing path.
   *
   * there's no "inflight" bundles then.
   *
   * that's why we wmlock the entire thing (except for scheduling consumers)
   */
  void evaluate(TransformT* trans, EvaluationBundleContext *c) {

    PValue* in1 = trans->getFirstInput();
    assert(in1);
    auto out = trans->getFirstOutput();
    assert(out);  /* even sink has valid @out */

    ASSERT_VALID_NUMA_NODE(this->_node);

    unique_lock<mutex> wmlock(trans->mtx_watermark);

    auto bundle_ptr = in1->getOneBundle(this->_node);
    assert(bundle_ptr);

		auto input_bundle = dynamic_pointer_cast<InputBundleT>(bundle_ptr);

		/* ---- punc path ----
		 * since we don't have inflight bundles, by the time we see punc,
		 * all previously enqueued bundles are consumed.
		 *
		 * also need to make sure no two puncs are handled concurrently, which is
		 * done by holding wm lock.
		 * */
		if (!input_bundle) {
			auto punc = dynamic_pointer_cast<Punc>(bundle_ptr);
			assert(punc && "neither bundle or punc, don't know what to do");

			process_punc(trans, c, punc, out);

			/* wm lock held. this deposit output bundles but don't spawn consumers */
			int num_consumers = fireAsRecordBundle(punc->min_ts, trans, c, true, true);

			wmlock.unlock();

			/* one wm and (optionally) one data bundle */
			assert(num_consumers == 1 || num_consumers == 2);

			for (int i = 0; i < num_consumers; i++) {
				c->SpawnConsumer(out, this->_node);
			}

			return;
		}

		/* ---- data bundle path ---- */

     /* we lock AddWindowsBundle() below: this ensures that data bundles
      * prior to a punc are procssed before a punc.
      * Instead of maintaining a "inflight" set of bundles and unlocking,
      * we lock the whole transform based on the assumption that AddWindowsBundle()
      * is fast. */

    trans->AddWindowsBundle(input_bundle);
    wmlock.unlock();

    I("after AddWindowsBundle(): transform internal #records = %lu: "
      "#windows = %lu",
        trans->GetRecordCount(), trans->GetWindowCount());
  }
};

#endif // WINDOWED_SUM_EVAL_H
