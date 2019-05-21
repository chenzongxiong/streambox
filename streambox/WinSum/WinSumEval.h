/* * WindowedSumEvalulator1.h
 *
 *  A generic eval for all WinSum_XXX transforms.
 *  The version that does partial aggregation and later combine.
 *  Rebase the implementation on SingleInputTransformEvaluator. So that
 *  the wm/punc handling is unified.
 */

#ifndef WINDOWEDSUMEVALULATOR1_H_
#define WINDOWEDSUMEVALULATOR1_H_

#include "core/EvaluationBundleContext.h"
#include "core/SingleInputTransformEvaluator.h"
#include "WinSum/WinSumBase.h"

template <class InputT, class OutputT, class TransformT>
class WinSumEval
		:public SingleInputTransformEvaluator<TransformT,
		  WindowsBundle<InputT>, RecordBundle<OutputT>>
{

    using InputBundleT = WindowsBundle<InputT>;
    using OutputBundleT = RecordBundle<OutputT>;

public:
  WinSumEval(int node)
		: SingleInputTransformEvaluator<TransformT,InputBundleT,OutputBundleT>(node) { }

#ifndef NDEBUG
  ptime last_punc_ts;
  bool once= true;
#endif

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
#if 0 /* do aggregation for each record. bursty workload */
  void flushState(TransformT* trans, const ptime up_wm,
  		vector<shared_ptr<OutputBundleT>>* output_bundles,
  	 EvaluationBundleContext* c, bool conlocked = true,
  	 bool purge = true) override
  {
  	 assert(output_bundles);

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

		ptime min_ts_flight = max_date_time;
		{
			if (!conlocked) {
				unique_lock < mutex > lock(trans->mtx_watermark);
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

		typename TransformT::windows_map winmap(
				trans->RetrieveWindows(purge, partial_min_ts));

		if (winmap.size() == 0) {  // this happens when no window is being closed.
			goto update_wm;
		}

		{
			// each output record also must carry a ts.
			auto output_bundle = make_shared<OutputBundleT>(winmap.size() * 4, /* just a guess of the # of records in the output bundle */
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
				output_bundle->add_record(Record<OutputT>(sum, it->first.window_end()));
			}

			output_bundles->push_back(output_bundle);

			//			out->depositOneBundle(output_bundle, this->_node);
			//    c->SpawnConsumer(out, this->_node);
		}

update_wm:
		/* update local wm (punc will be emitted later) */
		assert(conlocked && "XXX must lock here");
		trans->watermark = min(v->min_ts, min(up_wm, trans->min_ts));
  }
#endif

#if 0
	/* Pack all records from @winmap into a bundle, which is put to @output_bundles
	 * factor this logic out so it can be overriden
	 * */
	virtual void packIntoBundle( typename TransformT::AggResultT const & winmap,
					vector<shared_ptr<OutputBundleT>>* output_bundles, int node)
	{
		// each output record also must carry a ts.
		auto output_bundle = make_shared<OutputBundleT>(
						winmap.size() * 4, /* just a guess of the # of records in the output bundle */
						this->_node);

		// Iterate windows, each of which is associated with an aggregation
		// result OutputT.
		for (auto it = winmap.begin(); it != winmap.end(); it++) {
// 				VV("add one aggregated outpout (over %lu x vec)", it->second->size());
			output_bundle->add_record(Record<OutputT>(
							it->second,
							it->first.window_end())
			);
		}
		output_bundles->push_back(output_bundle);
	}

	virtual Record<OutputT> const & makeRecord(
					typename TransformT::AggResultT::const_iterator it)
	{
		return Record<OutputT>(it->second, /* the value */
		                       it->first.window_end() /* the ts */
		);
	}
#endif

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
  	 xzl_assert(!trans->force_emit_tvpair); /* must call specialized version */

  	 assert(output_bundles);

#if 0
 		// Compute a "partial watermark" for closing windows.
 		// We only consider upstream wm (snapshot), input, and pending
 		// (excluding internal state)
 		PValue* v = trans->getFirstInput(); // XXX deal with multiple ins
 		//assert(v);
		if(!v){
			std::cout << "ERROR: PValue* v = trans->getFirstInput()" << std::endl;
			abort();
		}
 		auto out = trans->getFirstOutput();
 		out = out; /* no warning */
 		assert(out);
 		assert(out->consumer);
#endif

#if 0 /* we don't need this */
 		ptime min_ts_flight = max_date_time;
 		{
 			if (!conlocked) {
 				unique_lock < mutex > lock(trans->mtx_watermark);
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
 				min(v->min_ts, min(up_wm, min_ts_flight));
#endif

// 		ptime start, end;
// 		start = boost::posix_time::microsec_clock::local_time();

 		typename TransformT::AggResultT winmap;
		ptime in_min_ts = trans->RetrieveState(&winmap, purge, up_wm);

		// no window is being closed.
 		if (winmap.size() == 0) {
 			/* internal state empty, and we just observed @up_wm: the transform's local
 			 * wm should be updated to @up_wm.*/
 			if (in_min_ts == max_date_time)
 				return up_wm;
 			else /* nothing flushed but internal state exists */
 				return min(in_min_ts, up_wm);
 		}

//		end = boost::posix_time::microsec_clock::local_time();
//		cout << "XXXX retrieve takes " << (end - start).total_milliseconds() << " ms" << endl;
//		cout << "XXXX total " << winmap.size() << " windows";

//		start = boost::posix_time::microsec_clock::local_time();
		/*  some windows are closed */
 		{
 			// each output record also must carry a ts.
 			auto output_bundle = make_shared<OutputBundleT>(
 					winmap.size() * 4, /* just a guess of the # of records in the output bundle */
 					this->_node);

 			// Iterate windows, each of which is associated with an aggregation
 			// result OutputT.
 			for (auto it = winmap.begin(); it != winmap.end(); it++) {
// 				VV("add one aggregated outpout (over %lu x vec)", it->second->size());
// 				output_bundle->add_record(Record<OutputT>(
// 							it->second,
// 							it->first.window_end())
// 						);
			  output_bundle->add_record(makeRecord(it->second, it->first));
 			}

 			output_bundles->push_back(output_bundle);
 		}

//		end = boost::posix_time::microsec_clock::local_time();
//		cout << "XXXX add records takes " << (end - start).total_milliseconds() << " ms" << endl;
//		cout << "total " << winmap.size() << " windows";

 		return min(in_min_ts, up_wm);

// 		if (in_min_ts > up_wm)
// 			return up_wm;
// 		else
// 			return in_min_ts;

#if 0
 update_wm:
 		/* update local wm (punc will be emitted later) */
 		assert(conlocked && "XXX must lock here");
 		trans->watermark = min(v->min_ts, min(up_wm, trans->min_ts));
#endif
  }


private:
  ptime flushStateSliding(TransformT* trans, const ptime up_wm,
  		vector<shared_ptr<OutputBundleT>>* output_bundles,
  	 EvaluationBundleContext* c)
  {
  	assert(output_bundles);
		assert(trans->multi > 1);
		xzl_assert(!trans->force_emit_tvpair); /* must call specialized version */

 		typename TransformT::AggResultT winmap;
		ptime in_min_ts = trans->RetrieveStateSliding(&winmap, up_wm, trans->multi);

		// no window is being closed.
 		if (winmap.size() == 0) {
 			/* internal state empty, and we just observed @up_wm: the transform's local
 			 * wm should be updated to @up_wm.*/
 			if (in_min_ts == max_date_time)
 				return up_wm;
 			else /* nothing flushed but internal state exists */
 				return min(in_min_ts, up_wm);
 		}

 		/*  some windows are closed & returned.
 		 *
 		 *  XXX handle the cases where some windows are missing at the end of the
 		 *  returned window range (can't miss in the middle due to growth on
 		 *  demand.
 		 * */
 		{
 			// each output record also must carry a ts.
 			auto output_bundle = make_shared<OutputBundleT>(
 					winmap.size() * 4, /* just a guess of the # of records in the output bundle */
 					this->_node);

 			int n = winmap.size() - trans->multi;
 			assert(n > 0); /* XXX what no sufficient deltas are returned? */

 			for (int i = 0; i < n; i++) {
				// Iterate sliding windows, each of which encompasses @trans->multi deltas;
				// each delta has a window state (result OutputT).
				OutputT out;
				TransformT::aggregate_init(&out);

				/* We can be smarter: caching intermediate results etc */
				for (auto it = std::next(winmap.begin(), i);
									it != std::next(winmap.begin(), i + trans->multi); it++) {
					TransformT::combine(out, it->second);
				}

//				output_bundle->add_record(Record<OutputT>(
//							out,
//							/* the ts for the sliding window's result == the end of the
//							 * last delta in the sliding win.
//							 */
//							std::next(winmap.begin(), i + trans->multi - 1)->first.window_end()
//							)
//				);

			  output_bundle->add_record(
				  makeRecord(out, std::next(winmap.begin(), i + trans->multi - 1)->first)
			  );
 			}
 			output_bundles->push_back(output_bundle);
 		}

 		return min(in_min_ts, up_wm);
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
#if 0
  int fireAsRecordBundle(ptime up_wm, TransformT* trans,
      EvaluationBundleContext *c,
      bool purge = false, bool conlocked = false)
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
    	if (!conlocked) {
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
				TransformT::aggregate_init(&sum);
				for (auto&& vec : it->second) {
						for (auto&& v : vec->vals) {
							TransformT::aggregate(&sum, v.data);
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
    assert(conlocked && "XXX must lock here");
//    trans->watermark = min(v->min_ts, min(up_wm, trans->min_ts));
    trans->SetWatermark(min(v->min_ts, min(up_wm, trans->min_ts)));

    auto new_punc = make_shared<Punc>(trans->watermark, this->_node);
		out->depositOneBundle(new_punc, this->_node);
		consumers ++;

		return consumers;
  }
#endif

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

#if 0
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
#endif

#if 0 /* this will override parent's evaluate() */
  /* lightweight: only "saving" incoming bundles into internal state. Leave heavy
   * aggregation work to the punc processing path.
   *
   * there's no "inflight" bundles then.
   *
   * that's why we conlock the entire thing (except for scheduling consumers)
   */
  void evaluate(TransformT* trans, EvaluationBundleContext *c) {

    PValue* in1 = trans->getFirstInput();
    assert(in1);
    auto out = trans->getFirstOutput();
    assert(out);  /* even sink has valid @out */

    ASSERT_VALID_NUMA_NODE(this->_node);

    unique_lock<mutex> conlock(trans->mtx_watermark);

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

			conlock.unlock();

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
    conlock.unlock();

    I("after AddWindowsBundle(): transform internal #records = %lu: "
      "#windows = %lu",
        trans->GetRecordCount(), trans->GetWindowCount());
  }
#endif

  /* do local aggregation and then combine to the transform.
   * XXX the use of AggResultT (MT safe) may be an overkill, since we don't need
   * MT safe here. */
	bool evaluateSingleInput (TransformT* trans,
      shared_ptr<InputBundleT> input_bundle,
      shared_ptr<OutputBundleT> output_bundle) override {

//		using AggResultT = map<Window, OutputT, Window>;
			using AggResultT = typename TransformT::AggResultT;

		AggResultT aggres;

		/* -- heavy lifting is done w/o locking the transform -- */

		for (auto it = input_bundle->vals.begin(); it != input_bundle->vals.end();
				it ++)
		{
			auto && win = it->first;
			auto && pfrag = it->second;

			OutputT sum;
			TransformT::aggregate_init(&sum);

			for (auto && rec : pfrag->vals) {
				TransformT::aggregate(&sum, rec.data);
			}

			assert (aggres.count(win) == 0 && "bug? duplicate windows in same bundle");

			aggres[win] = sum;
		}

		/* -- combine -- */
    trans->AddAggregatedResult(aggres);
//	  TransformT::AddAggregatedResult(aggres);

		return false; // no output bundle
	}

};

//#include "WindowedSumEvaluator1-tvpair.h"

#endif /* WINDOWEDSUMEVALULATOR1_H_ */
