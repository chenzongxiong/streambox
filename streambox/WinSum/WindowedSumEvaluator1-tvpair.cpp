/*
 * WindowedSumEvaluator1-tvpair.cpp
 *
 *  will produce ts,value pair, where ts is generated based on window and
 *  value is window result.
 *
 *  workaround: as WinSumEval's template specialization for tweet
 */

#ifndef WINDOWEDSUMEVALUATOR1_TVPAIR_H_
#define WINDOWEDSUMEVALUATOR1_TVPAIR_H_

#include "WinSumEval.h"

using OutputT = creek::tvpair;

template<>
ptime WinSumEval<long, creek::tvpair>::flushStateFixed(TransformT* trans, const ptime up_wm,
		vector<shared_ptr<OutputBundleT>>* output_bundles,
	 EvaluationBundleContext* c, bool purge)
{
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
				/* it->first is window; it->second is data (tvpair).
				 * we fill the 1st of tvpair with ts */
				it->second.first = it->first.start.ticks() + it->first.duration.ticks();
				output_bundle->add_record(Record<OutputT>(
							it->second,
							it->first.window_end())
						);
			}

			output_bundles->push_back(output_bundle);

			//			out->depositOneBundle(output_bundle, this->_node);
			//    c->SpawnConsumer(out, this->_node);
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


template<>
ptime WinSumEval<long, creek::tvpair>::flushStateSliding(TransformT* trans, const ptime up_wm,
		vector<shared_ptr<OutputBundleT>>* output_bundles,
	 EvaluationBundleContext* c)
{
	assert(output_bundles);
	assert(trans->multi > 1);

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

				auto & w = std::next(winmap.begin(), i + trans->multi - 1)->first;
				out.first = w.start.ticks() + w.duration.ticks(); /* encode the window info */
				output_bundle->add_record(Record<OutputT>(
							out,
							/* the ts for the sliding window's result == the end of the
							 * last delta in the sliding win.
							 */
							w.window_end()
							)
				);


//				output_bundle->add_record(Record<OutputT>(
//							out,
//							/* the ts for the sliding window's result == the end of the
//							 * last delta in the sliding win.
//							 */
//							std::next(winmap.begin(), i + trans->multi - 1)->first.window_end()
//							)
//				);


			}
			output_bundles->push_back(output_bundle);
		}

		return min(in_min_ts, up_wm);
}

#endif /* WINDOWEDSUMEVALUATOR1_TVPAIR_H_ */
