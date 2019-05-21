#ifndef SESSIONWINDOWINTOEVALUATOR_H
#define SESSIONWINDOWINTOEVALUATOR_H

/* @T: element type */
template <class T>
class SessionWindowIntoEvaluator
: public SingleInputTransformEvaluator<SessionWindowInto<T>,
	RecordBundle<T>, WindowsBundle<T>> {

		using RecordT = Record<T>;
		using TransformT = SessionWindowInto<T>;
		//  using InputBundleT = RecordBitmapBundle<T>;
		using InputBundleT = RecordBundle<T>;
		using OutputBundleT = WindowsBundle<T>;
		using SessionWindowT = SessionWindow<T>;

		// session windows produced by the evaluator that will be merged into
		// the transform's windows.
		set<SessionWindowT> windows;
		TransformT *trans = nullptr;

		/* windowing a record into the evaluator's *local* session window(s);
		 * update/merge windows are needed. Note we don't modify the transform's
		 * state -- no locking needed.
		 *
		 * XXX pretty much from SessionWindowInto::try_add_record()
		 *
		 * @return: true if any changes made to the windows
		 */
		bool try_add_record(RecordT const * rec,
				shared_ptr<InputBundleT> const bundle) {

			assert(rec);
			assert(trans);

			/*
			   since windows are sorted by their start time and they have no overlap:

			   get the first window that is *not* strictly smaller than rec's window.
			   then keep merging:
			   window0' <- window0 + rec's window
			   window0' <- window0' + window1
			   window0' <- window0' + window2
			   ...
			   stop when merge fails
			 */

			auto && w0 = std::lower_bound(windows.begin(), windows.end(), *rec);

			if (w0 == windows.end()) {
				// the new record is strictly larger (i.e. no overlap) than all
				// existing windows. create a new window for the record.

				// NB: okay to emplace() a session window and then add record to the
				// window, since try_add_record_notbefore() will not modify the set
				// key.
				auto newit = windows.emplace_hint(windows.end(),
						rec->ts, trans->record_duration, trans->record_duration);

				auto must_succ = newit->try_add_record_notbefore(rec, bundle);
				assert(must_succ == SessionWindowT::RETCODE::OK);

				return true;
			}

			auto ret = w0->try_add_record_notbefore(rec, bundle);
			if (ret == SessionWindowT::RETCODE::OK) {     // fast path
				/* window0's start unmodified.
				   its rank in the set does not change. */
				auto next = std::next(w0, 1);
				bool ret = false;
				while (next != windows.cend()) {
					auto && res = w0->merge_notbefore(*next);
					if (res == SessionWindowT::RETCODE::OK) {
						// merge okay. keep going. be careful we erase while iterating
						next = windows.erase(next);
						ret = true;
					} else { // merge failed
						// must be because that the next window is too far away
						assert (res == SessionWindowT::RETCODE::TOO_LATE);
						break;
					}
				}
				return ret;
			} else if (ret == SessionWindowT::RETCODE::TOO_EARLY) {
				// the new record is strictly smaller (i.e. no overlap) than all
				// existing windows. create a new window for it.
				auto ret = windows.emplace_hint(windows.begin(),
						rec->ts, trans->record_duration, trans->record_duration);
				auto must_succ = ret->try_add_record_notbefore(rec, bundle);
				assert(must_succ == SessionWindowT::RETCODE::OK);
				return true;
			}  else if (ret == SessionWindowT::RETCODE::MODIFY_START) {
				// slow path. have to extract window0 out of the set, modify it, and
				// insert it back. this is because set is immutable: we cannot modify
				// window0 in place.(emplace)
				SessionWindowT win0 = *w0;
				auto next = windows.erase(w0);
				auto ret = win0.try_add_record(rec, bundle);
				// since there's overlap, merge must succeed
				assert(ret);

				// keep merging/erasing until we cannot proceed
				while (next != windows.cend()) {
					auto res = w0->merge_notbefore(*next);
					if (res == SessionWindowT::RETCODE::OK) {
						// merge okay. keep going. be careful that we erase while iterating
						next = windows.erase(next);
					} else { // merge failed.
						// must be because that the next window is too far away
						assert (res == SessionWindowT::RETCODE::TOO_LATE);
						break;
					}
				}

				// then insert back win0 to the set. XXX this does copy
				windows.insert(windows.begin(), win0);

				return true;
			}

			/* (ret == SessionWindowT::RETCODE::TOO_LATE)
			 * since we find lower bound, rec can't be strictly later than the window
			 */
			E("bug?"); abort();
			return false;
		}

		bool evaluateSingleInput(TransformT* trans,
				shared_ptr<InputBundleT> input_bundle,
				shared_ptr<OutputBundleT> output_bundle) override {

			this->trans = trans;

			for (auto && it = input_bundle->begin();
					it != input_bundle->end(); ++it) {

				// since this transform owns the input bundle, and the input bundle
				// is immutable now, we assume the addresses of records inside
				// the input bundle will not change. we can safely take the addresses.
				this->try_add_record(&(*it), input_bundle);
			}

			trans->MergeWindowsSet(windows);

			return false;  // do not schedule downstream consumer yet
		}

		/*
		   emit the output bundle(s) by
		   materialize (and optionally purge) the transform's internal window state
		   into the actual window fragments.
		 */
		void fire(ptime up_wm, TransformT *trans, EvaluationBundleContext *c,
				bool purge = false) {

			/* Compute a "partial watermark" for closing windows.
			 * We only consider upstream, input, and pending
			 * (excluding internal state)
			 */
			PValue* v = trans->getFirstInput(); // XXX deal with multiple ins
			assert(v);
#if 0
			PTransform *upstream = v->producer;
			assert(upstream);
#endif

			ptime min_ts_flight = max_date_time;
			{
				unique_lock<mutex> lock(trans->mtx_watermark);
				for (auto && b : trans->inflight_bundles) {
					if (b->min_ts < min_ts_flight)
						min_ts_flight = b->min_ts;
				}
			}

			ptime partial_min_ts =
#if 0
				min(v->min_ts, min(upstream->watermark, min_ts_flight));
#endif
			min(v->min_ts, min(up_wm, min_ts_flight));

			/* materialize */
			typename TransformT::windows_map \
				winmap (trans->RetrieveWindows(purge, partial_min_ts));

			// this happens when no window is being closed or no records returned
			// from the transform
			if (winmap.empty())
				return;

			auto output_bundle = make_shared<OutputBundleT>();

			for (auto && w : winmap) {
				auto && fragptr = w.second;
				output_bundle->add(fragptr);
			}

			auto out = trans->getFirstOutput();
			assert(out);
			out->depositOneBundle(output_bundle);
			c->SpawnConsumer(out);
		}

		public:
		//  void OnNewUpstreamWatermark(ptime up_wm, TransformT* trans,
		//      EvaluationBundleContext* c) override {
		//    fire(up_wm, trans, c, true);
		//  }

		SessionWindowIntoEvaluator(int node)
			: SingleInputTransformEvaluator<TransformT,
			InputBundleT, OutputBundleT>(node) { }
	};

#endif /* SESSIONWINDOWINTOEVALUATOR_H */
