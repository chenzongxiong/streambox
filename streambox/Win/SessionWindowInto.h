#ifndef SESSIONWINDOWINTO_H
#define SESSIONWINDOWINTO_H

/* stateful.
 *
 * parallelization idea: each evaluator instance builds a local windows_map
 * (they touch their respective input bundles). no lock needed.
 * then each evaluator merges local windows_map with the transform's
 * windows_map. O(N)
 * this requires locking, but should be fast since no bundle
 * is touched.
 *
 */
template <typename InputT>
class SessionWindowInto : public PTransform {
	public:
		using time_duration = boost::posix_time::time_duration;
		using RecordT = Record<InputT>;
		// cannot be bitmapbundle: must be contig in mem
		using SessionWindowT = SessionWindow<InputT>;
		using WindowFragmentT = WindowFragment<InputT>;

		// intermediate maps before the final output, e.g. the maps passed by
		// evaluator (worker).
		//  using session_windows_map = map<SessionWindow, shared_ptr<WindowFragmentT>, Window>;

		// for final output, we don't need SessionWindow for saving ranges/bundles/etc
		using windows_map = map<Window, shared_ptr<WindowFragmentT>, Window>;

		using InputBundleT = RecordBundle<InputT>;
		using OutputBundleT = WindowsBundle<InputT>; /* contains WindowFragment */

	protected:
		// the age of the oldest internal record (the start of the oldest window)
		// *not* the watermark.
		ptime min_ts = max_date_time;

	public:
		const time_duration record_duration;

		// ordered windows, each containing pointers to input records
		// windows have *no* overlap
		set<SessionWindowT> windows;
		mutex _windows_mutex;

		// we hold pointers to the input bundles so they won't go away
		//  vector<shared_ptr<BundleBase>> pending_bundles;

		SessionWindowInto(string name, time_duration record_duration)
			: PTransform(name), record_duration(record_duration) { }


		/* debugging */
		void dump() {
			I("total %lu windows", this->windows.size());
			for (auto && w = this->windows.begin(); w!= this->windows.end(); w++) {
				I("window %s -- %s, %lu ranges",
						to_simplest_string(w->window_start()).c_str(),
						to_simplest_string(w->window_end()).c_str(),
						w->ranges.size());
				for (auto && range : w->ranges) {
					for (auto ptr = range.start; ptr <= range.end; ptr++)
						cout << (ptr->data) << " ";
					cout << endl;
				}
			}
		}

		void print(ostream & os) const {
			cout << "SessionWindowInto total windows: " << this->windows.size() << endl;
			for (auto && w = this->windows.begin(); w!= this->windows.end(); w++) {
				printf("\t\t window %s -- %s, %lu ranges: ",
						to_simplest_string(w->window_start()).c_str(),
						to_simplest_string(w->window_end()).c_str(),
						w->ranges.size());
				for (auto && range : w->ranges) {
					cout << " { ";
					for (auto ptr = range.start; ptr <= range.end; ptr++)
						cout << (ptr->data) << " ";
					cout << " } ";
				}
				cout << endl;
			}
		}

		/* accumulate a (task-local) bundle onto the transform-internal state.
		 * not as performant as MergeWindowsSet().
		 */
		void MergeRecordBundle(shared_ptr<InputBundleT> const bundle) {

			for (auto && recit : bundle->content) {
				try_add_record(&(recit), bundle); /* locked */
			}

			unique_lock<mutex> lock(mtx_watermark);

			if (!windows.empty()) {
				auto && ts = windows.begin()->window_start(); // the new min_ts after addition
				assert(watermark <= ts);

				// update the transform's internal min_ts (NOT the watermark)
				if (min_ts > ts) {
					VV("min_ts %s -> %s",
							to_simple_string(min_ts).c_str(),
							to_simple_string(ts).c_str());

					min_ts = ts;
				}
			}
		}

		/*
		 * merge an incoming set<SessionWindowT> to the transform's internal one
		 * This is preferred as multiple evaluators can build their own sets locally,
		 * and then call this func to merge their sets into the transform.
		 *
		 * potential ways to improve performance:
		 * 1. now we iterate over the incoming windows, so it is O(nlgn)...
		 *    since both sets are sorted, we can merge sort them in O(n)
		 * 2. N-way merge can be made even more efficient
		 */

		static struct WindowStrictComp<InputT> comp;

		void MergeWindowsSet(set<SessionWindowT> const & incoming) {

			unique_lock<mutex> lock(_windows_mutex);

			//    W("incoming windows %lu ...", incoming.size());

			for (auto it = incoming.begin(); it != incoming.end(); it++) {
				/* find the insertion/merge point.
				 *
				 * w0 is the 1st window that is strictly smaller than the incoming
				 * window *it (i.e. w0's end point is smaller than or equal
				 * to the start of @it). We can use lower_bound() on this because
				 * windows have no overlap.
				 *
				 * Note that we cannot use the default way to compare windows (using
				 * their start points).
				 */
				//        auto && w0 = std::lower_bound(windows.begin(), windows.end(), *it);
				auto && w0 = std::lower_bound(windows.begin(), windows.end(), *it, comp);

				if (w0 == windows.end()
						|| (it->window_end() <= w0->window_start())) {
					/* no window overlaps with @it; insert @it to the internal map
					 * directly. */
					windows.insert(*it); // copy
					continue;
				}

				// need merging. we extract (copy) w0, merge, and insert w0 back
				SessionWindowT win0 = *w0;
				//        W("w0 is"); cout << *w0;
				auto next = windows.erase(w0);
				auto ret = win0.merge(*it);
				assert(ret);

				// keep merging/erasing until we cannot proceed
				// this is necessary since the incoming window may "bridge" multiple
				// existing windows into one.
				while (next != windows.cend()) {
					//            W("win0 is"); cout << *w0;
					//            W("next is"); cout << *next;
					auto res = win0.merge_notbefore(*next);
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
			}

			assert(!windows.empty());
			// easy to find min_ts among the remaining windows as they're sorted
			unique_lock<mutex> wlock(mtx_watermark);
			min_ts = windows.begin()->window_start();

			//    W("after: windows %lu", windows.size());

			return;
		}

		/* DEPRECATED. okay for unit test purpose. see below.
		 *
		 * add a record to the session window(s) belonging to this transform.
		 update/merge windows are needed.
		 @return: true if any changes made to the windows

		 * (NB: it does not maintain min_ts)
		 */
		bool try_add_record(RecordT const * rec,
				shared_ptr<InputBundleT> const bundle) {
			assert(rec);

			/* we can be slow since we need to take the transform-wide lock.
			 * this may be improved, by, e.g. finer-grained locking, but adding
			 * individual records to transform is still not a good idea.
			 */

			unique_lock<mutex> lock(_windows_mutex);

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
			//    auto && w0 = std::lower_bound(windows.begin(), windows.end(), *rec,
			//        WindowRecComp<InputT>);

			auto && w0 = std::lower_bound(windows.begin(), windows.end(), *rec);

			if (w0 == windows.end()) {
				// the new record is strictly larger (i.e. no overlap) than all
				// existing windows. create a new window for the record.

#if 0
				// NB: we can hardly do emplace() and then add record to the
				// SessionWindow, since set is supposed to be immutable.
				SessionWindowT w (rec->ts, record_duration, record_duration);
				auto must_succ = w.try_add_record(rec);
				assert(must_succ);
				windows.insert(windows.end(), w);
#endif
				auto newit = windows.emplace_hint(windows.end(),
						rec->ts, record_duration, record_duration);

				auto must_succ = newit->try_add_record_notbefore(rec, bundle);
				assert(must_succ == SessionWindowT::RETCODE::OK);

				return true;
			}

			auto ret = w0->try_add_record_notbefore(rec, bundle);
			if (ret == SessionWindowT::RETCODE::OK) {     // fast path
				// window0's start unmodified.
				// its rank in the set does not change.
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
						rec->ts, record_duration, record_duration);
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
				//      assert(ret == SessionWindowT::RETCODE::OK);

				// keep merging/erasing until we cannot proceed
				while (next != windows.cend()) {
					auto res = win0.merge_notbefore(*next);
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

			/* (ret == SessionWindowT::RETCODE::TOO_LATE) */
			// bug?? since we find lower bound, rec can't be strictly later than
			// the window
			E("bug?"); abort();
			return false;
		}

		/* it's okay we return by copy, since @windows_map only contains shared_ptr
		 * to WindowFragment (not the actual value).
		 * note that @windows_map has different type than other transforms.
		 */
		windows_map RetrieveWindows(bool purge = false,
				ptime wm = max_date_time, int n = -1) {

			unique_lock<mutex> lock(_windows_mutex);

			windows_map ret;

			/*
			 * get the 1st window whose *end* time is no less than @wm.
			 * since @windows are sorted and not overlapping, all prior windows' end
			 * time must be smaller than @wm (any new records after @wm won't add to
			 * these windows) -- they need to be closed.
			 *
			 * when we close a window, we copy all the enclosed records from the
			 * the input bundle to the resultant window fragment. this can be intensive.
			 * so perhaps need to parallelize.
			 *
			 * see SessionWindow::operator<()
			 */

#if 0
			if (windows.empty())
				return ret;  // empty one
#endif
			W("dump all existing windows");
			for (auto w : windows) {
				//cout << w; //hym: we should override <<
			}
			W("--- end ---");

			auto endit = lower_bound(windows.begin(), windows.end(), wm);
			/* XXX: parallelize */
			for (auto it = windows.begin(); it != endit; it++) {
				Window w;  // {.start = it->start, .duration = it->duration};
				w.start = it->start;
				w.duration = it->duration;
				auto frag = make_shared<WindowFragmentT>(w);
				// materialize records in all ranges to the output win fragment
				for (auto && range : it->ranges) {
					for (auto ptr = range.start; ptr <= range.end; ptr++) {
						frag->add(*ptr); // copy records
					}
				}
				ret[w] = frag;
			}

			if (purge) {
				auto it = windows.begin();
				while (it != endit) {
					// will destroy session window, which will in turn destroy the
					// encompassed shared_ptr to input bundle, and therefore the
					// input bundle itself.
					it = windows.erase(it);
				}

				if (windows.size() == 0)
					min_ts = max_date_time;
				else // easy to find min_ts among the remaining windows as they're sorted
					min_ts = windows.begin()->window_start();
			}

			return ret; // by copy
		}
};


#endif /* SESSIONWINDOWINTO_H */
