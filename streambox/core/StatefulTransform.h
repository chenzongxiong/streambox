/*
 * StatefulTransform.h
 *
 *  Created on: Nov 26, 2016
 *      Author: xzl
 *
 *
 * An evaluator can aggregate its local input which is not yet submitted to
 * the transform.
 * later, after locking the transform, it can combine the local aggregation
 * result
 * into the transform's partial result, which we assume is fast.
 *
 * To close windows, partial results become complete and get emitted.
 * (the window info may be discarded)
 *
 * Common patterns:
 *  - has a internal state that is organized into fixed-length windows
 *  - eval path: add to the internal state
 *  - flush path: retrieve multiple windows of internal state; (optional)
 *  					    flush them; and output bundles
 */

#ifndef STATEFULTRANSFORM_H_
#define STATEFULTRANSFORM_H_

#define K2_NO_MEASUREMENT 1

#include "config.h"

extern "C" {
#include "measure.h"
}

#include "Transforms.h"
#include "Values.h"

#include <fstream>
#include <map>
extern std::map<Window, ptime, Window> window_keeper;
extern ofstream latency_file;
/* WindowResultT: one window's agg result. Not exposing to downstream
 * transforms.
 *
 * NB: WindowResultT will be *copied* back/forth in submitting/retrieving
 * states, so it'd better to be pointers.
 *
 */
template <typename TransformT, typename InputT, typename WindowResultT = InputT,
          typename LocalWindowResultT = WindowResultT>
class StatefulTransform : public PTransform {

  public:
    /* (partial) agg results for a set of windows.
     *
     * this is submitted by evaluator to trans, and also retrieved by evaluator
     * from trans.
     *
     * ordered but need lock protection.
     *
     * XXX Individual window may provide its own way (more precise) to track
     * its min ts. XXX
     * */
    using AggResultT = map<Window, WindowResultT, Window>;
    //	using AggResultT = unordered_map<Window, WindowResultT, Window>; /* need
    //a hash func */
    using LocalAggResultT = map<Window, LocalWindowResultT, Window>;
    using time_duration = boost::posix_time::time_duration;

  private:
    /* Internal agg results (state) organizes as a seq of windows.
     *
     * For this, we tried tbb::concurrent_unordered_map or
     * tbb::concurrent_vector
     * and abandon them. The reason is that the concurrent growth is limited.
     * See comments inline. */
    using SavedResultT = vector<WindowResultT>;

    //  atomic<long> window_count_;   // not very useful.

    /* take reader lock for aggregation (update a specific window's res) and
     * grow
     * (add a new window).
     * take writer lock for GC
     */
    boost::shared_mutex mtx_partial_results_;
    Window start_win_; /* window of the 0th item in the partial results */
    std::atomic<long> unflushed_; /* the lowest index of unflushed items */
    SavedResultT partial_results_;

    /* the min #items to grow the internal state on demand.
     * be careful: aggressive growth can be very slow if ResultT is a
     * sophisticated
     * d/s. */
    //  const long partial_results_grow_ = 4 * 1024 * 1024 / sizeof
    //  (WindowResultT);
    // const long partial_results_grow_ = 32;
    const long partial_results_grow_ = 32;
    // const long partial_results_grow_ = 1;

    /* the threshold for GC internal stale items */
    const long partial_results_gc_ = 10 * partial_results_grow_;
    // const long partial_results_gc_ = partial_results_grow_;

    bool is_start_win_assigned(void) {
        return (this->start_win_.duration !=
                boost::posix_time::time_duration(0, 0, 0));
    }

    void clear_start_win(void) {
        /* clean both start & duration */
        this->start_win_.start = boost::posix_time::time_duration(0, 0, 0);
        this->start_win_.duration = boost::posix_time::time_duration(0, 0, 0);
    }

    /* The duration of a window.
     * Unlike start_win, once assigned, these won't be changed/invalidated.
     * So no lock needed to read  */
    long win_duration_ms_ = -1;
    time_duration win_duration_;

    void set_win_duration(time_duration const &d) {
        if (win_duration_ms_ == -1) { /* unassigned */
            win_duration_ = d;
            win_duration_ms_ = d.total_milliseconds();
        } else {
            xzl_assert(win_duration_ms_ == d.total_milliseconds());
        }
    }

    /* Given a window, return the corresponding index in the internal state.
     * caller must hold rlock of the internal state
     * Will do extensive check.
     * */
    long win_to_index(Window const &win) {

        //		xzl_assert(win.duration == this->win_duration_);
        //
        //		long index = (win.start -
        //this->start_win_.start).total_milliseconds()
        //							/
        //this->win_duration_ms_;
        //
        //		xzl_assert(index >= 0);
        //
        //		xzl_assert(win.start.total_milliseconds()
        //				== index * this->win_duration_ms_
        //						+
        //this->start_win_.start.total_milliseconds());
        //		return index;

        auto start = this->start_win_.start;

        xzl_assert(win.duration == this->win_duration_);

        long index =
            (win.start - start).total_milliseconds() / this->win_duration_ms_;

        xzl_assert(index >= 0);

        if (!(win.start.total_milliseconds() ==
              index * this->win_duration_ms_ /* 2000 */
                  + start.total_milliseconds())) {
            EE("bug: start_win.start %lu win.start.total_milliseconds() %lu "
               "diff %lu index %lu duration %lu",
               start.total_milliseconds(), win.start.total_milliseconds(),
               (win.start - start).total_milliseconds(), index,
               this->win_duration_ms_);
        }

        xzl_assert(win.start.total_milliseconds() ==
                   index * this->win_duration_ms_ + start.total_milliseconds());
        return index;
    }

    /* Given @wm, return the highest index of the internal window
     * whose end <= @wm.
     *
     * return -1 if such window does not exist
     *
     * NB: @wm does not necessarily align to window boundary.
     *
     * Caller must hold rlock of the internal state
     *
     */
    long wm_to_index(ptime const &wm, bool must_align = false) {

        long span = (wm - this->start_win_.window_start()).total_milliseconds();
        long &winspan = this->win_duration_ms_;
        // std::cout << "wm_to_index: span = "
        //           << span << ", winspan = "
        //           << winspan << ", span % winspan = "
        //           <<  (span % winspan)
        //           << std::endl;
        xzl_assert(span >= 0);
        xzl_assert(winspan && winspan != -1); /* must have a valid start win */

        if (must_align && (span % winspan != 0)) {
            bug("does not align");
            return -1;
        }

        /* e.g. span==winspan: index = 0;
         * span == 1/2 winspan, index = -1;
         * span = 1.5 winspan, index = 0
         */
        long index = span / winspan - 1;

        return index;
    }

    /* round down the @wm to the window boundary.
     * caller needs rlock since we are reading this->win_duration_ms_
     * XXX align to microsecond?
     */
    ptime wm_rounddown(ptime const &wm) {

        long span = (wm - this->start_win_.window_start()).total_milliseconds();
        long &winspan = this->win_duration_ms_;

        xzl_assert(span >= 0);
        xzl_assert(winspan && winspan != -1); /* must have a valid start win */

        long offset = (span % winspan);

        return (wm - milliseconds(offset));
    }

    /* no lock needed -- aligned to microseconds.
     * nanoseconds() not available? */
    ptime wm_rounddown_epoch(ptime const &wm,
                             boost::posix_time::time_duration const &duration) {

        int64_t span = (wm - Window::epoch).total_microseconds();
        int64_t winspan = duration.total_microseconds();

        //		cout << to_simple_string(wm) << endl;

        xzl_assert(span >= 0);
        xzl_assert(winspan && winspan != -1);

        int64_t offset = (span % winspan);

        return (wm - microseconds(offset));
    }

#if 0
	bool is_wm_aligned(ptime const & wm) {

		long span = (wm - this->start_win_.start).total_milliseconds();
		long winspan = this->start_win_.duration.total_milliseconds();

		xzl_assert(winspan); /* must have a valid start win */

		return (span % winspan == 0) ? true : false;
	}
#endif

  public:
    StatefulTransform(string name) : PTransform(name), unflushed_(0) {
        /* so that we know start_win_ is unassigned */
        start_win_.duration = boost::posix_time::time_duration(0, 0, 0);
    }

    /* merge an incoming partial aggregation. Assume caller holds no _mutex.
     *
     * used by an evaluator to combine its local state (results) to the trans's
     * internal state.
     *
     * NB: HOT -- all evaluators call this.
     *
     * We expect to take the rlock all the time. When necessary (e.g. adding
     * new windows to the internal state), temporarily upgrade rlock to wlock
     * and downgrade after done.
     */
    void AddAggregatedResult(LocalAggResultT const &in) {
        /* reader lock.
         * see comment for mtx_partial_res. */

        boost::shared_lock<boost::shared_mutex> rlock(
            this->mtx_partial_results_);

        /* partial_res is a vector. need to compute its index */
        for (auto &win_res : in) {
            auto &win = win_res.first;
            auto &res = win_res.second;

        restart:
            xzl_assert(rlock.owns_lock());

            if (!is_start_win_assigned()) {
                /* no internal state. start_win_ should be garbage.
                 * set one.
                 * Need to upgrade rlock to wlock. */
                rlock.unlock();
                {
                    boost::unique_lock<boost::shared_mutex> wlock(
                        this->mtx_partial_results_);

                    /* Some concurrent writers may have slipped in & set the
                     * start
                     * window. check again */
                    if (!is_start_win_assigned()) {

                        xzl_assert(unflushed_ == 0);

                        /* This evaluator comes with @win, but a later eval may
                         * come w/
                         * even earlier window. Reset the start window
                         * (backwards) is a
                         * pain.
                         *
                         * Since this @win should not be too far away from the
                         * earliest window, we set @start_win_ to be a bit ahead
                         * of this @win.
                         *
                         * This also helps dealing with the initial a few
                         * sliding windows.
                         */

                        //						start_win_.start =
                        //win.start - seconds(10);
                        //						start_win_.duration
                        //= win.duration;

                        start_win_ = Window(
                            wm_rounddown_epoch(win.window_start() - seconds(10),
                                               win.duration),
                            win.duration);
                        // TODO: GC cannot happen if we force time_t to be 0
                        // if not set to 0, assertion in line 353 maybe fail
                        // start_win_ = Window(from_time_t(0), win.duration);
                        set_win_duration(win.duration);

                        // EE("start win unassigned. incoming win %s, set start
                        // win %s",
                        // 		to_simple_string(win.window_start()).c_str(),
                        // 		to_simple_string(start_win_.window_start()).c_str()
                        // );

                        /* XXXX adjust @unflushed_?? XXXX */
                    }
                }
                /* Let the wlock and downgrade to rlock.
                 * There's a chance of GC which cleans up start_win_ again. have
                 * to
                 * restart.
                 */
                rlock.lock();
                goto restart;
            }

            xzl_assert(win.start >= this->start_win_.start &&
                       "reset the start window? see comment above.");

            long index = win_to_index(win);

            //            EE("index is %ld", index);

            /* Need to grow internal state?
             *
             * Concurrent vector supports concurrent growth. but the problem
             * seems that elements under construction may be exposed to
             * concurrent
             * readers --> segfault
             *
             * With concurrent vector, we may overgrow (concurrent) but that's
             * okay. */
            if (index >= (long)partial_results_.size()) {
                /* Upgrade rlock to wlick. Will temporarily let go rlock.
                 * need to detect whether GC happens in between.
                 *
                 * if GC happens, start_win_ will advance.
                 */
                auto winstart = this->start_win_.start;
                bool need_restart = false;

                rlock.unlock();
                /* a concurrent growth/ GC may happen here. check size again. */
                {
                    // EE("try to grow");

                    boost::unique_lock<boost::shared_mutex> wlock(
                        this->mtx_partial_results_);

                    if (winstart == this->start_win_.start &&
                        is_start_win_assigned()) {
                        /* GC didn't happen. but concurrent growth may have
                         * happened.
                         * Need to exclusively get the state size again.
                         */
                        auto s = (long)partial_results_.size();
                        if (index >= s) {
                            /* grow by a reasonable delta */
                            //							WindowResultT
                            //acc;
                            //							TransformT::aggregate_init(&acc);
                            long delta =
                                max(index - s + 1, partial_results_grow_);
                            VV("grow internal state, delta=%lu", delta);
                            //							partial_results_.grow_by(delta,
                            //acc); /* copy ctor for the new items */
                            //							partial_results_.resize(s + delta,
                            //acc);
                            partial_results_.resize(s + delta);
                            for (int k = s; k < s + delta;
                                 k++) { /* init each window result explicitly */
                                TransformT::aggregate_init(
                                    &partial_results_[k]);
                            }
                        }
                    } else {
                        /* GC happened and invalidated @index. need restart */
                        need_restart = true;
                    }
                } /* Downgrade wlock -> rlock */
                /* a concurrent GC may happen here and invalidate @index. */
                rlock.lock();
                if (need_restart)
                    goto restart;
                if (winstart != this->start_win_.start) /* GC happened */
                    goto restart;
            }

            xzl_assert(rlock.owns_lock());

            /* finally we locate the window in the internal state */

            TransformT::combine(partial_results_[index],
                                res); /* this has to be MT safe */
        } /* next window in the incoming AggResult */
    }

    /* the core routine for state management.
     * Supports different ts for flush & retrieve
     *
     * @out: return the state of multiple windows
     *
     * @wm_ret: any window whose end <= @wm_ret will be returned
     * @wm_flush: any window whose end <= @wm_flush will be flushed (ie dropped)
     *           must be window aligned.
     *
     * @return: the min ts among the remaining internal state, or max_date_time
              if no internal state


          newer                  older
               3      2    1   0
          <--+---+---+---+---+---+
                           |         |     |
                           wmret     |     unflushed_=0
                                                          wm_flush

     ret: win0/1/2  flush: win0/1

           Expect to take & hold rlock. Upgrade to wlock only when necessary
           we don't need the lock unless we do GC: updating the index unflushed_
     can
           be atomic.
           In fact, we do not need to care for concurrent flushers. A later punc
     is only
           retrieved *after* the prior one is consumed.
     */
    ptime RetrieveState(AggResultT *out, ptime const &wm_ret,
                        ptime const &wm_flush) {

        auto &ret = *out;

        /* only take reader lock. start_win, partial_results_ won't change */
        boost::shared_lock<boost::shared_mutex> rlock(
            this->mtx_partial_results_);

        long s = partial_results_.size();
        // EE("paritial results size: %lu", s);
    restart:
        auto unflushed = unflushed_.load();

        /* internal results empty, iff unflushed == size */
        if (unflushed == s)
            return max_date_time;

        /* Now internal results are non empty.
         *
         * Calculate ranges.
         * retrieval range [unflushed, ret_index].
         * flush range [unflushed, flush_index].
         *
         * @wm_to_index() makes sure returned indices >= 0
         */

        long ret_index = wm_to_index(wm_ret);
        long flush_index = wm_to_index(wm_flush, true);
        long max_index = s - 1;
        xzl_assert(max_index >= 0);

        /* ret/flush ranges are bound by the actual # windows */
        ret_index = min(ret_index, max_index);
        flush_index = min(flush_index, max_index);
        // EE("flush %ld -- %ld, ret %ld -- %ld", unflushed, flush_index, unflushed,
        //   ret_index);


        if (ret_index < flush_index) {
            bug("bug? although this is supported");
        }

        if (ret_index < unflushed || flush_index < unflushed) {
            bug("bug? try to access GC'd state");
        }


        /* commit the new value of unflushed_ */
        if (!unflushed_.compare_exchange_strong(unflushed, flush_index + 1))
            goto restart; /* unflushed_ has changed. restart */

        /* we have advanced unflushed_ to flush_index + 1. Now no other thread
         * will touch the range to be flushed. */

        ptime ret_min_ts;

        I("flush %ld -- %ld, ret %ld -- %ld", unflushed, flush_index, unflushed,
          ret_index);

        k2_measure("retindex_computed");

        for (int i = unflushed; i <= ret_index; i++) {
            /* compute the window for each returned value. */
            Window w(start_win_.window_start() +
                     milliseconds(i * win_duration_ms_),
                     start_win_.duration);
            //  		ret[w] = partial_results_[i]; /* copy */
            ret.emplace(w, partial_results_[i]);
        }

        k2_measure("ret assembled");

        unflushed = flush_index + 1; /* could equal size(). see below. */

        /* determine the min ts for the remaining internal state */
        if (unflushed == s) {
            /* no window left in state */
            ret_min_ts = max_date_time;
        } else if (unflushed < s) {
            /* XXX okay? this window may contain zero result.
             * XXX If each individual window tracks its internal min ts, we
             * should
             * use that min ts, which will be tighter.
             * */
            ret_min_ts = start_win_.window_start() +
                         milliseconds(unflushed * win_duration_ms_);
        } else {
            xzl_assert(0 && "bug: unflushed exceeds end of latest window?");
        }

        /* GC: too many stale items in the partial results.
         * Since we return copies in flushing (or smart pointers), GC here is
         * safe.
         *
         * After GC, @unflushed_ becomes the 0-th item in the internal state.
         */
        if (unflushed > partial_results_gc_) {
            /* upgrade rlock -> wlock */
            rlock.unlock();
            {
                boost::unique_lock<boost::shared_mutex> wlock(
                    this->mtx_partial_results_);

                if (unflushed_ > partial_results_gc_) {
                    VV("GC start. #items=%ld", unflushed_.load());

                    if (unflushed_ == s) { /* empty internal state */

                        ptime curr_max_ts = min_date_time;
                        // cout << "window_keeper.size: " <<
                        // window_keeper.size() << endl;
                        boost::posix_time::ptime now1 =
                            boost::posix_time::microsec_clock::local_time();
                        for (auto it = window_keeper.begin();
                             it != window_keeper.end(); it++) {
                            auto win = it->first;
                            auto max_ts = it->second;
                            if (curr_max_ts < max_ts) {
                                curr_max_ts = max_ts;
                            }
                        }
                        boost::posix_time::ptime now =
                            boost::posix_time::microsec_clock::local_time();
                        // cout << "latency: " << (now -
                        // curr_max_ts).total_microseconds() / 1000.0 << " ms"
                        // << endl;
                        // latency_file << "latency: " << (now -
                        // curr_max_ts).total_microseconds() / 1000.0 << " ms"
                        // << endl;

                        clear_start_win();
                        window_keeper.clear();
                    } else { /* advance start_win_ */
                        start_win_.start +=
                            milliseconds(unflushed_ * win_duration_ms_);
                    }

                    partial_results_.erase(partial_results_.begin(),
                                           partial_results_.begin() +
                                               unflushed_);
                    unflushed_ = 0;
                    VV("GC done. #items=%ld", unflushed_.load());
                }
            }
        } /* GC done */

        return ret_min_ts;
    }

    /* To support sliding windows.
     *
     * @out: return the state of a superwindow, of which the end boundary is
     * determined by @wm.
     *
     * @multi: len of superwindow = multi * len of window
     *
     * as a side effect of retrieval, the 1st delta of the superwindow is
     * flushed.
     *
     * Limitation:
     * - Delta size should be ms-aligned.
     * - For code simplicity, we require a sliding window size is a multiply of
     * deltas.
     */
    ptime RetrieveStateSliding(AggResultT *out, ptime wm_up, int multi) {
        ptime t = wm_rounddown(wm_up); /* end of the last retrieved window */

        /* start of the latest sliding window. anything earlier than that won't
         * show up in later sliding windows. */
        ptime flush = t - milliseconds(multi * this->win_duration_ms_);
        xzl_assert(t > ptime(min_date_time) +
                           milliseconds(multi * this->win_duration_ms_) &&
                   "ptime value underflow.");

        return RetrieveState(out, wm_up, flush);
    }

    /* Old impl of window retrieval. Consider using the other RetrieveState()
     *
     * Get the recent N "windows" and purge them from the internal state
     if n unspecified, get all existing windows.

     @ purge: whether to purge the internal window states.

     @return: the min ts among the remaining internal state, or max_date_time
     if no internal state

     Expect to take & hold rlock. Upgrade to wlock only when necessary
     we don't need the lock unless we do GC: updating the index unflushed_ can
     be atomic.
     In fact, we do not need to care for concurrent flushers. A later punc is
     only
     retrieved *after* the prior one is consumed.

     */

    ptime RetrieveState(AggResultT *out, bool purge = false,
                        ptime wm = max_date_time, int n = -1) {
        auto &ret = *out;

        /* only take reader lock. start_win, partial_results_ won't change */
        boost::shared_lock<boost::shared_mutex> rlock(
            this->mtx_partial_results_);

        xzl_assert(purge && "not implemented");

    restart:
        auto unflushed = unflushed_.load();

        /* internal results empty, iff unflushed == size */
        if (unflushed == (long)partial_results_.size())
            return max_date_time;

        /* Now internal results are non empty.
         *
         * Calculate the flush range.
         * flush range [unflushed, high_index]
         */
        int high_index;
        int cnt;
        if (wm == max_date_time)
            cnt =
                partial_results_.size() - unflushed; /* all items. can be 0. */
        else {
            cnt = (wm - start_win_.window_start()).total_milliseconds() /
                      start_win_.duration.total_milliseconds() -
                  unflushed;
            /* e.g.
             * wm delta = 0.5x window: 	cnt = 0, high_index = -1
             * wm delta = 1x window: 		cnt = 1, high_index = 0
             * wm delta = 1.5x windows: cnt = 1, high_index = 0
             * */
        }

        //  	cout << "XXXX  cnt: " << cnt << "win duration: "
        //  start_win_.duration.total_milliseconds() << endl;
        //  	cout << "XXXX offset from window start"

        if (n != -1) { /* check with argument @n */
            cnt = std::min(cnt, n);
        }
        high_index =
            min(unflushed + cnt - 1, (long)(partial_results_.size() - 1));
        xzl_assert(high_index >= unflushed && "bug? flush GC'd range");

        /* commit the new value of unflushed_ */
        if (!unflushed_.compare_exchange_strong(unflushed, high_index + 1))
            goto restart; /* unflushed_ has changed. restart */

        /* we have advanced unflushed_ to high_index + 1. Now no other thread
         * will touch the range to be flushed. */

        ptime ret_min_ts;

        I("flush %ld -- %d", unflushed, high_index);

        for (int i = unflushed; i <= high_index; i++) {
            /* compute the window for each returned value. */
            Window w(
                start_win_.window_start() +
                    milliseconds(i * start_win_.duration.total_milliseconds()),
                start_win_.duration);
            //  		ret[w] = partial_results_[i]; /* copy */
            ret.emplace(w, partial_results_[i]);
        }

        unflushed = high_index + 1; /* could equal size(). see below. */

        /* determine the min ts for the remaining internal state */
        long s = partial_results_.size();
        if (unflushed == s) {
            /* no window left in state */
            ret_min_ts = max_date_time;
        } else if (unflushed < s) {
            /* XXX okay? the window may contain zero result.
             * XXX If each individual window tracks its internal min ts, we
             * should
             * use that min ts, which will be tighter.
             * */
            ret_min_ts = start_win_.window_start() +
                         milliseconds(unflushed *
                                      start_win_.duration.total_milliseconds());
        } else {
            xzl_assert(0 && "bug: unflushed too large?");
        }

        /* GC: too many stale items in the partial results.
         * Since we return copies in flushing (or smart pointers), GC here is
         * safe.
         *
         * After GC, @unflushed_ becomes the 0-th item in the internal state.
         */
        if (unflushed > partial_results_gc_) {
            /* upgrade rlock -> wlock */
            rlock.unlock();
            {
                boost::unique_lock<boost::shared_mutex> wlock(
                    this->mtx_partial_results_);

                if (unflushed_ > partial_results_gc_) {
                    EE("GC start. #items=%ld", unflushed_.load());

                    if (unflushed_ == s) { /* empty internal state */
                        clear_start_win();
                    } else { /* advance start_win_ */
                        start_win_.start += milliseconds(
                            unflushed_ *
                            start_win_.duration.total_milliseconds());
                    }

                    partial_results_.erase(partial_results_.begin(),
                                           partial_results_.begin() +
                                               unflushed_);
                    unflushed_ = 0;
                    EE("GC done. #items=%ld", unflushed_.load());
                }
            }
        } /* GC done */

        return ret_min_ts;
    }

    /* locking seems okay as long as this function is not hot. */
    unsigned long GetWindowCount() {
        xzl_assert(false && "todo");
        return -1;
    }

    /* #if 0 */
    /*   static WindowResultT const & aggregate_init(WindowResultT * acc); */

    /*   /\* aggregating *individual* input elements. to be specialized *\/ */
    /*   static WindowResultT const & aggregate(WindowResultT * acc, InputT
     * const & in); */

    /*   /\* combine the evaluator's (partial) aggregation results for a
     * particular window */
    /*    * to the tran's internal state. */
    /*    *\/ */
    /*   static WindowResultT const & combine(WindowResultT & mine,
     * WindowResultT const & others); */
    /* #endif */
};

#ifdef K2_NO_MEASUREMENT
#undef K2_NO_MEASUREMENT
#endif

#endif /* STATEFULTRANSFORM_H_ */
