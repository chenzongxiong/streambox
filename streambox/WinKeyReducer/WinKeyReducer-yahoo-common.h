//
// Created by manuelrenz on 12.04.18.
//

#ifndef CREEK_WINKEYREDUCER_YAHOO_COMMON_H
#define CREEK_WINKEYREDUCER_YAHOO_COMMON_H

/* ----------------------- specialization for yahoo --------------------------------- */

/* for a particular window */
/* not in use: since eval do local reduce and then combine.
 * And it requires to iterate HT, not supported by creek::map */
template <>
WindowResultPtr const & MyWinKeyReducer::aggregate
        (WindowResultPtr * acc, InputT const & in)
{
#if 0
    assert(acc);
	assert(WindowEqual()((*acc)->w, in.w));

	for (auto && kvs : in.vals) {
		auto & key = kvs.first;
		auto & v_container = kvs.second;
		(*acc)->add_vcontainer_safe(key, v_container);
	}
#endif
    xzl_bug("not impl");
    return *acc;
}

/* XXX if @mine is empty, should we just swap pointers? */
template <>
WindowResultPtr const & MyWinKeyReducer::combine
        (WindowResultPtr & mine, WindowResultPtr const & others)
{
    xzl_assert(mine && others);
    xzl_assert(WindowEqual()(mine->w, others->w));
    for (auto && kvs : others->vals) {
        auto & key = kvs.first;
        auto & v_container = kvs.second;
        mine->add_vcontainer_safe(key, v_container);
    }

    return mine;
}

template<>
std::pair<uint64_t, uint64_t> MyWinKeyReducer::do_reduce
        (uint64_t const & key, InternalValueContainerT const & vcontainer) {

    auto end = vcontainer.cend(); /* avoid calling it in each iteration */

    long sum = 0;
    for (auto it = vcontainer.cbegin(); it != end; ++it) {
        this->record_counter_.fetch_add(1, std::memory_order_relaxed);
        sum += *it;
    }

    return make_pair(key, sum);
}

template<>
std::pair<uint64_t, uint64_t> MyWinKeyReducer::do_reduce_unsafe
        (uint64_t const & key, InputValueContainerT const & vcontainer) {
    auto end = vcontainer.cend(); /* avoid calling it in each iteration */

    long sum = 0;
    for (auto it = vcontainer.cbegin(); it != end; ++it) {
        this->record_counter_.fetch_add(1, std::memory_order_relaxed);
        // if (*it != 1) {
        //     // assert(0 && "must be 1");
        // }
        sum += *it;
    }
    return make_pair(key, sum);
}

template<>
bool MyWinKeyReducer::ReportStatistics(PTransform::Statstics *stat) {
    /* internal accounting */
    static unsigned long total_records = 0, total_bytes = 0;
    /* last time we report */
    static unsigned long last_bytes = 0, last_records = 0;
    static ptime last_check, start_time;
    static int once = 1;

    /* only care about records */
    total_records = this->record_counter_.load(std::memory_order_relaxed);

    ptime now = boost::posix_time::microsec_clock::local_time();

    if (once) {
        once = 0;
        last_check = now;
        start_time = now;
        last_records = total_records;
        return false;
    }

    boost::posix_time::time_duration diff = now - last_check;

    {
        double interval_sec = (double) diff.total_milliseconds() / 1000;
        double total_sec = (double) (now - start_time).total_milliseconds() / 1000;

        stat->name = this->name.c_str();
        stat->mbps = (double) total_bytes / total_sec;
        stat->mrps = (double) total_records / total_sec;

        stat->lmbps = (double) (total_bytes - last_bytes) / interval_sec;
        stat->lmrps = (double) (total_records - last_records) / interval_sec;

        last_check = now;
        last_bytes = total_bytes;
        last_records = total_records;
    }

    return true;
}
#endif //CREEK_WINKEYREDUCER_YAHOO_COMMON_H
