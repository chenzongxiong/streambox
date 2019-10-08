#ifndef NEXMARK_WINKEY_REDUCER_COMMON_H
#define NEXMARK_WINKEY_REDUCER_COMMON_H

/* XXX if @mine is empty, should we just swap pointers? */
template <>
WindowResultPtr const & NexmarkWinKeyReducer::combine(WindowResultPtr & mine, WindowResultPtr const & others)
{
    // std::cout << "[DBG] NexmarkWinKeyReducer Combine" << std::endl;

    xzl_assert(mine && others);
    xzl_assert(WindowEqual()(mine->w, others->w));

    for (auto && kvs : others->vals) {
        auto & key = kvs.first;
        auto & v_container = kvs.second;
        // std::cout << "key: " << key << ", size: " << v_container.size() << std::endl;
        // if (v_container.size() != 1 || v_container.vals[0]->size() != 1) {
        // std::cout << "Error" << std::endl;
        //     abort();
        // }
        // for (auto && v: v_container) {
        //     std::cout << "v: " << v << std::endl;
        // }
        // auto & v = v_container.vals[0];
        // std::cout << "v: " << v->front() << std::endl;
        // std::cout << "v: " << v_container.vals[0] << std::endl;
        mine->add_vcontainer_safe(key, v_container);
    }

    // std::cout << "mine.size: " << mine->vals.size() << ", other.size: " << others->vals.size() << std::endl;
    // uint64_t count = 0;
    // uint64_t key;
    // std::map<uint64_t, uint64_t> keeper;
    // for (auto && kvs : mine->vals) {
    //     auto & key = kvs.first;
    //     auto & v_container = kvs.second;
    //     // std::cout << "key: " << key << ", size: " << v_container.size() << std::endl;
    //     if (v_container.size() > 0) {
    //         for (auto && vecPtr: v_container.vals) {
    //             // count += v;
    //             // std::cout << "vecPtr->front(): " << vecPtr->front() << ", vecPtr.size: " << vecPtr->size() << std::endl;
    //             count += vecPtr->front();
    //         }
    //         if (count > 0) {
    //             keeper.insert(std::pair<uint64_t, uint64_t>(key, count));
    //         }
    //     }
    //     // std::cout << "count: " << count << std::endl;
    //     count = 0;
    //     // std::cout << "================================================================================" << std::endl;
    // }

    // std::cout << "keeper.size: "  << keeper.size() << std::endl;
    // uint64_t max = 0;
    // uint64_t key = 0;
    // for (auto &&kv : keeper) {
    //     if (max < kv.second) {
    //         key = kv.first;
    //         max = kv.second;
    //     }
    // }
    // std::cout << "key: " << key << ", value: " << max << std::endl;

    // mine->vals.clear();

    return mine;
}

template<>
std::pair<uint64_t, uint64_t> NexmarkWinKeyReducer::do_reduce(uint64_t const & key, InternalValueContainerT const & vcontainer) {

    auto end = vcontainer.cend(); /* avoid calling it in each iteration */

    long sum = 0;
    for (auto it = vcontainer.cbegin(); it != end; ++it) {
        this->record_counter_.fetch_add(1, std::memory_order_relaxed);
        sum += *it;
    }
    return make_pair(key, sum);
    // return make_pair(key, 99);
}

template<>
std::pair<uint64_t, uint64_t> NexmarkWinKeyReducer::do_reduce_unsafe(uint64_t const & key, InputValueContainerT const & vcontainer) {
    auto end = vcontainer.cend(); /* avoid calling it in each iteration */

    long sum = 0;
    for (auto it = vcontainer.cbegin(); it != end; ++it) {
        this->record_counter_.fetch_add(1, std::memory_order_relaxed);
        sum += *it;
    }
    return make_pair(key, sum);
    // return make_pair(key, 99);
}

template<>
bool NexmarkWinKeyReducer::ReportStatistics(PTransform::Statstics *stat) {
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

#endif //NEXMARK_WINKEY_REDUCER_COMMON_H
