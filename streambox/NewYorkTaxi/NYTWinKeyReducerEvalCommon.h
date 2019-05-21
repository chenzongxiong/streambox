#ifndef NYT_WINKEY_REDUCER_EVAL_COMMON_H
#define NYT_WINKEY_REDUCER_EVAL_COMMON_H


template<>
void NYTWinKeyReduerEval::ReduceTopKParallel(TransformT* trans,
                                             typename TransformT::AggResultT const & winmap,
                                             vector<shared_ptr<OutputBundleT>>* output_bundles,
                                             EvaluationBundleContext* c, bool do_topK)
{
    // std::cout << "Hello from reduce top k parallel" << std::endl;
    // return ;
    assert(c);

    /* parallel execution.
     *
     * we have several *complete* windows and all the encompassed kvcontainer
     * pairs. Now we need to reduce on per window, per key.
     * the kvcontainers across windows may be unbalanced.
     */

    int total_keys = 0; /* may out of multiple windows */
    for (auto && w : winmap) {
        auto && frag = w.second;
        total_keys += frag->vals.size();
    }

    /* ad hoc values */
    const int min_keys_per_task = 50;
//	const int max_tasks = std::thread::hardware_concurrency();
    const int max_tasks = c->num_workers;

    int keys_per_task;

    if (total_keys / max_tasks > min_keys_per_task)
        keys_per_task = total_keys / max_tasks;
    else
        keys_per_task = min_keys_per_task;

    const int total_tasks = total_keys / keys_per_task;

    int task_key_count = 0;
    atomic<int> flag2(0);

//		EE("total_keys %lu", total_keys);

#if 0 // dbg  -- examine the data to reduce
    EE("will spawn %d tasks for %lu windows, each taking %d keys, tasks in queues: %d", total_tasks,
       winmap.size(), keys_per_task, c->executor_.q.size());

    for (auto & w : winmap) {
		auto & win = w.first;
        auto && frag = w.second;

        for (auto &kv : frag->vals) {
            auto k = kv.first;
            auto vcontainer = kv.second;
            long sum = 0;
            for (auto it = vcontainer.cbegin(); it != vcontainer.cend(); ++it) {
                EE("win.start: %s, frag.size: %lu. k: %lu, v: %lu", to_simplest_string(win.window_start()).c_str(),
                   frag->vals.size(),
                   k, *it);
            }
        }

		EE("window: %s #items = %lu",
           to_simplest_string(win.window_start()).c_str(), w.second->vals.size());
	}
#endif

    /* a worker's result */
//		using worker_res_t = map<Window, shared_ptr<vector<KVPair>>, Window>;
    using worker_res_t = unordered_map<Window,
                                       shared_ptr<vector<KVPair>>,
                                       WindowHash,
                                       WindowEqual>;

    /* the global result after merging */
//		using global_res_t = map<Window, vector<shared_ptr<vector<KVPair>>>, Window>;
    using global_res_t = unordered_map<Window,
                                       vector<shared_ptr<vector<KVPair>>>,
                                       WindowHash, WindowEqual>;

    vector < std::future < worker_res_t > > futures;

    /* decide the start/end iterators of the windows_map */
    xzl_assert(winmap.begin()->second);
    window_map_iterator2<KVPair> it_kvcontainer(&winmap, winmap.begin(),
                                                winmap.begin()->second->vals.begin());
    auto lastbutone = winmap.end();
    lastbutone--;
    const window_map_iterator2<KVPair> it_end(&winmap, lastbutone,
                                              lastbutone->second->vals.end());

    /* Distribute reduce work among multiple workers. Each worker gets a
     * contig range of keys (and their corresponding v containers).
     * A key does not span multi workers.
     *
     * Below: start from 1st key in the 1st win;
     * pass the iterator to a worker.
     * Advance the iterator and count keys.
     * When having enough keys, pass the iterator to the next worker.
     *
     * Carefully guard at @it_end.
     *
     * XXX The main thread will block at the join() point. So we shoulda assign
     * a portion of work to *this* thread too.
     */
    k2_measure("reduce_start");
    while (it_kvcontainer != it_end) {

        if (task_key_count == 0) { /* start a new task */

//				W("this start win: %s",
//							to_simplest_string(it_kvcontainer.get_current_window().window_start()).c_str());

#if 1
            /* NB: pass it_kvcontainer as copy */
            futures.push_back( // store a future
                c->executor_.push(
                    /* @it_kvcontainer: the start of this key range.
                     * @it_end: the end of the global key range.
                     * @keys_per_task: how many keys to processed by this task.
                     */
                    [trans, it_kvcontainer, &it_end, &keys_per_task, &flag2](int id)
                    /* ----- lambda  ---- */
                    {
                        int expected = 0;
                        int first = 0;
                        if (flag2.compare_exchange_strong(expected, 1)) {
                            k2_measure("1st_exec");
                            first = 1;
                        }

                        worker_res_t wres;

                        VV("worker: my start win: %s",
                           to_simplest_string(it_kvcontainer.get_current_window().window_start()).c_str());

                        /* advance iterator and increment the key counter, terminate
                         * whichever reaches its end.
                         *
                         * write reduce results a local window map */
#if 1   /* toggle this to debug: no reduce at all, wres empty */
                        int keys_processed = 0;
                        for (auto it = it_kvcontainer; it != it_end; ++it) {
                            auto && k = (*it).first;
                            auto && vcontainer = (*it).second;
                            auto win = it.get_current_window();

                            auto kvpair = trans->do_reduce(k, vcontainer);
                            // for dbg: zero reduce workload, change do_reduce

                            if (wres.count(win) == 0) {
                                VV("worker: save res to win: %s",
                                   to_simplest_string(win.window_start()).c_str());
                                wres[win] = make_shared<vector<KVPair>>();
                                wres[win]->reserve(4096); /* reserve capacity in local res */
                            }
                            wres[win]->push_back(kvpair);

                            if (++keys_processed == keys_per_task)
                                break;
                        }
#endif
                        if (first)
                            k2_measure("1st_done");

                        return wres;
                    }
                    /* ---- end of lambda ---- */
                    )
                );
#endif
        }

        ++it_kvcontainer;
        VV("move, current win: %s",
           to_simplest_string(it_kvcontainer.get_current_window().window_start()).c_str());

        if (++task_key_count == keys_per_task)
            task_key_count = 0;
    }

    k2_measure("task_cr_done");

    /* join all workers and combine their results (moving shared_ptrs) */

    global_res_t gres;

    // EE("total %lu futures...", futures.size());

    /* XXX the main thread will block here XXX */
    for (auto && f : futures) {
        auto && wres = f.get();
        for (auto && w : wres) {
            auto && win = w.first;
            auto && vec_ptr = w.second;
            /* if win does not exist in @gres, this will create a new vec
             * in place */
            gres[win].push_back(vec_ptr);
        }
    }
    k2_measure("join_done");
    futures.clear();

//		k2_measure_flush();

    /* NB: @res is unsorted, e.g. by the value (count). if we want it to be
     * sorted, we may apply extra transforms (like merge sort) */
    k2_measure("reduce_end");
//		k2_measure_flush();

    /* consuming global_res_t */
    auto output_bundle = make_shared<OutputBundleT>(
        winmap.size() * 4, /* just a guess of the # of records in the output bundle */
        this->_node);

#if 1
    do_topK = false;
    if (do_topK) { /* topK: single core using minheap. */
        const unsigned long k = 10;
        priority_queue<KVPair, vector<KVPair>, ReducedKVPairCompLess<KVPair>> minheap;

        /* the output may contain multiple windows. output topK for each of
         * the window.
         */
        VV("global res has %lu windows", gres.size());
        for (auto && r : gres) {
            auto && win = r.first;
            auto && vec_of_vecptr = r.second;

            /* @vecptr is a shared ptr pointing to a vector that contains KV pairs.
             * it is a window's worth of reduce result produced by one worker.
             */
            for (auto & vecptr : vec_of_vecptr) {
                for (auto && kv : *vecptr) {
                    if (minheap.size() < k) {
                        minheap.push(kv);
                        continue;
                    } else {
                        assert(minheap.size() == k);

                        /* the incoming element > the smallest element in the heap */
                        if (kv.second > minheap.top().second) {
                            minheap.pop();
                            minheap.push(kv);
                        }
                    }
                }
            }

            /* dump & empty the minheap (XXX should do it in descending order) */
#ifdef DEBUG
            cout << "dump minheap: ";
#endif
            while (minheap.size() > 0) {
                auto && kv = minheap.top();
//				output_bundle->add_record(win, kv);
                output_bundle->add_record(win, makeRecord(kv, win));
#ifdef DEBUG
                cout << "k:[" << kv.first << "] count: " << kv.second << endl;
#endif
                minheap.pop();
            }
//				cout << " -------------- \n";
//				abort();
        }
    } else {  /* simply output all counts */
        VV("global res has %lu windows", gres.size());
        for (auto && r : gres) {
            auto && win = r.first;
            auto && vec_of_vecptr = r.second;

            /* @vecptr is a shared ptr pointing to a vector that contains KV pairs.
             * it is a window's worth of reduce result produced by one worker.
             */
            for (auto & vecptr : vec_of_vecptr) {
                for (auto && kv : *vecptr) {
                // EE("win.start: %s, k: %lu, v: %lu", to_simplest_string(win.window_start()).c_str(),
                //    // kv.data.first, kv.data.second);
                //    kv.first, kv.second);

// //					output_bundle->add_record(win, kv);
//                     std::cout << "first: " << kv.first
//                               << ", second: " << kv.second
//                               << std::endl;
                    output_bundle->add_record(win, makeRecord(kv, win));
                }
            }
        }
    }
    k2_measure("combine_end");
#endif
    k2_measure_flush();

    output_bundles->push_back(output_bundle);

    W("reduce done. global res has %lu windows", gres.size());
    // EE("reduce done. global res has %lu windows", gres.size());
}

#endif //NYT_WINKEY_REDUCER_EVAL_COMMON_H
