#include "WinKeyReducerEval.h"

/* -----------------------------------------
 * specialized for netmon -- using partitioned hash table
 * ----------------------------------------- */

/* relies on partitioned HT. the way we distribute work among workers differs from
 * @ReduceTopKParallel*/
template<>
void
WinKeyReducerEval<
		std::pair<creek::ippair, long>,  /* kv in */
		WinKeyFragLocal_Simple,
		WinKeyFrag_SimpleMT,
		std::pair<creek::ippair, long>, /* kv out */
		WindowsBundle /* out bundle */
	>::
    ReduceTopKParallel(TransformT* trans,
					typename TransformT::AggResultT const & winmap,
					vector<shared_ptr<OutputBundleT>>* output_bundles,
					EvaluationBundleContext* c, bool do_topK)
{
	assert(c);

	/*
	 * we have several *complete* windows and all the encompassed kvcontainer
	 * pairs. Now we need to reduce on per window, per key.
	 *
	 * throw the whole @winmap to workers, let them decide which part to work on.
	 *
	 */

//	const int max_tasks = std::thread::hardware_concurrency();
	const int max_tasks = c->num_workers;
	int total_tasks = max_tasks;

	atomic<int> flag2(0);

	/* totally for debugging info. should be cheap */
	int total_buckets = 0;
	for (auto && w : winmap) {
		auto && frag = w.second;
		total_buckets += frag->vals.num_maps;
	}

	// EE("will spawn %d tasks for %lu windows, total %d maps", total_tasks,
	// 		winmap.size(), total_buckets);

	/* a worker's result */
	using worker_res_t = unordered_map<Window,
																		shared_ptr<vector<KVPair>>,
																		WindowHash,
																		WindowEqual>;

	/* the global result after merging */
	using global_res_t = unordered_map<Window,
																		 vector<shared_ptr<vector<KVPair>>>,
																		 WindowHash, WindowEqual>;

	vector <std::future < worker_res_t >> futures;

	/* decide the start/end iterators of the windows_map */

	/* Distribute reduce work among multiple workers. Each worker gets a
	 * contig range of buckets (of k/vcontainers)
	 * A key does not span multi workers.
	 *
	 * XXX The main thread will block at the join() point. So we shoulda assign
	 * a portion of work to *this* thread too.
	 */
	k2_measure("reduce_start");
	for (int task_id = 0; task_id < total_tasks; task_id++) {
			/* NB: pass it_kvcontainer as copy */
			futures.push_back( // store a future
					c->executor_.push(
							/*
							 * @winmap: the whole workloads to reduce
							 */
							[trans, &winmap, &total_tasks, task_id, &flag2](int id)
							/* ----- lambda  ---- */
							{
								int expected = 0;
								int first = 0;
								if (flag2.compare_exchange_strong(expected, 1)) {
									k2_measure("1st_exec");
									first = 1;
								}

								worker_res_t wres;

#if 1   /* toggle this to debug: no reduce at all, wres empty */

								for (auto && w : winmap) { /* for each window ... */
									auto && win = w.first;
									auto && frag = w.second;
									auto && creekmap = frag->vals; /* creek::maps */

									auto range = get_range(creekmap.num_maps, total_tasks, task_id);

									VV("worker %d: bucket range %d -- %d",
													task_id, range.first, range.first + range.second);

									for (int i = range.first; i < range.first + range.second; i++) {
										/* in a map, go through all keys */
//													for (auto it = maps[i].begin(); it != maps[i].end(); ++it) {
										for (auto && kvs : creekmap.get_map(i)) {
//														auto && k = it->first;
//														auto && vcontainer = it->second;

											auto && k = kvs.first;
											auto && vcontainer = kvs.second;

											auto kvpair = trans->do_reduce(k, vcontainer);

											if (wres.count(win) == 0) {
												VV("worker: save res to win: %s",
														to_simplest_string(win.window_start()).c_str());
												wres[win] = make_shared<vector<KVPair>>();
												wres[win]->reserve(4096); /* reserve capacity in local res */
											}
											wres[win]->push_back(kvpair);
										}

										/* wres[win] can be a ptr... */
//													if (wres.count(win))
//														EE("worker %d: handle %lu keys",
//																					task_id, wres[win]->size());
									}
								}
#endif
								if (first)
									k2_measure("1st_done");

								return wres;
							}
						/* ---- end of lambda ---- */
				)
			);
	}

	k2_measure("task_cr_done");

	/* join all workers and combine their results (moving shared_ptrs) */

	global_res_t gres;

	VV("total %lu futures...", futures.size());

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
				output_bundle->add_record(win, kv);
#ifdef DEBUG
				cout << "k:[" << kv.first << "] count: " << kv.second << endl;
#endif
				minheap.pop();
			}
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
					output_bundle->add_record(win, kv);
				}
			}
		}
	}
	k2_measure("combine_end");
#endif
	k2_measure_flush();

	output_bundles->push_back(output_bundle);

	W("reduce done. global res has %lu windows", gres.size());
}

#if 0 /* no longer specialized */
template<>
ptime WinKeyReducerEval<std::pair<creek::ippair, long>,
								WinKeyFragLocal_Simple, WinKeyFrag_SimpleMT>::
flushState(TransformT* trans, const ptime up_wm,
		vector<shared_ptr<OutputBundleT>>* output_bundles,
		EvaluationBundleContext* c, bool purge)
{
	assert(output_bundles);

	k2_measure("flush_state");

	//hym: XXX how to know the min_ts??? XXX
	//hym: XXX shall we just flush state according to up_wm? I think so.

	// all closed windows are copied to winmap
	typename TransformT::AggResultT winmap;
//		ptime in_min_ts = trans->RetrieveState(&winmap, purge, up_wm);  // old impl.
	ptime in_min_ts = trans->RetrieveState(&winmap, up_wm, up_wm);

	// no window is being closed.
	if (winmap.size() == 0) {
		k2_measure_flush();
		if (in_min_ts == max_date_time) /* internal state empty */
			return up_wm;
		else /* nothing flushed but internal state exists */
			return min(in_min_ts, up_wm);
	}

	/* -- some windows are closed -- */
// 		ReduceSerial(trans, winmap, output_bundles);
// 		ReduceTopKSerial(trans, winmap, output_bundles);

// 		ReduceTopKParallelHT(trans, winmap, output_bundles, c,
//// 				true /* execute topK */);
// 				false /* don't execute topK */);

	ReduceTopKParallelMT(trans, winmap, output_bundles, c,
// 				true /* execute topK */);
			false /* don't execute topK */);

	return min(in_min_ts, up_wm);

}//end flushState
#endif
