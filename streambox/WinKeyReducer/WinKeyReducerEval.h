#ifndef WINDOWED_KEY_RED_EVAL_H
#define WINDOWED_KEY_RED_EVAL_H

#define K2_NO_MEASUREMENT 1

extern "C" {
#include "measure.h"
}

#include "core/EvaluationBundleContext.h"
//#include "TransformEvaluator.h"
#include "core/SingleInputTransformEvaluator.h"
#include "WinKeyReducer/WinKeyReducer.h"
//#include "WinKeyReducerEval-yahoo-common.h"

/* support different output bundle formats. if @RecordBundle, the Window ts is encoded in the 1st item of the
 * value pair
 */


template <class KVIn,
					template<class> class InputWinKeyFragT,
					/* can be specialized based on key/val distribution */
					template<class> class InternalWinKeyFragT,
	        class KVOut,
	        template<class> class OutputBundleT_   /* WindowsBundle or RecordBundle */
					>
class WinKeyReducerEval
	: public SingleInputTransformEvaluator<
					WinKeyReducer<KVIn, InputWinKeyFragT, InternalWinKeyFragT, KVOut, OutputBundleT_>,
				WindowsKeyedBundle<KVIn, InputWinKeyFragT>,
	      OutputBundleT_<KVOut>			/* reduce results. often small */
		>
{

	using KVPair = KVIn;
	using KVPairOut = KVOut;
	using TransformT = WinKeyReducer<KVPair, InputWinKeyFragT, InternalWinKeyFragT, KVOut, OutputBundleT_>;
	using InputBundleT = WindowsKeyedBundle<KVPair, InputWinKeyFragT>;
	// the output should also be indexed by window.
	using OutputBundleT = OutputBundleT_<KVOut>;

public:
	WinKeyReducerEval(int node)
		: SingleInputTransformEvaluator<TransformT,InputBundleT,OutputBundleT>(node) { }

#ifndef NDEBUG
	ptime last_punc_ts;
	bool once= true;
#endif

private:
	/* converter: result record --> output record.
	 * default: result reord == output record
	 * can be specialized based on the need for output (e.g. assemble a record where first == window, etc.)*/
	Record<KVOut> const makeRecord(KVPair const & value, Window const & win) {
			return Record<KVOut>(value, win.window_end());
	}

	void ReduceSerial(TransformT* trans,
			typename TransformT::AggResultT const & winmap,
			vector<shared_ptr<OutputBundleT>>* output_bundles)
	{
 		/* reduction: serial execution
 		 *
 		 * hym: since there are not too many windows need to be closed
			 so serial execution should be fine, and we don't need
			 to use NUMA thread pool
 		 */

		k2_measure("ser_reduce_start");

		auto output_bundle = make_shared<OutputBundleT>(
						winmap.size() * 4, /* just a guess of the # of records in the output bundle */
				 		this->_node);

		for (auto && w : winmap) {
			auto && win = w.first;
			auto && pfrag = w.second;
			for (auto && kvcontainer : pfrag->vals) {
				 auto && k = kvcontainer.first;
				 auto && vcontainer = kvcontainer.second;
				 auto kvpair = trans->do_reduce(k, vcontainer);
//				 output_bundle->add_record(win, kvpair);
					output_bundle->add_record(win, makeRecord(kvpair, win));

				 //debugging
//				 cout << kvpair.first << ": " << kvpair.second << endl;
				 //res[win].push_back(kvpair);
				 //sum += kvpair.second;
			}
		}
		output_bundles->push_back(output_bundle);

		k2_measure("ser_reduce_end");
		k2_measure_flush();
	}

	void ReduceTopKSerial(TransformT* trans,
			typename TransformT::AggResultT const & winmap,
			vector<shared_ptr<OutputBundleT>>* output_bundles)
	{
		/* do reduction and stage the resultant kvpairs locally;
		 * run topK over the local kvpairs and produce the output bundle.
		 */

		/* we need separate results based on windows: across windows there can
		 * be duplicate keys.
		 */
		using MinHeap = priority_queue<KVPair,
																	vector<KVPair>,
																	ReducedKVPairCompLess<KVPair>>;

		MinHeap minheap; /* for one window. */

		const unsigned long K = 10;

		auto output_bundle = make_shared<OutputBundleT>(
						winmap.size() * 4, /* just a guess of the # of records in the output bundle */
				 		this->_node);

		/* Reduce per window. Store the results into the minheap and output to bundle */
		for (auto && w : winmap) {
			auto && win = w.first;
			auto && frag = w.second;

			for (auto && kvcontainer : frag.vals) {
				 auto && k = kvcontainer.first;
				 auto && vcontainer = kvcontainer.second;

				 auto kv = trans->do_reduce(k, vcontainer);

				 if (minheap.size() < K) {
					 /* store the kvpair to corresponding minheap */
					 minheap.push(kv);
				 } else {
						assert(minheap.size() == K);
						/* the new element larger than the smallest element in the heap */
						if (kv.second > minheap.top().second) {
							minheap.pop();
							minheap.push(kv);
						}
				 }
			}

			/* dump & empty the minheap */
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
		}
		cout << " -------------- \n";

		output_bundles->push_back(output_bundle);
	}

private:
	/* task_id: zero based.
	 * return: <start, cnt> */
	static pair<int,int> get_range(int num_items, int num_tasks, int task_id) {

		/* not impl yet */
		xzl_bug_on(num_items == 0);

		xzl_bug_on(task_id > num_tasks - 1);

		int items_per_task  = num_items / num_tasks;

		/* give first @num_items each 1 item */
		if (num_items < num_tasks) {
			if (task_id <= num_items - 1)
				return make_pair(task_id, 1);
			else
				return make_pair(0, 0);
		}

		/* task 0..n-2 has items_per_task items. */
		if (task_id != num_tasks - 1)
			return make_pair(task_id * items_per_task, items_per_task);

		if (task_id == num_tasks - 1) {
			int nitems = num_items - task_id * items_per_task;
			xzl_bug_on(nitems < 0);
			return make_pair(task_id * items_per_task, nitems);
		}

		xzl_bug("bug. never reach here");
		return make_pair(-1, -1);
	}

	void ReduceTopKParallel(TransformT* trans,
				typename TransformT::AggResultT const & winmap,
				vector<shared_ptr<OutputBundleT>>* output_bundles,
				EvaluationBundleContext* c, bool do_topK = true);
    void ReduceTopKParallel2(TransformT* trans,
                            typename TransformT::AggResultT const & winmap,
                            vector<shared_ptr<OutputBundleT>>* output_bundles,
                            EvaluationBundleContext* c, bool do_topK = true);
	/* parallel execution with hashtable partitioning
	 * using map's bucket interface. this speeds up task creation (no need iterate hashmap),
	 * but bucket interface seems to make things much much slower.
	 *
	 * this actuall works but it fails compilation because -MT HT has not bucket interface
	 * */
#if 0
	void ReduceTopKParallelHT(TransformT* trans,
					typename TransformT::AggResultT const & winmap,
					vector<shared_ptr<OutputBundleT>>* output_bundles,
					EvaluationBundleContext* c, bool do_topK = true)
		{
			assert(c);

			/*
			 * we have several *complete* windows and all the encompassed kvcontainer
			 * pairs. Now we need to reduce on per window, per key.
			 *
			 * throw the whole @winmap to workers, let them decide which part to work on.
			 *
			 */

			const int max_tasks = std::thread::hardware_concurrency();
			int total_tasks = max_tasks;

			atomic<int> flag2(0);

			/* totally for debugging info. should be cheap */
			int total_buckets = 0;
			for (auto && w : winmap) {
				auto && frag = w.second;
				total_buckets += frag->vals.unsafe_bucket_count();
			}

			// EE("will spawn %d tasks for %lu windows, total %d buckets", total_tasks,
			// 		winmap.size(), total_buckets);

			/* a worker's result */
			using worker_res_t = unordered_map<Window,
																				shared_ptr<vector<KVPair>>,
																				WindowHash,
																				WindowEqual>;

		/* Distribute reduce work among multiple workers. Each worker gets a
		 * contig range of keys (and their corresponding v containers).
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

									/* advance iterator and increment the key counter, terminate
									 * whichever reaches its end.
									 *
									 * write reduce results a local window map */
									int keys_processed = 0;
									for (auto it = it_kvcontainer; it != it_end; ++it) {
										auto && k = (*it).first;
										auto && vcontainer = (*it).second;
										auto win = it.get_current_window();

										auto kvpair = trans->do_reduce(k, vcontainer);
//										auto kvpair = make_pair("hello", 1);  // dbg

											VV("worker %d: bucket range %d -- %d",
															task_id, range.first, range.first + range.second);

											for (int i = range.first; i < range.first + range.second; i++) {
												/* in a bucket, go through all keys */
												for (auto it = m.unsafe_cbegin(i); it != m.unsafe_cend(i); it++) {
													auto && k = it->first;
													auto && vcontainer = it->second;

													auto kvpair = trans->do_reduce(k, vcontainer);

													if (wres.count(win) == 0) {
														VV("worker: save res to win: %s",
																to_simplest_string(win.window_start()).c_str());
														wres[win] = make_shared<vector<KVPair>>();
														wres[win]->reserve(4096); /* reserve capacity in local res */
													}
													wres[win]->push_back(kvpair);
												}
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
#endif

	/* relies on partitioned HT -- to be specialized. */
	void ReduceTopKParallelMT(TransformT* trans,
						typename TransformT::AggResultT const & winmap,
						vector<shared_ptr<OutputBundleT>>* output_bundles,
						EvaluationBundleContext* c, bool do_topK = true);
#if 0
			{
				assert(c);

				/*
				 * we have several *complete* windows and all the encompassed kvcontainer
				 * pairs. Now we need to reduce on per window, per key.
				 *
				 * throw the whole @winmap to workers, let them decide which part to work on.
				 *
				 */

				const int max_tasks = std::thread::hardware_concurrency();
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
													for (auto && kvs : creekmap.maps[i]) {
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
#endif

	/* Flush internal state (partial aggregation results) given the external
	 * wm.
	 *
	 * @up_wm: the external wm.
	 * return: the local wm to be updated after the flushing.
	 *
	 * Ref: WindowedSumEvaluator1::flushState()
	 *
	 * To be specialized
	 */
public:
	ptime flushState(TransformT* trans, const ptime up_wm,
			vector<shared_ptr<OutputBundleT>>* output_bundles,
			EvaluationBundleContext* c, bool purge = true) override
	{
        // std::cout << "flush state is called" << std::endl;
		assert(output_bundles);

		k2_measure("flush_state");

		//hym: XXX how to know the min_ts??? XXX
		//hym: XXX shall we just flush state according to up_wm? I think so.

		// all closed windows are copied to winmap
		typename TransformT::AggResultT winmap;
//		ptime in_min_ts = trans->RetrieveState(&winmap, purge, up_wm);  // old impl.
		ptime in_min_ts = trans->RetrieveState(&winmap, up_wm, up_wm);
        // std::cout << "winmap.size: " << winmap.size() << std::endl;
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
        ReduceTopKParallel(trans, winmap, output_bundles, c, true);
        // std::cout << "finshed flush state" << std::endl;
		return min(in_min_ts, up_wm);

	}//end flushState

	/* Open the bundle and do transform with InputT.
	 * Ref: the old AddWindowsKeyedBundle()
	 */
	bool evaluateSingleInputSimple (TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle)
	{

		/* Type of @input_bundle->vals happens to be the same
		 * as AggResultT. So just feed it.
		 * todo: local reduction of @frag */
		trans->AddAggregatedResult(input_bundle->vals);

#if 0
		for (auto && windowed : input_bundle->vals) {
			auto & win = windowed.first;
			auto & frag = windowed.second;

			/* todo: local reduction of @frag */
			trans->AddAggregatedResult(frag);
		}
#endif

		return false;
	}

	bool evaluateSingleInputReduce (TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle)
	{
		/* NB: local and intput_bundle only have to have same vcontainer */
		typename TransformT::LocalAggResultT local;

		/* single thread reduction */
		for (auto && w : input_bundle->vals) {
			auto && win = w.first;
			auto && frag = w.second;
			for (auto && kvcontainer : frag.vals) {
				 auto && k = kvcontainer.first;
				 auto && vcontainer = kvcontainer.second;
				 auto kvpair = trans->do_reduce_unsafe(k, vcontainer);
				 if (local.find(win) == local.end()) { /* shared ptr was not init */
					 TransformT::local_aggregate_init(&local[win]);
				 }
                 local[win]->add_kv_unsafe(kvpair, vcontainer.min_ts);
                 // local[win]->add_kv_unsafe(kvpair, vcontainer.min_ts, vcontainer.max_ts);
			}
		}

#if 0

        ptime frag_min_ts = max_date_time;
        ptime frag_max_ts = min_date_time;
        size_t counter = 0;

        for (auto & win_res : local) {
            auto & win = win_res.first;
            auto & res = win_res.second;

            if (frag_min_ts > res->min_ts) {
                frag_min_ts = res->min_ts;
            }
            if (frag_max_ts < res->max_ts) {
                frag_max_ts = res->max_ts;
            }
            counter ++;
        }
        // long frag_diff1 = (min_ts - frag_min_ts).total_milliseconds();
        // long frag_diff2 = (max_ts - frag_max_ts).total_milliseconds();
        // cout << "frag_diff1: " << frag_diff1
        //      << ", frag_diff2: " << frag_diff2
        //      << ", counter: " << counter
        //      << ", max_ts - min_ts: " << (frag_max_ts - frag_min_ts).total_milliseconds()
        //      << endl;

        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        long diff1 = (now - max_ts).total_milliseconds();
        if (diff1 < 0) {
            assert(0 && "diff1 is less than 0");
        }
        if (min_ts > max_ts) {
            if (max_ts != min_date_time && min_ts != max_date_time) {
                assert(0 && "min_ts > max_ts");
            }
        }
#endif
		trans->AddAggregatedResult(local);

#if 0
		/* debugging */
		int total_keys = 0; /* may out of multiple windows */
		for (auto && w : local) {   // partial_results_
            auto win = w.first;
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

			total_keys += frag->vals.size();
		}
		EE("total_keys %d", total_keys);
#endif
        // return false;
        // cout << "evaluate single input reduce : input_bundle.size: " << input_bundle->vals.size() << std::endl;
        // cout << "evaluate single input reduce : output_bundle.size: " << output_bundle->vals.size() << std::endl;
        return true;
	}

	/* the data bundle path */
	bool evaluateSingleInput (TransformT* trans,
				shared_ptr<InputBundleT> input_bundle,
				shared_ptr<OutputBundleT> output_bundle)  override
	{
		return evaluateSingleInputReduce(trans, input_bundle, output_bundle);
	}

}; //end WinKeyReducerEval

#endif // WINDOWED_KEY_RED_EVAL_H
