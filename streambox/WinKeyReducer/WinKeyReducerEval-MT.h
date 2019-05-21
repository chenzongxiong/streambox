//
// Created by xzl on 2/14/17.
//
// a partial specialization of WinKeyReducerEval.

/* -----------------------------------------
 * specialized for netmon -- using partitioned hash table
 * ----------------------------------------- */

/* relies on partitioned HT. the way we distribute work among workers differs from
 * @ReduceTopKParallel*/

#ifndef CREEK_WINKEYREDUCERMTEVAL_H
#define CREEK_WINKEYREDUCERMTEVAL_H

#include "WinKeyReducerEval.h"

template <class KVIn,
					class KVOut,
					template<class> class OutputBundleT_   /* WindowsBundle or RecordBundle */
				>
class WinKeyReducerEval<KVIn,
				WinKeyFragLocal_Simple, WinKeyFrag_SimpleMT,
				KVOut, OutputBundleT_>
	: public SingleInputTransformEvaluator<WinKeyReducer<KVIn,
					WinKeyFragLocal_Simple, WinKeyFrag_SimpleMT>,
					WindowsKeyedBundle<KVIn, WinKeyFragLocal_Simple>, /* input bundle */
					OutputBundleT_<KVOut>			/* reduce results. often small */
			>
{
	using InputWinKeyFragT = WinKeyFragLocal_Simple;
	using InternalWinKeyFragT = WinKeyFrag_SimpleMT;
	using KVPair = KVIn;
	using TransformT = WinKeyReducer<KVPair, InputWinKeyFragT, InternalWinKeyFragT>;
	using InputBundleT = WindowsKeyedBundle<KVPair, InputWinKeyFragT>;
	// the output should also be indexed by window.
	using OutputBundleT = OutputBundleT_<KVOut>;

public:
		WinKeyReducerEval(int node)
			: SingleInputTransformEvaluator<TransformT,InputBundleT,OutputBundleT>(node) { }

	void ReduceTopKParallelMT(TransformT* trans,
	                          typename TransformT::AggResultT const & winmap,
	                          vector<shared_ptr<OutputBundleT>>* output_bundles,
	                          EvaluationBundleContext* c, bool do_topK = true)
	{




	}


public:
	ptime flushState(TransformT* trans, const ptime up_wm,
	                 vector<shared_ptr<OutputBundleT>>* output_bundles,
	                 EvaluationBundleContext* c, bool purge = true) override
	{

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
			}
		}

		trans->AddAggregatedResult(local);

#if 0
		/* debugging */
		int total_keys = 0; /* may out of multiple windows */
		for (auto && w : local) {   // partial_results_
			auto && frag = w.second;
			total_keys += frag->vals.size();
		}
		EE("total_keys %d", total_keys);
#endif

		return false;
	}

	/* the data bundle path */
	bool evaluateSingleInput (TransformT* trans,
	                          shared_ptr<InputBundleT> input_bundle,
	                          shared_ptr<OutputBundleT> output_bundle)  override
	{
		return evaluateSingleInputReduce(trans, input_bundle, output_bundle);
	}
};
#endif //CREEK_WINKEYREDUCERMTEVAL_H
