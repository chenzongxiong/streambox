#ifndef SINGLE_INPUT_T_EVAL_H
#define SINGLE_INPUT_T_EVAL_H

#include "EvaluationBundleContext.h"
#include "TransformEvaluator.h"

template<class TransformT, class InputBundleT, class OutputBundleT>
class SingleInputTransformEvaluator: public TransformEvaulator<TransformT> {

public:
	SingleInputTransformEvaluator(int node) {
		ASSERT_VALID_NUMA_NODE(node);
		this->_node = node;
	}

	// process one incoming bundle.
	// to be implemented by concrete evaluator.
	// @output_bundle: an allocated output bundle to deposit records into
	// @return: true if a consumer should be scheduled to consume
	// the output bundle

	/* Concrete transform should override at least one of the following */

	/* legacy */
	virtual bool evaluateSingleInput(TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle) = 0;

	virtual bool evaluateSingleInput(TransformT* trans,
			shared_ptr<InputBundleT> input_bundle,
			shared_ptr<OutputBundleT> output_bundle,
			EvaluationBundleContext* c /* need this to access tp */) {

		return evaluateSingleInput(trans, input_bundle, output_bundle);

	}


	/*
	 * Given an "external" wm, attempt to flush the internal state.
	 *
	 * caller holds statelock.
	 *
	 * @up_wm: the wm considering received punc and inflight bundles.
	 * the function should attempt to flush based on @up_wm.
	 *
	 * @conlocked: whether the caller holds the conlock
	 * Often this func does not need this to update local wm
	 *
	 * @output_bundles: a local container to deposite output bundles
	 *
	 * return: the updated trans wm after flush.
	 * NB: this func should not update trans wm directly.
	 *
	 * for stateless transform, no flush is
	 * done and @up_wm is returned. for stateful transform, if the internal
	 * state cannot be flushed "up to" @up_wm, an earlier wm should be returned.
	 * This happens, e.g. when a window cannot be closed as of now.
	 *
	 * But the return @wm should not be later than @up_wm per definition.
	 *
	 */

	/* may grab the transform's internal lock */
	virtual ptime flushState(TransformT* trans, const ptime up_wm,
				vector<shared_ptr<OutputBundleT>>* output_bundles,
				EvaluationBundleContext* c, bool purge = true) {
        // std::cout << "flush state in single input transform evaluator" << std::endl;
			return up_wm;
	}

private:
#ifdef DEBUG
	ptime last_punc_ts;
	bool once = true;
#endif

	/* Process a punc. Flush the internal state and assign a punc to downstream
	 * as needed. NB the punc may not be immediately retrieved by downstream.
	 */
	bool process_punc(TransformT* trans, EvaluationBundleContext* c,
			shared_ptr<Punc> punc, PValue *out)
	{
		bool out_punc = false;

#ifdef DEBUG /* extra check: punc must be monotonic */
		if (once) {
			once = false;
		} else {
			assert(last_punc_ts < punc->min_ts && "bug: punc regression");
		}
		last_punc_ts = punc->min_ts;
#endif

		W("%s (node %d) get punc: %s", trans->name.c_str(), this->_node,
				to_simple_string(punc->min_ts).c_str());
		//assert(punc->min_ts > trans->GetWatermarkSafe());

#if 0
		punc->mark("retrieved by: " + trans->name);
#endif

		/* We retrieved a punc. This only happens after all bundles of the enclosing
		 * containers are consumed.
		 */

/* #if 0 */
/* 		/\* Deposit staged bundles to downstream before processing this punc. */
/* 		 * @staged_bundles is stored in the container, which won't go away before */
/* 		 * we decrement punc_refcnt. */
/* 		 *\/ */
/* 		assert(punc->staged_bundles); */
/* 		if (out->consumer) { */
/* 			for (auto & b : *(punc->staged_bundles)) { */
/* 				assert(b); */
/* 				out->consumer->depositOneBundle(b, b->node); */
/* 			} */

/* 			if (!punc->staged_bundles->empty()) { */
/* 				I("deposit %lu staged bundles to %s", punc->staged_bundles->size(), */
/* 						out->consumer->name.c_str()); */
/* 				c->SpawnConsumer(); /\* once *\/ */
/* 			} */
/* 		} */
/* #endif */

		/* Ready to flush internal state. */
		vector<shared_ptr<OutputBundleT>> flushed_bundles;
		ptime new_wm = flushState(trans, punc->min_ts, &flushed_bundles, c, true);

		/* XXX ugly XXX */
		vector<shared_ptr<BundleBase>> flushed_bundles_(flushed_bundles.begin(),
				flushed_bundles.end());

		if (out->consumer) {
/* #if 0 */
/* 			for (auto & b : flushed_bundles) { */
/* 				assert(b); */
/* 				out->consumer->depositOneBundle(b, b->node); */
/* #ifdef DEBUG_CONCURRENCY */
/* 				std::cout << trans->getName() << " depositOneBundle to: " << out->consumer->getName() << std::endl; */
/* #endif */
/* 			} */
/* #endif */
			trans->depositBundlesDownstream(out->consumer, punc, flushed_bundles_);
		}

		/* --- assign punc after data bundles --- */
//		shared_ptr<Punc> new_punc = nullptr;

		if (trans->SetWatermarkSafe(new_wm)) { /* local wm changed */
			if (out->consumer) {
				shared_ptr<Punc> new_punc = make_shared<Punc>(new_wm, punc->node);
/* #if 0 */
/* 				new_punc->inherit_markers(*punc); */
/* #endif */
//			out->consumer->depositOnePunc(new_punc, new_punc->node);
				trans->depositOnePuncDownstream(out->consumer, punc, new_punc, new_punc->node);
#ifdef DEBUG_CONCURRENCY
				std::cout << trans->getName() << " depositOnePunc to: " << out->consumer->getName() << std::endl;
#endif
				out_punc = true;
			} else {  /* no downstreawm (sink?) print out something */
				W("		%s passthrough a punc: %s", trans->name.c_str(),
						to_simple_string(trans->GetWatermarkSafe()).c_str());
	#if 0
				punc->	dump_markers();
	#endif
				c->OnSinkGetWm(punc->min_ts);
			}
		} else { /* local wm unchanged. */
			/*
			 * if the incoming wm is dropped (does not trigger wm output), we do
			 * not trace it for latency.
			 */
			if (out->consumer) {
				trans->cancelPuncDownstream(out->consumer, punc);
			} else {
				c->OnSinkGetWm(punc->min_ts);
				//bug("will this happen?");
			}
		}

		/* now the received punc is consumed and can be destroyed. This allows
		 * the enclosing container to be closed and the next punc to be retrieved
		 * by this tran.
		 */
//		auto t = punc->refcnt->fetch_sub(1);
//		t = t; assert(t == 1);

		long expected = PUNCREF_RETRIEVED;
		if (!punc->refcnt->compare_exchange_strong(expected, PUNCREF_CONSUMED)) {
			bug("bad punc's refcnt?");
		}

		/* --- spawning ---- */

		for (auto & b : flushed_bundles) {
			b = b; assert(b);
//			c->SpawnConsumer(out, b->node);
			c->SpawnConsumer();
		}

		if (out_punc)
//			c->SpawnConsumer(out, new_punc->node);
			c->SpawnConsumer();

		return out_punc;
	}

	bool process_bundle(TransformT* trans, EvaluationBundleContext* c,
			shared_ptr<InputBundleT> input_bundle, PValue *out)
	{

		bool ret = false;

		assert(out);

		vector<shared_ptr<OutputBundleT>> output_bundles;
		/* contains the output bundle */
		shared_ptr<OutputBundleT> output_bundle = nullptr;

		int output_size = (input_bundle->size() >= 0 ? input_bundle->size() : 128);

		output_bundle = make_shared<OutputBundleT>(output_size,
				input_bundle->node); /* same node as the input bundle */

		if (evaluateSingleInput(trans, input_bundle, output_bundle, c)) {
			output_bundles.push_back(output_bundle);
		}

		/* No extra action needed here. If the container owning the input bundle
		 * becomes empty because of this, next getOneBundle() will return its punc.
		 * Then this trans can flush the state & emit bundles them.
		 */

#if 0
		bool is_staged = true;

		if (!output_bundles.empty()) {
			ret = true;
			/* where to write the bundle? downstream or the staging area? */

			/* have to copy-construct a vector of basebundles shared ptr as the
			 * transform's tryStageBundles() takes BundleBase.
			 */
			vector<shared_ptr<BundleBase>> out_basebundles (output_bundles.begin(),
					output_bundles.end());

//			if (!trans->tryStageBundle(input_bundle, output_bundle)) {
			if (!trans->tryStageBundles(input_bundle, out_basebundles)) {
				is_staged = false;
				for (auto & b : output_bundles) {
					assert(b);
#ifdef DEBUG_CONCURRENCY
					std::cout << trans->getName() << " depositOneBundle to: " << out->consumer->getName() << std::endl;
#endif
					out->consumer->depositOneBundle(b, b->node);
				}
			}
		}
#endif

		if (!output_bundles.empty()) {
			/* XXX ugly XXX */
			vector<shared_ptr<BundleBase>> output_bundles_ (output_bundles.begin(),
					output_bundles.end());
			ret = true;
			trans->depositBundlesDownstream(out->consumer, input_bundle, output_bundles_);
		}

		auto oldref = input_bundle->refcnt->fetch_sub(1);  /* input bundle consumed */
		if(oldref <= 0){
			EE("bug: %s: oldref is %ld <= 0. container %lx",
					trans->name.c_str(), oldref, (unsigned long)(input_bundle->container));
			trans->dump_containers("bug");
			abort();
		}
		assert(oldref > 0);

		for (auto & b : output_bundles) {
			b = b; assert(b); /* must not be null */
			c->SpawnConsumer();
		}

#if 0
		if (!is_staged) { /* if staged, output bundles are not visible to downstream yet */
			for (auto & b : output_bundles) {
				b = b; assert(b); /* must not be null */
	//			c->SpawnConsumer(out, b->node);
				c->SpawnConsumer();
			}
		}
#endif

		return ret;
	}

public:
#if 0
	void evaluateOneBundle(TransformT* trans,
			EvaluationBundleContext* c, shared_ptr<BundleBase> bundle_ptr) override
	{

		assert(typeid(*trans) == typeid(TransformT));
		assert(bundle_ptr);

		PValue* in1 = trans->getFirstInput();
		//assert(in1);
		if(!in1){
			std::cout << "ERROR: in1 is NULL" << std::endl;
			abort();
		}
		auto out = trans->getFirstOutput();
		//assert(out);
		if(!out){
			std::cout << "ERROR: out is NULL" << std::endl;
			abort();
		}

		/* NB: Even sinks may have valid output, but that out has null consumer.
		 * See apply1() */

		/* at this time, the numa node of this eval must be determined. */
		ASSERT_VALID_NUMA_NODE(this->_node);

		vector<shared_ptr<OutputBundleT>> flushed_bundles;

#ifndef INPUT_ALWAYS_ON_NODE0
		/* this will fail if we don't have separate input queues and
		 * let evaluators to freely grab bundles from the queue...
		 */
		assert(this->_node == bundle_ptr->node);
#endif

		/* the bundle_ptr could be a data bundle or punc */
		auto input_bundle = dynamic_pointer_cast<InputBundleT>(bundle_ptr);

		if (!input_bundle) { /* punc path */
			auto punc = dynamic_pointer_cast<Punc>(bundle_ptr);
			assert(punc); /* neither bundle or punc, we don't know what to do */
			process_punc(trans, c, punc, out);
		} else { /* bundle path */
			process_bundle(trans, c, input_bundle, out);
		}
	}
#endif

	void evaluate(TransformT* trans, EvaluationBundleContext* c,
			shared_ptr<BundleBase> bundle_ptr = nullptr) override {

		assert(bundle_ptr && "must specify bundleptr. is old interface used?");

		assert(typeid(*trans) == typeid(TransformT));

		PValue* in1 = trans->getFirstInput();
		//assert(in1);
		if(!in1){
			std::cout << "ERROR: in1 is NULL" << std::endl;
			abort();
		}
		auto out = trans->getFirstOutput();
		//assert(out);
		if(!out){
			std::cout << "ERROR: out is NULL" << std::endl;
			abort();
		}
		/* NB: Even sinks may have valid output, but that out has null consumer.
		 * See apply1() */

		/* at this time, the numa node of this eval must be determined. */
		ASSERT_VALID_NUMA_NODE(this->_node);

		vector<shared_ptr<OutputBundleT>> flushed_bundles;

#if 0
		if (!bundle_ptr) {
			bundle_ptr = trans->getOneBundle(this->_node); /* from trans */
			if (!bundle_ptr) {
				/* this is possible when downstream has nothing to do? */
				W("warning: no bundle returned");
				return;
			}
		}
#endif

#ifndef INPUT_ALWAYS_ON_NODE0
		/* this will fail if we don't have separate input queues and
		 * let evaluators to freely grab bundles from the queue...
		 */
		assert(this->_node == bundle_ptr->node);
#endif

		/* the bundle_ptr could be a data bundle or punc */
		auto input_bundle = dynamic_pointer_cast<InputBundleT>(bundle_ptr);

		if (!input_bundle) { /* punc path */
			auto punc = dynamic_pointer_cast<Punc>(bundle_ptr);
			assert(punc); /* neither bundle or punc, we don't know what to do */
/*
#ifdef DEBUG
			std::cout << __FILE__ << __LINE__ << ": process a punc in " << trans->getName() << std::endl;
#endif
*/
			process_punc(trans, c, punc, out);
		} else { /* bundle path */
/*
#ifdef DEBUG
			std::cout << __FILE__ << __LINE__ << ": process a bundle in " << trans->getName() << std::endl;
#endif
*/
			process_bundle(trans, c, input_bundle, out);
		}

	}

};

#endif // SINGLE_INPUT_T_EVAL_H
