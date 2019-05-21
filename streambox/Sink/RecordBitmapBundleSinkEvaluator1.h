#ifndef RECORDBITMAPBUNDLESINKEVALUATOR1_H
#define RECORDBITMAPBUNDLESINKEVALUATOR1_H 

/* owned by hym
 * this is for sink of type SIDE_INFO_JD
 */

extern "C" {
#include "measure.h"
}

#include "core/EvaluationBundleContext.h"
#include "core/TransformEvaluator.h"
//#include "SingleInputTransformEvaluator.h"
#include "Join/Join.h"
#include "Sink/Sink.h"
#include "Values.h"

template <typename InputT>
class RecordBitmapBundleSinkEvaluator1: public TransformEvaulator<RecordBitmapBundleSink<InputT>> {
	using TransformT = RecordBitmapBundleSink<InputT>;
	using InputBundleT = RecordBitmapBundle<pair<long,vector<long>>>;
	using OutputBundleT = RecordBitmapBundle<pair<long,vector<long>>>;
	/*
	using TransformT = RecordBitmapBundleSink<InputT>;
	using InputBundleT = RecordBitmapBundle<InputT>;
	using OutputBundleT = RecordBitmapBundle<InputT>;
	*/
public:
	RecordBitmapBundleSinkEvaluator1(int node){
		ASSERT_VALID_NUMA_NODE(node);
		this->_node = node;
	}

	bool process_bundle(TransformT* trans, EvaluationBundleContext* c,
			shared_ptr<InputBundleT> input_bundle, PValue *out){

//		EE("do nothing but return");
//		return;

		bool ret = false;
		ret = ret; //fix warning
		shared_ptr<OutputBundleT> output_bundle = nullptr;
		vector<shared_ptr<OutputBundleT>> output_bundles;
		
		int output_size = (input_bundle->size() >= 0 ? input_bundle->size() : 128);
		output_bundle = make_shared<OutputBundleT>(output_size,
			input_bundle->node); /* same node as the input bundle */

		for (auto it = input_bundle->begin(); it != input_bundle->end(); ++it) {
			output_bundle->add_record(Record<pair<long,vector<long>>>
				( make_pair((*it).data.first,
				            (*it).data.second
					   ),
				  (*it).ts
				)
				
			);
		}
		
		output_bundles.push_back(output_bundle);

		vector<shared_ptr<BundleBase>> output_bundles_ (
			output_bundles.begin(),
			output_bundles.end());

		ret = true;
		trans->depositBundleDownstream_4(out->consumer, input_bundle, output_bundles_);
		auto oldref = input_bundle->refcnt->fetch_sub(1);
		oldref = oldref; // fix warning
		assert(oldref > 0);
		return true;
	}


	bool process_punc(TransformT* trans, EvaluationBundleContext* c,
			shared_ptr<Punc> punc, PValue *out){
		//hym: update trans' left_wm or right_wm
		if(punc->get_side_info() == 1){
			trans->left_wm = punc->min_ts;
			std::cout << "First sink receive a punc from left---------------" << std::endl;
		}else if(punc->get_side_info() == 2){
			trans->right_wm = punc->min_ts;
			std::cout << "First sink receive a pun from right ++++++++++++++" << std::endl;
		}else{
			assert(false && "Punc has wrong side info in Join's downstream");
		}

		//bundle_container *upcon = trans->localBundleToContainer(punc);
#ifdef MEASURE_LATENCY
		punc->mark("retrieved by: " + trans->name);
#endif
		ptime nw_wm = punc->min_ts;
		shared_ptr<Punc> new_punc = make_shared<Punc>(nw_wm, punc->node);
#ifdef MEASURE_LATENCY
		new_punc->inherit_markers(*punc);
#endif
		trans->depositPuncDownstream_4(out->consumer, punc, new_punc, new_punc->node);
		long expected = PUNCREF_RETRIEVED;
		if (!punc->refcnt->compare_exchange_strong(expected, PUNCREF_CONSUMED)) {
			bug("bad punc's refcnt?");
		}

		trans->SetWatermark(nw_wm);
		return true;
	}

	void evaluate(TransformT* trans, EvaluationBundleContext* c,
		shared_ptr<BundleBase> bundle_ptr = nullptr) override{
		
		assert(bundle_ptr && "must specify bundleptr. is old interface used?");
		assert(typeid(*trans) == typeid(TransformT));
	
//		EE("--- called ---- ");

		auto out = trans->getFirstOutput();
		assert(out);

		/* at this time, the numa node of this eval must be determined. */
		ASSERT_VALID_NUMA_NODE(this->_node);

		/* the bundle_ptr could be a data bundle or punc */
		auto input_bundle = dynamic_pointer_cast<InputBundleT>(bundle_ptr);
		
		if (!input_bundle) { /* punc path */
			//std::cout << __FILE__ << ": " <<  __LINE__ << " trans is " << trans->getName() <<std::endl;	
			auto punc = dynamic_pointer_cast<Punc>(bundle_ptr);
			assert(punc); /* neither bundle or punc, we don't know what to do */
			process_punc(trans, c, punc, out);
#ifdef DEBUG
			std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
#endif
		} else { //bundle path
			//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
			process_bundle(trans, c, input_bundle, out);
#ifdef DEBUG
			std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
#endif
		}

	}
};


#endif /* RECORDBITMAPBUNDLESINKEVALUATOR1_H */
