#ifndef UNBOUNDEDINMEMEVALUATOR_JOIN_H
#define UNBOUNDEDINMEMEVALUATOR_JOIN_H
//hym: new source for join
//     get ride of SimpleMapper
//     Source deposit bundles to Join directly
#include "UnboundedInMemEvaluator.h"
#include "Source/Unbounded_Join.h"
#include "UnboundedInMemEvaluatorBase.h"
#if 0
template<class T, template<class> class BundleT> class UnboundedInMemEvaluator_Join;  /* no generic impl */

template<template<class> class OutputBundleT>
class UnboundedInMemEvaluator_Join<long, OutputBundleT>
	: public UnboundedInMemEvaluatorBase<long, OutputBundleT> {

	using T = long;
	using TransformT = typename UnboundedInMemEvaluatorBase<T, OutputBundleT>::TransformT;
	//using OutputBundleT = RecordBitmapBundle<T>;
  	using BaseT = UnboundedInMemEvaluatorBase<T, OutputBundleT>;
	using KVPair = pair<long, long>;
	using OutputBundleT  = RecordBitmapBundle<pair<long, long>>;
	//hym: maybe should use the parameter.....
#endif
/*
template<template<class> class BundleT>
class UnboundedInMemEvaluator_Join<long, BundleT>
	: public UnboundedInMemEvaluatorBase<long, BundleT> {

	using T = long;
	using TransformT = typename UnboundedInMemEvaluatorBase<T, BundleT>::TransformT;
	//using OutputBundleT = RecordBitmapBundle<T>;
  	using BaseT = UnboundedInMemEvaluatorBase<T, BundleT>;
	using KVPair = pair<long, long>;
	using OutputBundleT  = RecordBitmapBundle<pair<long, long>>;
	//hym: maybe should use the parameter.....
*/
class UnboundedInMemEvaluator_Join
	: public TransformEvaulator<UnboundedInMem_Join>{
	using T = long;
	using OutputBundleT  = RecordBitmapBundle<pair<long, long>>;
	using TransformT = UnboundedInMem_Join;

public:
	ptime current_ts;

	/* we are not bound to any node */
	UnboundedInMemEvaluator_Join(int node) : current_ts(boost::gregorian::date(2016, Jan, 1)) { }

	void evaluate(TransformT* t, EvaluationBundleContext *c,
			shared_ptr<BundleBase> bundle_ptr = nullptr) override {

		auto out = t->getFirstOutput();
		assert(out);

		//# of bundles between two puncs
		const uint64_t bundle_per_interval = 1 * numa_num_configured_cpus();

		boost::posix_time::time_duration punc_interval =
			milliseconds(t->punc_interval_ms);

		boost::posix_time::time_duration delta =
			milliseconds(t->punc_interval_ms) / bundle_per_interval;

		const uint64_t records_per_bundle =
			t->records_per_interval / bundle_per_interval;

		EE(" ---- punc internal is %d sec (ev time) --- ",
			t->punc_interval_ms / 1000);

		//XXX WARNING: comment this temporarilly. Remember to restore this later!!!
		//const int num_nodes = numa_max_node() + 1;
		const int num_nodes = 1;

		uint64_t us_per_iteration =
			1e6 * t->records_per_interval * 2 / t->target_tput; /* the target us */
		//uint64_t offset1 = 0;
		//uint64_t offset2 = 0;
		uint64_t offset[10] = {0};//support 10 input streams at most
		int num_outputs = 2; //deposit 2 streams to Join's left and right

		while(true){
			boost::posix_time::ptime start_tick =
				boost::posix_time::microsec_clock::local_time();
			for(unsigned int i = 0; i < bundle_per_interval; i++){
				int nodeid = i % num_nodes;

				for(int oid = 0; oid < num_outputs; oid++){
					//OutputBundleT is RecordBitmapBundle<pair<long, long>>
					shared_ptr<OutputBundleT>
						bundle(make_shared<OutputBundleT>(
							records_per_bundle,
							nodeid
						));
					for(unsigned int j = 0; j < records_per_bundle; j++, offset[oid]++){
						if(offset[oid] == t->record_num){
							offset[oid] = 0; //wrap around
						}
						t->record_buffers[nodeid][offset[oid]].ts = current_ts + delta * i;
						bundle->add_record(
							//pair<long, long>
							Record<pair<long, long>>(
								pair<long, long>(
									t->record_buffers[nodeid][offset[oid]].data,
									oid //set to 0 or 1
								),
								current_ts + delta * i
							)//Record(data, ts)
						);//add_record
					}

					//XXX need pass side_info to this func
					//out[oid]->consumer->depositOneBundle(bundle, 0); //XXX only on node 0
					//std::cout << "deposit one bundle ++++++" << std::endl;
					//XXX out->consumer is Join
					if(oid == 0){ //to Join's left
						out->consumer->depositOneBundleToJoin_L(bundle, 0); //XXX only on node 0
					}else if(oid == 1){ //to Join's right
						out->consumer->depositOneBundleToJoin_R(bundle, 0); //XXX only on node 0
					}else{
						assert(false && "wrong oid???");
					}

					c->SpawnConsumer();
				}
			}//end for: bundle_per_interval

			t->byte_counter_.fetch_add(t->records_per_interval * t->record_len * 2,
					std::memory_order_relaxed);

			t->record_counter_.fetch_add(t->records_per_interval * 2,
					std::memory_order_relaxed);

			boost::posix_time::ptime end_tick =
				boost::posix_time::microsec_clock::local_time();
			auto elapsed_us = (end_tick - start_tick).total_microseconds();
			assert(elapsed_us > 0);

			if ((unsigned long)elapsed_us < us_per_iteration) {
				usleep(us_per_iteration - elapsed_us);
				//I("source pauses for %lu us", us_per_iteration - elapsed_us);
			}

			current_ts += punc_interval;
			current_ts += milliseconds(t->session_gap_ms);

			c->UpdateSourceWatermark(current_ts);

			/* Useful before the sink sees the 1st watermark */
			if (c->GetTargetWm() == max_date_time) { /* unassigned */
				c->SetTargetWm(current_ts);
			}

			static int wm_node = 0;
			for(int oid = 0; oid < num_outputs; oid++){
				/*
				out[oid]->consumer->depositOnePunc(
					make_shared<Punc>(BaseT::current_ts, wm_node), wm_node);
					//std::cout << "deposit one punc --------" << std::endl;
				*/
				if(oid  == 0){
					out->consumer->depositOnePuncToJoin_L(
						make_shared<Punc>(current_ts, wm_node), wm_node);
				}else if(oid == 1){
					out->consumer->depositOnePuncToJoin_R(
						make_shared<Punc>(current_ts, wm_node), wm_node);
				}else{
					assert(false && "wrong oid ??");
				}

				c->SpawnConsumer();
				if (++wm_node == numa_max_node()){
					wm_node = 0;
				}
			}

		}//end while
	}//end evaluate
};

#endif /* UNBOUNDEDINMEMEVALUATOR_JOIN_H */
