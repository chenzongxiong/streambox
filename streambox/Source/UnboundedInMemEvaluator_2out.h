#ifndef TWITTEREVALUATE_H
#define TWITTEREVALUATE_H

//hym: this is for UnboundedInMemTweets
//     emit 2 same output streams
//     include this file in UnboundedInMemEvaluator.h
void evaluate_2outputs(TransformT* t, EvaluationBundleContext *c,
	shared_ptr<BundleBase> bundle_ptr = nullptr) //override
{

	PValue* out[] = { t->getFirstOutput(), t->getSecondOutput() };
	assert(out[0]);

	int num_outputs = 1;
	if (out[1]){
		num_outputs = 2;
	}


	//# of bundles between two puncs
	const uint64_t bundle_per_interval = 1 * numa_num_configured_cpus(); /*# cores*/
	boost::posix_time::time_duration punc_interval =
		milliseconds(t->punc_interval_ms);

	boost::posix_time::time_duration delta =
		milliseconds(t->punc_interval_ms) / bundle_per_interval;

	const uint64_t records_per_bundle =
		t->records_per_interval / bundle_per_interval;

	EE(" ---- punc internal is %d sec (ev time) --- ",
			t->punc_interval_ms / 1000);

	//bool is_active = false;

	/* an infi loop that emit bundles to all NUMA nodes periodically.
	 * spawn downstream eval to consume the bundles.
	 * NB: we can only sleep in the main thread (which is not
	 * managed by the NUMA thread pool).
	 */

	const int num_nodes = numa_max_node() + 1;

	uint64_t us_per_iteration =
		1e6 * t->records_per_interval * 2 / t->target_tput; /* the target us */

//version 1
//hym: create bundle pointer for 2 streams
//     two streams will not share the bundle pointer
#if 1
	uint64_t offset[2] = {0};//support 10 input streams now
	while(true){
		boost::posix_time::ptime start_tick =
			boost::posix_time::microsec_clock::local_time();

		for(unsigned int i = 0; i < bundle_per_interval; i++){
			/* construct the bundles by reading NUMA buffers round-robin */
			int nodeid = i % num_nodes;

			for(int oid = 0; oid < num_outputs; oid++){
				/* Assumble a bundle by drawing records from the corresponding
				* NUMA buffer. */
				shared_ptr<BundleT<T>>
					bundle(make_shared<BundleT<T>>(
						records_per_bundle,  /* reserved capacity */
						nodeid));
				xzl_assert(bundle);
				VV("pack records in bundle ts %s:",
					to_simple_string(current_ts + delta * i).c_str());

				for(unsigned int j = 0; j < records_per_bundle; j++, offset[oid]++){

					if(offset[oid] == t->buffer_size_records){
						offset[oid] = 0; //wrap around
					}
					t->record_buffers[nodeid][offset[oid]].ts = BaseT::current_ts + delta * i;
					bundle->add_record(t->record_buffers[nodeid][offset[oid]]);
				}//end for: records_per_bundle
				out[oid]->consumer->depositOneBundle(bundle, nodeid);
				c->SpawnConsumer();
			}//end for: num_outputs
		}//end for: bundle_per_interval
#endif


//version 2
//hym: create one bundle pointer for two streams
//     the bundle pointer will be shared by two streams
#if 0
	uint64_t offset = 0;
	while(true){
		boost::posix_time::ptime start_tick =
			boost::posix_time::microsec_clock::local_time();

		for(unsigned int i = 0; i < bundle_per_interval; i++){
			/* construct the bundles by reading NUMA buffers round-robin */
			int nodeid = i % num_nodes;

			/* Assumble a bundle by drawing records from the corrsponding
			 * NUMA buffer. */
			shared_ptr<BundleT<T>>
				bundle(make_shared<BundleT<T>>(
							records_per_bundle,  /* reserved capacity */
							nodeid));
			xzl_assert(bundle);
			VV("pack records in bundle ts %s:",
					to_simple_string(current_ts + delta * i).c_str());

			for(unsigned int j = 0; j < records_per_bundle; j++, offset++){

				if(offset == t->buffer_size_records){
					offset = 0; //wrap around
				}
				t->record_buffers[nodeid][offset].ts = BaseT::current_ts + delta * i;
				bundle->add_record(t->record_buffers[nodeid][offset]);
			}//end for: records_per_bundle

			//deposit same bundles to two streams
			for(int oid = 0; oid < num_outputs; oid++){
				out[oid]->consumer->depositOneBundle(bundle, nodeid);
				c->SpawnConsumer();
			}//end for: num_outputs
		}//end for: bundle_per_interval
#endif

		t->byte_counter_.fetch_add(t->records_per_interval * t->string_len * 2,
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

		BaseT::current_ts += punc_interval;
		BaseT::current_ts += milliseconds(t->session_gap_ms);

		/* Make sure all data have left the source before
		 * advancing the watermark. otherwise, downstream transforms may see
		 * watermark advance before they see the (old) records.
		 */

		 /*
		  * update source watermark immediately, but propagate the watermark
		  * to downstream asynchronously.
		  */
		c->UpdateSourceWatermark(BaseT::current_ts);

		/* Useful before the sink sees the 1st watermark */
		if (c->GetTargetWm() == max_date_time) { /* unassigned */
			c->SetTargetWm(BaseT::current_ts);
		}

		/* where we process wm.
		 * we spread wm among numa nodes
		 * NB this->_node may be -1
		 */
		static int wm_node = 0;
		for(int oid = 0; oid < num_outputs; oid++){
			out[oid]->consumer->depositOnePunc(
					make_shared<Punc>(BaseT::current_ts, wm_node), wm_node);
			//std::cout << "deposit one punc --------" << std::endl;
			c->SpawnConsumer();
			if (++wm_node == numa_max_node()){
				wm_node = 0;
			}
		}

	}//end while
}//end evaluate
#endif /* TWITTEREVALUATE_H */
