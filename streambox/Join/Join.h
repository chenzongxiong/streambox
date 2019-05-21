#ifndef JOIN_H
#define JOIN_H
#include "core/Transforms.h"
#include "Values.h"

extern "C" {
#include "measure.h"
}


/* left/right inputs are streams of record<kvpair> */
template <class KVPair,
	template<class> class InputBundleT_ = RecordBitmapBundle,
	template<class> class OutputBundleT_ = RecordBitmapBundle>
class Join: public PTransform {
	using K = decltype(KVPair::first);
	using V = decltype(KVPair::second);
	using RecordKV = Record<KVPair>;
	using time_duration = boost::posix_time::time_duration;
	public:
	//hym: internal watermark for left and right sides
	//ptime left_wm;
	//ptime right_wm;
	// These two wm have been moved to PTransform

	/* 
	// This has beem moved to PTransform
	//hym: apply containr to Join transform
	// Join has two input containers streams, and two ouput container streams
	// Upstream will deposit bundles to left_containers_in or right_containers_in
	// Join evaluator will get bundles from those input containers, and 
	// deposit bundles to left_containers_out or right_containers_out
	// when Join process a punc to seal a container, it will move some bundles from
	// left_container_out or right_containerP_out to containers_, which is ordered 
	// and owned by Join't downstream
	// XXX TODO downstream grab bundles from left_containers_in or right_containers_in 
	// directly

	list<bundle_container> left_containers_in; //left input stream
	list<bundle_container> left_containers_out; //left outpout stream
	list<bundle_container> right_containers_in; //right input stream
	list<bundle_container> right_containers_out; //right output stream
	 */

#if 1 //hym: without uing Radix Join
	// containing kv pairs whose timestamps are bounded in the given range
	// (between two watermark advances)
	struct KVPairContainer {
		ptime start;
		ptime max;  // max ts carried by records
		//tbb::concurrent_unordered_map<K, RecordKV> vals;
		creek::concurrent_unordered_map<CONFIG_JOIN_HT_PARTITIONS, K, RecordKV> vals;

		KVPairContainer(ptime start)
			: start(start), max(start) { }

		void insert(RecordKV const & rec) {
			/*
			   if (vals.count(rec.data.first) != 0) {
			   printf("WARNING -- duplicate keys. existing %s: %ld %ld. "
			   "new %s: %ld %ld ",
			//"container size %lu",
			to_simplest_string(rec.ts).c_str(),
			rec.data.first,
			rec.data.second,
			to_simplest_string(vals[rec.data.first].ts).c_str(),
			vals[rec.data.first].data.first,
			vals[rec.data.first].data.second);
			//vals.size());
			//        abort();
			//      assert(vals.count(rec.data.first) == 0); // key must be unique?
			}
			 */
			vals[rec.data.first] = rec;

			if (rec.ts > max) { max = rec.ts; }
		}
	};
#endif


#if 0 //hum: using Radix Join
	//hym: now we only split the vals into 2 sets, can change this later
	//struct KVPairContainer_radix {
	struct KVPairContainer{
		ptime   start;
		ptime   max; //max ts carried by records
		//unordered_map<K, RecordKV> val_0; //group 0
		//unordered_map<K, RecordKV> val_1; //group 1
		tbb::concurrent_unordered_map<K, RecordKV> val_0; //group 0
		tbb::concurrent_unordered_map<K, RecordKV> val_1; //group 1

		//KVPairContainer_radix(ptime start)
		KVPairContainer(ptime start)
			: start(start), max(start){}

		//hym: limitation, the key only can be a NUMBER.
		void insert(RecordKV const & rec) {
			//hym: shall we handle duplicated key????
			//hym: we don't take care of this now, but later?
			if(rec.data.first % 2 == 0){
				val_0[rec.data.first] = rec;
				//                        cout << "+++++ insert into val_0" << endl;
			}else{
				val_1[rec.data.first] = rec;
				//                        cout << "+++++ insert into val_1" << endl;
			}

			if(rec.ts > max){
				max = rec.ts;
			}
		}
	};
#endif

#if 1
	// use KVPaireContainer's pointer in map
	struct WindowsKVPairContainer {
		// containers are ordered by their start time (watermark updates).
		// the records are not necessarily ordered within the same container or
		// across multiple containers.
		//map<ptime, KVPairContainer> containers;
		map<ptime, shared_ptr<KVPairContainer>> containers;
		boost::shared_mutex mtx_ctn;
		// start a new container with the given watermark as "start"
		// the new container will receive all records until the next watermark
		// update.
		// @return: the watermark of the "smallest" container (which is the lower
		// bound of the ts of all records)
		ptime watermark_update(ptime const & ts) {
			//containers.emplace(ts, KVPairContainer(ts));

			boost::unique_lock<boost::shared_mutex> writer_lock(mtx_ctn);

			containers.emplace(ts, make_shared<KVPairContainer>(ts));
			writer_lock.unlock();
			return containers.begin()->first;
		}

		bool add_record(RecordKV const & rec) {
			assert(containers.size());  // at least has one container.

			//XXX unsafe here
			// but if we get writer lock here, there will be deadlock.
			// because it has already gotten reader lock before add_record()
			/* Don't need this anymore. Purge will open new a new container. 
			   if(containers.empty()){
			   containers.emplace(min_date_time, make_shared<KVPairContainer>(min_date_time));
			   }
			 */

			//boost::shared_lock<boost::shared_mutex> reader_lock(mtx_ctn);

			// insert to the current (latest) container.
			//containers.rbegin()->second.insert(rec);
			containers.rbegin()->second->insert(rec);
			//reader_lock.unlock();
			return true;
		}

		// drop containers whose records are earlier than or equal to @cutoff
		// (i.e. only keep those who have records strictly later than @cutoff)
		// @return: the min_ts among the remaining containers.
		// NB: this is latency sensitive since it takes writer lock.
		// So we only actually destroy hashmaps after releasing writer lock.
		ptime purge_containers(ptime cutoff) {
			int ret = 0;
			assert(containers.size());

			/* temporarily hold the HT to be destroyed, but don't do so until
			 * we release the writer lock
			 */
			std::vector<shared_ptr<KVPairContainer>> vct;

			boost::unique_lock<boost::shared_mutex> writer_lock(mtx_ctn);
			//k2_measure("begin-purge");

			/* Don't need this any more, since purge will open a new container if there is no container.
			   assert(containers.size());
			   if(containers.empty()){
			   ptime now = boost::posix_time::microsec_clock::local_time(); 
			   return now;
			   }
			 */
			/* Note that this cannot be done by just examining the start watermark,
			 * although containers are ordered by their start watermarks.
			 e.g. ([] means start watermark)
			 [2] 3 15
			 [5] 7 6 12
			 [8] 9 10 13
			 when cutoff=8, we cannot drop the first two containers despite their
			 small start watermarks. Instead, we need to examine all containers.
			 */

			//ptime last_container_max = containers.rbegin()->second.max;
			//ptime last_container_max = containers.rbegin()->second->max;
			ptime last_container_max = containers.rbegin()->second->max;

			// note that we erase containers while iterating - need some special treatment
			for (auto && it = containers.begin(); it != containers.end();
					/* no increment */ ) {
				//if (it->second.max <= cutoff) {
				if (it->second->max <= cutoff) {
					//k2_measure("begin-erase");
					vct.push_back(it->second);
					containers.erase(it++);
					//k2_measure("end-erase");
					//k2_measure_flush();
					//it ++; // debugging
					ret ++;
					//				  std::cout << __FILE__ << __LINE__ << ": containers.erase(it)" << std::endl;
				}
				else
					++ it;
			}

			/* Note we even drop the last container if @cutoff exceeds its
			 * max ts. However, we remember this max: when we open a new container,
			 * we use the max as its "watermark". This is because that we will use
			 * the new container's watermark to approximate the min_ts of the
			 * left/right hash tables.
			 */

			//hym: have to unlock here, otherwise there will be deadlock since watermark_update() will get a lock again
			//writer_lock.unlock(); 

			if (containers.begin() == containers.end()) {
				//hym: have to unlock here, otherwise there will be deadlock since watermark_update() will get a lock again
				//writer_lock.unlock(); 
				//writer_lock.unlock(); 
				containers.emplace(last_container_max, make_shared<KVPairContainer>(last_container_max));
				//watermark_update(last_container_max);
				// so that the internal min_ts won't hold back transform's watermark

				//k2_measure("end-purge-open");
				//k2_measure_flush();
				writer_lock.unlock();
				//return max_date_time;
				return last_container_max;
			} else{
				// return the watermark of the "smallest" container, as an approximation
				// of the smallest ts of all the remaining records
				//k2_measure("end-purge-return");
				//k2_measure_flush();
				writer_lock.unlock();
				return containers.begin()->first;
			}

			xzl_bug("never here");
			}

#if 1 //hym: NOT using Radix Join
			// look for key match for records in time range [start, end)
			// return false if no found.
			// we use @rec to copy out found record in fear of concurrent purge.
			// XXX mutliple key matches??
			bool search(ptime start, ptime end, K const & key,
					RecordKV * rec) {

				//boost::shared_lock<boost::shared_mutex> reader_lock(mtx_ctn);

				/* the containers beyond @last all have start watermarks later than
				 * @end: they cannot contain records in range [start, end)
				 */
				auto && last = containers.lower_bound(end);
				for (auto && it = containers.begin(); it != last; it++) {
					//if (start > it->second.max) {
					if (start > it->second->max) {
						/* @start is larger than the max ts of the container. skip
						 * this container. */
						continue;
					}
					//auto && it_record = it->second.vals.find(key);
					auto && it_record = it->second->vals.find(key);
					//if (it_record == it->second.vals.end())
					if (it_record == it->second->vals.end())
						continue;
					// if there's a key match record out of the specified
					// time range, we skip to the next container.
					auto && record_ts = it_record->second.ts;
					if (!(record_ts >= start && record_ts < end))
						continue;
					*rec = it_record->second;  // found, copy.
					// XXX continue to check multi key matches?
					return true;
				}
				return false;
				}
#endif

				WindowsKVPairContainer() {
					watermark_update(min_date_time);
				}
			};
#endif //end WindowKVPairContainer

			//old design
#if 0
			struct WindowsKVPairContainer {
				// containers are ordered by their start time (watermark updates).
				// the records are not necessarily ordered within the same container or
				// across multiple containers.
				map<ptime, KVPairContainer> containers;

				// start a new container with the given watermark as "start"
				// the new container will receive all records until the next watermark
				// update.
				// @return: the watermark of the "smallest" container (which is the lower
				// bound of the ts of all records)
				ptime watermark_update(ptime const & ts) {
					containers.emplace(ts, KVPairContainer(ts));
					return containers.begin()->first;
				}

				bool add_record(RecordKV const & rec) {
					assert(containers.size());  // at least has one container.
					// insert to the current (latest) container.
					containers.rbegin()->second.insert(rec);

					return true;
				}

				// drop containers whose records are earlier than or equal to @cutoff
				// (i.e. only keep those who have records strictly later than @cutoff)
				// @return: the min_ts among the remaining containers.
				ptime purge_containers(ptime cutoff) {
					int ret = 0;
					assert(containers.size());

					/* Note that this cannot be done by just examining the start watermark,
					 * although containers are ordered by their start watermarks.
					 e.g. ([] means start watermark)
					 [2] 3 15
					 [5] 7 6 12
					 [8] 9 10 13
					 when cutoff=8, we cannot drop the first two containers despite their
					 small start watermarks. Instead, we need to examine all containers.
					 */

					ptime last_container_max = containers.rbegin()->second.max;

					// note that we erase containers while iterating - need some special treatment
					for (auto && it = containers.begin(); it != containers.end();
							/* no increment */ ) {
						if (it->second.max <= cutoff) {
							containers.erase(it++);
							ret ++;
						}
						else
							++ it;
					}

					/* Note we even drop the last container if @cutoff exceeds its
					 * max ts. However, we remember this max: when we open a new container,
					 * we use the max as its "watermark". This is because that we will use
					 * the new container's watermark to approximate the min_ts of the
					 * left/right hash tables.
					 */
					if (containers.begin() == containers.end()) {
						watermark_update(last_container_max);
						// so that the internal min_ts won't hold back transform's watermark
						return max_date_time;
					} else
						// return the watermark of the "smallest" container, as an approximation
						// of the smallest ts of all the remaining records
						return containers.begin()->first;
				}

#if 1 //hym: NOT using Radix Join
				// look for key match for records in time range [start, end)
				// return false if no found.
				// we use @rec to copy out found record in fear of concurrent purge.
				// XXX mutliple key matches??
				bool search(ptime start, ptime end, K const & key,
						RecordKV * rec) {
					/* the containers beyond @last all have start watermarks later than
					 * @end: they cannot contain records in range [start, end)
					 */
					auto && last = containers.lower_bound(end);
					for (auto && it = containers.begin(); it != last; it++) {
						if (start > it->second.max) {
							/* @start is larger than the max ts of the container. skip
							 * this container. */
							continue;
						}
						auto && it_record = it->second.vals.find(key);
						if (it_record == it->second.vals.end())
							continue;
						// if there's a key match record out of the specified
						// time range, we skip to the next container.
						auto && record_ts = it_record->second.ts;
						if (!(record_ts >= start && record_ts < end))
							continue;
						*rec = it_record->second;  // found, copy.
						// XXX continue to check multi key matches?
						return true;
					}
					return false;
				}
#endif

#if 0 //hym: Using Radix Join
				//hym: search using radix, KVPairContainer_radix
				//bool search_radix(ptime start, ptime end, K const & key, RecordKV * rec){
				bool search(ptime start, ptime end, K const & key, RecordKV * rec){
					auto && last = containers.lower_bound(end);
					for(auto && it = containers.begin(); it != last; it++){
						if(start > it->second.max){
							//start is larger than the max ts of the container, skip this container
							continue;
						}

						if(key % 2 == 0){
							//			cout << "search on val_0" << endl;
							auto && it_record = it->second.val_0.find(key);
							if(it_record == it->second.val_0.end()){
								//no match in this container
								continue;
							}
							auto && record_ts = it_record->second.ts;
							if(!(record_ts >= start && record_ts < end)){
								// ts out of time range
								continue;
							}
							//find record, and copy
							*rec = it_record->second;
							return true;
						}else{ //key %2 == 1
							//			cout << "search on val_1" << endl;
							auto && it_record = it->second.val_1.find(key);
							if(it_record == it->second.val_1.end()){
								//not match in this container
								continue;
							}
							auto && record_ts = it_record->second.ts;
							if(!(record_ts >= start && record_ts < end)){
								//out of time range
								continue;
							}
							*rec = it_record->second;
							return true;
						}
					}
					return false;
				}

#endif

				WindowsKVPairContainer() {
					watermark_update(min_date_time);
				}
			};
#endif //end WindowKVPairContainer

			//////////////////////////////////////////////////////////////

			// state of the entire transform. left, right.
				public:
			WindowsKVPairContainer win_containers[2];
			time_duration window_size;

				private:
			//public:
			//boost::shared_mutex _mutex[2];

				public:
			// caller must hold writer lock
			bool add_record(RecordKV const & rec, int i) {
				bool ret;
				//    boost::shared_lock<boost::shared_mutex> reader_lock(_mutex[i], boost::defer_lock);
				//boost::shared_lock<boost::shared_mutex> reader_lock(_mutex[i]);

				if (rec.ts < min_ts)  // XXX need a transform-wide lock??
					min_ts = rec.ts;

				//return win_containers[i].add_record(rec);
				//    reader_lock.lock();
				ret = win_containers[i].add_record(rec);
				//    reader_lock.unlock();

				return ret;
			}

			// caller must hold reader lock
			bool search(RecordKV const & in, RecordKV * out, int i) {
				//return win_containers[i].search(in.ts - window_size,
				//    in.ts + window_size, in.data.first, out);
				bool ret;
				//    boost::shared_lock<boost::shared_mutex> reader_lock(_mutex[i], boost::defer_lock);
				//boost::shared_lock<boost::shared_mutex> reader_lock(_mutex[i]);
				//    reader_lock.lock();
				ret = win_containers[i].search(in.ts - window_size, in.ts + window_size, in.data.first, out);
				//    reader_lock.unlock();
				return ret;
			}

			//hym: split this functio into watermark_update() and flush_state
#if 0
			// on receiving a new watermark on either of the two input streams
			// xzl: XXX bug? should we purge_containers(ts - window_size)?
			void watermark_update(ptime const & ts, int i) {
				ptime this_min_ts;

				// update containers of this stream
				//boost::unique_lock<boost::shared_mutex> this_writer_lock(_mutex[i]);
				this_min_ts = win_containers[i].watermark_update(ts);
				//this_writer_lock.unlock();

#if 1
				// and purge the state of the other stream
				//boost::unique_lock<boost::shared_mutex> other_writer_lock(_mutex[1-i]);
				ptime other_min_ts = win_containers[1-i].purge_containers(ts);
				(void)sizeof(other_min_ts); // prevent warning
				//other_writer_lock.unlock();
#endif

				//unique_lock<mutex> watermark_lock(mtx_watermark);
				//min_ts = min(this_min_ts, other_min_ts);
				//watermark_lock.unlock();
				// hym: we are not using min_ts and other_min_ts in new design
				/*
				   EE("new min_ts %s ts %s (%d: %s %d: %s)",
				   to_simplest_string(min_ts).c_str(),
				   to_simplest_string(ts).c_str(),
				   i, to_simplest_string(this_min_ts).c_str(),
				   1-i, to_simplest_string(other_min_ts).c_str()
				   );
				 */
			}
#endif

			//hym: these 2 funcs are from watermark_update() above
			void watermark_update(ptime const & ts, int i) {
				win_containers[i].watermark_update(ts);
			}

			/* this will internally lock join state */
			void flush_state(ptime const & ts, int i) {
				win_containers[1-i].purge_containers(ts);
			}

#if 0
			// called by the global evaluation context
			virtual ptime RefreshWatermark() override  {
				PValue* inputs[2] = {getFirstInput(), getSecondInput()};
				assert(inputs[0] && inputs[1]);

				PTransform* upstreams[2] = {inputs[0]->producer, inputs[1]->producer};
				assert(upstreams[0] && upstreams[1]);

				// consider inflights bundles for {left, right}
				ptime min_ts_flight[2];

				for (int i = 0; i < 2; i++) {
					min_ts_flight[i] = max_date_time;
					for (auto && b : dual_inflight_bundles[i]) {
						if (b->min_ts < min_ts_flight[i])
							min_ts_flight[i] = b->min_ts;
					}
				}

				// compute the watermarks of two streams (according to upstream, inputs,
				// and inflight). each stream will never see records earlier than the
				// updated watermark
				ptime min_ts_stream[2];
				for (int i = 0; i < 2; i++) {
					min_ts_stream[i] = min(upstreams[i]->watermark,
							min(min_ts_flight[i], inputs[i]->min_ts));
					watermark_update(min_ts_stream[i], i);
				}

				ptime pt = min(min_ts, min(min_ts_stream[0], min_ts_stream[1]));

				if (watermark > pt) {
					E("%s: watermark regress: current %s (new %s "
							"streams %s %s \t\t internal %s)",
							name.c_str(),
							to_simplest_string(watermark).c_str(),
							to_simplest_string(pt).c_str(),
							to_simplest_string(min_ts_stream[0]).c_str(),
							to_simplest_string(min_ts_stream[1]).c_str(),
							to_simplest_string(min_ts).c_str());
				}

				assert(watermark <= pt);

				if (watermark < pt) {
					watermark = pt;
					E("%s: update new watermark: %s\n"
							"\t\t (streams %s %s \t\t internal %s)",
							name.c_str(),
							to_simplest_string(watermark).c_str(),
							to_simplest_string(min_ts_stream[0]).c_str(),
							to_simplest_string(min_ts_stream[1]).c_str(),
							to_simplest_string(min_ts).c_str());
				} else {
					E("%s: watermark unchanged: current %s (new %s "
							"streams %s %s \t\t internal %s)",
							name.c_str(),
							to_simplest_string(watermark).c_str(),
							to_simplest_string(pt).c_str(),
							to_simplest_string(min_ts_stream[0]).c_str(),
							to_simplest_string(min_ts_stream[1]).c_str(),
							to_simplest_string(min_ts).c_str());
				}

				return pt;
			}
#endif

			Join(string name, time_duration window_size)
				: PTransform(name), window_size(window_size) {
					//XXX hym: should we init to current time???? 
					left_wm = boost::posix_time::microsec_clock::local_time();
					right_wm = boost::posix_time::microsec_clock::local_time();
				}

			void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
					shared_ptr<BundleBase> bundle_ptr) override;

				protected:
			// the age of the oldest internal record. *not* the watermark
			ptime min_ts = max_date_time;



				public:

			//hym: override Transform class' getOneBundle
			//get one bundle from left
			shared_ptr<BundleBase> getOneBundleOlderThan_left(ptime const & t,
					bool* has_pending_punc, int node = -1) {

				*has_pending_punc = false;

				/* fast path: if local wm more recent (>t), no work.
				 * (and any unemitted punc must >t) */

				if (this->GetWatermark() > t) {
					return nullptr;
				}


				/* slowpath: there may be some work in this trans. lock containers */

				shared_ptr<BundleBase> ret;

				//unique_lock<mutex> conlock(mtx_container);
				//unique_lock<mutex> conlock(left_mtx_container_in);
				//XXX rw
				boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
				boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
				unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
					= nullptr; /* will std::move later */
start:
				if (left_containers_in.empty()) {
					//			dump_containers("no containers");
					return nullptr;
				}

				/* Start from oldest container (at back()) and look for one available
				 * bundle.
				 * For each container, we first check available data bundles and then
				 * any ready punc. If no luck, move to a newer container.
				 *
				 * When encountering an empty container, check its refcnt.
				 */

				auto it = left_containers_in.rbegin();
				/* to remember the 1st container that 1) has invalidated punc; 2) has
				 * some bundles
				 */
				auto itt = left_containers_in.rend();

				while (it != left_containers_in.rend()) {

					if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
						assert(it->verifyPuncSafe() && "bug?");
						return nullptr;
					}

					/* punc has been canceled, if we see a later container assigned with
					 * a punc, we may return a bundle from *this* container.
					 */
					if (it->punc_refcnt == PUNCREF_CANCELED) {
						if (it->refcnt == 0) { /* dead container, clean it up. */
							/* because punc only canceled after bundles are deposited.*/

							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							assert(it->bundles.empty());
							it ++;
							it = list<bundle_container>::reverse_iterator(left_containers_in.erase(it.base()));
#ifdef DEBUG		
							std::cout << __FILE__ << __LINE__ << "left_containers_in.erase" << std::endl;
#endif
						}
						else if (!it->bundles.empty() && itt == left_containers_in.rend()) {
							/* bundles available + oldest container w/ punc canceled */
							itt = it;
							it ++;
						} else {
							/* else: could be no bundles in containers and some bundles are
							   outstanding... */
							it ++;
						}
						continue;
					}

					/* punc has *ever* been assigned */

					//if (it->punc->min_ts > t) { /* punc newer than wm. */
					if (it->getPuncSafe()->min_ts > t) { /* pfe()unc newer than wm. */
						return nullptr;
					}

					/* now we see a pending punc <=wm. */

					*has_pending_punc = true;

					/* grab from the oldest container w/ canceled punc */
					if (itt != left_containers_in.rend()) {
						//ret = itt->getBundleUnsafe();
						ret = itt->getBundleSafe();
						assert(ret);
						this->DecBundleCounter();
						//std::cout << __FILE__ << ": "<<__LINE__ << " getOneBundleOlderThan_left one" << std::endl; 
						return ret;
					}

					/* non-empty container. grab & go */
					//if ((ret = it->getBundleUnsafe())) {
					if ((ret = it->getBundleSafe())) {
						this->DecBundleCounter();
						//std::cout << __FILE__ << ": "<<__LINE__ << " getOneBundleOlderThan_left one" << std::endl;	
						return ret;
					}

					/* an empty container... */

					/* some data bundles outstanding.
					 * don't wait. move to a newer container. */
					if (it->refcnt != 0) {
						//				if (!(it->punc_refcnt != PUNCREF_RETRIEVED
						//										&& it->punc_refcnt != PUNCREF_CONSUMED)) {
						//					auto refcnt = it->punc_refcnt.load();
						//					E("bug: punc refcnt is %ld", refcnt);
						//				}
						assert(it->punc_refcnt != PUNCREF_RETRIEVED
								&& it->punc_refcnt != PUNCREF_CONSUMED);
						it ++;
						continue;
					}

					/* an empty container, punc has been assigned (CANCELED case is handled
					 * above). all bundles already consumed.
					 *
					 * we take a peek of punc_refcnt and decide.
					 * */

					auto punc_refcnt = it->punc_refcnt.load();

					/* first, clean any dead container so that it won't block
					 * later punc.
					 */
					if (punc_refcnt == PUNCREF_CONSUMED) { /* bundles+punc consumed */
						/* has to be the oldest container since we can't have more
						 * than one outstanding punc (the second oldest container, whose
						 * punc is assigned, must still holds its punc).
						 */

						//XXX rw
						if (!upgrade_locks(&rlock, &ulock, &pwlock)){
							goto start;
						}

						assert(it == left_containers_in.rbegin());
						//				it = containers_.erase(it);   /* erase w/ reverse iterator */
						it ++;
						it = list<bundle_container>::reverse_iterator(left_containers_in.erase(it.base()));
#ifdef DEBUG	
						//				std::cout << __FILE__ << __LINE__ << "left_containers_in.erase" << std::endl;
#endif			
						continue;
					} else if (punc_refcnt == PUNCREF_RETRIEVED) {
						/* punc outstanding. it may become 0 soon but we don't wait. */
						it ++;
						continue;
					} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
						if (it == left_containers_in.rbegin()) { /* oldest container. okay to emit */
							long expected = PUNCREF_ASSIGNED;
							if (!it->punc_refcnt.compare_exchange_strong(expected,
										PUNCREF_RETRIEVED)) {
								//bug("bad punc refcnt");
								it ++;
								continue;
							}
							//					auto r = it->punc_refcnt.fetch_sub(1);
							//					assert(r == 2); r = r;

							/* XXX: opt: we may check the next newer container and
							 * coalesce punc as needed. benefit unclear though.
							 */
							//return it->punc;
							return it->getPuncSafe();
						} else {
							/* not oldest container. can't emit punc before older ones. */
#if 0 /* not a good idea */
							/* opt: combine unemitted puncs.
							 *
							 * this is dangerous, since punc_refcnt may change under us.
							 * idea: temporarily decrements punc_refcnt to 1 so that
							 * other evaluators believe the punc is outstanding. after
							 * combing, increment to 2.  */
							auto & older = *(std::prev(it));
							if (older->punc_refcnt == 2) {
								older->wm = current->wm;
								//						it = containers_.erase(it); /* erase & continue */
								it ++;
								it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
							} else
#endif
								it ++;
							continue;
						}
					} else {
						EE("punc refcnt: %ld", punc_refcnt);
						assert(false && "illegal punc_perfcnt?");
					}

					assert(false && "bug?");
				} // while

				//		dump_containers("no bundles");

				return nullptr;
				}

				// get one bundle older than wm from rignt input stream
				shared_ptr<BundleBase> getOneBundleOlderThan_right(ptime const & t,
						bool* has_pending_punc, int node = -1) {

					*has_pending_punc = false;

					/* fast path: if local wm more recent (>t), no work.
					 * (and any unemitted punc must >t) */

					if (this->GetWatermark() > t) {
						return nullptr;
					}


					/* slowpath: there may be some work in this trans. lock containers */

					shared_ptr<BundleBase> ret;
					//unique_lock<mutex> conlock(mtx_container);
					//unique_lock<mutex> conlock(right_mtx_container_in);

					//XXX rw
					boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
					boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
					unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
						= nullptr; /* will std::move later */
start:
					if (right_containers_in.empty()) {
						//			dump_containers("no containers");
						return nullptr;
					}

					/* Start from oldest container (at back()) and look for one available
					 * bundle.
					 * For each container, we first check available data bundles and then
					 * any ready punc. If no luck, move to a newer container.
					 *
					 * When encountering an empty container, check its refcnt.
					 */

					auto it = right_containers_in.rbegin();
					/* to remember the 1st container that 1) has invalidated punc; 2) has
					 * some bundles
					 */
					auto itt = right_containers_in.rend();

					while (it != right_containers_in.rend()) {

						if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
							assert(it->verifyPuncSafe() && "bug?");
							return nullptr;
						}

						/* punc has been canceled, if we see a later container assigned with
						 * a punc, we may return a bundle from *this* container.
						 */
						if (it->punc_refcnt == PUNCREF_CANCELED) {
							if (it->refcnt == 0) { /* dead container, clean it up. */
								/* because punc only canceled after bundles are deposited.*/

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								assert(it->bundles.empty());
								it ++;
								it = list<bundle_container>::reverse_iterator(right_containers_in.erase(it.base()));
#ifdef DEBUG	
								std::cout  << __FILE__ << __LINE__ << "right_containers_in.erase" << std::endl;
#endif			
							} else if (!it->bundles.empty() && itt == right_containers_in.rend()) {
								/* bundles available + oldest container w/ punc canceled */
								itt = it;
								it ++;
							} else {
								/* else: could be no bundles in containers and some bundles are
								   outstanding... */
								it ++;
							}
							continue;
						}

						/* punc has *ever* been assigned */

						//if (it->punc->min_ts > t) { /* punc newer than wm. */
						if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
							return nullptr;
						}

						/* now we see a pending punc <=wm. */

						*has_pending_punc = true;

						/* grab from the oldest container w/ canceled punc */
						if (itt != right_containers_in.rend()) {
							//ret = itt->getBundleUnsafe();
							ret = itt->getBundleSafe();
							assert(ret);
							this->DecBundleCounter();
							//std::cout << __FILE__ << ": "<<__LINE__ << " getOneBundleOlderThan_right one" << std::endl;
							return ret;
						}

						/* non-empty container. grab & go */
						//if ((ret = it->getBundleUnsafe())) {
						if ((ret = it->getBundleSafe())) {
							this->DecBundleCounter();
							//std::cout << __FILE__ << ": "<<__LINE__ << " getOneBundleOlderThan_right one" << std::endl;
							return ret;
						}

						/* an empty container... */

						/* some data bundles outstanding.
						 * don't wait. move to a newer container. */
						if (it->refcnt != 0) {
							//				if (!(it->punc_refcnt != PUNCREF_RETRIEVED
							//										&& it->punc_refcnt != PUNCREF_CONSUMED)) {
							//					auto refcnt = it->punc_refcnt.load();
							//					E("bug: punc refcnt is %ld", refcnt);
							//				}
							assert(it->punc_refcnt != PUNCREF_RETRIEVED
									&& it->punc_refcnt != PUNCREF_CONSUMED);
							it ++;
							continue;
						}

						/* an empty container, punc has been assigned (CANCELED case is handled
						 * above). all bundles already consumed.
						 *
						 * we take a peek of punc_refcnt and decide.
						 * */

						auto punc_refcnt = it->punc_refcnt.load();

						/* first, clean any dead container so that it won't block
						 * later punc.
						 */
						if (punc_refcnt == PUNCREF_CONSUMED) { /* bundles+punc consumed */
							/* has to be the oldest container since we can't have more
							 * than one outstanding punc (the second oldest container, whose
							 * punc is assigned, must still holds its punc).
							 */

							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							assert(it == right_containers_in.rbegin());
							//				it = containers_.erase(it);   /* erase w/ reverse iterator */
							it ++;
							it = list<bundle_container>::reverse_iterator(right_containers_in.erase(it.base()));
#ifdef DEBUG
							//				std::cout << __FILE__ << __LINE__ << "right_containers_in.erase" << std::endl;
#endif			
							continue;
						} else if (punc_refcnt == PUNCREF_RETRIEVED) {
							/* punc outstanding. it may become 0 soon but we don't wait. */
							it ++;
							continue;
						} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
							if (it == right_containers_in.rbegin()) { /* oldest container. okay to emit */
								long expected = PUNCREF_ASSIGNED;
								if (!it->punc_refcnt.compare_exchange_strong(expected,
											PUNCREF_RETRIEVED)) {
									//bug("bad punc refcnt");
									it ++;
									continue;
								}
								//					auto r = it->punc_refcnt.fetch_sub(1);
								//					assert(r == 2); r = r;

								/* XXX: opt: we may check the next newer container and
								 * coalesce punc as needed. benefit unclear though.
								 */
								//return it->punc;
								return it->getPuncSafe();
							} else {
								/* not oldest container. can't emit punc before older ones. */
#if 0 /* not a good idea */
								/* opt: combine unemitted puncs.
								 *
								 * this is dangerous, since punc_refcnt may change under us.
								 * idea: temporarily decrements punc_refcnt to 1 so that
								 * other evaluators believe the punc is outstanding. after
								 * combing, increment to 2.  */
								auto & older = *(std::prev(it));
								if (older->punc_refcnt == 2) {
									older->wm = current->wm;
									//						it = containers_.erase(it); /* erase & continue */
									it ++;
									it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
								} else
#endif
									it ++;
								continue;
							}
						} else {
							EE("punc refcnt: %ld", punc_refcnt);
							assert(false && "illegal punc_perfcnt?");
						}

						assert(false && "bug?");
					} // while

					//		dump_containers("no bundles");

					return nullptr;
					}


					//call getOneBundleOlderThan_left or getOneBundleOlderThan_right
					shared_ptr<BundleBase> getOneBundleOlderThan(ptime const & t,
							bool* has_pending_punc, int node = -1) override {
						//hym: TODO
						//1. decide which side should grab bundle from
						//2. call getOneBundleOlderThan_left or getOneBundleOlderThan_right
						//XXX bad strategy now: grab from one side whose wm is small
						/*	shared_ptr<BundleBase> ret;
							if(this->left_wm =< this->right_wm){
							ret = getOneBundleOlderThan_left(t, has_pending_punc, node);
							if(ret != nullptr){
							return ret;
							}else{
							ret = getOneBundleOlderThan_right(t, has_pending_punc, node);
							return ret;
							}
							}else if(this->left_wm > this->right_wm){
							ret = getOneBundleOlderThan_right(t, has_pending_punc, node);
							if(ret != nullptr){
							return ret;
							}else{
							ret = getOneBundleOlderThan_left(t, has_pending_punc, node);
							return ret;
							}
							}else{
							assert(false && "Bug: dont's know left_wm and right_wm");	
							}
						 */      //this->dump_containers("getOneBundleOlderThan in Join");
						//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
						shared_ptr<BundleBase> ret;

						//test

						if(this->left_wm < this->right_wm){
							ret = getOneBundleOlderThan_left(t, has_pending_punc, node);
							if(ret != nullptr){
								//std::cout << __FILE__ << ": " << __LINE__ << " join get one from left in(older than)" << std::endl;
								return ret;
							}else{
								ret = getOneBundleOlderThan_right(t, has_pending_punc, node);
								if(ret != nullptr){
									//std::cout << __FILE__ << ": " << __LINE__ << " join get one from right in(older than)" << std::endl;
								}
								return ret;
							}
						}else{
							ret = getOneBundleOlderThan_right(t, has_pending_punc, node);
							if(ret != nullptr){
								//std::cout << __FILE__ << ": " << __LINE__ << " join get one from right in(older than)" << std::endl;
								return ret;
							}else{
								ret = getOneBundleOlderThan_left(t, has_pending_punc, node);
								if(ret != nullptr){
									//std::cout << __FILE__ << ": " << __LINE__ << " join get one from left in(older than)" << std::endl;
								}
								return ret;
							}
						}

					}


					shared_ptr<BundleBase> getOneBundle_left(int node = -1) {

						shared_ptr<BundleBase> ret;

						//XXX hym: here we whould use mtx_container_left and mutx_container_right

						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> conlock(left_mtx_container_in);

						//XXX rw
						boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
						boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
						unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
						if (left_containers_in.empty()) {
							//  		dump_containers("no containers");
							return nullptr;
						}

						/* Start from oldest container (at back()) and look for one available
						 * bundle.
						 * For each container, we first check available data bundles and then
						 * any ready-to-retrieve punc. If no luck, move to a newer container.
						 *
						 * when encountering an empty container (no bundles), check its refcnt.
						 */

						auto it = left_containers_in.rbegin();

						while (it != left_containers_in.rend()) {

							/* non-empty container. grab & go */
							if ((ret = it->getBundleSafe())) {
								VV("got one bundle from container %lx", (unsigned long)(&(*it)));
								this->DecBundleCounter();
								//std::cout << __FILE__ << ": " << __LINE__ << " getOneBundle_left one" << std::endl;
								return ret;
							}

							/* Can't get any bundle: an empty container... */

							if (it->punc_refcnt == PUNCREF_UNDECIDED) {
								/* punc to be assigned. since we allow multi "open" containers, we
								 * should continue to the next container... */
								assert(it->verifyPuncSafe());
								it ++;
								continue;
							}

							/* some data bundles are outstanding.
							 * don't wait. move to a newer container. */
							if (it->refcnt != 0) {
								assert(it->punc_refcnt != PUNCREF_CONSUMED
										&& it->punc_refcnt != PUNCREF_RETRIEVED);
								it ++;
								continue;
							}

							/* an empty container, punc has been determined.
							 * all bundles already consumed (refcnt = 0).
							 *
							 * take a peek of punc_refcnt ...
							 * */

							auto punc_refcnt = it->punc_refcnt.load();

							/* first, clean any dead container so that it won't block the next punc.
							 */

							if (punc_refcnt == PUNCREF_CONSUMED) {
								/* has to be the oldest container since we can't have more
								 * than one outstanding punc (the second oldest container, whose
								 * punc is assigned, must still holds its punc).
								 */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								assert(it == left_containers_in.rbegin());
								//				it = containers_.erase(it);   /* erase w/ reverse iterator */
								it ++;
								it = list<bundle_container>::reverse_iterator(left_containers_in.erase(it.base()));
								continue;
							} else if (punc_refcnt == PUNCREF_CANCELED) {
								/* erase */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								it ++;
								it = list<bundle_container>::reverse_iterator(left_containers_in.erase(it.base()));
							} else if (punc_refcnt == PUNCREF_RETRIEVED) {
								/* punc outstanding. it may become CONSUMED soon but we don't wait. */
								it ++;
								continue;
							} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
								if (it == left_containers_in.rbegin()) {
									/* oldest container. (older dead containers should have been
									 * cleaned up now). okay to emit */
									//					auto r = it->punc_refcnt.fetch_sub(1);
									//					assert(r == 2); r = r;

									long expected = PUNCREF_ASSIGNED;
									if (!it->punc_refcnt.compare_exchange_strong(expected,
												PUNCREF_RETRIEVED)) {
										//bug("bad punc refcnt");
										it ++;
										continue;
									}
									/* XXX: opt: we may check the next newer container and
									 * coalesce punc as needed. benefit unclear though.
									 */
									//return it->punc;
									return it->getPuncSafe();
								} else {
									/* not oldest container. can't emit punc before older ones doing so. */
#if 0 /* not a good idea */
									/* opt: combine unemitted puncs.
									 *
									 * this is dangerous, since punc_refcnt may change under us.
									 * idea: temporarily decrements punc_refcnt to 1 so that
									 * other evaluators believe the punc is outstanding. after
									 * combing, increment to 2.  */
									auto & older = *(std::prev(it));
									if (older->punc_refcnt == 2) {
										older->wm = current->wm;
										//						it = containers_.erase(it); /* erase & continue */
										it ++;
										it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
									} else
#endif
										it ++;
									continue;
								}
							} else {
								EE("punc refcnt: %ld", punc_refcnt);
								assert(false && "illegal punc_perfcnt?");
							}

							assert(false && "bug?");

						} // while

						//  	dump_containers("no bundles");

						return nullptr;
					}


					shared_ptr<BundleBase> getOneBundle_right(int node = -1) {

						shared_ptr<BundleBase> ret;

						//hym: XXX here we should use left_mtx_container and right_mtx_container

						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> conlock(right_mtx_container_in);

						//XXX rw
						boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
						boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
						unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
						if (right_containers_in.empty()) {
							//  		dump_containers("no containers");
							return nullptr;
						}

						/* Start from oldest container (at back()) and look for one available
						 * bundle.
						 * For each container, we first check available data bundles and then
						 * any ready-to-retrieve punc. If no luck, move to a newer container.
						 *
						 * when encountering an empty container (no bundles), check its refcnt.
						 */

						auto it = right_containers_in.rbegin();

						while (it != right_containers_in.rend()) {

							/* non-empty container. grab & go */
							if ((ret = it->getBundleSafe())) {
								VV("got one bundle from container %lx", (unsigned long)(&(*it)));
								this->DecBundleCounter();
								//std::cout << __FILE__ << ": " << __LINE__ << " getOneBUndle_right one" << std::endl;
								return ret;
							}

							/* Can't get any bundle: an empty container... */

							if (it->punc_refcnt == PUNCREF_UNDECIDED) {
								/* punc to be assigned. since we allow multi "open" containers, we
								 * should continue to the next container... */
								assert(it->verifyPuncSafe());
								it ++;
								continue;
							}

							/* some data bundles are outstanding.
							 * don't wait. move to a newer container. */
							if (it->refcnt != 0) {
#if 0				/* can no longer verify the following due to concurrent readers. unless we lock the container */
								assert(it->punc_refcnt != PUNCREF_CONSUMED
										&& it->punc_refcnt != PUNCREF_RETRIEVED);
#endif
								it ++;
								continue;
							}

							/* an empty container, punc has been determined.
							 * all bundles already consumed (refcnt = 0).
							 *
							 * take a peek of punc_refcnt ...
							 * */

							auto punc_refcnt = it->punc_refcnt.load();

							/* first, clean any dead container so that it won't block the next punc.
							 */

							if (punc_refcnt == PUNCREF_CONSUMED) {
								/* has to be the oldest container since we can't have more
								 * than one outstanding punc (the second oldest container, whose
								 * punc is assigned, must still holds its punc).
								 */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								assert(it == right_containers_in.rbegin());
								//				it = containers_.erase(it);   /* erase w/ reverse iterator */
								it ++;
								it = list<bundle_container>::reverse_iterator(right_containers_in.erase(it.base()));
								continue;
							} else if (punc_refcnt == PUNCREF_CANCELED) {
								/* erase */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								it ++;
								it = list<bundle_container>::reverse_iterator(right_containers_in.erase(it.base()));
							} else if (punc_refcnt == PUNCREF_RETRIEVED) {
								/* punc outstanding. it may become CONSUMED soon but we don't wait. */
								it ++;
								continue;
							} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
								if (it == right_containers_in.rbegin()) {
									/* oldest container. (older dead containers should have been
									 * cleaned up now). okay to emit */
									//					auto r = it->punc_refcnt.fetch_sub(1);
									//					assert(r == 2); r = r;

									long expected = PUNCREF_ASSIGNED;
									if (!it->punc_refcnt.compare_exchange_strong(expected,
												PUNCREF_RETRIEVED)) {
										//bug("bad punc refcnt");
										it ++;
										continue;
									}
									/* XXX: opt: we may check the next newer container and
									 * coalesce punc as needed. benefit unclear though.
									 */
									//return it->punc;
									return it->getPuncSafe();
								} else {
									/* not oldest container. can't emit punc before older ones doing so. */
#if 0 /* not a good idea */
									/* opt: combine unemitted puncs.
									 *
									 * this is dangerous, since punc_refcnt may change under us.
									 * idea: temporarily decrements punc_refcnt to 1 so that
									 * other evaluators believe the punc is outstanding. after
									 * combing, increment to 2.  */
									auto & older = *(std::prev(it));
									if (older->punc_refcnt == 2) {
										older->wm = current->wm;
										//						it = containers_.erase(it); /* erase & continue */
										it ++;
										it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
									} else
#endif
										it ++;
									continue;
								}
							} else {
								EE("punc refcnt: %ld", punc_refcnt);
								assert(false && "illegal punc_perfcnt?");
							}

							assert(false && "bug?");

						} // while

						//  	dump_containers("no bundles");

						return nullptr;
					}

					//hym: override PTransform's getOneBundle
					//XXX: bad strategy now, only grab a bundle from one side whose wm is smaller
					shared_ptr<BundleBase> getOneBundle(int node = -1) {
						/*	shared_ptr<BundleBase> ret;
							if(this->left_wm =< this->right_wm){
							ret = getOneBundle_left(node);
							if(ret != nullptr){
							return ret;
							}else{
							ret = getOneBundle_right(node);
							return ret;
							}
							}else if(this->left_wm > this->right_wm){
							ret = getOneBundle_right(node);
							if(ret != nullptr){
							return ret;
							}else{
							ret = getOneBundle_left(node);
							return ret;
							}
							}else{
							assert(false && "Bug: dont's know left_wm and right_wm");	
							}
						 */
						//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
						shared_ptr<BundleBase> ret;
						//test
						//ret = getOneBundle_left(node);
						//return ret;

						if(this->left_wm < this->right_wm){
							ret = getOneBundle_left(node);
							if(ret != nullptr){
								//std::cout << __FILE__ << ": " << __LINE__ << " join get one from left in" << std::endl;
								return ret;
							}else{
								ret = getOneBundle_right(node);
								if(ret != nullptr){
									//std::cout << __FILE__ << ": " << __LINE__ << " join get one from right in" << std::endl;
								}
								return ret;
							}
						}else {
							ret = getOneBundle_right(node);
							if(ret != nullptr){
								//std::cout << __FILE__ << ": " << __LINE__ << " join get one from right in" << std::endl;
								return ret;
							}else{
								ret = getOneBundle_left(node);
								if(ret != nullptr){
									//std::cout << __FILE__ << ": " << __LINE__ << " join get one from left in" << std::endl;
								}
								return ret;
							}
						}

					}
					/////////////////////////////////////////////////////////////
					//      new design: deposit bundle to downstream's left_unorderd_containers
					//	              or right_unordered_containers, instead of Join's left_unordered_
					//                  containers and right_unordered_containers
					/////////////////////////////////////////////////////////////

					//hym: open a new container on Join's downstream's left_unordered_containers
					int openJoinDownstreamContainer_left(bundle_container * const upcon,
							PTransform *downt){
						//hym TODO:
#ifdef DEBUG	
						//	std::cout << __FILE__ << __LINE__ << "openJoinDownstreamContainer_left" << std::endl;
#endif
						assert(downt && upcon);

						//  	/* fast path w/o locking */
						//  	if (upcon->downstream)
						//  		return 0;

						/* slow path
						 *
						 * in *this* trans, walk from the oldest container until @upcon,
						 * open a downstream container if there isn't one.
						 *
						 * (X = no downstream container;  V = has downstream container)
						 * possible:
						 *
						 * X X V V (oldest)
						 *
						 * impossible:
						 * X V X V (oldest)
						 *
						 */
						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> conlock(left_mtx_container_in);
						//unique_lock<mutex> conlock(left_mtx_container_in);

						//depositBundleDwonstream has get a big lock outside this func
#if 0
						unique_lock<mutex> conlock(mtx_container);
						unique_lock<mutex> down_conlock(downt->mtx_container);
#endif
						if (upcon->downstream) /* check again */
							return 0;

						int cnt = 0;
						bool miss = false; (void)sizeof(miss);
						for (auto it = this->left_containers_in.rbegin();
								it != this->left_containers_in.rend(); it ++) {

							/* also lock all downstream containers. good idea? */
							//  		unique_lock<mutex> down_conlock(downt->mtx_container);

							if (!it->downstream) {
								miss = true;
								/* alloc & link downstream container. NB: it->downstream is an atomic
								 * pointer, so it implies a fence here. */
								downt->left_unordered_containers.emplace_front();
								it->downstream = &(downt->left_unordered_containers.front());

								//hym: XXX 
								//assigne the side_info to this new container
								//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
								//it->downstream.load()->set_side_info(this->get_side_info());
								it->downstream.load()->set_side_info(it->get_side_info());
								assert(it->downstream.load()->get_side_info() == 1);	//left		

								cnt ++;
							} else { /* downstream link exists */
								assert (!miss && "bug? an older container misses downstream link.");
#ifdef DEBUG
								bool found = false;
								for (auto & con : downt->left_unordered_containers) {
									if (&con == it->downstream) {
										found = true;
										break;
									}
								}
								/* NB: downstream container, once received & processed the punc
								 * from the upstream, may be cleaned up prior to the upstream
								 * container. In this case, a dangling downstream pointer is possible.
								 */
								if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
									EE("%s: downstream (%s) container does not exist",
											this->name.c_str(), downt->name.c_str());
									cout << " -------------------------------\n";
									dump_containers("bug");
									downt->dump_containers("bug");
									cout << " -------------------------------\n\n";
									abort();
								}
#endif
							}
						}

#if 0 /* debug */
						cout << " -------------------------------\n";
						dump_containers("link containers");
						downt->dump_containers("link containers");
						cout << " -------------------------------\n\n";
#endif

						return cnt;
					}

					//open one containers in Join's downstream's right_unordered_containers
					int openJoinDownstreamContainer_right(bundle_container * const upcon,
							PTransform *downt){
						//hym TODO:

						assert(downt && upcon);

						//  	/* fast path w/o locking */
						//  	if (upcon->downstream)
						//  		return 0;

						/* slow path
						 *
						 * in *this* trans, walk from the oldest container until @upcon,
						 * open a downstream container if there isn't one.
						 *
						 * (X = no downstream container;  V = has downstream container)
						 * possible:
						 *
						 * X X V V (oldest)
						 *
						 * impossible:
						 * X V X V (oldest)
						 *
						 */
						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> conlock(right_mtx_container_in);
						//unique_lock<mutex> down_conlock(downt->right_mtx_unordered_containers);

						//depositBundleDownstream has gotten big locks outside this func
#if 0
						unique_lock<mutex> conlock(mtx_container);
						unique_lock<mutex> down_conlock(downt->mtx_container);
#endif
						if (upcon->downstream) /* check again */
							return 0;

						int cnt = 0;
						bool miss = false; (void)sizeof(miss);
						for (auto it = this->right_containers_in.rbegin();
								it != this->right_containers_in.rend(); it ++) {

							/* also lock all downstream containers. good idea? */
							//  		unique_lock<mutex> down_conlock(downt->mtx_container);

							if (!it->downstream) {
								miss = true;
								/* alloc & link downstream container. NB: it->downstream is an atomic
								 * pointer, so it implies a fence here. */
								downt->right_unordered_containers.emplace_front();
								it->downstream = &(downt->right_unordered_containers.front());

								//hym: XXX 
								//assigne the side_info to this new container
								//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
								//it->downstream.load()->set_side_info(this->get_side_info());
								it->downstream.load()->set_side_info(it->get_side_info());
								assert(it->downstream.load()->get_side_info() == 2);//right

								cnt ++;
							} else { /* downstream link exists */
								assert (!miss && "bug? an older container misses downstream link.");
#ifdef DEBUG
								bool found = false;
								for (auto & con : downt->right_unordered_containers) {
									if (&con == it->downstream) {
										found = true;
										break;
									}
								}
								/* NB: downstream container, once received & processed the punc
								 * from the upstream, may be cleaned up prior to the upstream
								 * container. In this case, a dangling downstream pointer is possible.
								 */
								if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
									EE("%s: downstream (%s) container does not exist",
											this->name.c_str(), downt->name.c_str());
									cout << " -------------------------------\n";
									dump_containers("bug");
									downt->dump_containers("bug");
									cout << " -------------------------------\n\n";
									abort();
								}
#endif
							}
						}

#if 0 /* debug */
						cout << " -------------------------------\n";
						dump_containers("link containers");
						downt->dump_containers("link containers");
						cout << " -------------------------------\n\n";
#endif

						return cnt;
					}


					//depositBundleDownstreamLR: left/right_unordered_containers
					//downt should has left_unordere_containers, right_unordered_containers, and
					//ordered_containers. 
					//downt's type should be 4
					void  depositBundleDownstreamLR(PTransform *downt,
							shared_ptr<BundleBase> input_bundle,
							vector<shared_ptr<BundleBase>> const & output_bundles, int node = -1){
						//TODO: deposit Left or Right
						assert(downt->get_side_info() == 4);

						bool try_open = false; (void)sizeof(try_open); /* have we tried open downstream containers? */
#ifdef DEBUG
						//	std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
						//std::cout << "trans is: " << downt->getName() << std::endl;
#endif
						assert(downt);
						bundle_container *upcon = localBundleToContainer(input_bundle);
						assert(upcon); /* enclosing container */

						//XXX rw
						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> down_conlock(downt->mtx_container);

						/* reader lock for *this* trans */
						boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
						/* writer lock for downstream trans */ 
						boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);

deposit:

						{		//{
							//		unique_lock<mutex> conlock(mtx_container);
							//		unique_lock<mutex> down_conlock(downt->mtx_container);

							if (upcon->downstream) {
								/* fast path. NB @downstream is atomic w/ seq consistency so it's
								 * safe. */
								downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
								downt->IncBundleCounter(output_bundles.size());
#ifdef DEBUG
								//				std::cout << __FILE__ << ": " << __LINE__ << "join deposits bundle to its out: " << upcon->get_side_info() << std::endl;
#endif
								return;
							}
							//}

#if 0		
							if(upcon->get_side_info() == 1){ //deposit type4's left_unorderd_containers
								unique_lock<mutex> conlock(downt->left_mtx_unordered_containers);
								if (upcon->downstream) {
									downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
									downt->IncBundleCounter(output_bundles.size());
									return;
								}
							}else if(upcon->get_side_info() == 2){ //deposit type4's right_unordered_containers
								unique_lock<mutex> conlock(downt->right_mtx_unordered_containers);
								if (upcon->downstream) {
									downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
									downt->IncBundleCounter(output_bundles.size());
									return;
								}
							}else{
								assert(false  && "Unknow side info");
							}
#endif
							assert(!try_open && "bug: already tried to open");

							/* it's possible somebody slipped in & opened the downstream container
							 * already. So ret == 0 then.
							 */
							/*
							//  	int ret =
							openDownstreamContainers(upcon, downt);
							//  	assert(ret && "no downstream container opened?");
							 */

							/*
							   if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
							   openDownstreamContainers(upcon, downt);
							   }else if(this->get_side_info() == 1){ //Simplemapper 1, left
							   openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in		
							   }else if(this->get_side_info() == 2){ //Simplemapper 2, right
							   openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
							   }else{
							   assert("Bug: Unknown know side info");
							   }
							 */
							if(upcon->get_side_info() == 1){
								//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;	
								//openJoinOutContainer_left(upcon, downt);
								openJoinDownstreamContainer_left(upcon, downt);
								//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
							}else if(upcon->get_side_info() == 2){
								//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
								//openJoinOutContainer_right(upcon, downt);
								openJoinDownstreamContainer_right(upcon, downt);
								//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
							}else{
								assert("Join: wrong side info");
							}
						}//end lock
						assert(upcon->downstream);
						try_open = true;
						goto deposit;
					}

					//hym: Join deposits a punc to a container in its downstream's 
					//     left_unordered_containers or right_unordered_containers
					void depositPuncDownstreamLR(PTransform* downt,
							shared_ptr<Punc> const & input_punc, shared_ptr<Punc> punc, int node = -1)
					{
						bool try_open = false; (void)sizeof(try_open);

						assert(downt);
						bundle_container *upcon = localBundleToContainer(input_punc);
						assert(upcon); /* enclosing container */

						//XXX rw
						//unique_lock<mutex> conlock(mtx_container);
						//unique_lock<mutex> down_conlock(downt->mtx_container);

						/* reader lock for *this* trans */
						boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
						/* writer lock for downstream trans */ 
						boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);

deposit:
						{
							//{
							//		unique_lock<mutex> conlock(mtx_container);
							//		unique_lock<mutex> down_conlock(downt->mtx_container);
							if (upcon->downstream) { /* fast path. NB @downstream is atomic, implying fence */
								downt->depositOnePuncToContainer(upcon->downstream, punc, node);
								//upcon->has_punc = 1; //the upcon has a container already
								upcon->downstream.load()->has_punc = 1;
								//std::cout << __FILE__ << ": " <<  __LINE__ << " " << this->getName() << " deposit a punc to " << downt->getName() << " to side: " << upcon->get_side_info()  << std::endl;
								return;
							}
							//}
#if 0
							if(upcon->get_side_info() == 1){ //deposit type4's left_unorderd_containers
								unique_lock<mutex> conlock(downt->left_mtx_unordered_containers);
								if (upcon->downstream) {
									downt->depositOnePuncToContainer(upcon->downstream, punc, node);
									upcon->downstream.load()->has_punc = 1;
									return;
								}
							}else if(upcon->get_side_info() == 2){ //deposit type4's right_unordered_containers
								unique_lock<mutex> conlock(downt->right_mtx_unordered_containers);
								if (upcon->downstream) {
									downt->depositOnePuncToContainer(upcon->downstream, punc, node);
									upcon->downstream.load()->has_punc = 1;
									return;
								}
							}else{
								assert(false  && "Unknow side info");
							}
#endif

							assert(!try_open && "bug: already tried to open");

							//openDownstreamContainers(upcon, downt);
							//hym:XXX
							/*	if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
								openDownstreamContainers(upcon, downt);
								}else if(this->get_side_info() == 1){ //Simplemapper 1, left
								openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in		
								}else if(this->get_side_info() == 2){ //Simplemapper 2, right
								openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
								}else{
								assert("Bug: Unknown know side info");
								}
							 */

							if(upcon->get_side_info() == 1){
								//openJoinOutContainer_left(upcon, downt);
								openJoinDownstreamContainer_left(upcon, downt);
							}else if(upcon->get_side_info() == 2){
								//openJoinOutContainer_right(upcon, downt);
								openJoinDownstreamContainer_right(upcon, downt);
							}else{
								assert("Join: wrong side info");
							}
						}//end lock
						assert(upcon->downstream);
						try_open = true;
						goto deposit;
					}

#if 0
					////////////////////////////////////////////////////
					// hym: Following funcs are for Join's new source UnboundedInMemEvaluator_Join
					///////////////////////////////////////////////////

					//hym: depositOneBundleToJoin_L and depositOneBundleToJoin_R will be called by Source
					//hym: Source deposit one bundle to Join's left_containers_in
					void depositOneBundleToJoin_L(shared_ptr<BundleBase> bundle, int node = -1) {

						/* lock all containers until reliably getting the target container.
						 * (in case concurrent additions of containers) */
						//XXX rw
						//unique_lock<mutex> conlock(mtx_container);
						boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
						boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
						unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
						if (left_containers_in.empty()) {
							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}
							left_containers_in.emplace_front(); /* can't copy */
						}

						auto current_con = left_containers_in.begin();
						//current_con->set_side_info(this->get_side_info());
						current_con->set_side_info(1); //left

						//if (current_con->punc) { /* the latest container already has a punc */
						if (current_con->getPuncSafe()) { /* the latest container already has a punc */
							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							assert(current_con->punc_refcnt != PUNCREF_UNDECIDED && "invalid punc refcnt");

							/* The latest container is already dead: clean up.
							 * NB: even if we see punc_refcnt == 1 here, it may be dead
							 * (decremented to 0) beneath us any time. It will be cleaned by future
							 * calls.
							 */
							if (current_con->punc_refcnt == PUNCREF_CONSUMED) {
								/* extensive sanity check... */

								//assert(current_con->punc && "dead container but no valid punc?");
								assert(current_con->getPuncSafe() && "dead container but no valid punc?");
								assert(current_con->bundles.size() == 0);
								assert(current_con->refcnt == 0);

								current_con = left_containers_in.erase(current_con);
								//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
								if (left_containers_in.empty()) {
									/* after cleanup, no container left. start a new one. */
									left_containers_in.emplace_front();
									current_con = left_containers_in.begin();
									goto open_container;
								} else {
									/* current_con now points to the 2nd most recent container.
									 * sanity check it */
									assert(current_con->punc_refcnt != PUNCREF_CONSUMED
											&& "bug: can't have two dead containers");
								}
							}

							/* careful when reading from an atomic/volatile value ... */
							auto prefcnt = current_con->punc_refcnt.load();
							prefcnt = prefcnt;
							assert( prefcnt == PUNCREF_CONSUMED   /* latest container dies just now */
									|| prefcnt == PUNCREF_RETRIEVED		/* outstanding */
									|| prefcnt == PUNCREF_ASSIGNED		/* unretrieved yet */
							      );

							left_containers_in.emplace_front();
							current_con --;
						}

						/* Now @current_con points to an open container. Insert the bundle.
						 *
						 * @current_con won't go away since the punc won't be emitted
						 * concurrently.
						 */
open_container:
						/* should have any type of lock */
						//XXX rw
						assert(rlock.owns_lock() || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

						//current_con->putBundleUnsafe(bundle);
						current_con->putBundleSafe(bundle);
						this->IncBundleCounter();

					} //end depositOneBundleToJoin_L


					//hym: Source deposit one bundle to Join's right_containers_in
					void depositOneBundleToJoin_R(shared_ptr<BundleBase> bundle, int node = -1) {

						/* lock all containers until reliably getting the target container.
						 * (in case concurrent additions of containers) */
						//XXX rw
						//unique_lock<mutex> conlock(mtx_container);
						boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
						boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
						unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
						if (right_containers_in.empty()) {
							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}
							right_containers_in.emplace_front(); /* can't copy */
						}

						auto current_con = right_containers_in.begin();
						//current_con->set_side_info(this->get_side_info());
						current_con->set_side_info(2); //right
						//if (current_con->punc) { /* the latest container already has a punc */
						if (current_con->getPuncSafe()) { /* the latest container already has a punc */
							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							assert(current_con->punc_refcnt != PUNCREF_UNDECIDED && "invalid punc refcnt");

							/* The latest container is already dead: clean up.
							 * NB: even if we see punc_refcnt == 1 here, it may be dead
							 * (decremented to 0) beneath us any time. It will be cleaned by future
							 * calls.
							 */
							if (current_con->punc_refcnt == PUNCREF_CONSUMED) {
								/* extensive sanity check... */

								//assert(current_con->punc && "dead container but no valid punc?");
								assert(current_con->getPuncSafe() && "dead container but no valid punc?");
								assert(current_con->bundles.size() == 0);
								assert(current_con->refcnt == 0);

								current_con = right_containers_in.erase(current_con);
								//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
								if (right_containers_in.empty()) {
									/* after cleanup, no container left. start a new one. */
									right_containers_in.emplace_front();
									current_con = right_containers_in.begin();
									goto open_container;
								} else {
									/* current_con now points to the 2nd most recent container.
									 * sanity check it */
									assert(current_con->punc_refcnt != PUNCREF_CONSUMED
											&& "bug: can't have two dead containers");
								}
							}

							/* careful when reading from an atomic/volatile value ... */
							auto prefcnt = current_con->punc_refcnt.load();
							prefcnt = prefcnt;
							assert( prefcnt == PUNCREF_CONSUMED   /* latest container dies just now */
									|| prefcnt == PUNCREF_RETRIEVED		/* outstanding */
									|| prefcnt == PUNCREF_ASSIGNED		/* unretrieved yet */
							      );

							right_containers_in.emplace_front();
							current_con --;
						}

						/* Now @current_con points to an open container. Insert the bundle.
						 *
						 * @current_con won't go away since the punc won't be emitted
						 * concurrently.
						 */
open_container:
						/* should have any type of lock */
						//XXX rw
						assert(rlock.owns_lock() || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

						//current_con->putBundleUnsafe(bundle);
						current_con->putBundleSafe(bundle);
						this->IncBundleCounter();

					} //end depositOneBundleToJoin_R

					// hym: this func is only called by source
					//      deposit a punc to Join's left_containers_in
					void depositOnePuncToJoin_L(shared_ptr<Punc> punc, int node = -1) {

						//std::cout << "depositOnePunc() in trans: " << this->getName() << std::endl;
#ifdef MEASURE_LATENCY

						//  	char buf[50];
						//  	long cnt = this->getNumBundles();
						//  	snprintf(buf, 50, "deposit to: %s (%ld bundles) ",
						//  			this->name.c_str(), cnt);
						//  	punc->mark(buf);
						/*
						   ostringstream oss;
						   vector<long> cnts;
						   this->getNumBundles(&cnts);
						   oss << "deposit to: " << this->name << " (bundles: ";
						   for (auto & cnt : cnts) {
						   oss << cnt << " ";
						   }
						   oss << ")";
						   punc->mark(oss.str());
						 */
						ostringstream oss;
						vector<cont_info> cnts;
						this->getNumBundlesTs(&cnts);
						oss << "deposit to: " << this->name <<  " " << boost::posix_time::microsec_clock::local_time() <<" " << "\n (bundles: ";
						for (auto & cnt : cnts) {
							if (cnt.bundle >= 0)
								oss << cnt.bundle << " ";
							else {
								auto r = cnt.bundle + PUNCREF_MAX;
								oss << puncref_key(r) << " ";
							}
							oss << cnt.side_info << " ";
							//  			if (r >= 0) { /* punc ever assigned. show punc ts */
							if (cnt.punc != max_date_time) {
								oss << to_simple_string(cnt.punc) << "\t";
							}
						}
						oss << ")";
						punc->mark(oss.str());

						//  	punc->mark("deposit to: " + this->name);
#endif
						//XXX rw
						//unique_lock<mutex> conlock(mtx_container);
						boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
						boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
						unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
						if (left_containers_in.empty()) {
							//std::cout << __FILE__ << ":" << __LINE__ << " containers_.empty() wrong???" << std::endl;
							/* Whole transform drained. we don't even have pending flushed
							 * bundles (otherwise there will be containers).
							 * Create a new container and seal it with @wm.
							 */
							//   		containers_.push_front(bundle_container());

							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}
							left_containers_in.emplace_front();
							//containers_.begin()->setPunc(punc);
							left_containers_in.begin()->setPuncSafe(punc);
							return;
						}
						//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
						/* assign or update the wm of the most recent container */
						auto current = left_containers_in.begin();
						if (current->punc_refcnt == PUNCREF_RETRIEVED) {
							/* punc already emitted. start a new container. */

							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							assert(current->punc && "must have a valid punc ptr");
							assert(!current->refcnt && "all bundles must be consumed");
							//   		containers_.push_front(bundle_container());
							left_containers_in.emplace_front();
							//containers_.begin()->setPunc(punc);
							left_containers_in.begin()->setPuncSafe(punc);
							return;
						}
						//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
						/* clean up any dead container. */
						if (current->punc_refcnt == PUNCREF_CONSUMED) {
							/* extensive sanity check... */

							//XXX rw
							if (!upgrade_locks(&rlock, &ulock, &pwlock)){
								goto start;
							}

							//assert(current->punc && "dead container but no valid punc?");
							assert(current->getPuncSafe() && "dead container but no valid punc?");
							assert(current->bundles.size() == 0);
							assert(current->refcnt == 0);

							current = left_containers_in.erase(current);
							//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
							if (left_containers_in.empty()) {
								/* after cleanup, no container left. start a new one. */
								left_containers_in.emplace_front();
								//containers_.begin()->setPunc(punc);
								left_containers_in.begin()->setPuncSafe(punc);
								return;
							} else {
								/* current now points to the 2nd most recent container. */
								assert(current->punc_refcnt != 0
										&& "bug: can't have two dead containers");
							}
						}
						//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
						//if (!current->punc || current->punc_refcnt == PUNCREF_ASSIGNED) {
						if (!current->getPuncSafe() || current->punc_refcnt == PUNCREF_ASSIGNED) {
							//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
							//current->setPunc(punc); /* assign or overwrite punc */
							current->setPuncSafe(punc); /* assign or overwrite punc */
							//current->set_side_info(this->get_side_info());
						} else {
							assert(false && "bug?");
						}
						}


						// hym: this func is only called by source
						//      deposit a punc to Join's left_containers_in
						void depositOnePuncToJoin_R(shared_ptr<Punc> punc, int node = -1) {

							//std::cout << "depositOnePunc() in trans: " << this->getName() << std::endl;
#ifdef MEASURE_LATENCY

							//  	char buf[50];
							//  	long cnt = this->getNumBundles();
							//  	snprintf(buf, 50, "deposit to: %s (%ld bundles) ",
							//  			this->name.c_str(), cnt);
							//  	punc->mark(buf);
							/*
							   ostringstream oss;
							   vector<long> cnts;
							   this->getNumBundles(&cnts);
							   oss << "deposit to: " << this->name << " (bundles: ";
							   for (auto & cnt : cnts) {
							   oss << cnt << " ";
							   }
							   oss << ")";
							   punc->mark(oss.str());
							 */
							ostringstream oss;
							vector<cont_info> cnts;
							this->getNumBundlesTs(&cnts);
							oss << "deposit to: " << this->name <<  " " << boost::posix_time::microsec_clock::local_time() <<" " << "\n (bundles: ";
							for (auto & cnt : cnts) {
								if (cnt.bundle >= 0)
									oss << cnt.bundle << " ";
								else {
									auto r = cnt.bundle + PUNCREF_MAX;
									oss << puncref_key(r) << " ";
								}
								oss << cnt.side_info << " ";
								//  			if (r >= 0) { /* punc ever assigned. show punc ts */
								if (cnt.punc != max_date_time) {
									oss << to_simple_string(cnt.punc) << "\t";
								}
							}
							oss << ")";
							punc->mark(oss.str());

							//  	punc->mark("deposit to: " + this->name);
#endif
							//XXX rw
							//unique_lock<mutex> conlock(mtx_container);
							boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
							boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
							unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
								= nullptr; /* will std::move later */
start:
							if (right_containers_in.empty()) {
								//std::cout << __FILE__ << ":" << __LINE__ << " containers_.empty() wrong???" << std::endl;
								/* Whole transform drained. we don't even have pending flushed
								 * bundles (otherwise there will be containers).
								 * Create a new container and seal it with @wm.
								 */
								//   		containers_.push_front(bundle_container());

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}
								right_containers_in.emplace_front();
								//containers_.begin()->setPunc(punc);
								right_containers_in.begin()->setPuncSafe(punc);
								return;
							}
							//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
							/* assign or update the wm of the most recent container */
							auto current = right_containers_in.begin();
							if (current->punc_refcnt == PUNCREF_RETRIEVED) {
								/* punc already emitted. start a new container. */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								assert(current->punc && "must have a valid punc ptr");
								assert(!current->refcnt && "all bundles must be consumed");
								//   		containers_.push_front(bundle_container());
								right_containers_in.emplace_front();
								//containers_.begin()->setPunc(punc);
								right_containers_in.begin()->setPuncSafe(punc);
								return;
							}
							//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
							/* clean up any dead container. */
							if (current->punc_refcnt == PUNCREF_CONSUMED) {
								/* extensive sanity check... */

								//XXX rw
								if (!upgrade_locks(&rlock, &ulock, &pwlock)){
									goto start;
								}

								//assert(current->punc && "dead container but no valid punc?");
								assert(current->getPuncSafe() && "dead container but no valid punc?");
								assert(current->bundles.size() == 0);
								assert(current->refcnt == 0);

								current = right_containers_in.erase(current);
								//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
								if (right_containers_in.empty()) {
									/* after cleanup, no container left. start a new one. */
									right_containers_in.emplace_front();
									//containers_.begin()->setPunc(punc);
									right_containers_in.begin()->setPuncSafe(punc);
									return;
								} else {
									/* current now points to the 2nd most recent container. */
									assert(current->punc_refcnt != 0
											&& "bug: can't have two dead containers");
								}
							}
							//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
							//if (!current->punc || current->punc_refcnt == PUNCREF_ASSIGNED) {
							if (!current->getPuncSafe() || current->punc_refcnt == PUNCREF_ASSIGNED) {
								//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
								//current->setPunc(punc); /* assign or overwrite punc */
								current->setPuncSafe(punc); /* assign or overwrite punc */
								//current->set_side_info(this->get_side_info());
							} else {
								assert(false && "bug?");
							}
							}//end depositOnePuncToJoin_R

#endif

							};
#endif /* JOIN_H */
