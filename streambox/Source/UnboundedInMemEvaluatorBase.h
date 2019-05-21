/*
 * An evaluator for the unbounded source; has internal parallelism
 * output: bundles of strings. each string may contain spaces.
 *
 * (later one string can go through mapper and become kv pairs.)
 *
 * @T: element type, e.g. long
 */
#ifndef UNBOUNDED_INMEM_EVALBASE_H
#define UNBOUNDED_INMEM_EVALBASE_H

#include <numa.h>
//#include "utilities/threading.hh"
//using namespace Kaskade;

#include "Source/Unbounded.h"
#include "core/TransformEvaluator.h"
#include "core/EvaluationBundleContext.h"

template<class T, template<class> class BundleT>
class UnboundedInMemEvaluatorBase
: public TransformEvaulator<UnboundedInMem<T, BundleT>> {

	public:
		using TransformT = UnboundedInMem<T, BundleT>;

		/* obsoleted: by bundle_container */
		bool is_pressure_high_pvalue(EvaluationBundleContext *c) {
			assert(!c->_values.empty());

			/* an ad hoc way to calculate pressure. can be improved.
			 * keep in mind that we have, e.g. 40 cpu cores
			 */
			const int local_threshold = 160;
			const int global_threshold = 400;

			int total_pending = 0;
			for (auto && v : c->_values) {
				/* just peek. we don't grab lock */
				int pending = v->size();
				if (pending > local_threshold)
					return true;
				total_pending += pending;
			}

			if (total_pending > global_threshold)
				return true;
			else
				return false;
		}

		/* the "bundle container" version. bundles no longer saved at pvalues */
		bool is_pressure_high_bundle_containers(EvaluationBundleContext *c) {
			assert(!c->_transforms.empty());

			/* an ad hoc way to calculate pressure. can be improved.
			 * keep in mind that we have, e.g. 40 cpu cores
			 *
			 * idea: should check the *trend* of the bundle/container growth.
			 * also only react if the wm latency is high.
			 *
			 * when wm lat is high:
			 * 	if there're backlog for bundles/containers --> reduce the source rate
			 * 	otherwise, we shouldn't (reducing rate further inc delay)
			 *
			 */

			const int local_threshold = 160;
			const int global_threshold = 400;

			//    const int local_threshold = 40;
			//    const int global_threshold = 100;

			//        const int local_threshold = 40;
			//        const int global_threshold = 80;

			/* container threshold?
			 * in fact, we should check *empty* containers */
			//    const int local_container_threshold = 5;
			//    const int global_container_threshold = 20;

			const int local_container_threshold = 999;
			const int global_container_threshold = 999;

			long total_pending = 0, total_con = 0;
			for (auto && t : c->_transforms) {
				long pending = t->getNumBundles();
				if (pending > local_threshold)
					return true;
				total_pending += pending;

				long con = t->containers_.size();
				if (con > local_container_threshold)
					return true;
				total_con += con;
			}

			if (total_pending > global_threshold
					|| total_con > global_container_threshold)
				return true;
			else
				return false;
		}

		/* Sensing all downstream pressure. if any value is clogged, the source
		 * simply backs off.
		 *
		 * This should be pretty cheap; the source can afford the check from time
		 * to time.
		 */
		bool is_pressure_high(EvaluationBundleContext *c) {
			return is_pressure_high_bundle_containers(c);
		}


		/* for debugging */
		void dump_all_containers(EvaluationBundleContext *c) {
			for (auto && t : c->_transforms) {
				t->dump_containers();
			}
		}

	protected:
		/* @bytes, @records: the "net" bytes and records that have just been
		 * processed.
		 * return: true if we print out anything */
#if 0
		bool report_progress(uint64_t bytes, uint64_t records) {
			uint64_t total_bytes

				total_bytes += bytes;
			total_records += records;

			ptime now = boost::posix_time::microsec_clock::local_time();

			if (once) {
				once = 0;
				last_check = now;
				start_time = now;
				last_bytes = bytes;
				last_records = records;
				return false;
			}

			boost::posix_time::time_duration diff = now - last_check;

			if (diff.total_milliseconds() > 3000) { /* report interval */
				double interval_sec = (double) diff.total_milliseconds() / 1000;
				double total_sec = (double) (now - start_time).total_milliseconds() / 1000;

				double mbps = (double) total_bytes / (1024 * 1024) / total_sec;
				double mrps = (double) total_records / total_sec;

				double lmbps = (double) (total_bytes - last_bytes) / (1024 * 1024)
					/ interval_sec;
				double lmrps = (double) (total_records - last_records)
					/ interval_sec;

				//  		E("", mbps, mrps);
				E("recent: %.2f MB/s %.2f records/s    avg: %.2f MB/s %.2f records/s",
						lmbps, lmrps, mbps, mrps);

				last_check = now;
				last_bytes = total_bytes;
				last_records = total_records;

				return true;
			}

			return false;
		}
#endif

#if 0 /* moved to UnboundedInMem */
		/* internal accounting */
		uint64_t total_bytes = 0, total_records = 0;
		/* last time we report */
		uint64_t last_bytes = 0, last_records = 0;
		ptime last_check, start_time;
		int once = 1;
#endif

		/* the eval pauses between issuing two consecutive "waves" of bundles.
		 * if 0, does not pause but just sched_yield() */

		const int pause_ms;

		/* starting point of the ts that we generate. note that we shoudl avoid making
		 * it epoch because the way we calculate start win, see AddAggregatedResult */
		//  ptime current_ts = Window::epoch;
		ptime current_ts;

		UnboundedInMemEvaluatorBase(int pause_ms = 0)
            // : pause_ms (pause_ms), current_ts(boost::gregorian::date(2018, Oct, 23)) { }
            : pause_ms (pause_ms), current_ts(boost::posix_time::second_clock::local_time()) {}
		void pause_between_waves() {
			if (this->pause_ms == 0)
				sched_yield();
			else {
				W(" ----------- to sleep %d ms (XXX shouldn't sleep?) ", this->pause_ms);
				usleep(this->pause_ms * 1000);
			}
		}

};

#endif // UNBOUNDED_INMEM_EVALBASE_H
