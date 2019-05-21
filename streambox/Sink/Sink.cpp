#define K2_NO_DEBUG 1

#include "Sink.h"
#include "RecordBitmapBundleSinkEvaluator.h"
#include "RecordBitmapBundleSinkEvaluator1.h"

#ifdef USE_TBB_DS
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#endif

#ifdef USE_FOLLY_HASHMAP
#include "folly/AtomicHashMap.h"
#endif

#include "Source/UnboundedInMemEvaluator.h"

#ifdef MEASURE_LATENCY
//static int first_1 = 1;
boost::posix_time::ptime start;
#endif

/* ----  specialization: different ways to dump bundle contents  --- */

template<>
void RecordBitmapBundleSink<pair<long, vector<long>>>::printBundle
	(RecordBitmapBundle<pair<long,vector<long>>> & input_bundle) {

	I("got one bundle");

	for (auto it = input_bundle.begin(); it != input_bundle.end(); ++it) {
		auto & k = (*it).data.first;
		auto & vec = (*it).data.second;

		cout << "K: " << k << endl;
		for (auto & v : vec) {
			cout << v << " ";
		}
		cout << " -------------- " << endl;
	}

	return;
}

template<>
void RecordBitmapBundleSink<creek::string>::printBundle
	(RecordBitmapBundle<creek::string> & input_bundle) {

  	static atomic<long> count (0);

    //W("---- got one bundle: %ld ---- ", count);

#ifdef MEASURE_LATENCY
#if 0
    for (auto it = input_bundle.begin(); it != input_bundle.end(); ++it) {
			std::string str((*it).data);
			std::string spec("http://www.hongyumiao.com");
			if(str.compare(spec) == 0){
				if(first_1){ //only find the first match
					std::cout << "++++++find: http://www.hongyumiao.com+++++" << std::endl;
					first_1 = 0;
					boost::posix_time::ptime end = boost::posix_time::second_clock::local_time();
					boost::posix_time::time_duration diff = end - start;
					double interval_sec = (double) diff.total_milliseconds() / 1000 ;
					std::cout << "start time: " << boost::posix_time::to_simple_string(start) << std::endl;
					std::cout << "End time: " << boost::posix_time::to_simple_string(end) << std::endl;
					std::cout << "Latency: " << interval_sec << std::endl;
				}
			}
    }
#endif

#if 0
    /* sample bundles to dump their propagation path */
    static long const skip = 1024;
    static mutex mtx;
    if ((count % skip) == 0) {
    	unique_lock<mutex> lck(mtx);
    	input_bundle.dump_markers();
    }
#endif
    
#endif

#ifdef DEBUG
//#if 0
    for (auto it = input_bundle.begin(); it != input_bundle.end(); ++it) {
    	stringstream ss; /* xzl: this can be expensive per vtune */
    	ss << (*it).data;
    	I("%s", ss.str().c_str());
//    	cout << (*it).data << endl;
    }
    I("--------------------------");
#endif

    count ++;
}

template<>
bool RecordBitmapBundleSink<creek::string>::report_progress
	(RecordBitmapBundle<creek::string> & input_bundle){
       
  	//total_bytes += bytes;
  	total_bytes += input_bundle.size();
  	//total_records += records;
  	total_records += 0; // XXX don't care about records now.

  	ptime now = boost::posix_time::microsec_clock::local_time();

  	if (once) {
  		once = 0;
  		last_check = now;
  		start_time = now;
  		//last_bytes = bytes;
		last_bytes = input_bundle.size();
  		//last_records = records;
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

        	//std::cout << "In RecordBitmapBundleSink" << std::endl; 
//  		EE("", mbps, mrps);
  		EE("--Output-- recent: %.2f MB/s %.2f records/s    avg: %.2f MB/s %.2f records/s",
  				lmbps, lmrps, mbps, mrps);

  		last_check = now;
  		last_bytes = total_bytes;
  		last_records = total_records;

  		return true;
  	}

  	return false;
}


/* -------------------- ExecEvaluator --------------------
 * out of line to avoid circular dependency
 */

template<class T>
void RecordBitmapBundleSink<T>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
/*
	RecordBitmapBundleSinkEvaluator<T> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
*/
	if(this->get_side_info() == SIDE_INFO_JD){
#ifdef WORKAROUND_JOIN_JDD
		RecordBitmapBundleSinkEvaluator<T> eval(nodeid); /* construct a normal eval */
#else
		RecordBitmapBundleSinkEvaluator1<T> eval(nodeid);
#endif

		//RecordBitmapBundleSinkEvaluator1<pair<long,vector<long>>> eval(nodeid);
		eval.evaluate(this, c, bundle_ptr);
	}else { /* this is a default type sink */
		RecordBitmapBundleSinkEvaluator<T> eval(nodeid);
		eval.evaluate(this, c, bundle_ptr);
	}

}

/* ---- instantiation for concrete types --- */

template
void RecordBitmapBundleSink<pair<long, vector<long>>>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void RecordBitmapBundleSink<creek::string>::ExecEvaluator(
		int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

