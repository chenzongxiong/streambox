#ifndef _SINK_H
#define _SINK_H

#include <iostream>

#include "config.h"

extern "C" {
#include "log.h"
#include "measure.h"
}

#include "core/Transforms.h"

using namespace std;

///////////////////////////////////////////////////////

/* A bunch of records. This is good for that each window outputs a single value
 * @T: element type, e.g. long */
template<class T>
class RecordBundleSink : public PTransform {
	using InputBundleT = RecordBundle<T>;

	public:
	RecordBundleSink(string name) : PTransform(name) { }

	static void printBundle(const InputBundleT & input_bundle) {
		I("got one bundle");
	}

	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
			shared_ptr<BundleBase> bundle_ptr) override;

};

///////////////////////////////////////////////////////

/* @T: element type, e.g. long */
template<class T>
class RecordBitmapBundleSink : public PTransform {
	using InputBundleT = RecordBitmapBundle<T>;

	public:
	RecordBitmapBundleSink(string name) : PTransform(name) { }

	/* default behavior.
	 * can't do const. see below */
	static void printBundle(InputBundleT & input_bundle) {
		I("got one bundle");
	}

	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
			shared_ptr<BundleBase> bundle_ptr) override;

	static bool report_progress(InputBundleT & input_bundle){ return false;} //hym: defined in Sink.cpp

	private:
	//hym: for measuring throughput(report_progress())
	static uint64_t total_bytes;
	static uint64_t total_records;
	static uint64_t last_bytes;
	static uint64_t last_records;
	static ptime last_check;
	static ptime start_time;
	static int once;
#if 0
	/* internal accounting */
	uint64_t total_bytes = 0, total_records = 0;
	/* last time we report */
	uint64_t last_bytes = 0, last_records = 0;
	ptime last_check, start_time;
	int once = 1;
#endif
};
template<class T>
uint64_t RecordBitmapBundleSink<T>::total_bytes = 0;
template<class T>
uint64_t RecordBitmapBundleSink<T>::total_records = 0;
template<class T>
uint64_t RecordBitmapBundleSink<T>::last_bytes = 0;
template<class T>
uint64_t RecordBitmapBundleSink<T>::last_records = 0;

template<class T>
ptime RecordBitmapBundleSink<T>::last_check = boost::posix_time::microsec_clock::local_time();
template<class T>
ptime RecordBitmapBundleSink<T>::start_time = boost::posix_time::microsec_clock::local_time();;

template<class T>
int RecordBitmapBundleSink<T>::once = 1;
///////////////////////////////////////////////////////

template<class T>
class WindowsBundleSink : public PTransform {
	using InputBundleT = WindowsBundle<T>;

	public:
WindowsBundleSink(string name)  : PTransform(name), record_counter_(0) { }

	static void printBundle(const InputBundleT & input_bundle) {
		I("got one bundle");
	}

	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
			shared_ptr<BundleBase> bundle_ptr) override;

    atomic<unsigned long> record_counter_;
// #ifndef LSDS_STREAMBOX
    bool ReportStatistics(Statstics* stat) override {
        /* internal accounting */
        static unsigned long total_records = 0, total_bytes = 0;
        /* last time we report */
        static unsigned long last_bytes = 0, last_records = 0;
        static ptime last_check, start_time;
        static int once = 1;

        /* only care about records */
        total_records = this->record_counter_.load(std::memory_order_relaxed);

        ptime now = boost::posix_time::microsec_clock::local_time();

        if (once) {
            once = 0;
            last_check = now;
            start_time = now;
            last_records = total_records;
            return false;
        }

        boost::posix_time::time_duration diff = now - last_check;

        {
            double interval_sec = (double) diff.total_milliseconds() / 1000;
            double total_sec = (double) (now - start_time).total_milliseconds() / 1000;

            stat->name = this->name.c_str();
            stat->mbps = (double) total_bytes / total_sec;
            stat->mrps = (double) total_records / total_sec;

            stat->lmbps = (double) (total_bytes - last_bytes) / interval_sec;
            stat->lmrps = (double) (total_records - last_records) / interval_sec;

#if 0
            EE("recent: %.2f MB/s %.2f K records/s    avg: %.2f MB/s %.2f K records/s",
              stat->lmbps, stat->lmrps/1000, stat->mbps, stat->mrps/1000);
#endif

            last_check = now;
            last_bytes = total_bytes;
            last_records = total_records;
        }
        return true;
    }
// #endif
};

#endif // _SINK_H
