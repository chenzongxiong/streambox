#ifndef YAHOOPARSER_H
#define YAHOOPARSER_H
#include <re2/re2.h>
#include "Mapper/Mapper.h"
#include "Yahoo/YahooRecord.h"
#include "Values.h"


template <typename InputT = string_range, typename OutputT = YahooRecord, template<class> class BundleT_ = RecordBundle>
class YahooParser : public Mapper<InputT> {
	using OutputBundleT = BundleT_<OutputT>;
public:
YahooParser(string name = "yahoo_parser") : Mapper<InputT>(name), record_counter_(0), byte_counter_(0) { }

	Record<OutputT> do_map(Record<InputT> const & in) {
        // deep copy
        ptime ts = boost::posix_time::microsec_clock::local_time(); //  has negative impact on throughput
        return Record<OutputT>((OutputT)in.data, ts);
	}

  	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
    	shared_ptr<BundleBase> bundle_ptr = nullptr) override;

    atomic<unsigned long> record_counter_;
    atomic<unsigned long> byte_counter_;

    bool ReportStatistics(PTransform::Statstics* stat) override {
        /* internal accounting */
        static unsigned long total_records = 0, total_bytes = 0;
        /* last time we report */
        static unsigned long last_bytes = 0, last_records = 0;
        static ptime last_check, start_time;
        static int once = 1;

        /* only care about records */
        total_records = this->record_counter_.load(std::memory_order_relaxed);
        total_bytes = this->byte_counter_.load(std::memory_order_relaxed);

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

            last_check = now;
            last_bytes = total_bytes;
            last_records = total_records;
        }

//  	E("bundle counter %lu", this->bundle_counter_.load());

        return true;
    }
};

#endif /* YAHOOPARSER_H */
