#ifndef NEXMARKAGGREGATION_HPP
#define NEXMARKAGGREGATION_HPP

#include "Values.h"
#include "Mapper/Mapper.h"

using namespace std;

template <class InputT,
          class OutputT,
          template<class> class BundleT
          >
class NexmarkAggregation : public Mapper<InputT> {

    using InputBundleT = BundleT<InputT>;
    using OutputBundleT = BundleT<OutputT>;
    using TransformT = NexmarkAggregation<InputT, OutputT, BundleT>;

public:

    NexmarkAggregation(string name="nexmark_mapper") : Mapper<InputT>(name), record_counter_(0) {}

    /* Based on @mode, this is specialized into different ops
     *
     * note that we pass in @output_bundle, since one input Record (string_range)
     * may result in multiple output records (words)
     * @return: the # of records emitted
     */
    uint64_t do_map(Record<InputT> const & in,
                    shared_ptr<OutputBundleT> output_bundle);

    void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                       shared_ptr<BundleBase> bundle = nullptr) override;

    atomic<unsigned long> record_counter_;
    bool ReportStatistics(PTransform::Statstics* stat) override {
        /* internal accounting */
        static unsigned long total_records = 0, total_bytes = 0;
        /* last time we report */
        static unsigned long last_bytes = 0, last_records = 0;
        static ptime last_check, start_time;
        static int once = 1;

        /* only care about records */
        total_records = TransformT::record_counter_.load(std::memory_order_relaxed);

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
            E("recent: %.2f MB/s %.2f K records/s    avg: %.2f MB/s %.2f K records/s",
              lmbps, lmrps/1000, mbps, mrps/1000);
#endif
            last_check = now;
            last_bytes = total_bytes;
            last_records = total_records;
        }
//  	E("bundle counter %lu", this->bundle_counter_.load());
        return true;
    }
};


#endif //NEXMARKAGGREGATION_HPP
