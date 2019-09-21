#ifndef NEXMARK_FILTER_H
#define NEXMARK_FILTER_H

#include "Mapper/Mapper.h"
#include "Nexmark/NexmarkRecord.hpp"
#include "Values.h"

template <typename InputT=NexmarkRecord,
          typename OutputT=NexmarkRecord,
          template<class> class BundleT_=RecordBundle
          >
class NexmarkFilter : public Mapper<InputT> {
    using OutputBundleT = BundleT_<OutputT>;

public:
    NexmarkFilter(string name="nexmark_filter") : Mapper<InputT>(name), record_counter_(0) {}

    bool do_map(Record<InputT> const & in) {
        // pick up a event from {"view", "click", "purchase"}
        // if the event is the same as the given event, return true
        // else, return false.
        // return true, means the input record will flow into the next
        // pipeline.
        // const char* event_type = in.data.event_type;
        // if (strcmp(event_type, "view") == 0) {
        //     // std::cout << "true: event_type: " << event_type << std::endl;
        //     return true;
        // } else {
        //     // std::cout << "false: event_type: " << event_type << std::endl;
        //     return false;
        // }

        uint64_t auction = in.data.auction;
        uint64_t price = in.data.price;
        // std::cout << "auction: " << auction << ", "
        //           << "price: " << price << std::endl;

        // auto tuple = in.data;
        // tuple.print();
        return true;

    }

    void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                       shared_ptr<BundleBase> bundle_ptr = nullptr) override;

    atomic<unsigned long> record_counter_;

    bool ReportStatistics(PTransform::Statstics* stat) {
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

            last_check = now;
            last_bytes = total_bytes;
            last_records = total_records;
        }
        return true;
    }
};


#endif //NEXMARK_FILTER_H
