#ifndef NYTQUERY_H
#define NYTQUERY_H

#include <re2/re2.h>
#include "Mapper/Mapper.h"
#include "NewYorkTaxi/NYTRecord.h"
#include "Values.h"

template <typename InputT,
          typename OutputT,
          // typename OutputT=pair<uint64_t, uint64_t>
          template<class> class BundleT_>

class NYTQuery : public Mapper<InputT> {
    using OutputBundleT = BundleT_<OutputT>;

public:
    NYTQuery(string name = "nyt_query") : Mapper<InputT>(name), record_counter_(0) { }

    const double DISTANCE = 5.0;
    // "Select #trips, avg distance where distance > 5 miles group by region (over the last 2 sec)"
    // bool do_map(Record<InputT> const & in) {
    bool do_map(Record<InputT> &in) {
        // in.data.print();
        double trip_distance = in.data.trip_distance;
        char* vendor_id = in.data.vendorID;
        if (trip_distance > DISTANCE && vendor_id[0] == 'V') {
            // std::cout << "trip_distance: " << trip_distance
            //           << ", vendor_id: " << vendor_id << std::endl;
            return true;
        } else {
            return false;
        }

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


#endif //NYTQUERY_H
