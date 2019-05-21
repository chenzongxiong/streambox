#ifndef LRB_TOLL_NOTIFICATION_H
#define LRB_TOLL_NOTIFICATION_H

#include <string>
#include <map>

#include "Values.h"
#include "Mapper/Mapper.h"
#include "LinearRoad/LinearRoadRecord.h"
#include "LinearRoad/LinearRoadTollRecord.h"


template <class InputT = LinearRoadRecord,
          class OutputT = LinearRoadTollRecord,
          template<class> class BundleT_ = RecordBundle
          >

class LinearRoadTollNotification : public Mapper<InputT> {

    using OutputBundleT = BundleT_<OutputT>;

public:
LinearRoadTollNotification(string name="linear_road_toll") : Mapper<InputT>(name), record_counter_(0), time_of_last_toll(0) {
    acc_dec = new AccidentDetection();
    avg_seg_speed_map.reserve(SEGMENT_NUM);
  }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                     shared_ptr<BundleBase> bundle = nullptr) override;

  atomic<unsigned long> record_counter_;

  size_t time_of_last_toll;
  std::map<int, StopMap> stop_map; /* vid, StopMap */
  // std::map<int, AvgSpeedMap> avg_seg_speed_map;  /* segment id, AvgSpeedMap */
  std::vector<AvgSpeedMap> avg_seg_speed_map;
  AccidentDetection* acc_dec;


  bool ReportStatistics(PTransform::Statstics* stat) override {
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

#endif  // LRB_TOLL_NOTIFICATION_H
