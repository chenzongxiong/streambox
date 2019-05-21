#ifndef LRB_PARSER_H
#define LRB_PARSER_H

#include "Values.h"
#include "Mapper/Mapper.h"
#include "LinearRoad/LinearRoadRecord.h"

using namespace std;

template <class InputT = string_range,
  class OutputT = LinearRoadRecord,
  template<class> class BundleT = RecordBundle
  >
  class LinearRoadParser : public Mapper<InputT> {

	using InputBundleT = BundleT<InputT>;
	using OutputBundleT = BundleT<OutputT>;
	using TransformT = LinearRoadParser<InputT,OutputT,BundleT>;

public:
  atomic<unsigned long> record_counter_;
    atomic<unsigned long> byte_counter_;
public:
    LinearRoadParser(string name = "linear_road_parser") : Mapper<InputT>(name), record_counter_(0), byte_counter_(0) { }

  Record<OutputT> do_map(Record<InputT> const & in) {
    return Record<OutputT>((OutputT)in.data, in.ts);
  }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                     shared_ptr<BundleBase> bundle = nullptr) override;

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

  	return true;
  }

  };

#endif // LRB_PARSER_H
