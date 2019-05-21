/*
 * extract words from string_range and pack them as records into bundles.
 */

#ifndef WC_MAPPER_H
#define WC_MAPPER_H

#include "Values.h"
#include "Mapper/Mapper.h"

using namespace std;

/* op to run (static as template argument) */
namespace wc_mapper{
	enum mode {
			WC, /* word split */
			LINE, /* line split */
	};
}

template <class InputT = string_range,
					class OutputT = pair<creek::string, long>,
//					class InputBundleT = RecordBitmapBundle<InputT>,
//					class OutputBundleT = RecordBitmapBundle<OutputT>,
					template<class> class BundleT = RecordBitmapBundle,
				  wc_mapper::mode mode = wc_mapper::WC /* which op to run */
				 >
class WordCountMapper : public Mapper<InputT> {

//  using OutputBundleT = RecordBitmapBundle<OutputT>;

	using InputBundleT = BundleT<InputT>;
	using OutputBundleT = BundleT<OutputT>;
	using TransformT = WordCountMapper<InputT,OutputT,BundleT,mode>;

private:
  static atomic<unsigned long> record_counter_;

public:

  WordCountMapper(string name = "wc_mapper") : Mapper<InputT>(name) { }

  /* Based on @mode, this is specialized into different ops
   *
   * note that we pass in @output_bundle, since one input Record (string_range)
   * may result in multiple output records (words)
   * @return: the # of records emitted
   */
  static uint64_t do_map(Record<InputT> const & in,
          shared_ptr<OutputBundleT> output_bundle);
#if 0
  uint64_t do_map(Record<InputT> const & in,
        shared_ptr<OutputBundleT> output_bundle) {
  	if (mode_ == WC)
  		return do_map_wc(in, output_bundle);
  	if (mode_ == LINE)
  		return do_map_line(in, output_bundle);
  	else {
  		bug("unknown mode");
  		return 0;
  	}
  }

  static uint64_t do_map_wc(Record<InputT> const & in,
      shared_ptr<OutputBundleT> output_bundle) {

    uint64_t i = 0, cnt = 0;

#if 0 /* this breaks our design that input bundles are const */
    /* rewrite the input string in place */
    for (uint64_t i = 0; i < in.data.len; i++)
      in.data.data[i] = toupper(in.data.data[i]);
#endif

    while (i < in.data.len) {
//      while (i < in.data.len && (in.data.data[i] < 'A' || in.data.data[i] > 'Z'))
      while (i < in.data.len && (!isalpha(in.data.data[i])))
        i++;
      uint64_t start = i;
//      while (i < in.data.len && ((in.data.data[i] >= 'A' && in.data.data[i] <= 'Z') || in.data.data[i] == '\''))
      while (i < in.data.len && (isalpha(in.data.data[i]) || in.data.data[i] == '\''))
        i++;
      if (i > start) {
          string str(in.data.data + start, i - start);
          std::transform(str.begin(), str.end(),str.begin(), ::toupper);

          output_bundle->add_record(Record<KVPair>(KVPair(str, 1), in.ts));
          cnt ++;
      }
    }

    return cnt;
  }

  static uint64_t do_map_line(Record<InputT> const & in,
      shared_ptr<OutputBundleT> output_bundle) {

    uint64_t i = 0, cnt = 0;

    while (i < in.data.len) {
      /* locate the next '\n' */
    	while (i < in.data.len && (in.data.data[i]) != '\n')
        i++;
      uint64_t start = i;
      /* locate the next next \n */
      while (i < in.data.len && (in.data.data[i]) != '\n')
        i++;
      if (i > start) {
				string str(in.data.data + start, i - start);
				output_bundle->add_record(Record<string>(str, in.ts));
				cnt ++;
      }
    }

    return cnt;
  }
#endif

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

#endif // WC_MAPPER_H
