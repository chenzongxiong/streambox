#ifndef UNBOUNDED_TWEETS_H_
#define UNBOUNDED_TWEETS_H_

template<template<class> class BundleT>   /* can't use default argument = RecordBitmapBundle<string_range> */
class UnboundedInMem<string_range, BundleT> : public PTransform {

	using T = string_range;

public:
  const char * input_fname;

  /* these are event time, irrelevant to the engine processing timing */
  const int punc_interval_ms = 1000;
  const unsigned long records_per_interval;  /* # records between two puncs */
  uint64_t string_len; /* the length covered by each record */
  const int session_gap_ms; /* gap between "bursts" of bundles. for testing session windows */
  const uint64_t target_tput; /* in record/sec */

  /* array of prefilled buffers of records, one for each NUMA node */
//  vector<vector<Record<T>>> records;
  vector<Record<T> *> record_buffers;
  uint64_t buffer_size_records = 0;

private:
  vector<char const *> buffers; /* buffers, one for each NUMA node */
  uint64_t buffer_size = 0;

public:
  UnboundedInMem (string name, const char *input_fname,
  		unsigned long rpi, uint64_t tt,
  		uint64_t record_size, int session_gap_ms)
    : PTransform(name), input_fname(input_fname),
      records_per_interval(rpi),
      string_len(record_size),
      session_gap_ms(session_gap_ms),
      target_tput(tt),
      byte_counter_(0), record_counter_(0)
  {
    int fd;
    struct stat finfo;
    int num_nodes = numa_num_configured_nodes();

    CHECK_ERROR((fd = open(input_fname, O_RDONLY)) < 0);
    CHECK_ERROR(fstat(fd, &finfo) < 0);

//    NumaThreadPool& pool = NumaThreadPool::instance();
    buffer_size = records_per_interval * string_len * 2;

    /* sanity check: file long enough for the buffer? */
    if ((int64_t)buffer_size > finfo.st_size) {
    	EE("input data not enough. need %.2f MB. has %.2f MB",
    			(double) buffer_size / 1024 / 1024, (double) finfo.st_size / 1024 / 1024);
    	abort();
    }

    buffer_size_records = buffer_size / string_len;
    xzl_assert(buffer_size_records > 0);

    printf("---- source configuration ---- \n");
    printf("source file: %s\n", input_fname);
    printf("source file size: %.2f MB\n", (double)finfo.st_size/1024/1024);
    printf("buffer size: %.2f MB\n", (double)buffer_size/1024/1024);
    printf("%10s %10s %10s %10s %10s %10s %10s %10s\n",
    		"#nodes:", "KRec", "MB", "epoch/ms", "KRec/epoch", "MB/epoc",
    		"target:KRec/S", "RecSize" );
    printf("%10d %10lu %10lu %10d %10lu %10lu %10lu %10lu\n",
    		num_nodes, buffer_size_records/1000, buffer_size/1024/1024,
    		punc_interval_ms, records_per_interval/1000,
    		records_per_interval * string_len /1024/1024,
    		target_tput/1000, string_len);

    /* get per-node buffers and fill the buffers with the file contents */
    for (int i = 0; i < num_nodes; i++) {

    	printf("node %d: loading buffer %.2f MB\n", i, (double)buffer_size/1024/1024);

    	boost::progress_display show_progress(buffer_size);

    	char *p = (char *)numa_alloc_onnode(buffer_size, i);
			assert(p);

			uint64_t r = 0;
			while (r < buffer_size) {
				auto ret = pread(fd, p + r, buffer_size, r);
				if (ret == 0 || ret == -1) {
					perror("read failure");
					abort();
				} else {
					r += ret;
					show_progress += ret;
				}
			}
			assert(r == buffer_size);

			buffers.push_back(p);
    }

    /* fill the buffers of records */

    for (int i = 0; i < num_nodes; i++) {
    	Record<T> * record_buffer =
    			(Record<T> *) numa_alloc_onnode(sizeof(Record<T>) * buffer_size_records,
    					i);

    	xzl_assert(record_buffer);

    	for (unsigned int j = 0; j < buffer_size_records; j++) {
    		record_buffer[j].data.data = buffers[i] + j * string_len;
    		record_buffer[j].data.len = string_len;
    		/* record_buffer[j].ts will be filled by eval */
    	}

    	record_buffers.push_back(record_buffer);
    }

//    EE("filled %lu buffers, each with %lu records,"
//    		"each covering size %lu MB from file %s.",
//        buffers.size(), buffer_size_records, buffer_size/1024/1024, input_fname);
  }

  // source, no-op
  virtual ptime RefreshWatermark(ptime wm) override {
    return wm;
  }

  /* internal accounting  -- to be updated by the evaluator*/
  atomic<unsigned long> byte_counter_, record_counter_;

  bool ReportStatistics(PTransform::Statstics* stat) override {

    /* internal accounting */
     unsigned long total_records =
    		 record_counter_.load(std::memory_order_relaxed);
     unsigned long total_bytes =
    		 byte_counter_.load(std::memory_order_relaxed);

     /* last time we report */
     static unsigned long last_bytes = 0, last_records = 0;
     static ptime last_check, start_time;
     static int once = 1;

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

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase>) override;

//  void ExecEvaluator(int nodeid, EvaluationBundleContext *c) override {
//  	/* instantiate an evaluator */
//  	UnboundedInMemEvaluator<string_range> eval(nodeid);
//		eval.evaluate(this, c);
//  }

};


#endif /* UNBOUNDED_TWEETS_H_ */
