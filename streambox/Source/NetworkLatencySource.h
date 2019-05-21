#ifndef NETWORK_MONITOR_SOURCE_H
#define NETWORK_MONITOR_SOURCE_H 
#include "Source/Unbounded.h"
#include "Values.h" //for struct srcdst_rtt

#include <arpa/inet.h>

//for network mointor app 
//template<>
class NetworkLatencySource: public PTransform{
public:
	int interval_ms;  /* the time difference between consecutive bundles */
	const char * input_fname;
	const int punc_interval_ms = 1000;
	const unsigned long records_per_interval; // # records between two puncs
	int record_len = 0; /* the length covered by each record */
	const int session_gap_ms; // for tesging session windows
	const uint64_t target_tput; // in record/sec
	
  	//vector<Record<long> *> record_buffers;
//  	vector<Record<struct srcdst_rtt> *> record_buffers;
	vector<Record<pair<creek::ippair, long>>*> record_buffers;

	uint64_t buffer_size_records = 0;

  	vector<struct srcdst_rtt *> buffers;
  	uint64_t buffer_size = 0;
  	long record_num = 0;
	
	/* internal accounting  -- to be updated by the evaluator*/
  	atomic<unsigned long> byte_counter_, record_counter_;
	
	//XXX only use one node for test
	//int num_nodes = numa_num_configured_nodes();
	int num_nodes = 1;
	NetworkLatencySource (string name, const char *input_fname,
				unsigned long rpi, uint64_t tt,
				uint64_t record_size, int session_gap_ms)
  		: PTransform(name), input_fname(input_fname),
				records_per_interval(rpi),
				record_len(record_size), //= sizeof(strut srcdst_rtt)
				session_gap_ms(session_gap_ms),
				target_tput(tt),
				byte_counter_(0),
				record_counter_(0){

		struct stat finfo;

		FILE *file = fopen(input_fname, "rb");
		if(!file){
			assert(false && "open file failed !!!!");
		}
		CHECK_ERROR(fstat(fileno(file), &finfo) < 0);
		
		if (record_len != sizeof(struct srcdst_rtt)) {
			xzl_bug("record_len must equal size of srcdst_rtt. abort.");
		}

		buffer_size = records_per_interval * record_len * 4;

    /* sanity check: file long enough for the buffer? */
    if ((int64_t)buffer_size > finfo.st_size) {
    	EE("input data not enough. need %.2f MB. has %.2f MB",
    			(double) buffer_size / 1024 / 1024, (double) finfo.st_size / 1024 / 1024);
    	abort();
    }

		struct srcdst_rtt sdr;

#if 0 /* xzl: why do we do this?! */
		// get record_num
		while(fread(&sdr, sizeof(struct srcdst_rtt), 1, file)){
			record_num ++;
			memset(&sdr, 0, sizeof(struct srcdst_rtt));
		}
#endif

		record_num = buffer_size / record_len;
		xzl_assert(record_num > 0);

//		printf("total size %.2f MB\n", (double) record_num * sizeof(struct srcdst_rtt) /1024/1024);
//		std::cout << "record_num is " << record_num << std::endl;

    printf("---- source configuration ---- \n");
    printf("source file: %s\n", input_fname);
    printf("source file size: %.2f MB\n", (double)finfo.st_size/1024/1024);
    printf("buffer size: %.2f MB\n", (double)buffer_size/1024/1024);
    printf("%10s %10s %10s %10s %10s %10s %10s %10s\n",
    		"#nodes:", "KRec", "MB", "epoch/ms", "KRec/epoch", "MB/epoc",
    		"target:KRec/S", "RecSize" );
    printf("%10d %10lu %10lu %10d %10lu %10lu %10lu %10d\n",
    		num_nodes, record_num/1000, buffer_size/1024/1024,
    		punc_interval_ms, records_per_interval/1000,
    		records_per_interval * record_len /1024/1024,
    		target_tput/1000, record_len);

		/*
		struct stat finfo;
		stat(input_fname, &finfo);
		std::cout << "size2 is " << finfo.st_size << std::endl;
		buffer_size = finfo.st_size;
		*/
		//fill buffer
		for(int i = 0; i < num_nodes; i++){
			struct srcdst_rtt *p = (struct srcdst_rtt *)numa_alloc_onnode(buffer_size, i);
			assert(p);
			fseek(file, 0, SEEK_SET); // seek back to beginning of file
			long j = 0;
			while(fread(&sdr, sizeof(struct srcdst_rtt), 1, file)){
				//p[j++] = sdr;
				//strncpy(sdr.sd_ip , sd_ip, 32);
				strncpy(p[j].sd_ip, sdr.sd_ip, 32);
				p[j].rtt = sdr.rtt; 
				//std::cout << p[j].sd_ip << std::endl;	
				j++;
				if (j == record_num)  /* don't read until file end */
					break;
				/*
				strncpy(p->sd_ip, sdr.sd_ip, 32);
				p->rtt = sdr.rtt;
				std::cout << p->sd_ip << std::endl;
				p++;
				*/
				memset(&sdr, 0, sizeof(struct srcdst_rtt));
			}
			buffers.push_back(p);
		}
		
		//file the buffers of records
		for(int i = 0; i < num_nodes; i++){
			Record<pair<creek::ippair, long>> * record_buffer =
				(Record<pair<creek::ippair, long>> *) numa_alloc_onnode(sizeof(Record<pair<creek::ippair, long>>) * record_num, i);
			assert(record_buffer);
			for(long j = 0; j < record_num; j++){
#if 0 /* hym's */
				//record_buffer[j].data = buffers[i][j];
				//std::cout << buffers[i][j].sd_ip << std::endl;
				strncpy(record_buffer[j].data.sd_ip, buffers[i][j].sd_ip, 32);
				record_buffer[j].data.rtt = buffers[i][j].rtt;
				//std::cout << record_buffer[j].data.sd_ip << std::endl; //Debugging
				//record_buffer[j].data.len = record_len; //sizeof(long)
				//record_buffer[j].ts will be filled by eval
#endif
				record_buffer[j].data = make_pair(convert_ip_pair(buffers[i][j].sd_ip), buffers[i][j].rtt);
			}
			record_buffers.push_back(record_buffer);
		}//end for
	
		fclose(file);

	}//end init 

	/* @str max 32 bytes. in the following format:
	 * 100.65.229.132-100.79.74.130
	 * 3x8(digits)+3x2(dots)+1(dash)=31bytes
	 */
#define SLEN 32
	uint64_t convert_ip_pair(const char * str) {
		char s[SLEN], *s1 = nullptr;

		memcpy(s, str, SLEN);

		if (s[SLEN-1] != '\0') {
			xzl_bug("not zero terminated");
		}

		/* split s */
		for (int i = 0; i < SLEN; i++) {
			if (s[i] == '-') {
				s[i] = '\0';
				s1 = s + i + 1;
			}
		}
		if (!s1)
			xzl_bug("bag string");

		/* now s and s1 point to two zero-term ip strings */
		uint64_t res;
		int r;
		r = inet_pton(AF_INET, s, &res);  /* 1st int */
		if (r != 1)
			xzl_bug("illegal ip addr");

		r = inet_pton(AF_INET, s1, ((char *)&res) + 4); /* 2nd int */
		if (r != 1)
			xzl_bug("illegal ip addr");

		return res;
	}

	// source, no-op. note that we shouldn't return the transform's wm
	virtual ptime RefreshWatermark(ptime wm) override {
		return wm;
	}

	bool ReportStatistics(PTransform::Statstics* stat) override {
		//TODO
		//return false;
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

}; //end UnboundedInMem
#endif /* NETWORK_MONITOR_SOURCE_H */
