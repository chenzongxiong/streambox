#ifndef UNBOUNDED_JOIN_H
#define UNBOUNDED_JOIN_H 
#include "Unbounded.h"
#include "Values.h"

//hym: new source for join
//     get ride of SimpleMapper
//     Source deposit bundles to Join directly 
#if 0
template<class T,  								/* record data type */
	template<class> class BundleT		/* bundle type */
	>
	class UnboundedInMem_Join; /* generic impl not available */


	template<template<class> class BundleT>
	class UnboundedInMem_Join<long, BundleT> : public PTransform {
		//class UnboundedInMem_Join : public PTransform {
		using T = long;
		using OutputBundleT = RecordBitmapBundle<T>;
#endif

		class UnboundedInMem_Join: public PTransform{
			using OutputBundleT = RecordBitmapBundle<pair<long, long>>;
			using T = long;
			public:
			int interval_ms;  /* the time difference between consecutive bundles */
			//const int session_gap_ms; /* gap between "bursts" of bundles. for testing session windows */

			/*new variables*/
			const char * input_fname;
			const int punc_interval_ms = 1000;
			const unsigned long records_per_interval; // # records between two puncs
			int record_len = 0; /* the length covered by each record */
			const int session_gap_ms; // for tesging session windows
			const uint64_t target_tput; // in record/sec
			//uint64_t string_len; /* the length covered by each record */

			//hym: here we don't read ???
			//array of prefilled buffers of record, one for each NUMA node
			//vector<Record<T> *> record_buffers;
			vector<Record<long> *> record_buffers;
			uint64_t buffer_size_records = 0;

			//vector<long const *> buffers; /*buffers, one for each NUMA node*/
			vector<long *> buffers;
			uint64_t buffer_size = 0;
			unsigned long record_num = 0;
#if 0
			UnboundedInMem (string name, int interval_ms = 100, int session_gap_ms = 0)
				: PTransform(name), interval_ms(interval_ms),
				session_gap_ms(session_gap_ms) { }
#endif
			UnboundedInMem_Join (string name, const char *input_fname,
					unsigned long rpi, uint64_t tt,
					uint64_t record_size, int session_gap_ms)
				: PTransform(name), input_fname(input_fname),
				records_per_interval(rpi),
				record_len(record_size), //= sizeof(long)
				session_gap_ms(session_gap_ms),
				target_tput(tt),
				byte_counter_(0),
				record_counter_(0){
					//TODO
					interval_ms = 50;
					//int fd;
					struct stat finfo;
					//XXX only use one node for test
					//	int num_nodes = numa_num_configured_nodes();
					int num_nodes = 1;

					std::cout << "WARNING: num_nodes is set to 1!!!!! Reset it in Unbounded.h and UnboundedInMemEvaluator.h!!!" << std::endl;

					//scan the file to see how many record(intige) the file has
					//int record_num = 0;
					//long i = 0;
					FILE * file = fopen(input_fname, "r");
					if(!file){
						assert(false && "open file failed!!!");		
					}

					std::cout << "WARNING: record_num is 282090931. Should set the value according to input file!!!!" << std::endl;
					record_num = 282090931;
					//XXX comment this temporarilly!!! Remember to restore this
#if 0
					// get # of records
					// http://stackoverflow.com/questions/4600797/read-int-values-from-a-text-file-in-c
					fscanf(file, "%ld", &i);
					while(!feof(file)){
						record_num ++;
						//std::cout << i << " ";
						fscanf(file, "%ld", &i);
					}
#endif
					// get file size
					// http://stackoverflow.com/questions/238603/how-can-i-get-a-files-size-in-c
					unsigned long fsize;
#if 0   // method 1
					fseek(file, 0, SEEK_END);// seek to end of file
					fsize = ftell(file); // get current file pointer
					fseek(file, 0, SEEK_SET); // seek back to beginning of file
#endif
					// method 2
					stat(input_fname, &finfo);
					fsize = finfo.st_size;
					fsize = fsize; //fix warning

					//	buffer_size = record_num * record_len; // record_num * sizeof(long)
					buffer_size = records_per_interval * record_len * 2 * 10;

					/* sanity check: file long enough for the buffer? */
					if ((int64_t)buffer_size > finfo.st_size) {
						EE("input data not enough. need %lu KB. has %lu KB",
								buffer_size / 1024, finfo.st_size / 1024);
						abort();
					}
					record_num = buffer_size / record_len;

					// print source config info
					printf("---- source configuration ---- \n");
					printf("source file: %s\n", input_fname);
					printf("source file size: %.2f MB\n", (double)finfo.st_size/1024/1024);
					printf("buffer size: %.2f MB\n", (double)buffer_size/1024/1024);
					printf("Number of Records: %ld \n", record_num);

					// allocate and fill buffers: should replace genRandomElement() function
					// XXX shall we allocate a buffer on each node??
					// This may be not a good idea. It's better that two input streams on the same NUMA node
					for(int i = 0; i < num_nodes; i++){
						long *p = (long *)numa_alloc_onnode(buffer_size, i);
						assert(p);
						fseek(file, 0, SEEK_SET); // seek back to beginning of file

						long rcd;
						unsigned long j = 0;
						int ret = fscanf(file, "%ld", &rcd);
						while(!feof(file)){
							if (j == record_num)
								break;
							p[j++] = rcd;	
							ret = fscanf(file, "%ld", &rcd);
						}
						ret = ret; //fix warning
						std::cout << "add " << j << " records(long) to a buffer" << std::endl;	
						buffers.push_back(p);
					}

					//fill the buffers of records
					for(int i = 0; i < num_nodes; i++){
						Record<T> * record_buffer =
							(Record<T> *) numa_alloc_onnode(sizeof(Record<T>) * record_num,i);
						assert(record_buffer);
						for(unsigned long j = 0; j < record_num; j++){
							//record_buffer[j].data.data = buffers[i][j];
							record_buffer[j].data = buffers[i][j];
							//std::cout << "record_buffer[j.data] is " << record_buffer[j].data << std::endl;
							//record_buffer[j].data.len = record_len; //sizeof(long)
							//record_buffer[j].ts will be filled by eval
						}
						record_buffers.push_back(record_buffer);
					}

					fclose(file);

				}//end UnboundedInMem init

			/* lockfree
			 *
			 * @ts: (hacking): the ts for each record in the bundle. to be removed when
			 * the source can generate ts.
			 * return: # records actually filled in
			 * */
			uint64_t FillBundle(int out_id, OutputBundleT & bundle,
					uint64_t capacity, ptime ts) {
				/*
				   uint64_t i;
				   for (i = 0; i < capacity; i++) {
				   bundle.add_record(Record<T>(genRandomElement(out_id), ts));
				   }
				   return i;
				 */
				return 0;
			}

			// source, no-op. note that we shouldn't return the transform's wm
			virtual ptime RefreshWatermark(ptime wm) override {
				return wm;
			}


			/* internal accounting  -- to be updated by the evaluator*/
			atomic<unsigned long> byte_counter_, record_counter_;
			bool ReportStatistics(PTransform::Statstics* stat) override {
				//TODO
				/*
				   std::cout << __FILE__ << __LINE__ << "  TODO: ReportStatistics is needed!!!!!!!!!!!!!!" << std::endl;
				//assert(false && "todo...");
				return false;
				 */
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

			private:

			T genRandomElement(int output_id) {

				assert(output_id >= 0 && output_id < 2);

				static atomic<T> cnt0 (0);
				static atomic<T> cnt1 (0);

				if (output_id == 0)
					return cnt0.fetch_add(1);
				else
					return cnt1.fetch_add(1);
  }

};


#endif /* UNBOUNDED_JOIN_H */
