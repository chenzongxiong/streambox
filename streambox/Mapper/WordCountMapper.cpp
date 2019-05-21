#define K2_NO_DEBUG 1

#include "Values.h"
#include "WordCountMapperEvaluator.h"
#include "WordCountMapper.h"

//template <class InputT, class OutputT, class InputBundleT, class OutputBundleT,
template <class InputT, class OutputT, template<class> class BundleT,
					wc_mapper::mode mode>
//void WordCountMapper<InputT, OutputT, InputBundleT, OutputBundleT, mode>::ExecEvaluator(int nodeid,
void WordCountMapper<InputT, OutputT, BundleT, mode>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
//	using InputBundleT = BundleT<InputT>;
//	using OutputBundleT = BundleT<OutputT>;

#ifndef NDEBUG // if evaluator get stuck ..
	static atomic<int> outstanding (0);
#endif

	/* instantiate an evaluator */
	WordCountMapperEvaluator<InputT, OutputT, BundleT, mode> eval(nodeid);

#ifndef NDEBUG
	outstanding ++;
#endif

	eval.evaluate(this, c, bundle_ptr);

#ifndef NDEBUG
	outstanding --;
	int i = outstanding;
	(void)i;
	I("end eval... outstanding = %d", i);
#endif
}

/* ------- template specialization -------
 *
 * NB: we only impl a subset of InputT/OutputT/mode combinations
 *
 * Method cannot do partial template specialization.
 *
 * (tedious...can be more concise?)
 */

//using InputT = string_range;
//using OutputT = pair<string, long>;
//using KVPair = pair<string, long>;
//using OutputBundleT = RecordBitmapBundle<OutputT>;

//template <class InputBundleT, class OutputBundleT>
//uint64_t WordCountMapper<string_range, pair<string, long>, InputBundleT, OutputBundleT, wc_mapper::WC>
//	::do_map(Record<string_range> const & in, shared_ptr<OutputBundleT> output_bundle)

template <>
uint64_t WordCountMapper<string_range, pair<creek::string, long>,
//												 RecordBundle<string_range>,
//												 RecordBundle<pair<string,long>>,
													RecordBundle,
												 wc_mapper::WC>
			::do_map(Record<string_range> const & in, shared_ptr<RecordBundle<pair<creek::string,long>>> output_bundle)
{

	using KVPair = pair<creek::string, long>;

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
        creek::string str(in.data.data + start, i - start);
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);

//        output_bundle->add_record(Record<KVPair>(KVPair(str, 1), in.ts));
        output_bundle->emplace_record(KVPair(str, 1), in.ts);
        cnt ++;
    }
  }

  record_counter_.fetch_add(cnt, std::memory_order_relaxed);
  return cnt;
}

//using InputT = string_range;
//using OutputT = string;
//using OutputBundleT = RecordBitmapBundle<OutputT>;
template<>
uint64_t WordCountMapper<string_range, creek::string,
//												RecordBundle<string_range>, RecordBundle<string>,
													RecordBundle,
													wc_mapper::LINE>
				::do_map(Record<string_range> const & in,
							   shared_ptr<RecordBundle<creek::string>> output_bundle)
{
  uint64_t i = 0, cnt = 0, start;

//  E("start map, len %lu", in.data.len);

  while (i < in.data.len) {
    /* locate the next '\n' */
  	while (i < in.data.len && (in.data.data[i] != '\n'))
      i++;
    i ++;
    start = i;
    /* locate the next next \n */
    while (i < in.data.len && (in.data.data[i] != '\n'))
      i++;
//    E("here i=%lu", i);
    if (i > start && i < in.data.len) {
			creek::string str(in.data.data + start, i - start);
			output_bundle->add_record(Record<creek::string>(str, in.ts));
			cnt ++;
//			cout << "record: " << str << endl;
    }
  }

  record_counter_.fetch_add(cnt, std::memory_order_relaxed);

  return cnt;
}


#if 0
template<typename InputBundleT, typename OutputBundleT>
atomic<unsigned long>
	WordCountMapper<string_range, pair<string, long>,
								 InputBundleT, OutputBundleT, wc_mapper::WC>::record_counter_(0);

template<typename InputBundleT, typename OutputBundleT>
atomic<unsigned long>
WordCountMapper<string_range, string,
								InputBundleT, OutputBundleT, wc_mapper::LINE>::record_counter_(0);
#endif

/* no template specilization */
template <class InputT,
					class OutputT,
//					class InputBundleT,
//					class OutputBundleT,
					template<class> class BundleT,
				  wc_mapper::mode Mode>
atomic<unsigned long>
//WordCountMapper<InputT, OutputT, InputBundleT, OutputBundleT, Mode>::record_counter_(0);
	WordCountMapper<InputT, OutputT, BundleT, Mode>::record_counter_(0);

/* -------instantiation concrete classes------- */

/* using record bundle for input/output */
template
void WordCountMapper<string_range, creek::string,
											RecordBundle, wc_mapper::LINE>::ExecEvaluator
		(int nodeid, EvaluationBundleContext *c,
				shared_ptr<BundleBase> bundle = nullptr);

template
void WordCountMapper<string_range, pair<creek::string, long>,
										RecordBundle, wc_mapper::WC>::ExecEvaluator
		(int nodeid, EvaluationBundleContext *c,
				shared_ptr<BundleBase> bundle = nullptr);

/* todo: using record bitmap bundle for input/output */
