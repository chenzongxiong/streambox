#define K2_NO_DEBUG 1

#include "Values.h"
#include "Yahoo/YahooRecord.h"
#include "Yahoo/YahooAggregationEvaluator.h"
#include "Yahoo/YahooAggregation.h"

template <class InputT, class OutputT, template<class> class BundleT
          >
//        wc_mapper::mode mode>
void YahooAggregation<InputT, OutputT, BundleT>::ExecEvaluator(int nodeid,
                                                                    EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{

#ifndef NDEBUG // if evaluator get stuck ..
    static atomic<int> outstanding (0);
#endif

    /* instantiate an evaluator */
//    YahooAggregationEvaluator<InputT, OutputT, BundleT, mode> eval(nodeid);
    YahooAggregationEvaluator<InputT, OutputT, BundleT> eval(nodeid);

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

template <>
uint64_t YahooAggregation<YahooRecord, pair<uint64_t, uint64_t>,
                          RecordBundle>
                          //wc_mapper::WC>
::do_map(Record<YahooRecord> const & in, shared_ptr<RecordBundle<pair<uint64_t, uint64_t>>> output_bundle)
{
    // using KVPair = pair<uint64_t, long>;
    using KVPair = pair<uint64_t, uint64_t>;
    uint64_t bucketPos = *((uint64_t*) in.data.campaign_id + 1);
    // std::cout << "in.data.campaign_id: " << in.data.campaign_id <<
#if 0
    std::cout << "bucketPos: " << bucketPos << std::endl;
#endif
    output_bundle->emplace_record(KVPair(bucketPos, 1), in.ts);
    return 1;
}

/* no template specilization */
template <class InputT,
        class OutputT,
//					class InputBundleT,
//					class OutputBundleT,
          template<class> class BundleT>

//        wc_mapper::mode Mode>
atomic<unsigned long>
//YahooAggregation<InputT, OutputT, InputBundleT, OutputBundleT, Mode>::record_counter_(0);
//        YahooAggregation<InputT, OutputT, BundleT, Mode>::record_counter_(0);
YahooAggregation<InputT, OutputT, BundleT>::record_counter_(0);

/* -------instantiation concrete classes------- */

/* using record bundle for input/output */

// template
// void YahooAggregation<YahooRecord, pair<uint64_t, long>, RecordBundle>::ExecEvaluator
//         (int nodeid, EvaluationBundleContext *c,
//          shared_ptr<BundleBase> bundle = nullptr);

template
void YahooAggregation<YahooRecord, pair<uint64_t, uint64_t>, RecordBundle>::ExecEvaluator
        (int nodeid, EvaluationBundleContext *c,
         shared_ptr<BundleBase> bundle = nullptr);

/* todo: using record bitmap bundle for input/output */
