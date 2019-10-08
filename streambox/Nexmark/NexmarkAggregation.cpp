#include "Values.h"
#include "Nexmark/NexmarkRecord.hpp"
#include "Nexmark/NexmarkAggregationEvaluator.hpp"
#include "Nexmark/NexmarkAggregation.hpp"


template <class InputT,
          class OutputT,
          template<class> class BundleT
          >
void NexmarkAggregation<InputT, OutputT, BundleT>::ExecEvaluator(int nodeid,
                                                                 EvaluationBundleContext *c,
                                                                 shared_ptr<BundleBase> bundle_ptr) {
    NexmarkAggregationEvaluator<InputT, OutputT, BundleT> eval(nodeid);
    eval.evaluate(this, c, bundle_ptr);
}

using KVPair = pair<uint64_t, uint64_t>;

template<>
uint64_t NexmarkAggregation<NexmarkRecord, KVPair, RecordBundle>::do_map(Record<NexmarkRecord> const& in,
                                                                         shared_ptr<RecordBundle<KVPair>> output_bundle) {
    ptime ts = boost::posix_time::microsec_clock::local_time(); //  has negative impact on throughput
    uint64_t auction = in.data.auction;
    // if ((in.ts - Window::epoch).total_microseconds() < 0) {
    //     std::cout << "name: " << this->name << ", ts: " << (in.ts - Window::epoch).total_microseconds() << std::endl;
    //     abort();
    // }
    output_bundle->emplace_record(KVPair(auction, 1), ts);
    // output_bundle->emplace_record(KVPair(auction, 1), in.ts);
    return 1;

    // uint64_t price = in.data.price;
    // // if ((in.ts - Window::epoch).total_microseconds() < 0) {
    // //     std::cout << "name: " << this->name << ", ts: " << (in.ts - Window::epoch).total_microseconds() << std::endl;
    // //     abort();
    // // }
    // if (maxPrice < price) {
    //     maxPrice = price;
    //     output_bundle->emplace_record(KVPair(0, price), ts);
    // }

    // // output_bundle->emplace_record(KVPair(auction, 1), ts);
    // // output_bundle->emplace_record(KVPair(auction, 1), in.ts);
    // return 1;

}

template<>
uint64_t NexmarkAggregation<KVPair, KVPair, RecordBundle>::do_map(Record<KVPair> const& in,
                                                                  shared_ptr<RecordBundle<KVPair>> output_bundle) {
    // if ((in.ts - Window::epoch).total_microseconds() < 0) {
    //     std::cout << "name: " << this->name << ", ts: " << (in.ts - Window::epoch).total_microseconds() << std::endl;
    //     abort();
    // }
    ptime ts = boost::posix_time::microsec_clock::local_time(); //  has negative impact on throughput
    output_bundle->emplace_record(in.data, ts);
    // output_bundle->emplace_record(in.data, in.ts);

    return 1;
}


template<>
uint64_t NexmarkAggregation<NexmarkRecord, NexmarkRecord, RecordBundle>::do_map(Record<NexmarkRecord> const& in,
                                                                                shared_ptr<RecordBundle<NexmarkRecord>> output_bundle) {
    uint64_t price = in.data.price * 89 / 100;
    NexmarkRecord tuple(in.data.auction, in.data.bidder, price, in.data.dateTime);
    output_bundle->emplace_record(tuple, in.ts);
    return 1;
}

template<>
uint64_t NexmarkAggregation<NexmarkRecord, NexmarkOutputRecord, RecordBundle>::do_map(Record<NexmarkRecord> const& in,
                                                                                      shared_ptr<RecordBundle<NexmarkOutputRecord>> output_bundle) {
    uint64_t price = in.data.price * 89 / 100;
    NexmarkOutputRecord tuple(in.data.auction, price);
    output_bundle->emplace_record(tuple, in.ts);
    return 1;
}


template
void NexmarkAggregation<NexmarkRecord, KVPair, RecordBundle>::ExecEvaluator(int nodeid,
                                                                            EvaluationBundleContext *c,
                                                                            shared_ptr<BundleBase> bundle=nullptr);

template
void NexmarkAggregation<KVPair, KVPair, RecordBundle>::ExecEvaluator(int nodeid,
                                                                     EvaluationBundleContext *c,
                                                                     shared_ptr<BundleBase> bundle=nullptr);

template
void NexmarkAggregation<NexmarkRecord, NexmarkRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                                   EvaluationBundleContext *c,
                                                                                   shared_ptr<BundleBase> bundle=nullptr);

template
void NexmarkAggregation<NexmarkRecord, NexmarkOutputRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                                         EvaluationBundleContext *c,
                                                                                         shared_ptr<BundleBase> bundle=nullptr);
