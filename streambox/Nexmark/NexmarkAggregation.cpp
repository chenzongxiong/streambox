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
    uint64_t auction = in.data.auction;
    output_bundle->emplace_record(KVPair(auction, 1), in.ts);
    return 1;
}

template<>
uint64_t NexmarkAggregation<KVPair, KVPair, RecordBundle>::do_map(Record<KVPair> const& in,
                                                                  shared_ptr<RecordBundle<KVPair>> output_bundle) {
    // uint64_t auction = in.data.auction;
    // output_bundle->emplace_record(KVPair(auction, 1), in.ts);
    // std::cout << "nexmark aggregation 2" << std::endl;
    output_bundle->emplace_record(in.data, in.ts);

    return 1;
}


template<>
uint64_t NexmarkAggregation<NexmarkRecord, NexmarkRecord, RecordBundle>::do_map(Record<NexmarkRecord> const& in,
                                                                                shared_ptr<RecordBundle<NexmarkRecord>> output_bundle) {
    uint64_t price = in.data.price * 89 / 100;
    NexmarkRecord tuple(in.data.auction, in.data.bidder, price, in.data.dateTime);
    output_bundle->emplace_record(tuple, in.ts);
    return 0;
}

template<>
uint64_t NexmarkAggregation<NexmarkRecord, NexmarkOutputRecord, RecordBundle>::do_map(Record<NexmarkRecord> const& in,
                                                                                      shared_ptr<RecordBundle<NexmarkOutputRecord>> output_bundle) {
    uint64_t price = in.data.price * 89 / 100;
    NexmarkOutputRecord tuple(in.data.auction, price);
    output_bundle->emplace_record(tuple, in.ts);
    return 0;
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
