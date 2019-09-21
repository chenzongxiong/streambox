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

template
void NexmarkAggregation<NexmarkRecord, KVPair, RecordBundle>::ExecEvaluator(int nodeid,
                                                                            EvaluationBundleContext *c,
                                                                            shared_ptr<BundleBase> bundle=nullptr);
