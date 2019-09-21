#include "Nexmark/NexmarkRecord.hpp"
#include "Nexmark/NexmarkFilter.hpp"
#include "Nexmark/NexmarkFilterEvaluator.hpp"

template<class InputT,
         class OutputT,
         template<class> class BundleT_>
void NexmarkFilter<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
                                                             EvaluationBundleContext *c,
                                                             shared_ptr<BundleBase> bundle_ptr) {
    /* instantiate an evaluator */
    NexmarkFilterEvaluator<InputT, OutputT, BundleT_> eval(nodeid);
    eval.evaluate(this, c, bundle_ptr);
}

template
void NexmarkFilter<NexmarkRecord, NexmarkRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                              EvaluationBundleContext *c,
                                                                              shared_ptr<BundleBase> bundle_ptr);
