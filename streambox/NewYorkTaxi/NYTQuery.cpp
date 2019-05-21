#include "NewYorkTaxi/NYTRecord.h"
#include "NewYorkTaxi/NYTQuery.h"
#include "NewYorkTaxi/NYTQueryEvaluator.h"


template<class InputT, class OutputT, template<class> class BundleT_>
void NYTQuery<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
                                                        EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr) {

    /* instantiate an evaluator */
    NYTQueryEvaluator<InputT, OutputT, BundleT_> eval(nodeid);

    eval.evaluate(this, c, bundle_ptr);
}

// template
// void NYTQuery<NYTRecord,
//               NYTRecord,
//               WindowsBundle>::ExecEvaluator(int nodeid,
//                                             EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);


// template
// void NYTQuery<NYTRecord,
//               NYTRecord,
//               RecordBundle>::ExecEvaluator(int nodeid,
//                                             EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);


template
void NYTQuery<NYTRecord,
              // region_id, trip_distance
              std::pair<uint64_t, uint64_t>,
              RecordBundle>::ExecEvaluator(int nodeid,
                                           EvaluationBundleContext *c,
                                           shared_ptr<BundleBase> bundle_ptr);
