#include "Values.h"
#include "LinearRoad/LinearRoadRecord.h"
#include "LinearRoad/LinearRoadTollRecord.h"
#include "LinearRoad/LinearRoadTollNotification.h"
#include "LinearRoad/LinearRoadTollNotificationEvaluator.h"


template <class InputT, class OutputT, template<class> class BundleT>
void LinearRoadTollNotification<InputT, OutputT, BundleT>::ExecEvaluator(int nodeid,
                                                                         EvaluationBundleContext *c,
                                                                         shared_ptr<BundleBase> bundle_ptr)
{
    /* instantiate an evaluator */
    LinearRoadTollNotificationEvaluator<InputT, OutputT, BundleT> eval(nodeid);

    eval.evaluate(this, c, bundle_ptr);
}

template
void LinearRoadTollNotification<LinearRoadRecord, LinearRoadTollRecord, RecordBundle>::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void LinearRoadTollNotification<LinearRoadRecord, pair<uint64_t, long>, RecordBundle>::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
