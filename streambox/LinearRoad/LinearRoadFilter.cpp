#include "LinearRoad/LinearRoadRecord.h"
#include "LinearRoad/LinearRoadFilter.h"
#include "LinearRoad/LinearRoadFilterEvaluator.h"

template<class InputT, class OutputT, template<class> class BundleT_>
void LinearRoadFilter<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
                                                                EvaluationBundleContext *c,
                                                                shared_ptr<BundleBase> bundle_ptr)
{
  /* instantiate an evaluator */
  LinearRoadFilterEvaluator<InputT, OutputT, BundleT_> eval(nodeid);

  eval.evaluate(this, c, bundle_ptr);
}

template
void LinearRoadFilter<LinearRoadRecord, LinearRoadRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                                      EvaluationBundleContext *c,
                                                                                      shared_ptr<BundleBase> bundle_ptr);
