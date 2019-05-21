#include "Values.h"
#include "NetworkLatencySource.h"
#include "NetworkLatencySourceEvaluator.h"

void NetworkLatencySource::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr){

	NetworkLatencySourceEvaluator eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}
