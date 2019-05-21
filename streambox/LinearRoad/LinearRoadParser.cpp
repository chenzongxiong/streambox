#include "Values.h"
#include "LinearRoad/LinearRoadParser.h"
#include "LinearRoad/LinearRoadParserEvaluator.h"


template <class InputT, class OutputT, template<class> class BundleT>
void LinearRoadParser<InputT, OutputT, BundleT>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* instantiate an evaluator */
	LinearRoadParserEvaluator<InputT, OutputT, BundleT> eval(nodeid);

	eval.evaluate(this, c, bundle_ptr);
}

template
void LinearRoadParser<string_range, LinearRoadRecord, RecordBundle>::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
