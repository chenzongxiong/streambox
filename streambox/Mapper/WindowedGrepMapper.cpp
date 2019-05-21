#include "WindowedGrepMapper.h"
#include "WindowedGrepMapperEvaluator.h"

template<>
void WindowedGrepMapper<>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle)
{
	/* instantiate an evaluator */
		WindowedGrepMapperEvaluator<string_range, creek::string> eval(nodeid);
		eval.evaluate(this, c, bundle);
}
