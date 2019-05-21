#define K2_NO_DEBUG 1
#include "GrepMapper.h"
#include "GrepMapperEvaluator.h"

template<class InputT, class OutputT, template<class> class BundleT_>
void GrepMapper<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
#ifndef NDEBUG /* if evaluators get stuck ...*/
	static atomic<int> outstanding (0);
#endif

	/* instantiate an evaluator */
	GrepMapperEvaluator<InputT, OutputT, BundleT_> eval(nodeid);

#ifndef NDEBUG		// for debug
//	I("begin eval...");
	outstanding ++;
#endif

	eval.evaluate(this, c, bundle_ptr);

#ifndef NDEBUG   // for debug
	outstanding --; int i = outstanding; i = i; I("end eval... outstanding = %d", i);
#endif
}

/* -------instantiation concrete classes------- */
template
void GrepMapper<string_range, creek::string, RecordBundle>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

