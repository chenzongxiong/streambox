#define K2_NO_DEBUG 1

#include "Yahoo/YahooRecord.h"
#include "Yahoo/YahooParser.h"
#include "Yahoo/YahooParserEvaluator.h"

template<class InputT, class OutputT, template<class> class BundleT_>
void YahooParser<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
	EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{

#ifndef NDEBUG /* if evaluators get stuck ...*/
	static atomic<int> outstanding (0);
#endif

	/* instantiate an evaluator */
	YahooParserEvaluator<InputT, OutputT, BundleT_> eval(nodeid);

#ifndef NDEBUG		// for debug
//	I("begin eval...");
	outstanding ++;
#endif

	eval.evaluate(this, c, bundle_ptr);
	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;

#ifndef NDEBUG   // for debug
	outstanding --; int i = outstanding; i = i; I("end eval... outstanding = %d", i);
#endif
}

template
void YahooParser<string_range, YahooRecord, RecordBundle>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
