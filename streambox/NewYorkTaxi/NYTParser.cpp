#define K2_NO_DEBUG 1

#include "NewYorkTaxi/NYTRecord.h"
#include "NewYorkTaxi/NYTParser.h"
#include "NewYorkTaxi/NYTParserEvaluator.h"

template<class InputT, class OutputT, template<class> class BundleT_>
void NYTParser<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
                                                         EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr) {

	/* instantiate an evaluator */
	NYTParserEvaluator<InputT, OutputT, BundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}

template
void NYTParser<string_range, NYTRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                     EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
