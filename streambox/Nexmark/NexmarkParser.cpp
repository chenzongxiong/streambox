#include "Nexmark/NexmarkParser.hpp"
#include "Nexmark/NexmarkParserEvaluator.hpp"

template<class InputT,
         class OutputT,
         template<class> class BundleT_>
void NexmarkParser<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
                                                             EvaluationBundleContext *c,
                                                             shared_ptr<BundleBase> bundle_ptr) {
	/* instantiate an evaluator */
	NexmarkParserEvaluator<InputT, OutputT, BundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}

template
void NexmarkParser<string_range, NexmarkRecord, RecordBundle>::ExecEvaluator(int nodeid,
                                                                             EvaluationBundleContext *c,
                                                                             shared_ptr<BundleBase> bundle_ptr);
