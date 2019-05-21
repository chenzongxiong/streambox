#include "config.h"

// #include "WinKeyReducer.h"
// #include "WinKeyReducerEval.h"

#include "NewYorkTaxi/NYTWindowedSum.h"
#include "NewYorkTaxi/NYTWindowedSumEvaluator.h"

template <typename InputT,
          typename OutputT>
void NYTWindowedSum<InputT, OutputT>::ExecEvaluator(int nodeid,
                                                    EvaluationBundleContext *c,
                                                    shared_ptr<BundleBase> bundle_ptr) {
	/* instantiate an evaluator */
	NYTWindowedSumEvaluator<InputT, OutputT, NYTWindowedSum<InputT, OutputT>> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}

template
void NYTWindowedSum<std::pair<uint64_t, uint64_t>,
                    std::pair<uint64_t, uint64_t>>::ExecEvaluator(int nodeid,
                                                                  EvaluationBundleContext *c,
                                                                  shared_ptr<BundleBase> bundle_ptr);
