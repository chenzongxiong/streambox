//
// to be included from WinSum_XXX.cpp

#ifndef CREEK_WINSUM_DISPATCH_EVAL_H
#define CREEK_WINSUM_DISPATCH_EVAL_H

#include "WinSumEval.h"

// ----------------------------------------------------- //
// dispatch to eval (who needs to see concrete types)
// ----------------------------------------------------- //
#ifndef MyWinSumEval  /* WinSum_XXX.cpp can optionally define this */
#define MyWinSumEval WinSumEval
#endif

template <class InputT, class OutputT>
void MyWinSum<InputT,OutputT>::ExecEvaluator(int nodeid,
                 EvaluationBundleContext *c, std::shared_ptr<BundleBase> bundle)
{
    /* instantiate an evaluator */
    MyWinSumEval<InputT,OutputT,MyWinSum<InputT,OutputT>> eval(nodeid);
    eval.evaluate(this, c, bundle);
}

#undef MyWinSumEval

#if 0
/* instantiate */
template
void MyWinSum<InputT,OutputT>::ExecEvaluator(int nodeid,
				EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);
#endif

#endif //CREEK_WINSUM_DISPATCH_EVAL_H
