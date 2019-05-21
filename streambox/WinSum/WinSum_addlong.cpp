// ----------------------------------------------------- //
// eval selector
// ----------------------------------------------------- //

//template <class InputT, class OutputT>
//void MyWinSum<InputT,OutputT>::ExecEvaluator(int nodeid,
//		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle)
//{
//	/* instantiate an evaluator */
//	WinSumEval<InputT,OutputT,MyWinSum> eval(nodeid);
//	eval.evaluate(this, c, bundle);
//}

// ----------------------------------------------------- //
// operators
// ----------------------------------------------------- //

#include "WinSum_addlong.h"

template <>
long const & WinSum_addlong<long, long>::aggregate(long * acc,
		long const & in) {
	*acc += in;
	return *acc;
}

template <>
long const & WinSum_addlong<long, long>::aggregate_init(long * acc) {
	*acc = 0;
	return *acc;
}

template <>
long const & WinSum_addlong<long, long>::combine
	(long & mine, long const & others) {
	mine += others;
	return mine;
}

/* must come after static member methods. type alias does not work?  */

//template <class InputT, class OutputT>
//using MyWinSum = WinSum_addlong<InputT,OutputT>;

#define MyWinSum WinSum_addlong
#include "WinSum-dispatch-eval.h"

/* TODO: instantiate concreate classes, otherwise won't link */

//template <class InputT, class OutputT>
//void MyWinSum<InputT,OutputT>::ExecEvaluator(int nodeid,
//                                             EvaluationBundleContext *c, shared_ptr<BundleBase> bundle)
//{
//	/* instantiate an evaluator */
//	WinSumEval<InputT,OutputT,MyWinSum<InputT,OutputT>> eval(nodeid);
//	eval.evaluate(this, c, bundle);
//}

/* instantiate concreate classes, otherwise won't link */

template
void
WinSum_addlong<long, long>
::ExecEvaluator(int nodeid,
                EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);
