#include "Values.h"
#include "Unbounded_Join.h"
#include "UnboundedInMemEvaluator_Join.h"

void UnboundedInMem_Join::ExecEvaluator(int nodeid, EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr){
	UnboundedInMemEvaluator_Join eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}
#if 0
//hym:
template<template<class> class BundleT>
void UnboundedInMem_Join<long, BundleT>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr) {
	/* instantiate an evaluator */
	//UnboundedInMemEvaluator_Join<long, BundleT> eval(nodeid);
	//UnboundedInMemEvaluator_Join<long, RecordBitmapBundle> eval(nodeid);
	UnboundedInMemEvaluator_Join<long, BundleT<pair<long, long>>> eval(nodeid);
	
	eval.evaluate(this, c, bundle_ptr);
}
//hym: instantiation
template
void UnboundedInMem_Join<long, RecordBitmapBundle>
			::ExecEvaluator(int nodeid,
					EvaluationBundleContext *c,
					shared_ptr<BundleBase> bundle_ptr);
#endif
