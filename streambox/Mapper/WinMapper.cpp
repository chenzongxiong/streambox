/*
 * WinMapper.cpp
 *
 *  Created on: Jan 23, 2017
 *      Author: xzl
 */

#include "WinMapper.h"
#include "WinMapperEvaluator.h"

template<class InputT, class OutputT, template<class> class BundleT_>
void WinMapper<InputT, OutputT, BundleT_>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	/* instantiate an evaluator */
	WinMapperEvaluator<InputT, OutputT, BundleT_> eval(nodeid);
	eval.evaluate(this, c, bundle_ptr);
}

/* -------instantiation concrete classes------- */
template
void WinMapper<string_range, long, RecordBundle>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
