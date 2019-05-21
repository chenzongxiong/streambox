#define K2_NO_DEBUG 1

#include "JoinEvaluator1.h"
#include "core/EvaluationBundleContext.h"

template<class KVPair,
template<class> class InputBundleT_,
template<class> class OutputBundleT_
>
void Join<KVPair, InputBundleT_, OutputBundleT_>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
#ifdef DEBUG
	std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
#endif
#ifndef NDEBUG /* if evaluators get stuck ...*/
	static atomic<int> outstanding (0);
#endif

	/* instantiate an evaluator */
	//GrepMapperEvaluator<string_range,string> eval(nodeid);
//	JoinEvaluator1<pair<long, long>, InputBundleT_, OutputBundleT_> eval(nodeid);
	JoinEvaluator1<KVPair, InputBundleT_, OutputBundleT_> eval(nodeid);
	//JoinEvaluator1<KVPair> eval(nodeid);

#ifndef NDEBUG		// for debug
//	I("begin eval...");
	outstanding ++;
#endif

	eval.evaluate(this, c, bundle_ptr);

#ifndef NDEBUG   // for debug
	outstanding --; int i = outstanding; i = i; I("end eval... outstanding = %d", i);
#endif
}

//template class Join<pair<long, long>>; //Fix error:  reference undefined ExecEvaluator
//http://stackoverflow.com/questions/115703/storing-c-template-function-definitions-in-a-cpp-file

template
void Join<std::pair<long, long>, RecordBitmapBundle, RecordBitmapBundle>::ExecEvaluator(int nodeid,
	EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);

template
void Join<std::pair<long, long>, RecordBundle, RecordBundle>::ExecEvaluator(int nodeid,
	EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr);
