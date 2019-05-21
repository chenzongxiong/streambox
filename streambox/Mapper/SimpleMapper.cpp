#ifndef SIMPLEMAPPER_CPP
#define SIMPLEMAPPER_CPP

#define K2_NO_DEBUG 1

#include "SimpleMapper.h"
#include "SimpleMapperEvaluator.h"

template<>
void SimpleMapper<>::ExecEvaluator(int nodeid,
	EvaluationBundleContext *c, shared_ptr<BundleBase> bundle_ptr)
{
	
#ifndef NDEBUG /* if evaluators get stuck ...*/
	static atomic<int> outstanding (0);
#endif
	//std::cout << "SimpleMapper<>::ExecEvaluator" << std::endl;
	
	/* instantiate an evaluator */
	SimpleMapperEvaluator<long, pair<long, long>> eval(nodeid);

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



#endif /* SIMPLEMAPPER_CPP */
