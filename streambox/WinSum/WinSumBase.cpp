
#include "WinSumBase.h"
//#include "WinSumEvaluator.h"
#include "WinSumEval.h"

//template <class InputT, class OutputT>
//void trans<InputT,OutputT>::ExecEvaluator(int nodeid,
//		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle)
//{
//	/* instantiate an evaluator */
//	WinSumEval<InputT,OutputT> eval(nodeid);
//	eval.evaluate(this, c, bundle);
//}

/* -------instantiation concrete classes-------
 * they might not be used. but that's fine.
 * */
#if 0 /* never need to instantiation this base class */
template
void trans<creek::string, creek::concurrent_vector_ptr<creek::string>>
	::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

template
void trans<creek::string, creek::concurrent_unordered_set<creek::string>>
        ::ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                        shared_ptr<BundleBase> bundle);

/* -- partitioned hashtable -- */

template
void trans<creek::string, creek_set_array::SetArrayPtr>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

/* -- tweet (sentiment) -- */

template
void trans<long,long>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

/* --- tweet (tvpair to ease join) -- */

template
void trans<long,creek::tvpair>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

template
void trans<string,vector<string>>::ExecEvaluator(int nodeid,
		EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);
#endif

#if 0
template <>
map<Window, vector<string>, Window> const &
	trans<string, vector<string>>::combine
	(map<Window, vector<string>, Window> & mine,
			map<Window, vector<string>, Window> const & others)
{
	/* if @mine does not have such a window, this will auto create one */
	for (auto & wvec : others) {
		mine[wvec.first].insert(mine[wvec.first].end(),
				wvec.second.begin(), wvec.second.end());
	}

	return mine;
}
#endif
