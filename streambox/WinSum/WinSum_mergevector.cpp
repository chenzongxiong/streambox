#include "WinSum_mergevector.h"

//using namespace creek; /* NB: we are not directly mapping to tbb etc. */
using crstring = creek::string;
using stringvec = creek::concurrent_vector<creek::string>;
using stringvec_ptr = std::shared_ptr<stringvec>;

// ----------------------------------------------------- //
// operators
// ----------------------------------------------------- //

/* --- pass around vector of strings. turns out to be expensive. --- */

template <>
stringvec const & WinSum_mergevector<creek::string, stringvec>::aggregate
		(stringvec * acc, creek::string const & in) {
	acc->push_back(in);
	return *acc;
}

template <>
stringvec const & WinSum_mergevector<creek::string, stringvec>::combine
	(stringvec & mine, stringvec const & others) {
	mine.grow_by(others.begin(), others.end());
	return mine;
}

template <>
stringvec const & WinSum_mergevector<creek::string, stringvec>::aggregate_init
		(stringvec* acc) {
	acc->clear();
	return *acc;
}

/* --- shared_ptr version --- */

template <>
stringvec_ptr const & WinSum_mergevector<creek::string, stringvec_ptr>::aggregate
		(stringvec_ptr * acc, creek::string const & in) {
	assert(*acc);
	(*acc)->push_back(in);
	return *acc;
}

template <>
stringvec_ptr const & WinSum_mergevector<creek::string, stringvec_ptr>::combine
	(stringvec_ptr & mine, stringvec_ptr const & others) {
	assert(mine && others);
	mine->grow_by(others->begin(), others->end());
	return mine;
}

template <>
stringvec_ptr const & WinSum_mergevector<creek::string, stringvec_ptr>::aggregate_init
		(stringvec_ptr* acc) {
	assert(acc);
	*acc = make_shared<stringvec>();
	(*acc)->clear();
	return *acc;
}

/* --- unlocked version --- */

template <>
vector<string> const & WinSum_mergevector<string, vector<string>>::aggregate_init
(vector<string> * acc) {
	acc->clear();
	return *acc;
}

template <>
vector<string> const & WinSum_mergevector<string, vector<string>>::aggregate
(vector<string> * acc, string const & in) {
	acc->push_back(in);
	return *acc;
}

/* caller should hold lock */
template <>
vector<string> const & WinSum_mergevector<string, vector<string>>::combine
(vector<string> & mine, vector<string> const & others) {
	mine.insert(mine.end(), others.begin(), others.end());
	return mine;
}

/* make life easier */
//template <class InputT, class OutpuT>
//using MyWinSum = WinSum_stringvector<InputT, OutpuT>;

#define MyWinSum WinSum_mergevector
#include "WinSum-dispatch-eval.h"

/* instantiate concreate classes, otherwise won't link */

template
void WinSum_mergevector<creek::string, stringvec_ptr>
::ExecEvaluator(int nodeid,
                EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);
