#include "WinSum_mergeset.h"
#include "WinSumEval.h"

// ----------------------------------------------------- //
// operators
// ----------------------------------------------------- //

/* -- for distinct count -- */

#ifdef USE_FOLLY_HASHMAP

/* (uint64_t as key, for folly::AtomicHashMap) */

using Set = folly::AtomicHashMap<uint64_t, uint64_t>;
using SetPtr = shared_ptr<Set>;

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate
		(SetPtr * acc, creek::string const & in) {
	assert(*acc);
//	(*acc)->insert( (long unsigned int)(tbb::tbb_hasher(in)) );  /* XXX if std::string, use default Hash */
//	uint64_t s = (uint64_t)(tbb::tbb_hasher(in));
	(*acc)->insert(make_pair(tbb::tbb_hasher(in), 1));
	return *acc;
}

/* AtomicHashMap has no range insertion */
template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::combine
	(SetPtr & mine, SetPtr const & others) {
	assert(mine && others);
	for (auto it = others->begin(); it != others->end(); ++it)
		mine->insert(*it);
	return mine;
}

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate_init
		(SetPtr* acc) {
	assert(acc);
	*acc = make_shared<Set>(8 * 1024); /* init hashset size; does not help much */
	return *acc;
}

#elif defined(USE_CUCKOO_HASHMAP)

/* --  (string as key) -- */

using Set = cuckoohash_map<creek::string, int, CityHasher<creek::string>>;
using SetPtr = shared_ptr<Set>;

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate
		(SetPtr * acc, creek::string const & in) {
	assert(*acc);
	(*acc)->insert(in, 0);
	return *acc;
}

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::combine
	(SetPtr & mine, SetPtr const & others) {
	assert(mine && others);
	const auto & lck = others->lock_table();
	for (auto it = lck.begin(); it != lck.end(); it++) {
		mine->insert(it->first, it->second);
	}
	return mine;
}

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate_init
		(SetPtr* acc) {
	assert(acc);
	*acc = make_shared<Set>();
	return *acc;
}

#elif defined(USE_TBB_HASHMAP)
/* --  (string as key) -- */
using Set = tbb::concurrent_unordered_set<creek::string>;
using SetPtr = shared_ptr<Set>;


template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate
		(SetPtr * acc, creek::string const & in) {
	assert(*acc);
	(*acc)->insert(in);
	return *acc;
}

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::combine
	(SetPtr & mine, SetPtr const & others) {
	assert(mine && others);
	mine->insert(others->begin(), others->end());
	return mine;
}

template <>
SetPtr const & WinSum_mergeset<creek::string, SetPtr>::aggregate_init
		(SetPtr* acc) {
	assert(acc);
	*acc = make_shared<Set>();
	return *acc;
}

/* using partitioned hash table */
template <>
creek_set_array::SetArrayPtr const & WinSum_mergeset<creek::string, creek_set_array::SetArrayPtr>::aggregate
		(creek_set_array::SetArrayPtr * acc, creek::string const & in) {
	assert(*acc);
	creek_set_array::insert(*acc, in);
	return *acc;
}

template <>
creek_set_array::SetArrayPtr const & WinSum_mergeset<creek::string, creek_set_array::SetArrayPtr>::combine
	(creek_set_array::SetArrayPtr & mine, creek_set_array::SetArrayPtr const & others) {
	assert(mine && others);
	creek_set_array::merge(mine, others);
	return mine;
}

template <>
creek_set_array::SetArrayPtr const & WinSum_mergeset<creek::string, creek_set_array::SetArrayPtr>::aggregate_init
		(creek_set_array::SetArrayPtr* acc) {
	assert(acc);
	*acc = make_shared<creek_set_array::SetArray>();
	return *acc;
}

//template <typename InputT, typename OutputT>
//using MyWinSum = WinSum_mergeset<InputT, OutputT>;
#define MyWinSum WinSum_mergeset
#include "WinSum-dispatch-eval.h"

/* instantiate concreate classes, otherwise won't link */

template
void WinSum_mergeset<creek::string, creek_set_array::SetArrayPtr>
::ExecEvaluator(int nodeid,
                EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

template
void WinSum_mergeset<creek::string, SetPtr>
::ExecEvaluator(int nodeid,
                EvaluationBundleContext *c, shared_ptr<BundleBase> bundle);

#endif // USE_TBB_HASHMAP
