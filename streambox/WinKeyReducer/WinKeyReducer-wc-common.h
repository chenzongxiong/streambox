/* to be included in WinKeyReducer-wc-XXX.cpp
 *
 * because differet bundle types share the same transform ops
 * */

#ifndef CREEK_WINKEYREDUCER_WC_OPS_H
#define CREEK_WINKEYREDUCER_WC_OPS_H


/* ----------------------- specialization for wordcount --------------------------------- */

/* for a particular window */
/* not in use: since eval do local reduce and then combine.
 * And it requires to iterate HT, not supported by creek::map */
template <>
WindowResultPtr const & MyWinKeyReducer::aggregate
(WindowResultPtr * acc, InputT const & in)
{
#if 0
assert(acc);
	assert(WindowEqual()((*acc)->w, in.w));

	for (auto && kvs : in.vals) {
		auto & key = kvs.first;
		auto & v_container = kvs.second;
		(*acc)->add_vcontainer_safe(key, v_container);
	}
#endif
	xzl_bug("not impl");
	return *acc;
}

/* XXX if @mine is empty, should we just swap pointers? */
template <>
WindowResultPtr const & MyWinKeyReducer::combine
(WindowResultPtr & mine, WindowResultPtr const & others)
{
	xzl_assert(mine && others);
	xzl_assert(WindowEqual()(mine->w, others->w));

	for (auto && kvs : others->vals) {
	auto & key = kvs.first;
	auto & v_container = kvs.second;
	mine->add_vcontainer_safe(key, v_container);
	}

	return mine;
}

template<>
std::pair<creek::string, long> MyWinKeyReducer::do_reduce
(creek::string const & key, InternalValueContainerT const & vcontainer) {

	auto end = vcontainer.cend(); /* avoid calling it in each iteration */

	long sum = 0;
	for (auto it = vcontainer.cbegin(); it != end; ++it)
	sum += *it;

	return make_pair(key, sum);
}

template<>
std::pair<creek::string, long> MyWinKeyReducer::do_reduce_unsafe
(creek::string const & key, InputValueContainerT const & vcontainer) {

	auto end = vcontainer.cend(); /* avoid calling it in each iteration */

	long sum = 0;
	for (auto it = vcontainer.cbegin(); it != end; ++it)
	sum += *it;

	return make_pair(key, sum);
}



#endif //CREEK_WINKEYREDUCER_WC_OPS_H

