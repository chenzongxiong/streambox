//
// Created by xzl on 2/15/17.
//

#ifndef CREEK_WINKEYREDUCEROPS_H
#define CREEK_WINKEYREDUCEROPS_H

//template <class InputT, class LocalWindowResultT, class WindowResultT,


template <
		typename KVPairIn,
		/* the input windowkeyfrag. its vcontainer must be the same as the
		 * upstream's output, e.g wingbk.
		 * the format can be really cheap, i.e. no concurrent support.
		 */
		template<class> class InputWinKeyFragT,
		/* internal state. concurrent d/s
		 * can be specialized based on key/val distribution.
		 * regarding the d/s for holding eval's local aggregation results,
		 * see comment below. */
		template<class> class InternalWinKeyFragT, // = WindowKeyedFragment
		>
class WinKeyReducerOps
{

	using K = decltype(KVPairIn::first);
	/* container can be unsafe */
	using InputValueContainerT = typename InputWinKeyFragT<KVPairIn>::ValueContainerT;
	/* container must be safe */
	using InternalValueContainerT = typename InternalWinKeyFragT<KVPairIn>::ValueContainerT;
	using InputT = InputWinKeyFragT<KVPairIn>;
	using WindowResultT = shared_ptr<InternalWinKeyFragT<KVPairIn>>;
	using LocalWindowResultT = WindowResultT;

	// reduce on one specific key and all its values (as in a value container)
	virtual KVPairIn do_reduce(K const &, InternalValueContainerT const &);
	virtual KVPairIn do_reduce_unsafe(K const &, InputValueContainerT const &);

	/* Take one InputT. to be specialized */
	static WindowResultT const & aggregate_init(WindowResultT * acc) {
		xzl_assert(acc);
		(*acc) = make_shared<typename WindowResultT::element_type>();
		return *acc;
	}

	static LocalWindowResultT const & local_aggregate_init(LocalWindowResultT * acc) {
		xzl_assert(acc);
		/* more generic -- don't have to change back/forth */
		(*acc) = make_shared<typename LocalWindowResultT::element_type>();
		return *acc;
	}

	static WindowResultT const & aggregate(WindowResultT * acc, InputT const & in);

	/* combine the evaluator's (partial) aggregation results (for a particular window)
	 * to the tran's internal state.
	 */
	static WindowResultT const & combine(WindowResultT & mine, LocalWindowResultT const & others);
};

#endif //CREEK_WINKEYREDUCEROPS_H
