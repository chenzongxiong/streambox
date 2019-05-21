#ifndef WINDOWED_KEY_REC_H
#define WINDOWED_KEY_REC_H

#include "Values.h"
#include "core/StatefulTransform.h"

/* Stateful reducer of windowed/keyed records.
 * Execute DoFn on each window / each key
 *
 * InputT : WindowKeyedFragment<KVPair>
 *
 */
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
	typename KVPairOut,
	template<class> class OutputBundleT_
	>
class WinKeyReducer
		: public StatefulTransform<
		  		WinKeyReducer<KVPairIn,InputWinKeyFragT,InternalWinKeyFragT,KVPairOut,OutputBundleT_>,
		  		InputWinKeyFragT<KVPairIn>,			/* InputT */
		  		shared_ptr<InternalWinKeyFragT<KVPairIn>>,		/* WindowResultT */
		  		shared_ptr<InternalWinKeyFragT<KVPairIn>>		/* LocalWindowResultT */
//		  		shared_ptr<InputWindowKeyedFragmentT<KVPair>>		/* LocalWindowResultT */
		  >
{

public:
	using K = decltype(KVPairIn::first);
	using V = decltype(KVPairIn::second);

	/* container can be unsafe */
	using InputValueContainerT = typename InputWinKeyFragT<KVPairIn>::ValueContainerT;
	/* container must be safe */
	using InternalValueContainerT = typename InternalWinKeyFragT<KVPairIn>::ValueContainerT;

	//  using windows_map = map<Window, WindowKeyedFragment<KVPair>, Window>;
	using TransformT = WinKeyReducer<KVPairIn, InputWinKeyFragT, InternalWinKeyFragT,
	      KVPairOut,OutputBundleT_>;
	using InputT = InputWinKeyFragT<KVPairIn>;

	/* Aggregation types
	 *
	 * localwindowresult: used by eval to hold local result
	 * windowresult: used for transform's internal state.
	 *
	 * the intention of separating these two is that localwindowresult does not have
	 * to be concurrent d/s and therefore can be cheap.
	 * However, when making localwindowresult as based on InputWindowKeyedFragmentT
	 * (e.g. std::unordered_map with SimpleValueContainerUnsafe), the performance
	 * is even worse than (localwindowresult == windowresult).
	 *
	 * This can be tested by toggling the comments below.
	 */
	using WindowResultT = shared_ptr<InternalWinKeyFragT<KVPairIn>>;
	using LocalWindowResultT = WindowResultT;
	//  using LocalWindowResultT = shared_ptr<InputWindowKeyedFragmentT<KVPair>>;

	/* AggResultT and LocalAggResultT are decl in StatefulTransform */
	//  using LocalAggResultT = map<Window, LocalWindowResultT, Window>;

public:
	WinKeyReducer(string name)
		: StatefulTransform<TransformT, InputT, WindowResultT, LocalWindowResultT>(name), record_counter_(0) { }

	////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////

	// reduce on one specific key and all its values (as in a value container)
	virtual KVPairIn do_reduce(K const &, InternalValueContainerT const &);
	virtual KVPairIn do_reduce_unsafe(K const &, InputValueContainerT const &);

	/* Take one InputT. to be specialized */
	static WindowResultT const & aggregate_init(WindowResultT * acc) {
		xzl_assert(acc);
		//  	(*acc) = make_shared<InternalWindowKeyedFragmentT<KVPair>>();
		(*acc) = make_shared<typename WindowResultT::element_type>();
		return *acc;
	}

	static LocalWindowResultT const & local_aggregate_init(LocalWindowResultT * acc) {
		xzl_assert(acc);
		//  	(*acc) = make_shared<InputWindowKeyedFragmentT<KVPair>>();
		//  	(*acc) = make_shared<InternalWindowKeyedFragmentT<KVPair>>();
		/* more generic -- don't have to change back/forth */
		(*acc) = make_shared<typename LocalWindowResultT::element_type>();
		return *acc;
	}

	static WindowResultT const & aggregate(WindowResultT * acc, InputT const & in);

	/* combine the evaluator's (partial) aggregation results (for a particular window)
	 * to the tran's internal state.
	 */
	static WindowResultT const & combine(WindowResultT & mine, LocalWindowResultT const & others);

	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
			shared_ptr<BundleBase> bundle_ptr) override;

    std::atomic<unsigned long> record_counter_;
    // zxchen: lsds yahoo benchmark doesn't implement this function, work-around method
    bool ReportStatistics(PTransform::Statstics* stat) override { return false;}
};


/* the default version that outputs WindowsBundle */
template <typename KVPairIn, template<class> class InputWinKeyFragT, template<class> class InternalWinKeyFragT /* = WindowKeyedFragment*/>
    using WinKeyReducer_winbundle = WinKeyReducer<KVPairIn, InputWinKeyFragT, InternalWinKeyFragT, KVPairIn, /* kv out */WindowsBundle>;

#endif // WINDOWED_KEY_REC_H
