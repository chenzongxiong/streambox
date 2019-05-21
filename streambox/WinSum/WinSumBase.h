/*
 * Windowed aggregation.
 *
 * Use the new interface of StatefulTransform. For old implementation
 * see old/
 *
 * Only a base, e.g. the common vars needed by concrete WinSum_XXX.
 * For static polymorphism, still need to pass the concrete class (which has the
 * intense ops) up to @StatefulTransform.
 *
 * Should never be directly invoked.
 */

#ifndef WINDOWED_SUM_H
#define WINDOWED_SUM_H

#include "config.h"

#include "Values.h"
#include "core/StatefulTransform.h"

/* Aggregation.
 *
 * - stateful.
 * - does not expect kvpair
 * - purge depending on watermark advance
 *
 * Limitation: internal windowstate is also OutputT, which may be overkill,
 * e.g. for distinct count
 */
template <
				template<class,class> class WinSumConcreteT, /* need to pass this to StatefulTransform */
				typename InputT,
			typename OutputT = InputT>
class WinSumBase
		: public StatefulTransform<WinSumConcreteT<InputT,OutputT>, InputT, OutputT> {

public:

	using TransformT = WinSumConcreteT<InputT,OutputT>;

	/* multi: =1 for fixed window;
	 * >1 for sliding window (x factor for the sliding window size)
	 */
  WinSumBase(string name, int multi = 1, bool force_emit_tvpair = false)
	: StatefulTransform<TransformT, InputT, OutputT>(name), multi(multi),
	  force_emit_tvpair(force_emit_tvpair) { }

#if 0	
  /* for aggregating one single input. to be specialized */
  static OutputT const & aggregate_init(OutputT * acc);
  static OutputT const & aggregate(OutputT * acc, InputT const & in);

  /* combine the evaluator's (partial) aggregation results (for a particular window)
   * to the tran's internal state.
   */
  static OutputT const & combine(OutputT & mine, OutputT const & others);
#endif
	
  /* multiplier for supporting sliding window.
   * = 1 for tumbling window.  */
  const int multi;

  /* force record data to be <long, record_data>, where long represents the
   * window's end ts. hate this.
   */
  const bool force_emit_tvpair;

	/* ExecEvaluator should be in concrete WinSum_XXX */
#if 0
  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> bundle) override {
	  TransformT::ExecEvaluator(nodeid, c, bundle);
  }
#endif

//  void OnNewUpstreamWatermark(int nodeid, EvaluationBundleContext *c,
//  		ptime wm) override;

};

#if 0  /* sometimes cause link failure, why */
template <>
inline long const & WinSum<long, long>::aggregate(long & acc,
		long const & in) {
	acc += in;
	return acc;
}

template <>
inline vector<string> const & WinSum<string, vector<string>>::aggregate
		(vector<string> & acc, string const & in) {
	acc.push_back(in);
	return acc;
}

template <>
inline long const & WinSum<long, long>::aggregate_init(long & acc) {
	acc = 0;
	return acc;
}

template <>
inline vector<string> const & WinSum<string, vector<string>>::aggregate_init
		(vector<string> & acc) {
	acc.clear();
	return acc;
}
#endif

#endif // WINDOWED_SUM_H
