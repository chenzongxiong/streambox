//#ifndef CREEK_WINSUM_TEMPLATE_H_H
//#define CREEK_WINSUM_TEMPLATE_H_H
#include "config.h"
#include "Values.h"
#include "core/StatefulTransform.h"

/* to be included in concrete WinSum_XXX.h */

/* This is
 *
 * Aggregation.
 *
 * - stateful.
 * - does not expect kvpair
 * - purge depending on watermark advance
 *
 * Limitation: internal windowstate is also OutputT, which may be overkill,
 * e.g. for distinct count
 */
template <typename InputT, typename OutputT>
class MyWinSum
				: public StatefulTransform<MyWinSum<InputT,OutputT>, InputT, OutputT> {

public:

		using TransformT = MyWinSum<InputT,OutputT>;

		/* multi: =1 for fixed window;
		 * >1 for sliding window (x factor for the sliding window size)
		 */
		MyWinSum(string name, int multi = 1, bool force_emit_tvpair = false)
						: StatefulTransform<TransformT, InputT, OutputT>(name), multi(multi),
						  force_emit_tvpair(force_emit_tvpair) { }

		/* for aggregating one single input. to be specialized */
    static OutputT const & aggregate_init(OutputT * acc);
    static OutputT const & aggregate(OutputT * acc, InputT const & in);

	  /* combine the evaluator's (partial) aggregation results (for a particular window)
	   * to the tran's internal state.
	   */
	  static OutputT const & combine(OutputT & mine, OutputT const & others);

		/* multiplier for supporting sliding window.
		 * = 1 for tumbling window.  */
		const int multi;

		/* force record data to be <long, record_data>, where long represents the
		 * window's end ts. hate this.
		 */
		const bool force_emit_tvpair;

		void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
		                   shared_ptr<BundleBase> bundle) override;

};

//#endif //CREEK_WINSUM_TEMPLATE_H_H
