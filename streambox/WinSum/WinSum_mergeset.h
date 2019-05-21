#ifndef WINSUM_MERGESET_H
#define WINSUM_MERGESET_H

#include "WinSumBase.h"

//#define MyWinSum WinSum_mergeset
//#include "trans-template.h"
//#undef MyWinSum

#if 1
template <typename InputT, typename OutputT>
class WinSum_mergeset : public WinSumBase<WinSum_mergeset, InputT, OutputT> {

public:

	WinSum_mergeset(string name, int multi = 1)
					: WinSumBase<WinSum_mergeset, InputT, OutputT>(name, multi) { }
	
  /* for aggregating one single input. to be specialized */
  static OutputT const & aggregate_init(OutputT * acc);
  static OutputT const & aggregate(OutputT * acc, InputT const & in);

  /* combine the evaluator's (partial) aggregation results (for a particular window)
   * to the tran's internal state.
   */
  static OutputT const & combine(OutputT & mine, OutputT const & others);

	void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		std::shared_ptr<BundleBase> bundle) override;
};
#endif

#endif /* WINSUM_MERGESET_H */

