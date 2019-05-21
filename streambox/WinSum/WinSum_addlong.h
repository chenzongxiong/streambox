
#ifndef WINDOWEDSUMADD_H
#define WINDOWEDSUMADD_H

#include "WinSumBase.h"

//#define trans WinSum_addlong
//#include "trans-template.h"

#if 1
/* only have to define static ops here */
template <typename InputT, typename OutputT>
class WinSum_addlong : public WinSumBase<WinSum_addlong, InputT, OutputT> {

public:

	WinSum_addlong(string name, int multi = 1)
		: WinSumBase<WinSum_addlong, InputT, OutputT>(name, multi) { }
	
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

#endif /* WINDOWEDSUMADD_H */

