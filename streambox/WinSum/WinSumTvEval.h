
/*
 * will produce record whose *value* is pair of <ts,value>, where ts
 * is generated based on window and value is window result.
 *
 * See  WindowedSumEvaluator1-tvpair.cpp for old design
 * As WinSumEval's specialization for tweet
*/

#ifndef CREEK_WINSUMTSEVAL_H
#define CREEK_WINSUMTSEVAL_H

#include "WinSumEval.h"

template <class InputT, class OutputT, class TransformT>
class WinSumTvEval : public WinSumEval<InputT, OutputT, TransformT>
{
public:
		WinSumTvEval(int node) : WinSumEval<InputT, OutputT, TransformT>(node) { }

private:
		/* @value must be a kvpair */
		Record<OutputT> const makeRecord(OutputT & value,
		                                   Window const & win) override {
			/* rewrite the record ts according to the window */
			value.first = win.start.ticks() + win.duration.ticks();
			return Record<OutputT>(value, win.window_end());
		}
};

#endif //CREEK_WINSUMTSEVAL_H
