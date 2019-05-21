/*
 * Stateless. a combo of windowing and grouping by key. since both transforms
 * are re-organizing the records, each will have to touch & move every record.
 * Thus, it's more reasonable to combine the two transforms.
 */
#ifndef WINGBK_H
#define WINGBK_H

#include "core/Transforms.h"

// input bundle: RecordBundle or RecordBitmapBUndle
// output bundle: WindowsKeyedBundle
template <class KVPair,
					template<class> class InputBundleT,
					template<class> class WindowKeyedFragmentT  /* don't have to be MT safe */
					>
class WinGBK : public PTransform {
  using InputT = Record<KVPair>;

public:
  const boost::posix_time::time_duration window_size;
  const ptime start; // the starting point of windowing.
  WinGBK(string name,
      boost::posix_time::time_duration window_size,
      ptime start = Window::epoch)
    : PTransform(name), window_size(window_size), start(start) { }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  			shared_ptr<BundleBase> bundle_ptr) override;

//  void ExecEvaluator(int nodeid, EvaluationBundleContext *c) override {
//  	/* instantiate an evaluator */
//  	WinGBKEvaluator<KVPair> eval(nodeid);
//		eval.evaluate(this, c);
//  }

};

#endif // WINGBK_H
