/*
 * The transform evaluator interface.
 * They are mostly transient and stateless.
 */
#ifndef TRANSFORM_EVAL_H
#define TRANSFORM_EVAL_H

class EvaluationBundleContext;

template<class TransformT>
class TransformEvaulator {
public:

  virtual void evaluate(TransformT* transform, EvaluationBundleContext* c,
  		shared_ptr<BundleBase> ptr)
    { assert(0); }

#if 0
	void evaluateOneBundle (TransformT* trans,
			EvaluationBundleContext* c, shared_ptr<BundleBase>) { assert(0); }
#endif

  virtual ~TransformEvaulator() { }

  // callback when a watermark change has been propagated from upstream, and
  // hit the transform
  // @watermark: the new watermark after the update.
//  virtual void OnNewUpstreamWatermark(ptime watermark,
//      TransformT* trans, EvaluationBundleContext* c) { }

  /* numa node. most, but not all transform eval uses this. e.g. one source
   * transform eval may output to multiple nodes.
   *
   * when node in use, a concrete evaluator should determine the node id.
   */
  int _node = -1;
};

#endif // TRANSFORM_EVAL_H
