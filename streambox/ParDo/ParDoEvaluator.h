#ifndef PARDOEVALUATOR_H
#define PARDOEVALUATOR_H

// DEPRECATED.
// see beam: ParDoEvaluator.java
// ParDo::evaluateHelper() for a trivial evaluator implementation
// arguments: the types of the input/output element types of the ParDo
//
// the version that supports bundle
template <class InputT, class OutputT>
class ParDoEvaluator : public TransformEvaulator<ParDo<InputT, OutputT>> {
public:

  void evaluate(ParDo<InputT, OutputT>* trans, EvaluationBundleContext* c) {

    PValue* in1 = trans->getFirstInput(); // XXX handle extra inputs
    assert(in1);

    // create Bundle-specific ctxs (which holds the output Bundle)
    PerBundleContext<InputT, OutputT> ctx;

    unique_lock<mutex> lock(trans->mtx_watermark);

    // get one pending bundle from the input
    auto bundle = \
        dynamic_pointer_cast<RecordBundle<InputT>>(in1->getOneBundle(this->_node));
    assert(bundle);

    assert(trans->inflight_bundles.count(bundle) == 0);
    trans->inflight_bundles.insert(bundle);
    lock.unlock();

    // XXX: start Bundle
    trans->getFn()->startBundle(&ctx);

    // Go over each element. Note that we pass in reference so that
    // they can be updated in place. Optimization: we may pointer-swing
    // the input bundle to the output bundle.
    for (auto&& v : bundle->content)
        trans->getFn()->processElementSimple(&ctx, v);

    trans->getFn()->finishBundle(&ctx);

    // since ParDo is 1:1 mapping, the bundle's min_ts remains unchanged
    ctx.outputBundle->min_ts = bundle->min_ts;

    // deposit the output Bundle (staged in the context)
    // to the output PValue
    auto out = trans->getFirstOutput();
    assert(out);

    lock.lock();
    out->depositOneBundle(ctx.outputBundle);

    // now the input bundle is gone and output bundle is commited.
    assert(trans->inflight_bundles.count(bundle) == 1);
    trans->inflight_bundles.erase(bundle);
    lock.unlock();

    c->SpawnConsumer(out, this->_node);
  }

  ParDoEvaluator(int node) { this->_node = node; }

  ~ParDoEvaluator() { }
};

#endif /* PARDOEVALUATOR_H */
