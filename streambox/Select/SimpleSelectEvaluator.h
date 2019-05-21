#ifndef SIMPLESELECTEVALUATOR_H
#define SIMPLESELECTEVALUATOR_H

/* mask records by marking the bitmap that comes with the bundle */
template <class InputT>
class SimpleSelectEvaluator
    : public TransformEvaulator<SimpleSelect<InputT>> {

  using TransformT = SimpleSelect<InputT>;
  using InputBundleT = RecordBitmapBundle<InputT>;
  using OutputBundleT = RecordBitmapBundle<InputT>;

public:
  SimpleSelectEvaluator(int node) {
  	ASSERT_VALID_NUMA_NODE(node);
  	this->_node = node;
  }

  void evaluate(TransformT* trans, EvaluationBundleContext* c) {

    PValue* in1 = trans->getFirstInput();
    assert(in1);

    // get one pending bundle from the input.
    // this will update input's min_ts. note that the bundle is still
    // a "pending" work in the current transform

    unique_lock<mutex> lock(trans->mtx_watermark);

    auto bundle = \
        dynamic_pointer_cast<InputBundleT>(in1->getOneBundle());
    assert(bundle);

    assert(trans->inflight_bundles.count(bundle) == 0);
    trans->inflight_bundles.insert(bundle);
    lock.unlock();

    // go through Records w/ iterator, which auto skips the "masked" records.
    for (auto && it = bundle->begin(); it != bundle->end(); ++it) {
        if (!(TransformT::do_select(*it)))
          it.mask();
//        else
//          E("--- leave one unmasked --- ");
    }

    // deposit the output Bundle to the output PValue
    auto out = trans->getFirstOutput();
    assert(out);

    lock.lock(); // protect against concurrent watermk refresh
    // output the same bundle, whose bitmap selector has been updated in place.
    out->depositOneBundle(bundle);

    // now the input bundle is gone and output bundle is commited.
    assert(trans->inflight_bundles.count(bundle) == 1);
    trans->inflight_bundles.erase(bundle);
    lock.unlock();

    c->SpawnConsumer(out);
  }
};

#endif /* SIMPLESELECTEVALUATOR_H */
