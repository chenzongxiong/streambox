#ifndef DUALINPUTTRANSFORMEVALUATOR_H
#define DUALINPUTTRANSFORMEVALUATOR_H

// this is pretty much for "join"
template <class TransformT,
  class LeftInputBundleT, class RightInputBundleT,
  class OutputBundleT>
class DualInputTransformEvaluator : public TransformEvaulator<TransformT> {

public:
  // to be implemented by subclass
  // @return: whether a consumer should be scheduled to consume
  // the output bundle
  virtual bool evaluateDualInput(TransformT* trans,
      shared_ptr<LeftInputBundleT> left_input_bundle,
      shared_ptr<RightInputBundleT> right_input_bundle,
      shared_ptr<OutputBundleT> output_bundle) = 0;

 /* get two input bundles, one from each input; output one bundle.
  *
  * note that either input may not always have a bundle. One reason is
  * that the prodcuer of one input schedules this evaluator, without
  * knowing if the other input has any deposited bundle.
  */
  void evaluate(TransformT* trans, EvaluationBundleContext* c) override {
    PValue* left_in = trans->getFirstInput();
    PValue* right_in = trans->getSecondInput();
    assert(left_in && right_in);

    // XXX use separate lock for each input stream?
    unique_lock<mutex> lock(trans->mtx_watermark);
    auto left_input_bundle = dynamic_pointer_cast<LeftInputBundleT>(
        left_in->getOneBundle(-2));
    auto right_input_bundle = dynamic_pointer_cast<RightInputBundleT>(
        right_in->getOneBundle(-2));

//    assert(left_input_bundle && right_input_bundle);

    W("left bundle %s; right bundle %s",
        left_input_bundle ? to_simplest_string(left_input_bundle->min_ts).c_str() : "(null)",
        right_input_bundle ? to_simplest_string(right_input_bundle->min_ts).c_str() : "(null)");

    if (left_input_bundle) {
      assert(trans->dual_inflight_bundles[0].count(left_input_bundle) == 0);
      trans->dual_inflight_bundles[0].insert(left_input_bundle);
    }

    if (right_input_bundle) {
      assert(trans->dual_inflight_bundles[1].count(right_input_bundle) == 0);
      trans->dual_inflight_bundles[1].insert(right_input_bundle);
    }

    lock.unlock();

    auto output_bundle = make_shared<OutputBundleT>(128, this->_node);

    if (evaluateDualInput(trans, left_input_bundle, right_input_bundle,
        output_bundle)) {
      auto out = trans->getFirstOutput();
      assert(out);

      lock.lock(); // protect against concurrent watermk refresh
      out->depositOneBundle(output_bundle, this->_node);

      // now the input bundle is gone and output bundle is commited.
      if (left_input_bundle) {
        assert(trans->dual_inflight_bundles[0].count(left_input_bundle) == 1);
        trans->dual_inflight_bundles[0].erase(left_input_bundle);
      }
      if (right_input_bundle) {
        assert(trans->dual_inflight_bundles[1].count(right_input_bundle) == 1);
        trans->dual_inflight_bundles[1].erase(right_input_bundle);
      }
      lock.unlock();
      c->SpawnConsumer(out, this->_node);  /* often the downstream is sink... */
    } else { // output bundle has no records
        lock.lock();
        if (left_input_bundle) {
          assert(trans->dual_inflight_bundles[0].count(left_input_bundle) == 1);
          trans->dual_inflight_bundles[0].erase(left_input_bundle);
        }
        if (right_input_bundle) {
          assert(trans->dual_inflight_bundles[1].count(right_input_bundle) == 1);
          trans->dual_inflight_bundles[1].erase(right_input_bundle);
        }
        lock.unlock();
    }
  }

  virtual ~DualInputTransformEvaluator() { }
};

#endif /* DUALINPUTTRANSFORMEVALUATOR_H */
