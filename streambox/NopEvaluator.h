#ifndef NOPEVALUATOR_H
#define NOPEVALUATOR_H

// simply move over all bundles of the ingress to the egress.
template<typename T>
class NopEvaluator : public TransformEvaulator<PTransformNop<T>> {
  using Transform = PTransformNop<T>;

public:
  NopEvaluator(int node) { this->_node = 0; /* don't care */ }

  void evaluate(Transform* trans, EvaluationBundleContext* c) override  {

    PValue* input = trans->getFirstInput();
    PValue* out = trans->getFirstOutput();
    assert(input && out);

    while (auto bundle = input->getOneBundle(this->_node)) {
        out->depositOneBundle(bundle, this->_node);
//        c->ScheduleConsumer(out);
        c->SpawnConsumer(out, this->_node);
    }
  }
};

#endif /* NOPEVALUATOR_H */
