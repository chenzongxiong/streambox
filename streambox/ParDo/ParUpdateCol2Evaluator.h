#ifndef PARUPDATECOL2EVALUATOR_H
#define PARUPDATECOL2EVALUATOR_H

// Update the bundle, as two columns, in place.
// does not produce new columns
template <class T1, class T2
//  // Fns on each column
//  class DoFn1 = DoFnSimple<T1>,
//  class DoFn2 = DoFnSimple<T2>,
//  // type of transform to be evaluated
//  class Transform = ParUpdateCol2<T1, T2>
>
//  // type of the columar store
//  class Header = Header2<T1, T2>,
//  int N = 2>
class ParUpdateCol2Evaluator
    : public TransformEvaulator<ParUpdateCol2<T1, T2>> {

  // Fns on each column
  using DoFn1 = DoFnSimple<T1>;
  using DoFn2 = DoFnSimple<T2>;
  // type of transform to be evaluated
  using Transform = ParUpdateCol2<T1, T2>;

public:
  ParUpdateCol2Evaluator(int node) { this->_node = node; }

  void evaluate(Transform* trans, EvaluationBundleContext *c) {
    PValue *input = trans->getFirstInput();
    PValue *output = trans->getFirstOutput();
    assert(input && output);

    auto bundle_in = \
        dynamic_pointer_cast<BundleCol2<T1,T2>>(input->getOneBundle());
    assert(bundle_in);

#if 0
    // go over column 1
    auto col = std::get<1>(store);
    for (auto&& v : *col) {
        std::get<0>(trans->_fns)(nullptr, v);
    }

    // go over column 2
    auto col = std::get<2>(store);
    for (auto&& v : *col) {
        std::get<1>(trans->_fns)(nullptr, v);
    }
#endif

    // std::for_each cannot iterate over two vectors
    // at the same time.
    // in processing each element, we pass in ptr to the whole column.

    // go over column 1 (allow update in place)
    {
      auto col = std::get<0>(bundle_in->cols);
      for (unsigned int i = 0; i < col->size(); i++) {
          if ((*bundle_in->selector)[i])
            trans->fn1->processElement(nullptr,
                col, bundle_in->selector, i);
      }
    }

    {
      // go over column 2 (allow update in place)
      auto col = std::get<1>(bundle_in->cols);
      for (unsigned int i = 0; i < col->size(); i++) {
          if ((*bundle_in->selector)[i])
            trans->fn2->processElement(nullptr,
                            col, bundle_in->selector, i);
      }
    }

    auto bundle_out = make_shared<BundleCol2<T1,T2>>();
    bundle_out->selector = bundle_in->selector;
    bundle_out->cols = bundle_in->cols; // just copy

    output->depositOneBundle(bundle_out);
    c->RunConsumer(output);
  }

};

#endif /* PARUPDATECOL2EVALUATOR_H */
