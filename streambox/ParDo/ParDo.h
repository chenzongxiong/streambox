// xzl: we simplify things. just one @ParDo instead of @Bound and @UnBound
// a very simple ParDo transform.

template <class InputT, class OutputT>
class ParDoEvaluator;

template <class InputT, class OutputT>
class ParDo : public PTransform {
//class ParDo : public PTransform <ParDo> {
private:
  DoFn<InputT, OutputT>* _fn;

public:
  ParDo(string name, DoFn<InputT, OutputT>* fn) :
    PTransform(name), _fn(fn) { }

  DoFn<InputT, OutputT>* getFn() {
    return _fn;
  }

  // XXX dirty: makes evaluator dispatch eaiser
  ParDoEvaluator<InputT, OutputT>* CreateEvaluator();
};
