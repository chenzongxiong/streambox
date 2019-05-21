#ifndef MINIMALPROCESSCONTEXT_H
#define MINIMALPROCESSCONTEXT_H

// this is for evaluating one particular element (overkill??)
template<class InputT, class OutputT>
class MinimalProcessContext : public DoFn<InputT, OutputT>::ProcessContext {
public:
  InputT _content;
  MinimalContext<InputT, OutputT>* _ctx; // also refer to the xform eval ctx

  MinimalProcessContext(InputT content,
      MinimalContext<InputT, OutputT> *ctx)
    : _content(content), _ctx(ctx) { }

  InputT element() override { return _content; }
  ptime timestamp() override { return min_date_time; }
  BoundedWindow* window() override { return nullptr; }
  void output(OutputT o) override {
    assert(_ctx);
    _ctx->output(o);
  }
};

#endif /* MINIMALPROCESSCONTEXT_H */
