#ifndef MINIMALCONTEXT_H
#define MINIMALCONTEXT_H

#include "ParDo/ParDo.h"
// Contexts for evaluating ParDo. They are all templated.
// XXX they should be part of the global context (for eval the pipeline)

template<class InputT, class OutputT>
class MinimalContext : public DoFn<InputT, OutputT>::Context {
public:
  vector<OutputT>* _buf;       // context-wide staging buffer
  void output (OutputT o) override {
    _buf->push_back(o);
  }

  MinimalContext() : _buf(new vector<OutputT>) { }
};

#endif /* MINIMALCONTEXT_H */
