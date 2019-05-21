#ifndef PERBUNDLECONTEXT_H
#define PERBUNDLECONTEXT_H

#include "Values.h"
//#include "ParDo.h"
// the state for evaluating one bundle. per-thread.
template<class InputT, class OutputT>
class PerBundleContext : public DoFn<InputT, OutputT>::Context {
public:
  // context-wide staging buffer for a bundle.
  // make it shared_ptr because it may be
  // potentially sent to multiple downstream PValues
  shared_ptr<RecordBundle<OutputT>> outputBundle;
  void output (OutputT o) override {
    outputBundle->content.push_back(o);
  }

  PerBundleContext() :
    outputBundle(make_shared<RecordBundle<OutputT>>()) { }
};

#endif /* PERBUNDLECONTEXT_H */
