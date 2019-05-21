#ifndef FIXEDWINDOWINTO_H
#define FIXEDWINDOWINTO_H

#include "core/Transforms.h"
#include "boost/date_time/posix_time/ptime.hpp"

template <typename InputT, template<class> class InputBundle>
class FixedWindowInto : public PTransform {
public:
  const boost::posix_time::time_duration window_size;
  const ptime start; // the starting point of windowing.

  FixedWindowInto(string name,
      boost::posix_time::time_duration window_size,
      ptime start = Window::epoch)
    : PTransform(name), window_size(window_size), start(start) { }

  void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> bundle_ptr = nullptr) override;

};

#endif /* FIXEDWINDOWINTO_H */
