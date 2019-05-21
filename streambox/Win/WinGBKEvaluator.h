#ifndef WINGBK_EVAL_H
#define WINGBK_EVAL_H

#include "core/SingleInputTransformEvaluator.h"
#include "Win/WinGBK.h"
#include "Values.h"
#include <map>

extern std::map<Window, ptime, Window> window_keeper;

/* FixedWindowInto + GroupByKey. Stateless */
template <class KVPair,
				template<class> class InputBundleT_,
				/* For output kvpairs, WinKeyFrag can be specialized based on key/val distribution */
				template<class> class WindowKeyedFragmentT
				>
class WinGBKEvaluator
    : public SingleInputTransformEvaluator<WinGBK<KVPair, InputBundleT_, WindowKeyedFragmentT>,
      	InputBundleT_<KVPair>, WindowsKeyedBundle<KVPair, WindowKeyedFragmentT>>
{
  using TransformT = WinGBK<KVPair, InputBundleT_, WindowKeyedFragmentT>;
  using InputBundleT = InputBundleT_<KVPair>;
//  using InputBundleT = RecordBitmapBundle<KVPair>;
  using OutputBundleT = WindowsKeyedBundle<KVPair,WindowKeyedFragmentT>; /* this is fixed */

public:

#if 0
  void evaluate(TransformT* trans, EvaluationBundleContext* c) {

    PValue* in1 = trans->getFirstInput();
    assert(in1);

    // get one pending bundle from the input.
    // this will update input's min_ts. note that the bundle is still
    // a "pending" work in the current transform

    unique_lock<mutex> lock(trans->mtx_watermark);

    auto input_bundle = \
        dynamic_pointer_cast<InputBundleT>(in1->getOneBundle());
    assert(input_bundle);

    assert(trans->inflight_bundles.count(input_bundle) == 0);
    trans->inflight_bundles.insert(input_bundle);
    lock.unlock();

    auto output_bundle = make_shared<OutputBundleT>();

    // go through Records in input bundle (the iterator automatically
    // skips "masked" Records.
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {

        // the time offset within a window
        long offset = ((*it).ts - trans->start).total_microseconds() \
            % (trans->window_size).total_microseconds();

//        auto rec = (*it);
//        assert(rec.data.first == 12 && rec.data.second == 1234);

        // add_value() will organize records by windows and by keys
        // the bundle's, and the underlying v container's min_ts
        // will also be updated.
        output_bundle->add_record(
              Window((*it).ts - microseconds(offset), trans->window_size),
              *it);
    }

    // deposit the output Bundle to the output PValue
    auto out = trans->getFirstOutput();
    assert(out);

    lock.lock(); // protect against concurrent watermk refresh
    out->depositOneBundle(output_bundle);

    // now the input bundle is gone and output bundle is commited.
    assert(trans->inflight_bundles.count(input_bundle) == 1);
    trans->inflight_bundles.erase(input_bundle);
    lock.unlock();

    c->SpawnConsumer(out);
  }
#endif

    bool evaluateSingleInput(TransformT* trans,
        shared_ptr<InputBundleT> input_bundle,
        shared_ptr<OutputBundleT> output_bundle) override
    {

      bool ret = false;


      /* Since we need to pass in a window in adding each record, we reuse
       * this one instead of constructing a new one each time.
       */
      Window ww;
      ww.duration = trans->window_size;

      /* Go through Records in input bundle (the iterator automatically
         skips "masked" Records.
         When records are fine-grained (e.g. words), this can be hot. */
      ptime max_ts = min_date_time;
      for (auto && it = input_bundle->begin();
          it != input_bundle->end(); ++it) {
      	/* when input_bundle is RecordBundle, @it is a vanilla vector iterator */

          // the time offset within a window, in ms.
          long offset_ms = ((*it).ts - trans->start).total_milliseconds()
              % (trans->window_size).total_milliseconds();
          if (max_ts < it->ts) {
              max_ts = it->ts;
          }
          // add_record() will organize records by windows and by keys
          // the bundle's, and the underlying v container's min_ts
          // will also be updated.
//          output_bundle->add_record(
//              	Window((*it).ts - milliseconds(offset), trans->window_size),
//              	*it);
          ww.start = (*it).ts - milliseconds(offset_ms) - Window::epoch;
          output_bundle->add_record(ww, *it); /* hot. WindowsKeyedBundle op inline? */

          ret = true;
      }
      auto it = window_keeper.find(ww);
      if (it == window_keeper.end()) {
          window_keeper[ww] = max_ts;
      } else {
          long diff_ = (it->second - max_ts).total_milliseconds();
          if (diff_ < 0) {
              it->second = max_ts;
          }
      }

      // boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
      // long diff = (now - max_ts).total_milliseconds();
      // if (diff < 0) {
      //     assert(false && "WinGBK ERROR");
      // }
      // return false;
      return ret;
    }

    WinGBKEvaluator(int node)
  	: SingleInputTransformEvaluator<TransformT,
  	  			InputBundleT, OutputBundleT>(node) { }
};

#endif // WINGBK_EVAL_H
