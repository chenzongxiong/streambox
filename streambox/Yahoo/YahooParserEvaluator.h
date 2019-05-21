#ifndef YAHOOPARSEREVALUATOR_H
#define YAHOOPARSEREVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Yahoo/YahooRecord.h"
#include "Yahoo/YahooParser.h"

/* convert a stream of records to a stream of <YahooRecord> records */
template <typename InputT, typename OutputT, template<class> class BundleT_>
class YahooParserEvaluator
    : public SingleInputTransformEvaluator<YahooParser<InputT,OutputT>,
      BundleT_<InputT>, BundleT_<OutputT>> {

  using TransformT = YahooParser<InputT,OutputT>;
  using InputBundleT = BundleT_<InputT>;
  using OutputBundleT = BundleT_<OutputT>;

public:
  bool evaluateSingleInput (TransformT* trans,
        shared_ptr<InputBundleT> input_bundle,
        shared_ptr<OutputBundleT> output_bundle) override {

      // simply return false to block all tuples from unbound-inmem to yahoo-parser
      // that means the rest parts of pipeline don't have any loads.
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        auto out_record = trans->do_map(*it);
        output_bundle->add_record(out_record);
    }
    trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
    trans->byte_counter_.fetch_add(78*input_bundle->size(), std::memory_order_relaxed);

    return true;
    // return false;
  }

  YahooParserEvaluator(int node)
  	: SingleInputTransformEvaluator<TransformT,
  	  			InputBundleT, OutputBundleT>(node) { }

};

#endif /* YAHOOPARSEREVALUATOR_H */
