#ifndef NYTPARSEREVALUATOR_H
#define NYTPARSEREVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "NewYorkTaxi/NYTRecord.h"
#include "NewYorkTaxi/NYTParser.h"

/* convert a stream of records to a stream of <NYTRecord> records */
template <typename InputT, typename OutputT, template<class> class BundleT_>
class NYTParserEvaluator
    : public SingleInputTransformEvaluator<NYTParser<InputT,OutputT, BundleT_>,
                                           BundleT_<InputT>, BundleT_<OutputT>> {

    using TransformT = NYTParser<InputT, OutputT, BundleT_>;
    using InputBundleT = BundleT_<InputT>;
    using OutputBundleT = BundleT_<OutputT>;

public:
    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {

        // simply return false to block all tuples from unbound-inmem to nyt-parser
        // that means the rest parts of pipeline don't have any loads.
        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            auto out_record = trans->do_map(*it);
            // out_record.data.print();
            output_bundle->add_record(out_record);
        }
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        trans->byte_counter_.fetch_add(sizeof(NYTRecord)*input_bundle->size(), std::memory_order_relaxed);

        return true;
    }

    NYTParserEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT,
                                        InputBundleT, OutputBundleT>(node) { }

};

#endif /* NYTPARSEREVALUATOR_H */
