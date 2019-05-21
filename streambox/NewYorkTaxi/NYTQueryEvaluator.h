#ifndef NYTQUERYEVALUATOR_H
#define NYTQUERYEVALUATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "NewYorkTaxi/NYTRecord.h"
#include "NewYorkTaxi/NYTQuery.h"

template <typename InputT, typename OutputT, template<class> class BundleT_>
class NYTQueryEvaluator
    : public SingleInputTransformEvaluator<NYTQuery<InputT, OutputT, BundleT_>,
                                           BundleT_<InputT>, BundleT_<OutputT>> {

    using TransformT = NYTQuery<InputT, OutputT, BundleT_>;
    using InputBundleT = BundleT_<InputT>;
    using OutputBundleT = BundleT_<OutputT>;

public:
    bool evaluateSingleInput (TransformT* trans,
                              shared_ptr<InputBundleT> input_bundle,
                              shared_ptr<OutputBundleT> output_bundle) override {
        //shared_ptr<RecordBundle<std::pair<uint64_t, uint64_t>>> output_bundle) override {
        // std::atomic<uint64_t> trip_distance;

        // uint64_t trip_distance = (uint64_t)(it->data.trip_distance);
        for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
            // auto win = it->first;
            // auto val = it->second;

            // std::cout << "win.start: " << win.window_start()
            //           << "\twin.end: " << win.window_end()
            //           << "\tval: " << val->vals.size()
            //           << std::endl;

            // for (auto v_it = val->vals.begin(); v_it != val->vals.end(); v_it ++) {
            //     if (trans->do_map(*v_it)) {
            //         output_bundle->add_record(win, *v_it);
            //     }
            // }
            if (trans->do_map(*it)) {
                uint64_t longitude = (uint64_t)((it->data.pickup_long+750)*100.0);
                uint64_t latitude = (uint64_t)((it->data.pickup_lang+1400)*100.0);
                uint64_t region = (longitude + 12345*latitude) % 100;
                // cout << "region: " << region << endl;
                output_bundle->emplace_record(std::pair<uint64_t, uint64_t>(region, (uint64_t)it->data.trip_distance), it->ts);
                // output_bundle->add_record(*it);
            }
        }
        // std::cout << "input_bundle.size: " << input_bundle->size()
        //           << ", output_bundle.size: " << output_bundle->size()
        //           << std::endl;
        trans->record_counter_.fetch_add(input_bundle->size(), std::memory_order_relaxed);
        // trans->byte_counter_.fetch_add(sizeof(NYTRecord)*input_bundle->size(), std::memory_order_relaxed);
        return true;

    }

    NYTQueryEvaluator(int node)
        : SingleInputTransformEvaluator<TransformT,
                                        InputBundleT, OutputBundleT>(node) { }

};

#endif //NYTQUERYEVALUATOR_H
