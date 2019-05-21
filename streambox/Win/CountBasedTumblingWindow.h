#ifndef COUNT_BASED_TUMBLING_WINDOW_H
#define COUNT_BASED_TUMBLING_WINDOW_H

#include "core/Transforms.h"

template <typename InputT, template<class> class InputBundle>
class CountBasedTumblingWindow
    : public PTransform {

public:
    // const uint64_t maxCount;
    const uint64_t MAX_COUNT;
    std::atomic<uint64_t> counter;
    std::atomic<uint64_t> win_id;
public:
    CountBasedTumblingWindow(const std::string name, const uint64_t maxCount)
        : PTransform(name), MAX_COUNT(maxCount), counter(0), win_id(10) {
    }

    void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
                       shared_ptr<BundleBase> bundle_ptr = nullptr) override;
};

#endif // COUNT_BASED_TUMBLING_WINDOW_H
