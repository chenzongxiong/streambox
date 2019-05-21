#ifndef SIMPLEMAPPEREVAULATOR_H
#define SIMPLEMAPPEREVAULATOR_H

#include "core/SingleInputTransformEvaluator.h"
#include "Mapper/SimpleMapper.h"

/* convert a stream of records to a stream of <KV> records */
template <typename InputT, typename KVPair>
class SimpleMapperEvaluator
    : public SingleInputTransformEvaluator<SimpleMapper<InputT,KVPair>,
      RecordBitmapBundle<InputT>, RecordBitmapBundle<KVPair>> {

  using TransformT = SimpleMapper<InputT,KVPair>;
  using InputBundleT = RecordBitmapBundle<InputT>;
  using OutputBundleT = RecordBitmapBundle<KVPair>;

public:
#if 0
  void evaluate(TransformT* trans, EvaluationBundleContext* c) {

    PValue* in1 = trans->getFirstInput();
    assert(in1);

    // get one pending bundle from the input.
    // this will update input's min_ts. note that the bundle is still
    // a "pending" work in the current transform

    unique_lock<mutex> lock(trans->mtx_watermark);

    auto bundle = \
        dynamic_pointer_cast<InputBundleT>(in1->getOneBundle());
    assert(bundle);

    assert(trans->inflight_bundles.count(bundle) == 0);
    trans->inflight_bundles.insert(bundle);
    lock.unlock();

    auto output_bundle = make_shared<OutputBundleT>();

    // go through Records in input bundle (the iterator automatically
    // skips "masked" Records.
    for (auto && it = bundle->begin(); it != bundle->end(); ++it) {
        auto out_record = trans->do_map(*it);
        output_bundle->add_record(out_record);
//        assert(out_record.data.first == 12 && out_record.data.second == 1234);
    }

    // deposit the output Bundle to the output PValue
    auto out = trans->getFirstOutput();
    assert(out);

    lock.lock(); // protect against concurrent watermk refresh
    out->depositOneBundle(output_bundle);

    // now the input bundle is gone and output bundle is commited.
    assert(trans->inflight_bundles.count(bundle) == 1);
    trans->inflight_bundles.erase(bundle);
    lock.unlock();

    c->SpawnConsumer(out);
  }
#endif

  bool evaluateSingleInput (TransformT* trans,
        shared_ptr<InputBundleT> input_bundle,
        shared_ptr<OutputBundleT> output_bundle) override {

    //std::cout << "SimpleMapperEvaluator.h line 68: evaluateSingleInput" << std::endl;
    //std::cout << "SimpleMapperEvaluator.h line 69: side_info is " << trans->side_info << std::endl;
    // go through Records in input bundle (the iterator automatically
    // skips "masked" Records.
    for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
        /*
	//std::cout << "1" << std::endl;
	auto out_record = trans->do_map(*it);
	//std::cout << "2" << std::endl;
        output_bundle->add_record(out_record);
	//std::cout << "3" << std::endl;
	*/
	//hym: to verify join result
	auto out_record = trans->do_map(*it, trans->get_side_info());
	//std::cout << "2" << std::endl;
        output_bundle->add_record(out_record);
    }

    return true;
  }

  SimpleMapperEvaluator(int node)
  	: SingleInputTransformEvaluator<TransformT,
  	  			InputBundleT, OutputBundleT>(node) { }
#if 0
  	ASSERT_VALID_NUMA_NODE(node);
  	this->_node = node;
  }
#endif

};

#endif /* SIMPLEMAPPEREVAULATOR_H */
