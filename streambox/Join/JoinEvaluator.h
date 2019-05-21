#ifndef JOINEVALUATOR_H
#define JOINEVALUATOR_H

template <class KVPair>
class JoinEvaluator : public DualInputTransformEvaluator<
		      Join<KVPair>,
		      RecordBitmapBundle<KVPair>, RecordBitmapBundle<KVPair>,
		      RecordBitmapBundle<pair<decltype(KVPair::first),
		      vector<decltype(KVPair::second)>>>
				    > {

					    using K = decltype(KVPair::first);
					    using V = decltype(KVPair::second);
					    using TransformT = Join<KVPair>;
					    using RecordKV = Record<KVPair>;
					    using RecordJoined = Record<pair<K, vector<V>>>;
					    using InputBundleT = RecordBitmapBundle<KVPair>;
					    using OutputBundleT = RecordBitmapBundle<pair<K, vector<V>>>;

					    public:
					    JoinEvaluator(int node) {
						    ASSERT_VALID_NUMA_NODE(node);
						    this->_node = node;
					    }

					    private:
					    // @i: which side of bundle we're consuming. left 0, right 1
					    void consumeOneBundle(TransformT* trans, shared_ptr<InputBundleT> bundle,
							    shared_ptr<OutputBundleT> output_bundle, int i) {
						    assert(trans);

						    //    boost::unique_lock<boost::shared_mutex> this_writer_lock(trans->_mutex[i]);
						    //    boost::shared_lock<boost::shared_mutex> other_reader_lock(trans->_mutex[1-i]);

#if 0 //hym: NOT use lock
						    boost::unique_lock<boost::shared_mutex> left_writer_lock(trans->_mutex[0], boost::defer_lock);
						    boost::unique_lock<boost::shared_mutex> right_writer_lock(trans->_mutex[1], boost::defer_lock);
						    boost::shared_lock<boost::shared_mutex> left_reader_lock(trans->_mutex[0], boost::defer_lock);
						    boost::shared_lock<boost::shared_mutex> right_reader_lock(trans->_mutex[1], boost::defer_lock);

						    if (i == 0) {
							    W("0:waiting for left writer lock...");
							    left_writer_lock.lock();
							    W("0:waiting for right reader lock...");
							    right_reader_lock.lock();
						    } else {
							    W("1:waiting for left reader lock...");
							    left_reader_lock.lock();
							    W("1:waiting for right writer lock...");
							    right_writer_lock.lock();
						    }
#endif //hym: Not use lock

						    W("got locks. bundle has %lu records", bundle->content.size());
						    RecordKV kv;
						    for (auto && it = bundle->begin(); it != bundle->end(); ++it) {
							    // search the other side. if match, emit the joined pair
							    if (trans->search(*it, &kv, 1-i)) {
								    output_bundle->add_record(RecordJoined(
											    make_pair(kv.data.first, vector<V>({(*it).data.second, kv.data.second})),
											    kv.ts
											    ));
							    }
							    // if there's match, shall we still add the record to this side hashtable?
							    trans->add_record(*it, i);
						    }

						    //    other_reader_lock.unlock();
						    //    this_writer_lock.unlock();

#if 0 //hym: NOT use lock
						    if (i == 0) {
							    right_reader_lock.unlock();
							    left_writer_lock.unlock();
						    } else {
							    right_writer_lock.unlock();
							    left_reader_lock.unlock();
						    }
#endif //hym: NOT use lock
					    }

					    bool evaluateDualInput(TransformT* trans,
							    shared_ptr<InputBundleT> left_input_bundle,
							    shared_ptr<InputBundleT> right_input_bundle,
							    shared_ptr<OutputBundleT> output_bundle) override {

						    assert(trans);
#if 0
						    // consume the left input bundle
						    for (auto && it = left_input_bundle->begin();
								    it != left_input_bundle->end(); ++it) {
							    trans->add_record(*it);
						    }

						    // then the right input bundle
						    RecordKV kv;
						    for (auto && it = right_input_bundle->begin();
								    it != right_input_bundle->end(); ++it) {
							    if (trans->search(*it, &kv)) {
								    output_bundle->add_record(RecordJoined(kv.ts,
											    make_pair(kv.data.first, vector<V>({(*it).data.second, kv.data.second}
													    ))));
							    }
						    }
#endif
						    VV("to do bundle 0");
						    if (left_input_bundle)
							    consumeOneBundle(trans, left_input_bundle, output_bundle, 0);

						    VV("done bundle 0; to do bundle 1");
						    if (right_input_bundle)
							    consumeOneBundle(trans, right_input_bundle, output_bundle, 1);
						    VV("done bundle 1");

						    if (output_bundle->content.size())
							    return true;
						    else
							    return false;
					    }
				    };

#endif /* JOINEVALUATOR_H */
