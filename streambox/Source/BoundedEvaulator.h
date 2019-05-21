#ifndef BOUNDEDEVAULATOR_H
#define BOUNDEDEVAULATOR_H

// T: the element type of the source
// only supports one NUMA node
template<class T>
class BoundedEvaulator: public TransformEvaulator<Bounded<T>> {
	public:

		const int numa_node_id = 0;

		void evaluate(Bounded<T>* t, EvaluationBundleContext* c) override {
			// we already have the source
			BoundedSource<T> *s = t->_source;

			// the output PValue -- XXX assume there's only one source
			auto out = t->getFirstOutput();
			assert(out);

			// create a new Bundle
			//    auto bundle = new Bundle<T>();
			shared_ptr<RecordBundle<T>> bundle(make_shared<RecordBundle<T>>(128, numa_node_id));
			assert(bundle);

			// construct a reader every time?
			auto reader = s->createReader(NULL);
			for (bool available = reader->start(); available; available =
					reader->advance()) {
				VV("%s: get value (type %s), ts %s",
						__func__,
						typeid(reader->getCurrent()).name(),
						to_simple_string(reader->getCurrentTimestamp()).c_str());

				// xzl: XXX very inefficient?
				bundle->content.push_back(reader->getCurrent());

				if (bundle->content.size() == max_bundle_size) {
					I("emit one bundle to the output PValue ...");
					out->depositOneBundle(bundle);
					c->SpawnConsumer(out, numa_node_id);
					// create & produce to a new bundle
					//          bundle = new Bundle<T>();
					bundle = make_shared<RecordBundle<T>>();
					assert(bundle);
				}
			}

			if (bundle->content.size()) {
				I("emit one bundle to the output PValue (node 0)...");
				out->depositOneBundle(bundle, numa_node_id);
				//        c->ScheduleConsumer(out);
				c->SpawnConsumer(out, numa_node_id);
			}
		}

		~BoundedEvaulator() {
		}

		BoundedEvaulator(unsigned int max_bundle_size = 64)
			: max_bundle_size (max_bundle_size) { }

	private:
		unsigned int max_bundle_size;
};


#endif /* BOUNDEDEVAULATOR_H */
