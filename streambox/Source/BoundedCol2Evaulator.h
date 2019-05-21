#ifndef BOUNDEDCOL2EVAULATOR_H
#define BOUNDEDCOL2EVAULATOR_H

// Read in the list of tuples; internally convert them to column format.
//
// T: the element type of the source
// reference design: Read.java registerDefaultTransformEvaluator()
// save the results as a list of columns
template<class T1, class T2    // types of all columns, in order
//  class Header = Header2<T1, T2>,
//  class RawTuple = std::tuple<T1,T2>
>
class BoundedCol2Evaulator
: public TransformEvaulator<Bounded<std::tuple<T1,T2>>> {

	using Header = Header2<T1, T2>;
	using RawTuple = std::tuple<T1,T2>;

	public:
	BoundedCol2Evaulator(int node) { this->_node = node; }

	void evaluate(Bounded<RawTuple>* t, EvaluationBundleContext* c) override {
		// we already have the source
		BoundedSource<RawTuple> *s = t->_source;

		// containers for columns
		auto col1 = make_shared<vector<T1>>();
		auto col2 = make_shared<vector<T2>>();

		// have a list internal to pipeline runner
		// construct a reader every time?
		auto reader = s->createReader(NULL);
		bool available;
		for (available = reader->start(); available; available =
				reader->advance()) {
			RawTuple t = reader->getCurrent();
			printf("%s: get value (%ld %.2f) ts %s\n",
					__func__,
					std::get<0>(t), std::get<1>(t),
					to_simple_string(reader->getCurrentTimestamp()).c_str());

			col1->push_back(std::get<0>(t));
			col2->push_back(std::get<1>(t));

			// also init the ts columns?
		}

		// hopefully this is a space-efficient impl.
		auto selector = make_shared<vector<bool>>(col1->size(), true);
		auto cols = make_tuple(col1, col2);

		//    auto h = new Header();
		//    h->cols = cols;  // copy
		//    h->selector = selector;

		auto bundle = make_shared<BundleCol2<T1,T2>>();

		bundle->cols = cols; // copy. each column pointer is also copied
		bundle->selector = selector;

		auto out = t->getFirstOutput();
		assert(out);

		out->depositOneBundle(bundle);
		c->ScheduleConsumer(out);
	}

	~BoundedCol2Evaulator() {
	}
};

#if 0
template<class T1, class T2,    // types of all columns, in order
	class Header = Header2<T1, T2>,
	class Transform = PrinterCol2<T1, T2>,
	class RawTuple = std::tuple<T1,T2>>
	class PrinterCol2Evaulator
	: public TransformEvaulator<Bounded<RawTuple>> {
		public:
			void evaluate(Transform* trans, EvaluationContext* c = NULL) override {
				PValue *input = c->getInput(trans);

				assert(input);

				// the underlying content
				auto store = (Header*) (c->getPValue(c->getInput(trans)));

				if (!store) {
					printf("cast contents (%lx) seemed to fail:%s\n",
							(unsigned long)(input),
							typeid(store).name());
				}
				assert(store);

				// XXX why can't we pass in *(store->selector)[i] as
				// reference?
				// also, std::for_each cannot iterate over two vectors
				// at the same time.

				// go over column 1 (allow update in place)
				{
					auto col = std::get<0>(store->cols);
					for (int i = 0; i < col->size(); i++) {
						cout << col[i];
					}
					cout << endl;
				}

				{
					// go over column 2 (allow update in place)
					auto col = std::get<1>(store->cols);
					for (int i = 0; i < col->size(); i++) {
						cout << col[i];
					}
					cout << endl;
				}
			}

			~PrinterCol2Evaulator() {
			}
	};
#endif

#endif /* BOUNDEDCOL2EVAULATOR_H */
