#ifndef DOFN_H_
#define DOFN_H_

#include "boost/date_time/posix_time/posix_time.hpp"

template <class InputT, class OutputT>
class DoFn {

	//// inner classes
	public:

		class Context {
			public:
				virtual void output(OutputT o) = 0;
				virtual ~Context() { }
		};

#if 1 // xzl: XXX get rid of this
		//////////////////////////////////////////////////////

		class ProcessContext : public Context {
			public:
				virtual InputT element() = 0;
				virtual ptime timestamp() = 0;
				virtual BoundedWindow* window() = 0;
				virtual ~ProcessContext() { }
		};
#endif

		//////////////////////////////////////////////////////

	public:
		DoFn() { }

		virtual void startBundle(Context* c) = 0;
		virtual void processElement(ProcessContext *c) = 0;
		virtual void finishBundle(Context* c) = 0;
		virtual void processElementSimple(Context *c, Record<InputT>& element) {
			assert(0); // must be implemented by subclasses
		}

		virtual ~DoFn() { }

};

// xzl -- a simple version: no context.
// allow update element in place.
template <class InputT>
class DoFnSimple {

	//// inner classes
	public:

		class Context {
			public:
				virtual ~Context() { }
		};

		//////////////////////////////////////////////////////

	public:
		DoFnSimple() { }

		virtual void startBundle(Context* c) = 0;

		// pass in reference: allow in-place update
		// since this is not temporal, we keep the timestamp, window... unchanged

		// to avoid the cost of maintaining shared_ptr in copy, we pass by
		// reference. however, we expect the function should never drop the reference
		// count (e.g. by calling reset()).
		virtual void processElement(Context *c,
				const shared_ptr<vector<InputT>>&,
				const shared_ptr<vector<bool>>& selector, int i) = 0;

#if 0
		using OutputBundleT = RecordBundle<InputT>;
		virtual unsigned int dofn(Record<InputT> const & in,
				OutputBundleT & ouptut) = 0;
#endif

		virtual void finishBundle(Context* c) = 0;

		virtual ~DoFnSimple() { }

};

/* even simpler. no context at all */
template <class InputT, class OutputT>
class DoFnSimpler {
	using OutputBundleT = RecordBundle<OutputT>;

	virtual unsigned int dofn(Record<InputT> const & in,
			OutputBundleT & ouptut) = 0;
};

/////////////////////////////////////////////////////////
// Trivial DoFns ...

// increments each long by one.
#include <vector>
using namespace std;
class DoFnPlus1: public DoFn<long, long> {
	private:
		vector<long> _elements;

	public:
		void startBundle(Context *c) override {
			_elements.clear();
		}
		void finishBundle(Context *c) override {
			for (long ele : _elements)
				c->output(ele);
		}

		void processElement(ProcessContext *c) override {
			long ele = c->element();
			printf("%s: element is %ld\n", __PRETTY_FUNCTION__, ele);
			_elements.push_back(ele + 1);
		}

		void processElementSimple(Context *c, Record<long>& t) override {
			c->output(t.data);
		}
};

// do that for each element of in a tuple
#include <tuple>
template <class T1, class T2,
	 class Tuple = std::tuple<T1,T2>>
	 class DoFnPlus1Tuple2 : public DoFn<Tuple,Tuple> {
		 private:
			 vector<Tuple> _elements;

		 public:
			 void startBundle(typename DoFn<Tuple,Tuple>::Context *c) override {
				 _elements.clear();
			 }
			 void finishBundle(typename DoFn<Tuple,Tuple>::Context *c) override {
				 for (Tuple ele : _elements)
					 c->output(ele);
			 }

			 void processElement(typename DoFn<Tuple,Tuple>::ProcessContext *c) override {
				 Tuple ele = c->element();
				 printf("%s: got one element \n", __PRETTY_FUNCTION__);
				 _elements.push_back(
						 make_tuple(get<0>(ele)+ 1, get<1>(ele)+1)
						 );
			 }

			 void processElementSimple(typename DoFn<Tuple,Tuple>::Context *c,
					 Record<Tuple>& t) override {
				 c->output(make_tuple(get<0>(t.data)+ 1, get<1>(t.data)+1));
			 }
	 };

template <class T>
class DoFnPrinter: public DoFn<T,T> {
	private:
	public:
		void startBundle(typename DoFn<T,T>::Context *c) override { }
		void finishBundle(typename DoFn<T,T>::Context *c) override { }

		void processElement(typename DoFn<T,T>::ProcessContext *c) override { }

		void processElementSimple(typename DoFn<T,T>::Context *c,
				Record<T>& t) override {
			cout << "===========" << t.ts << ": " << t.data << "=======" << endl;
		}
};

template <class T>
class DoFnSimplePlus1 : public DoFnSimple<T> {
	public:
		void startBundle(typename DoFnSimple<T>::Context *c) override { }

		void finishBundle(typename DoFnSimple<T>::Context *c) override { }

		void processElement(typename DoFnSimple<T>::Context *c,
				const shared_ptr<vector<T>>& element,
				const shared_ptr<vector<bool>>& mask, int i) override {
			(*element)[i] += 1;
		}
};

// for debugging
template <class T>
class DoFnSimplePrinter : public DoFnSimple<T> {
	public:
		void startBundle(typename DoFnSimple<T>::Context *c) override {
			cout << "startBundle" << endl;
		}

		void finishBundle(typename DoFnSimple<T>::Context *c) override {
			cout << "endBundle" << endl;
		}

		void processElement(typename DoFnSimple<T>::Context *c,
				const shared_ptr<vector<T>>& element,
				const shared_ptr<vector<bool>>& mask, int i) override {
			cout << "i: " << (*element)[i] << endl;
		}
};

template <class T>
class DoFnSimplerPlus1 : public DoFnSimpler<T, T> {
	using OutputBundleT = RecordBundle<T>;

	unsigned int dofn(Record<T> const & in,
			OutputBundleT & output_bundle) override {
		output_bundle.add_record(Record<T>(in.data + 1, in.ts));
		return 1;
	}
};

template <class T>
class DoFnSimplerPrinter : public DoFnSimpler<T, T> {
	using OutputBundleT = RecordBundle<T>;

	unsigned int dofn(Record<T> const & in,
			OutputBundleT & output_bundle) override {
		cout << "i: " << in.data << endl;
		return 0;
	}
};

#if 0
template <
class K, class V,
      class KV = std::pair<K,V>,
      class Record = Record<KV>
      >
      template <class T,
      class Record = Record<T>>
      class DoFnWindowInto : public DoFn<Record, Record> {
	      private:
		      WindowedBundle * bundle;
	      public:
		      void starBundle(typename DoFn<Record, Record>::Context *c)  override {
			      bundle = new WindowedBundle();
		      }
      };
#endif

#endif /* DOFN_H_ */
