#ifndef FIXEDWINDOWINTOEVALUATOR_H
#define FIXEDWINDOWINTOEVALUATOR_H

#include "core/EvaluationBundleContext.h"
#include "core/TransformEvaluator.h"
#include "core/SingleInputTransformEvaluator.h"
#include "Win/FixedWindowInto.h"

/* output bundle format is fixed */
template <class InputT, template<class> class InputBundleT_>
class FixedWindowIntoEvaluator
: public SingleInputTransformEvaluator<FixedWindowInto<InputT, InputBundleT_>,
	InputBundleT_<InputT>, WindowsBundle<InputT>> {

		using TransformT = FixedWindowInto<InputT, InputBundleT_>;
		//  using InputBundleT = RecordBitmapBundle<InputT>;
		using InputBundleT = InputBundleT_<InputT>;
		using OutputBundleT = WindowsBundle<InputT>;

		public:
		FixedWindowIntoEvaluator(int node)
			: SingleInputTransformEvaluator<TransformT, InputBundleT, OutputBundleT>(node) {
        }

		bool evaluateSingleInput (TransformT* trans,
				shared_ptr<InputBundleT> input_bundle,
				shared_ptr<OutputBundleT> output_bundle) override {

			/* go through Records w/ iterator, which deals with the "masked" records.
			 * let @WindowsBundle takes care of data layout.
			 */
            // processing time implementation
            boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();

            long diff = (now - trans->start).total_microseconds();
            long window_size = trans->window_size.total_microseconds();
            long offset = diff % window_size;

            // put all records inside this window
            for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
                output_bundle->add_record(Window(now - microseconds(offset), trans->window_size),
                                          *it);
            }

            // std::cout << "win.start: " << win.window_start()
            //           << "\twin.end: " << win.window_end()
            //           << "\tval: " << val->vals.size()
            //           << std::endl;

            // std::cout << "output_bundle->size: " << output_bundle->vals.size() << std::endl;

            // event time implementation
			// for (auto && it = input_bundle->begin(); it != input_bundle->end(); ++it) {
			// 	// the time offset within a window
			// 	long offset = ((*it).ts - trans->start).total_microseconds() \
			// 		      % (trans->window_size).total_microseconds();

			// 	output_bundle->add_value(
			// 			Window((*it).ts - microseconds(offset), trans->window_size),
			// 			*it);
			// }

            return true;

		}
	};


#if 0
// XXX refactor through inheriting from SingleInputTransformEvaluator
template <class InputT>
class FixedWindowIntoEvaluator
: public TransformEvaulator<FixedWindowInto<InputT>> {
	public:
		void evaluate(FixedWindowInto<InputT>* trans, EvaluationBundleContext* c) {

			PValue* in1 = trans->getFirstInput();
			assert(in1);

			// get one pending bundle from the input.
			// this will update input's min_ts. note that the bundle is still
			// a "pending" work in the current transform
			unique_lock<mutex> lock(trans->mtx_watermark);

			auto bundle = \
				      dynamic_pointer_cast<RecordBitmapBundle<InputT>>(in1->getOneBundle());
			assert(bundle);

			assert(trans->inflight_bundles.count(bundle) == 0);
			trans->inflight_bundles.insert(bundle);
			lock.unlock();

			auto output_bundle = make_shared<WindowsBundle<InputT>>();

			// go through Records w/ iterator, which deals with the "masked" records.
			for (auto && it = bundle->begin(); it != bundle->end(); ++it) {
				// the time offset within a window
				long offset = ((*it).ts - trans->start).total_microseconds() \
					      % (trans->window_size).total_microseconds();

				output_bundle->add_value(
						Window((*it).ts - microseconds(offset), trans->window_size),
						*it);
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

			// NB: the consumer may consume one bundle (which may contain data of
			// multi windows) at a time; or it may hold on to consume multi bundles
			// altogether.
			// The consumer needs to keep a series of per-window state internally.
		}
};
#endif

#endif /* FIXEDWINDOWINTOEVALUATOR_H */
