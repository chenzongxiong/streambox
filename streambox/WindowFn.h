#include <vector>
#include "boost/date_time/posix_time/posix_time.hpp"
using namespace boost::posix_time;
using namespace std;

/**
 * The argument to the {@link Window} transform used to assign elements into
 * windows and to determine how windows are merged.  See {@link Window} for more
 * information on how {@code WindowFn}s are used and for a library of
 * predefined {@code WindowFn}s.
 *
 * <p>Users will generally want to use the predefined
 * {@code WindowFn}s, but it is also possible to create new
 * subclasses.
 *
 * <p>To create a custom {@code WindowFn}, inherit from this class and override all required
 * methods.  If no merging is required, inherit from {@link NonMergingWindowFn}
 * instead.  If no merging is required and each element is assigned to a single window, inherit from
 * {@code PartitioningWindowFn}.  Inheriting from the most specific subclass will enable more
 * optimizations in the runner.
 *
 * @param <T> type of elements being windowed
 * @param <W> {@link BoundedWindow} subclass used to represent the
 *            windows used by this {@code WindowFn}
 */

template <class T, class W>
class WindowFn {  /* HasDisplayData NOT implemented*/
	public:
		/**
		 * Information available when running {@link #assignWindows}.
		 */
		virtual void mergeWindows(MergeContext *c) { };
		virtual vector<W*> assignWindows(AssignContext *c) { };
		virtual bool isNonMerging() { return false; }
		virtual bool assignsToSingleWindow() { return false; }
		virtual ~WindowFn() { }
		virtual ptime getOutputTime(ptime inputTimestamp, GlobalWindow* window) = 0;
		///////////////////////////////////////////////////////////
		/**
		 * Information available when running {@link #assignWindows}.
		 */
		class AssignContext {
			public:
				virtual T element() = 0;
				virtual ptime timestamp() = 0;
				/**
				 * Returns the windows the current element was in, prior to this
				 * {@code WindowFn} being called.
				 */
				vector<BoundedWindow*> windows() = 0;

				virtual ~AssignContext() { }
		};
		///////////////////////////////////////////////////////////
		class MergeContext {
			public:
				virtual vector<W*>*  windows() = 0;
				virtual void merge(vector<W *> toBeMerged, W* mergeResult) = 0;
				virtual ~MergeContext() { }
		};
};

extern vector<BoundedWindow *> GLOBAL_WINDOWS; // Values.cpp

template <class T, class W>
class GlobalWindows : public WindowFn {
	public:
		vector<BoundedWindow *>
			assignWindows(WindowFn<T,W>::AssignContext *c) override {
				return GLOBAL_WINDOWS; // by copy
			}
		bool isNonMerging() override { return true; }
		bool assignsToSingleWindow() override { return true; }
		ptime getOutputTime(ptime inputTimestamp, GlobalWindow* window) override {
			return inputTimestamp;
		}
};

