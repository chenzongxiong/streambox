/*
 * Source.h
 *
 *  Created on: Jun 21, 2016
 *      Author: xzl
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include <mutex>
#include <string>
#include <climits>
#include <algorithm>
#include <list>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "BoundedWindow.h"
#include "OffsetRangeTracker.h"

using namespace std;
using namespace boost::posix_time;

/////////////////////////////////////////////////////////////

// @param <T> Type of elements read by the source.
template<class T>
class Source {
public:
	virtual void validate() = 0;
	virtual ~Source() { }
};

// in the original design, @Reader is part of the @Source.
template <class T>
class Reader {
public:
	virtual bool start() = 0;
	virtual bool advance() = 0;
	virtual T getCurrent() = 0;
	virtual ptime getCurrentTimestamp() = 0;
	virtual void close() = 0;
	virtual Source<T>* getCurrentSource() = 0;
	virtual ~Reader() { }
};

/////////////////////////////////////////////////////////////

/**
 * A {@link Source} that reads a finite amount of input and, because of that, supports
 * some additional operations.
 *
 * <p>The operations are:
 * <ul>
 * <li>Splitting into bundles of given size: {@link #splitIntoBundles};
 * <li>Size estimation: {@link #getEstimatedSizeBytes};
 * <li>Telling whether or not this source produces key/value pairs in sorted order:
 * {@link #producesSortedKeys};
 * <li>The accompanying {@link BoundedReader reader} has additional functionality to enable runners
 * to dynamically adapt based on runtime conditions.
 *     <ul>
 *       <li>Progress estimation ({@link BoundedReader#getFractionConsumed})
 *       <li>Tracking of parallelism, to determine whether the current source can be split
 *        ({@link BoundedReader#getSplitPointsConsumed()} and
 *        {@link BoundedReader#getSplitPointsRemaining()}).
 *       <li>Dynamic splitting of the current source ({@link BoundedReader#splitAtFraction}).
 *     </ul>
 *     </li>
 * </ul>
 *
 * <p>To use this class for supporting your custom input type, derive your class
 * class from it, and override the abstract methods. For an example, see {@link DatastoreIO}.
 *
 * @param <T> Type of records read by the source.
 */

template<class T> class BoundedReader;

template<class T>
class BoundedSource: public Source<T> {
public:

	virtual list<BoundedSource<T>*>* splitIntoBundles(
			long desiredBundleSizeBytes, PipelineOptions* options) = 0;

	virtual long getEstimatedSizeBytes(PipelineOptions* options) = 0;

	virtual bool producesSortedKeys(PipelineOptions* options) = 0;

	virtual BoundedReader<T>* createReader(PipelineOptions* options) = 0;

	virtual ~BoundedSource() { }
};

// xzl -- note that there's only one template argument, for
// both the source and the enclosed reader.
template<class T>
class BoundedReader: public Reader<T> {
public:
	// return negative when cannot estimate
	double getFranctionConsumed() {
		return -1;
	}

	const static long SPLIT_POINTS_UNKNOWN = -1;

	long getSplitPointsConsumed() {
		return SPLIT_POINTS_UNKNOWN;
	}

	long getSplitPointsRemaining() {
		return SPLIT_POINTS_UNKNOWN;
	}

	virtual BoundedSource<T>* getCurrentSource() = 0;

	BoundedSource<T>* splitAtFraction(double fraction) {
		return NULL;
	}

	/**
	 * By default, returns the minimum possible timestamp.
	 */
	ptime getCurrentTimestamp() {
		return BoundedWindow::TIMESTAMP_MIN_VALUE;
	}

	virtual ~BoundedReader() { }
};

/////////////////////////////////////////////////////////////

/**
 * A {@link BoundedSource} that uses offsets to define starting and ending positions.
 *
 * <p>{@link OffsetBasedSource} is a common base class for all bounded sources where the input can
 * be represented as a single range, and an input can be efficiently processed in parallel by
 * splitting the range into a set of disjoint ranges whose union is the original range. This class
 * should be used for sources that can be cheaply read starting at any given offset.
 * {@link OffsetBasedSource} stores the range and implements splitting into bundles.
 *
 * <p>Extend {@link OffsetBasedSource} to implement your own offset-based custom source.
 * {@link FileBasedSource}, which is a subclass of this, adds additional functionality useful for
 * custom sources that are based on files. If possible implementors should start from
 * {@link FileBasedSource} instead of {@link OffsetBasedSource}.
 *
 * <p>Consult {@link RangeTracker} for important semantics common to all sources defined by a range
 * of positions of a certain type, including the semantics of split points
 * ({@link OffsetBasedReader#isAtSplitPoint}).
 *
 * @param <T> Type of records represented by the source.
 * @see BoundedSource
 * @see FileBasedSource
 * @see RangeTracker
 */

template <class T>
class OffsetBasedSource: public BoundedSource<T> {

private:
	long _startOffset;
	long _endOffset;
	long _minBundleSize;

public:

	virtual ~OffsetBasedSource() { }

	/**
	 * @param startOffset starting offset (inclusive) of the source. Must be non-negative.
	 *
	 * @param endOffset ending offset (exclusive) of the source. Use {@link Long#MAX_VALUE} to
	 *        indicate that the entire source after {@code startOffset} should be read. Must be
	 *        {@code > startOffset}.
	 *
	 * @param minBundleSize minimum bundle size in offset units that should be used when splitting the
	 *                      source into sub-sources. This value may not be respected if the total
	 *                      range of the source is smaller than the specified {@code minBundleSize}.
	 *                      Must be non-negative.
	 */

	OffsetBasedSource(long startOffset, long endOffset, long minBundleSize) :
			_startOffset(startOffset), _endOffset(endOffset), _minBundleSize(
					minBundleSize) {
	}

	long getStartOffset() {
		return _startOffset;
	}

	long getEndOffset() {
		return _endOffset;
	}

	long getMinBundleSize() {
		return _minBundleSize;
	}

	long getEstimatedSizeBytes(PipelineOptions* options) {
		long trueEndOffset =
				(_endOffset == LONG_MAX) ? getMaxEndOffset(options) : _endOffset;
		return getBytesPerOffset() * (trueEndOffset - getStartOffset());
	}

//	list<OffsetBasedSource<T>*>* splitIntoBundles(long desiredBundleSizeBytes,
	list<BoundedSource<T>*>* splitIntoBundles(long desiredBundleSizeBytes,
			PipelineOptions* options) {

		long desiredBundleSizeOffsetUnits = max(
				max(1L, desiredBundleSizeBytes / getBytesPerOffset()),
				_minBundleSize);

//		list<OffsetBasedSource<T>*>* subSources = new list<OffsetBasedSource<T>*>();
		list<BoundedSource<T>*>* subSources = new list<BoundedSource<T>*>();
		long start = _startOffset;
		long maxEnd = min(_endOffset, getMaxEndOffset(options));

		while (start < maxEnd) {
			long end = start + desiredBundleSizeOffsetUnits;
			end = min(end, maxEnd);

			// Avoid having a too small bundle at the end and ensure that we respect minBundleSize.
			long remaining = maxEnd - end;
			if ((remaining < desiredBundleSizeOffsetUnits / 4)
					|| (remaining < _minBundleSize)) {
				end = maxEnd;
			}
			subSources->push_back(createSourceForSubrange(start, end));

			start = end;
		}
		return subSources;
	}

	void validate() {
			assert(_startOffset >= 0 && _endOffset >=0 && _startOffset < _endOffset
					&& _minBundleSize >= 0);
		}

	string toString() {
			return "[" + to_string(_startOffset) + ", " + to_string(_endOffset) + ")";
		}

	/**
	 * Returns approximately how many bytes of data correspond to a single offset in this source.
	 * Used for translation between this source's range and methods defined in terms of bytes, such
	 * as {@link #getEstimatedSizeBytes} and {@link #splitIntoBundles}.
	 *
	 * <p>Defaults to {@code 1} byte, which is the common case for, e.g., file sources.
	 */
	long getBytesPerOffset() {
		return 1L;
	}

	/**
	 * Returns the actual ending offset of the current source. The value returned by this function
	 * will be used to clip the end of the range {@code [startOffset, endOffset)} such that the
	 * range used is {@code [startOffset, min(endOffset, maxEndOffset))}.
	 *
	 * <p>As an example in which {@link OffsetBasedSource} is used to implement a file source, suppose
	 * that this source was constructed with an {@code endOffset} of {@link Long#MAX_VALUE} to
	 * indicate that a file should be read to the end. Then this function should determine
	 * the actual, exact size of the file in bytes and return it.
	 */
	virtual long getMaxEndOffset(PipelineOptions* options) = 0;

	/**
	 * Returns an {@link OffsetBasedSource} for a subrange of the current source. The
	 * subrange {@code [start, end)} must be within the range {@code [startOffset, endOffset)} of
	 * the current source, i.e. {@code startOffset <= start < end <= endOffset}.
	 */
	virtual OffsetBasedSource<T>* createSourceForSubrange(long start,
			long end) = 0;

	/**
	 * Whether this source should allow dynamic splitting of the offset ranges.
	 *
	 * <p>True by default. Override this to return false if the source cannot
	 * support dynamic splitting correctly. If this returns false,
	 * {@link OffsetBasedSource.OffsetBasedReader#splitAtFraction} will refuse all split requests.
	 */
	bool allowsDynamicSplitting() {
		return true;
	}
};

/**
 * A {@link Source.Reader} that implements code common to readers of all
 * {@link OffsetBasedSource}s.
 *
 * <p>Subclasses have to implement:
 * <ul>
 *   <li>The methods {@link #startImpl} and {@link #advanceImpl} for reading the
 *   first or subsequent records.
 *   <li>The methods {@link #getCurrent}, {@link #getCurrentOffset}, and optionally
 *   {@link #isAtSplitPoint} and {@link #getCurrentTimestamp} to access properties of
 *   the last record successfully read by {@link #startImpl} or {@link #advanceImpl}.
 * </ul>
 */
template <class T>
class OffsetBasedReader: public BoundedReader<T> {
private:
	OffsetBasedSource<T>* _source;

	/** The {@link OffsetRangeTracker} managing the range and current position of the source. */
	OffsetRangeTracker _rangeTracker;

	/**
	 * Returns true if the last call to {@link #start} or {@link #advance} returned false.
	 */
	bool isDone() {
		return _rangeTracker.isDone();
	}

	/**
	 * Returns true if there has been a call to {@link #start}.
	 */
	bool isStarted() {
		return _rangeTracker.isStarted();
	}

public:

	virtual ~OffsetBasedReader() { }

	/**
	 * @param source the {@link OffsetBasedSource} to be read by the current reader.
	 */
	OffsetBasedReader(OffsetBasedSource<T>* source) :
			_source(source), _rangeTracker(
					OffsetRangeTracker(_source->getStartOffset(),
							_source->getEndOffset())) {
	}

	/**
	 * Returns the <i>starting</i> offset of the {@link Source.Reader#getCurrent current record},
	 * which has been read by the last successful {@link Source.Reader#start} or
	 * {@link Source.Reader#advance} call.
	 * <p>If no such call has been made yet, the return value is unspecified.
	 * <p>See {@link RangeTracker} for description of offset semantics.
	 */

	bool start() {
		return (startImpl()
				&& _rangeTracker.tryReturnRecordAt(isAtSplitPoint(),
						getCurrentOffset()))
				|| _rangeTracker.markDone();
	}

	bool advance() {
		return (advanceImpl()
				&& _rangeTracker.tryReturnRecordAt(isAtSplitPoint(),
						getCurrentOffset()))
				|| _rangeTracker.markDone();
	}

	// xzl: synchronized method in java
	OffsetBasedSource<T>* getCurrentSource() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);
		return _source;
	}

	double getFractionConsumed() {
		return _rangeTracker.getFractionConsumed();
	}

	long getSplitPointsConsumed() {
		return _rangeTracker.getSplitPointsProcessed();
	}

	long getSplitPointsRemaining() {
		if (isDone()) {
			return 0;
		} else if (!isStarted()) {
			// Note that even if the current source does not allow splitting, we don't know that
			// it's non-empty so we return UNKNOWN instead of 1.
			return BoundedSource<T>::BoundedReader.SPLIT_POINTS_UNKNOWN;
		} else if (!getCurrentSource()->allowsDynamicSplitting()) {
			// Started (so non-empty) and unsplittable, so only the current task.
			return 1;
		} else if (getCurrentOffset()
				>= _rangeTracker.getStopPosition() - 1) {
			// If this is true, the next element is outside the range. Note that even getCurrentOffset()
			// might be larger than the stop position when the current record is not a split point.
			return 1;
		} else {
			// Use the default.
			return BoundedSource<T>::BoundedReader.getSplitPointsRemaining();
		}
	}

	// TODO: improve tracing output
	OffsetBasedSource<T>* splitAtFraction(double fraction) {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		if (!getCurrentSource()->allowsDynamicSplitting()) {
			return NULL;
		}
		if (_rangeTracker.getStopPosition() == LONG_MAX) {
			printf(
					"Refusing to split unbounded OffsetBasedReader {} at fraction\n");
			return NULL;
		}
		long splitOffset = _rangeTracker.getPositionForFractionConsumed(
				fraction);
		printf(
				"Proposing to split OffsetBasedReader at fraction %.2f (offset %d)",
				fraction, splitOffset);
		if (!_rangeTracker.trySplitAtPosition(splitOffset)) {
			return NULL;
		}

		long start = _source->getStartOffset();
		long end = _source->getEndOffset();
		OffsetBasedSource<T>* primary = _source->createSourceForSubrange(
				start, splitOffset);
		OffsetBasedSource<T>* residual = _source->createSourceForSubrange(
				splitOffset, end);

		// xzl: XXX need to destroy the original _source ???
		this->_source = primary;
		return residual;
	}

protected:
	virtual long getCurrentOffset() = 0;

	/**
	 * Returns whether the current record is at a split point (i.e., whether the current record
	 * would be the first record to be read by a source with a specified start offset of
	 * {@link #getCurrentOffset}).
	 *
	 * <p>See detailed documentation about split points in {@link RangeTracker}.
	 */
	bool isAtSplitPoint() {
		return true;
	}

	/**
	 * Initializes the {@link OffsetBasedSource.OffsetBasedReader} and advances to the first record,
	 * returning {@code true} if there is a record available to be read. This method will be
	 * invoked exactly once and may perform expensive setup operations that are needed to
	 * initialize the reader.
	 *
	 * <p>This function is the {@code OffsetBasedReader} implementation of
	 * {@link BoundedReader#start}. The key difference is that the implementor can ignore the
	 * possibility that it should no longer produce the first record, either because it has exceeded
	 * the original {@code endOffset} assigned to the reader, or because a concurrent call to
	 * {@link #splitAtFraction} has changed the source to shrink the offset range being read.
	 *
	 * @see BoundedReader#start
	 */
	virtual bool startImpl() = 0;

	/**
	 * Advances to the next record and returns {@code true}, or returns false if there is no next
	 * record.
	 *
	 * <p>This function is the {@code OffsetBasedReader} implementation of
	 * {@link BoundedReader#advance}. The key difference is that the implementor can ignore the
	 * possibility that it should no longer produce the next record, either because it has exceeded
	 * the original {@code endOffset} assigned to the reader, or because a concurrent call to
	 * {@link #splitAtFraction} has changed the source to shrink the offset range being read.
	 *
	 * @see BoundedReader#advance
	 */
	virtual bool advanceImpl() = 0;

};
#endif /* SOURCE_H_ */
