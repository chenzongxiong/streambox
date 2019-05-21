/*
 * OffsetRangeTracker.h
 *
 *  Created on: Jun 22, 2016
 *      Author: xzl
 */

#ifndef OFFSETRANGETRACKER_H_
#define OFFSETRANGETRACKER_H_

#include <assert.h>
#include <math.h>
#include <mutex>
#include <string>
#include <algorithm>
#include <climits>

#include "boost/format.hpp"

using namespace std;

/**
 * A {@code RangeTracker} is a thread-safe helper object for implementing dynamic work rebalancing
 * in position-based {@link org.apache.beam.sdk.io.BoundedSource.BoundedReader}
 * subclasses.
 *
 * <h3>Usage of the RangeTracker class hierarchy</h3>
 * The abstract {@code RangeTracker} interface should not be used per se - all users should use its
 * subclasses directly. We declare it here because all subclasses have roughly the same interface
 * and the same properties, to centralize the documentation. Currently we provide one
 * implementation - {@link OffsetRangeTracker}.
 *
 * <h3>Position-based sources</h3>
 * A position-based source is one where the source can be described by a range of positions of
 * an ordered type and the records returned by the reader can be described by positions of the
 * same type.
 *
 * <p>In case a record occupies a range of positions in the source, the most important thing about
 * the record is the position where it starts.
 *
 * <p>Defining the semantics of positions for a source is entirely up to the source class, however
 * the chosen definitions have to obey certain properties in order to make it possible to correctly
 * split the source into parts, including dynamic splitting. Two main aspects need to be defined:
 * <ul>
 *   <li>How to assign starting positions to records.
 *   <li>Which records should be read by a source with a range {@code [A, B)}.
 * </ul>
 * Moreover, reading a range must be <i>efficient</i>, i.e., the performance of reading a range
 * should not significantly depend on the location of the range. For example, reading the range
 * {@code [A, B)} should not require reading all data before {@code A}.
 *
 * <p>The sections below explain exactly what properties these definitions must satisfy, and
 * how to use a {@code RangeTracker} with a properly defined source.
 *
 * <h3>Properties of position-based sources</h3>
 * The main requirement for position-based sources is <i>associativity</i>: reading records from
 * {@code [A, B)} and records from {@code [B, C)} should give the same records as reading from
 * {@code [A, C)}, where {@code A <= B <= C}. This property ensures that no matter how a range
 * of positions is split into arbitrarily many sub-ranges, the total set of records described by
 * them stays the same.
 *
 * <p>The other important property is how the source's range relates to positions of records in
 * the source. In many sources each record can be identified by a unique starting position.
 * In this case:
 * <ul>
 *   <li>All records returned by a source {@code [A, B)} must have starting positions
 *   in this range.
 *   <li>All but the last record should end within this range. The last record may or may not
 *   extend past the end of the range.
 *   <li>Records should not overlap.
 * </ul>
 * Such sources should define "read {@code [A, B)}" as "read from the first record starting at or
 * after A, up to but not including the first record starting at or after B".
 *
 * <p>Some examples of such sources include reading lines or CSV from a text file, reading keys and
 * values from a BigTable, etc.
 *
 * <p>The concept of <i>split points</i> allows to extend the definitions for dealing with sources
 * where some records cannot be identified by a unique starting position.
 *
 * <p>In all cases, all records returned by a source {@code [A, B)} must <i>start</i> at or after
 * {@code A}.
 *
 * <h3>Split points</h3>
 *
 * <p>Some sources may have records that are not directly addressable. For example, imagine a file
 * format consisting of a sequence of compressed blocks. Each block can be assigned an offset, but
 * records within the block cannot be directly addressed without decompressing the block. Let us
 * refer to this hypothetical format as <i>CBF (Compressed Blocks Format)</i>.
 *
 * <p>Many such formats can still satisfy the associativity property. For example, in CBF, reading
 * {@code [A, B)} can mean "read all the records in all blocks whose starting offset is in
 * {@code [A, B)}".
 *
 * <p>To support such complex formats, we introduce the notion of <i>split points</i>. We say that
 * a record is a split point if there exists a position {@code A} such that the record is the first
 * one to be returned when reading the range  {@code [A, infinity)}. In CBF, the only split points
 * would be the first records in each block.
 *
 * <p>Split points allow us to define the meaning of a record's position and a source's range
 * in all cases:
 * <ul>
 *   <li>For a record that is at a split point, its position is defined to be the largest
 *   {@code A} such that reading a source with the range {@code [A, infinity)} returns this record;
 *   <li>Positions of other records are only required to be non-decreasing;
 *   <li>Reading the source {@code [A, B)} must return records starting from the first split point
 *   at or after {@code A}, up to but not including the first split point at or after {@code B}.
 *   In particular, this means that the first record returned by a source MUST always be
 *   a split point.
 *   <li>Positions of split points must be unique.
 * </ul>
 * As a result, for any decomposition of the full range of the source into position ranges, the
 * total set of records will be the full set of records in the source, and each record
 * will be read exactly once.
 *
 * <h3>Consumed positions</h3>
 * As the source is being read, and records read from it are being passed to the downstream
 * transforms in the pipeline, we say that positions in the source are being <i>consumed</i>.
 * When a reader has read a record (or promised to a caller that a record will be returned),
 * positions up to and including the record's start position are considered <i>consumed</i>.
 *
 * <p>Dynamic splitting can happen only at <i>unconsumed</i> positions. If the reader just
 * returned a record at offset 42 in a file, dynamic splitting can happen only at offset 43 or
 * beyond, as otherwise that record could be read twice (by the current reader and by a reader
 * of the task starting at 43).
 *
 * <h3>Example</h3>
 * The following example uses an {@link OffsetRangeTracker} to support dynamically splitting
 * a source with integer positions (offsets).
 * <pre> {@code
 *   class MyReader implements BoundedReader<Foo> {
 *     private MySource currentSource;
 *     private final OffsetRangeTracker tracker = new OffsetRangeTracker();
 *     ...
 *     MyReader(MySource source) {
 *       this.currentSource = source;
 *       this.tracker = new MyRangeTracker<>(source.getStartOffset(), source.getEndOffset())
 *     }
 *     ...
 *     boolean start() {
 *       ... (general logic for locating the first record) ...
 *       if (!tracker.tryReturnRecordAt(true, recordStartOffset)) return false;
 *       ... (any logic that depends on the record being returned, e.g. counting returned records)
 *       return true;
 *     }
 *     boolean advance() {
 *       ... (general logic for locating the next record) ...
 *       if (!tracker.tryReturnRecordAt(isAtSplitPoint, recordStartOffset)) return false;
 *       ... (any logic that depends on the record being returned, e.g. counting returned records)
 *       return true;
 *     }
 *
 *     double getFractionConsumed() {
 *       return tracker.getFractionConsumed();
 *     }
 *   }
 * } </pre>
 *
 * <h3>Usage with different models of iteration</h3>
 * When using this class to protect a
 * {@link org.apache.beam.sdk.io.BoundedSource.BoundedReader}, follow the pattern
 * described above.
 *
 * <p>When using this class to protect iteration in the {@code hasNext()/next()}
 * model, consider the record consumed when {@code hasNext()} is about to return true, rather than
 * when {@code next()} is called, because {@code hasNext()} returning true is promising the caller
 * that {@code next()} will have an element to return - so {@link #trySplitAtPosition} must not
 * split the range in a way that would make the record promised by {@code hasNext()} belong to
 * a different range.
 *
 * <p>Also note that implementations of {@code hasNext()} need to ensure
 * that they call {@link #tryReturnRecordAt} only once even if {@code hasNext()} is called
 * repeatedly, due to the requirement on uniqueness of split point positions.
 *
 * @param <PositionT> Type of positions used by the source to define ranges and identify records.
 *
 * xzl: this is a "helper". it does not operate on the source, but just maintain the record positions.
 * [ the reader needs to take the position and do actual read. ]
 */
template<class PositionT>
class RangeTracker {
public:
	/**
	 * Returns the starting position of the current range, inclusive.
	 */
	virtual PositionT getStartPosition() = 0;

	/**
	 * Returns the ending position of the current range, exclusive.
	 */
	virtual PositionT getStopPosition() = 0;

	/**
	 * Atomically determines whether a record at the given position can be returned and updates
	 * internal state. In particular:
	 * <ul>
	 *   <li>If {@code isAtSplitPoint} is {@code true}, and {@code recordStart} is outside the current
	 *   range, returns {@code false};
	 *   <li>Otherwise, updates the last-consumed position to {@code recordStart} and returns
	 *   {@code true}.
	 * </ul>
	 * <p>This method MUST be called on all split point records. It may be called on every record.
	 */
	virtual bool tryReturnRecordAt(bool isAtSplitPoint,
			PositionT recordStart) = 0;

	/**
	 * Atomically splits the current range [{@link #getStartPosition}, {@link #getStopPosition})
	 * into a "primary" part [{@link #getStartPosition}, {@code splitPosition})
	 * and a "residual" part [{@code splitPosition}, {@link #getStopPosition}), assuming the current
	 * last-consumed position is within [{@link #getStartPosition}, splitPosition)
	 * (i.e., {@code splitPosition} has not been consumed yet).
	 *
	 * <p>Updates the current range to be the primary and returns {@code true}. This means that
	 * all further calls on the current object will interpret their arguments relative to the
	 * primary range.
	 *
	 * <p>If the split position has already been consumed, or if no {@link #tryReturnRecordAt} call
	 * was made yet, returns {@code false}. The second condition is to prevent dynamic splitting
	 * during reader start-up.
	 */
	virtual bool trySplitAtPosition(PositionT splitPosition) = 0;

	/**
	 * Returns the approximate fraction of positions in the source that have been consumed by
	 * successful {@link #tryReturnRecordAt} calls, or 0.0 if no such calls have happened.
	 */
	virtual double getFractionConsumed() = 0;

	virtual ~RangeTracker() { }
};

/**
 * A {@link RangeTracker} for non-negative positions of type {@code long}.
 */

class OffsetRangeTracker : public RangeTracker<long> {

private:
	long _startOffset;
	long _stopOffset;
	long _lastRecordStart = -1;
	long _offsetOfLastSplitPoint = -1;
	long _splitPointsSeen = 0;
	bool _done = false;

public:
	  /**
	   * Offset corresponding to infinity. This can only be used as the upper-bound of a range, and
	   * indicates reading all of the records until the end without specifying exactly what the end is.
	   *
	   * <p>Infinite ranges cannot be split because it is impossible to estimate progress within them.
	   */
	static const long OFFSET_INFINITY = LONG_MAX;

	OffsetRangeTracker(long startOffset, long stopOffset)
		: _startOffset(startOffset), _stopOffset(stopOffset) { }

	bool isStarted() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		return (_offsetOfLastSplitPoint != -1)|| _done;
	}

	bool isDone() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		return _done;
	}

	long getStartPosition() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		return _startOffset;
	}

	long getStopPosition() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		return _stopOffset;
	}

	bool tryReturnRecordAt(bool isAtSplitPoint, long recordStart) {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		if (!isStarted() && !isAtSplitPoint) {
			printf("err: The first record [starting at %ld] must be at a split point\n", recordStart);
			assert(0);
		}
		if (recordStart < _startOffset) {
			printf("err:Trying to return record [starting at %ld] which is before the start offset [%ld]",
					recordStart, _startOffset);
			assert(0);
		}

		_lastRecordStart = recordStart;

		if (isAtSplitPoint) {
			if (recordStart == _offsetOfLastSplitPoint) {
				printf("err: Record at a split point has same offset as the previous split point: "
					"previous split point at %ld, current record starts at %ld\n",
					_offsetOfLastSplitPoint, recordStart);
				assert(0);
			}
			if (recordStart >= _stopOffset) {
				_done = true;
				return false;
			}
			_offsetOfLastSplitPoint = recordStart;
			++ _splitPointsSeen;
		}

		return true;
	}

	bool trySplitAtPosition(long splitOffset) {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

		if (_stopOffset == OFFSET_INFINITY) {
			printf("Refusing to split {} at {}: stop position unspecified\n");
			return false;
		}
		if (!isStarted()) {
			printf("Refusing to split {} at {}: unstarted\n");
			return false;
		}

	    // Note: technically it is correct to split at any position after the last returned
	    // split point, not just the last returned record.
	    // TODO: Investigate whether in practice this is useful or, rather, confusing.
	    if (splitOffset <= _lastRecordStart) {
	      printf(
	          "Refusing to split {} at {}: already past proposed split position\n");
	      return false;
	    }
	    if (splitOffset < _startOffset || splitOffset >= _stopOffset) {
	      printf(
	          "Refusing to split {} at {}: proposed split position out of range\n");
	      return false;
	    }
	    printf("Agreeing to split {} at {}\n");
	    this->_stopOffset = splitOffset;
	    return true;
	}

	long getPositionForFractionConsumed(double fraction) {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

	    if (_stopOffset == OFFSET_INFINITY) {
	      printf(
	          "getPositionForFractionConsumed is not applicable to an unbounded range\n");
	      assert(0);
	    }

	    return (long) ceil(_startOffset + fraction * (_stopOffset - _startOffset));
	  }

	double getFractionConsumed() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);

	    if (!isStarted()) {
	      return 0.0;
	    } else if (isDone()) {
	      return 1.0;
	    } else if (_stopOffset == OFFSET_INFINITY) {
	      return 0.0;
	    } else if (_lastRecordStart >= _stopOffset) {
	      return 1.0;
	    } else {
	      // E.g., when reading [3, 6) and lastRecordStart is 4, that means we consumed 3,4 of 3,4,5
	      // which is (4 - 3 + 1) / (6 - 3) = 67%.
	      // Also, clamp to at most 1.0 because the last consumed position can extend past the
	      // stop position.
	      return min(1.0, 1.0 * (_lastRecordStart - _startOffset + 1) / (_stopOffset - _startOffset));
	    }
	}

	/**
	 * Returns the total number of split points that have been processed.
	 *
	 * <p>A split point at a particular offset has been seen if there has been a corresponding call
	 * to {@link #tryReturnRecordAt(boolean, long)} with {@code isAtSplitPoint} true. It has been
	 * processed if there has been a <em>subsequent</em> call to
	 * {@link #tryReturnRecordAt(boolean, long)} with {@code isAtSplitPoint} true and at a larger
	 * offset.
	 *
	 * <p>Note that for correctness when implementing {@link BoundedReader#getSplitPointsConsumed()},
	 * if a reader finishes before {@link #tryReturnRecordAt(boolean, long)} returns false,
	 * the reader should add an additional call to {@link #markDone()}. This will indicate that
	 * processing for the last seen split point has been finished.
	 *
	 * @see org.apache.beam.sdk.io.OffsetBasedSource for a {@link BoundedReader}
	 * implemented using {@link OffsetRangeTracker}.
	 */

	long getSplitPointsProcessed() {
		static mutex _mutex;
		unique_lock < mutex > lock(_mutex);
		if (!isStarted()) {
			return 0;
		} else if (isDone()) {
			return _splitPointsSeen;
		} else {
			// There is a current split point, and it has not finished processing.
			assert(_splitPointsSeen > 0);
			// "A started rangeTracker should have seen > 0 split points",
			return _splitPointsSeen - 1;
		}
	}

	  /**
	   * Marks this range tracker as being done. Specifically, this will mark the current split point,
	   * if one exists, as being finished.
	   *
	   * <p>Always returns false, so that it can be used in an implementation of
	   * {@link BoundedReader#start()} or {@link BoundedReader#advance()} as follows:
	   *
	   * <pre> {@code
	   * public boolean start() {
	   *   return startImpl() && rangeTracker.tryReturnRecordAt(isAtSplitPoint, position)
	   *       || rangeTracker.markDone();
	   * }} </pre>
	   */

		bool markDone() {
			static mutex _mutex;
			unique_lock < mutex > lock(_mutex);

			_done = true;
			return false;
		}

		string toString() {
		    string stopString = (_stopOffset == OFFSET_INFINITY)
		    		? "infinity" : to_string(_stopOffset);

		    if (_lastRecordStart >= 0) {
		      return str(boost::format(
		          "<at [starting at %ld] of offset range [%ld, %s)>")
		          % _lastRecordStart
		          % _startOffset
		          % stopString.c_str());
		    } else {
		      return str(boost::format("<unstarted in offset range [%d, %s)>")
		    		  % _startOffset % stopString.c_str());
		    }
		}

		// todo: copy()
};

#endif /* OFFSETRANGETRACKER_H_ */
