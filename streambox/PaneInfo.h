/*
 * PaneInfo.h
 *
 *  Created on: Jul 1, 2016
 *      Author: xzl
 */
#ifndef _PANE_INFO
#define _PANE_INFO

#include "assert.h"
#include <map>

using namespace std;

namespace pane_info {
  // scoped enum
  enum class Timing {
    EARLY,
    ON_TIME,
    LATE,
    UNKNOWN
  };

  class PaneInfo {
  public:

  private:
    bool _isFirst, _isLast;
    Timing _timing;
    long _index, _nonSpeculativeIndex;

  public:
    // xzl: unclear what @onTimeIndex is for.
    PaneInfo(bool isFirst, bool isLast, Timing timing,
        long index, long onTimeIndex)
      : _isFirst(isFirst), _isLast(isLast), _timing(timing),
        _index(index), _nonSpeculativeIndex(onTimeIndex) { }

    bool isUnknown() { return (_timing == Timing::UNKNOWN); }

    bool isFirst() { return _isFirst; }

    bool isLast() { return _isLast; }

    Timing getTiming() { return _timing; }

    long getIndex() { return _index; }

    long getNonSpeculativeIndex() { return _nonSpeculativeIndex; }

  };

  // xzl: apparently, (the last arg) onTimeIndex == -1 iff Timing::EARLY
  static PaneInfo BYTE_TO_PANE_INFO[] = {
      PaneInfo(true, true, Timing::EARLY, 0, -1),
      PaneInfo(true, false, Timing::ON_TIME, 0, 0),
      PaneInfo(false, true, Timing::LATE, 0, 0),
      PaneInfo(false, false, Timing::UNKNOWN, 0, 0),
  };

  static PaneInfo* createPane(bool isFirst, bool isLast, Timing timing,
      long index, long onTimeIndex) {
    if (isFirst || timing == Timing::UNKNOWN) {
      return BYTE_TO_PANE_INFO + static_cast<int>(timing);
    } else {
      return new PaneInfo(isFirst, isLast, timing, index, onTimeIndex);
    }
  }

#if 0 // xzl: not called at all?
  static PaneInfo* createPane(bool isFirst, bool isLast, Timing timing) {
    assert(isFirst);
    return createPane(isFirst, isLast, timing, 0,
        timing == Timing::EARLY ? -1 : 0);
  }
#endif

  // xzl: XXX separate the following to cpp files??
  /**
   * {@code PaneInfo} to use for elements on (and before) initial window assignemnt (including
   * elements read from sources) before they have passed through a {@link GroupByKey} and are
   * associated with a particular trigger firing.
   */
  static PaneInfo* NO_FIRING = \
      createPane(true, true, Timing::UNKNOWN, 0, 0);
}

#endif /* _PANE_INFO */
