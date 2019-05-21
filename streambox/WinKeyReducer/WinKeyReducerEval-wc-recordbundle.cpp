
/*
 * WindowKeyedReducerEval-wc.cpp
 *
 *  Created on: Jan 21, 2017
 *      Author: xzl
 *
 *  Specialization of WindowKeyedReducerEval for wordcount.
 *  Separate because if left in the header, have error
 *  "Specialization of member function template after instantiation error"
 */

#include "WinKeyReducerEval.h"

/* -----------------------------------------
 * specialized for wordcount -- using normal hashtable
 * NB: we have to do full specialization per C++
 * ----------------------------------------- */

using MyKVPairIn = std::pair<creek::string, long>;
using MyKVPairOut = creek::tvpair;

using MyWinKeyReduerEval = WinKeyReducerEval<
		MyKVPairIn, /* kv in */
		WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
		MyKVPairOut, /* kvout, note that this is different than in pair */
		RecordBundle /* output bundle */
>;

/* specialize the output record converter */
template<>
Record<MyKVPairOut> const MyWinKeyReduerEval::makeRecord(MyKVPairIn const & value, Window const & win) {
//	value.first = win.start.ticks() + win.duration.ticks();
//	return Record<KVOut>(value, win.window_end());

	/* avoid updating @value in place which may be r/o */
	return Record<MyKVPairOut>(make_pair(win.start.ticks() + win.duration.ticks(), value.second),
											 win.window_end()
	);
}

#include "WinKeyReducerEval-wc-common.h"
