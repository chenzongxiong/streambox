#include "WinKeyReducer/WinKeyReducerEval.h"

using NYTWinKeyReduerEval = WinKeyReducerEval<std::pair<uint64_t, uint64_t>, /* kv in */
                                              WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
                                              std::pair<uint64_t, uint64_t>, /* kvout */
                                              WindowsBundle /* output bundle */
                                              >;

#include "NYTWinKeyReducerEvalCommon.h"
