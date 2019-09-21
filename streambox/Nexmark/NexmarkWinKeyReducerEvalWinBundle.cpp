#include "WinKeyReducer/WinKeyReducerEval.h"

using NexmarkWinKeyReduerEval = WinKeyReducerEval<std::pair<uint64_t, uint64_t>, /* kv in */
                                                  WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
                                                  std::pair<uint64_t, uint64_t>, /* kvout */
                                                  WindowsBundle /* output bundle */
                                                  >;

#include "NexmarkWinKeyReducerEvalCommon.hpp"
