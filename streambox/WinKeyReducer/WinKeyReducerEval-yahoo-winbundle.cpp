//
// Created by manuelrenz on 12.04.18.
//


#include "WinKeyReducerEval.h"

/* -----------------------------------------
 * specialized for yahoo -- using normal hashtable
 * NB: we have to do full specialization per C++
 * ----------------------------------------- */

// using MyWinKeyReduerEval = WinKeyReducerEval<std::pair<uint64_t, long>, /* kv in */
//         WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
//         std::pair<uint64_t, long>, /* kvout */
//         WindowsBundle /* output bundle */
// >;

using MyWinKeyReduerEval = WinKeyReducerEval<std::pair<uint64_t, uint64_t>, /* kv in */
        WinKeyFragLocal_Std, WinKeyFrag_Std, /* internal map format */
        std::pair<uint64_t, uint64_t>, /* kvout */
        WindowsBundle /* output bundle */
>;

#include "WinKeyReducerEval-yahoo-common.h"
