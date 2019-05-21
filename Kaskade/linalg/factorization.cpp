/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>

#include "umfpack.h"

#include "linalg/factorization.hh"

namespace Kaskade
{

template<>
void Factorization<double,int>::tripletToCompressedColumn(int nRows, int nCols, size_t nNonZeros,
    std::vector<int> const &ridx, std::vector<int> const &cidx, std::vector<double> const &values,
    std::vector<int> &Ap, std::vector<int> &Ai, std::vector<double> &Az, std::vector<int> &map)
  {
    assert(cidx.size()==nNonZeros &&ridx.size()==nNonZeros && values.size()==nNonZeros);
    Ap.resize(nCols+1);
    Ai.resize(nNonZeros);
    Az.resize(nNonZeros);
    int status = umfpack_di_triplet_to_col(nRows, nCols, nNonZeros,
       &ridx[0], &cidx[0], &values[0], &Ap[0], &Ai[0], &Az[0], &map[0]);
    umfpack_di_report_status(0,status);
    assert(status==UMFPACK_OK);
  }

template<>
void Factorization<double,long>::tripletToCompressedColumn(long nRows, long nCols, size_t nNonZeros,
    std::vector<long> const &ridx, std::vector<long> const &cidx, std::vector<double> const &values,
    std::vector<long> &Ap, std::vector<long> &Ai, std::vector<double> &Az, std::vector<long> &map)
  {
    assert(cidx.size()==nNonZeros &&ridx.size()==nNonZeros && values.size()==nNonZeros);
    Ap.resize(nCols+1);
    Ai.resize(nNonZeros);
    Az.resize(nNonZeros);
    long status = umfpack_dl_triplet_to_col(nRows, nCols, nNonZeros,
       &ridx[0], &cidx[0], &values[0], &Ap[0], &Ai[0], &Az[0], &map[0]);
    umfpack_dl_report_status(0,status);
    assert(status==UMFPACK_OK);
  }
  
}  // namespace Kaskade