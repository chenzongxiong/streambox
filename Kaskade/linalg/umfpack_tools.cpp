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

#include "umfpack.h"
#include "linalg/umfpack_tools.hh"

void umfpack_triplet_to_col(
    int rows,
    int cols,
    std::vector<int> const& ridx,
    std::vector<int> const& cidx,
    std::vector<double> const& values,
    std::vector<int>& Ap,
    std::vector<int>& Ai,
    std::vector<double>& Az)
{
  Ap.resize(cols+1);
  Ai.resize(values.size());
  Az.resize(values.size());
  umfpack_di_triplet_to_col(rows,cols,ridx.size(),
      &ridx[0],&cidx[0],&values[0],
      &Ap[0],&Ai[0],&Az[0],
      0);
}

void umfpack_triplet_to_col(
    int rows,
    int cols,
    std::vector<long> const& ridx,
    std::vector<long> const& cidx,
    std::vector<double> const& values,
    std::vector<long>& Ap,
    std::vector<long>& Ai,
    std::vector<double>& Az)
{
  Ap.resize(cols+1);
  Ai.resize(values.size());
  Az.resize(values.size());
  umfpack_dl_triplet_to_col(rows,cols,ridx.size(),
      &ridx[0],&cidx[0],&values[0],
      &Ap[0],&Ai[0],&Az[0],
      0);
}

void umfpack_col_to_triplet(
    std::vector<int> const& Ap,
    std::vector<int>& Ti)
{
  Ti.resize(Ap.back());
  umfpack_di_col_to_triplet(Ap.size()-1,&Ap[0],&Ti[0]);
}

void umfpack_col_to_triplet(
    std::vector<long> const& Ap,
    std::vector<long>& Ti)
{
  Ti.resize(Ap.back());
  umfpack_dl_col_to_triplet(Ap.size()-1,&Ap[0],&Ti[0]);
}



