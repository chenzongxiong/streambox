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

/**
 * \file
 * \author Bodo Erdmann , Felix Lehmann
 */


#ifndef ICCPRECOND_HH
#define ICCPRECOND_HH

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include <boost/timer/timer.hpp>

#include <dune/istl/preconditioners.hh>

#include "protospp.h" // "Malloc" from ITSOL

#include "linalg/triplet.hh"

extern "C" {
#include "taucs.h"
}

namespace Kaskade
{

//---------------------------------------------------------------------
//---------------------------------------------------------------------

/**
 * \ingroup linalg
 *
 * PrecondType::ICC Incomplete Cholesky Preconditioner from TAUCS library
 *
 * 'modified' is defined in TAUCS, if is set true, the rows of the factorization L
 * will be scaled so that LL^T has same rowsum as A
 *
 * we expect to be passed only a lower triangle of A (upper = 0).
 * If this is not the case, nnz_of_one_triangle should be passed with the number of
 * nonzeros of the lower triangle of A (including nonzeroes on diagonal)
 *
 */
template <class Op>
class ICCPreconditioner: public Dune::Preconditioner<typename Op::range_type, typename Op::range_type>
{
  typedef typename Op::range_type range_type;
  typedef range_type               domain_type;
  typedef typename Op::field_type field_type;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<field_type,1,1> > BCRS_Matrix;
  typedef typename BCRS_Matrix::ConstColIterator ColIter;


public:
  static int const category = Dune::SolverCategory::sequential;

  // nnz_of_one_triangle has to be set if lower AND upper triangle
  // are passed. We will only touch the lower triangle
  explicit ICCPreconditioner(Op& op, double droptol=1.0e-2 , int modified = false , int nnz_of_one_triangle = -1 )
  {
    boost::timer::cpu_timer iluTimer;
    std::unique_ptr<BCRS_Matrix> const A(op.template getPointer<BCRS_Matrix>());
    int n = (*A).N() , nnz;
    if(nnz_of_one_triangle > 0) nnz = nnz_of_one_triangle;
    else nnz = (*A).nonzeroes();

    taucs_logfile((char*)"stdout");

    taucs_ccs_matrix *M = taucs_ccs_create(n,n,nnz,TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);

    std::vector<std::vector<int> > entries;
    entries.resize(n);
    for(int i=0;i<n;i++)
      for (ColIter cI=(*A)[i].begin(); cI!=(*A)[i].end() && cI.index()<=i; ++cI)
        entries[cI.index()].push_back(i);

    int idx_ = 0;
    for(int k=0;k<n;k++)
    {
      (M->colptr)[k] = idx_;
      for (std::vector<int>::iterator it = entries[k].begin(); it != entries[k].end(); it++)
      {
        (M->rowind)[idx_] = *it;
        (M->values).d[idx_] = (*A)[*it][k];
        idx_ += 1;
      }
    }
    (M->colptr)[n] = nnz;

    L = taucs_ccs_factor_llt(M, droptol, modified);

    rhs     = (double *)Malloc(n*sizeof(double), "main");
    sol     = (double *)Malloc(n*sizeof(double), "main");

    std::cout << "PrecondType::ICC time = " << (double)(iluTimer.elapsed().user)/1e9 << "s\n";
  }

  ~ICCPreconditioner()
  {
    taucs_ccs_free(L);
    free(sol);
    free(rhs);
  }

  virtual void pre (domain_type&, range_type&) {}
  virtual void post (domain_type&) {}

  virtual void apply (domain_type& x, range_type const& y) {
    y.write(rhs);
    taucs_ccs_solve_llt(L,sol,rhs);
    x.read(sol);
  }

private:
  taucs_ccs_matrix *L;
  double *data;
  double *sol, *rhs;
};

/* undef stupid macro names defined in taucs.h */
#undef max
#undef min
#undef realloc
#undef malloc
#undef calloc
#undef free

}  // namespace Kaskade
#endif
