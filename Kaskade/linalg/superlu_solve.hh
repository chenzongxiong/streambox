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

#ifndef SUPERLU_SOLVE_HH
#define SUPERLU_SOLVE_HH

#include <vector>
#include <iostream>
#include <memory>

#include "slu_ddefs.h"

#include "linalg/factorization.hh"

namespace Kaskade
{

/**
 * \ingroup linalg
 * \brief Factorization of sparse linear systems with mumps
 *
 */
template <class Scalar,class SparseIndexInt=int, class DIL=int>
class SUPERLUFactorization: public Factorization<Scalar,SparseIndexInt>
{
  using Factorization<Scalar,SparseIndexInt>::tripletToCompressedColumn;
public:

  ///Version of constructor keeping input data in triplet format (aka coordinate format) constant.
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  SUPERLUFactorization(SparseIndexInt n_,
                   SparseIndexInt dummy,
                   std::vector<SparseIndexInt> const& ridx,
                   std::vector<SparseIndexInt> const& cidx,
                   std::vector<Scalar> const& values)
  : N(ridx.size()), n(n_)
	{
	  assert(cidx.size()==N && values.size()==N);

      verbose = Factorization<Scalar>::getVerbose();

	  if (this->getVerbose()>=2)
	    {
	      std::cout << "SuperLU" << " solver, n=" << n << ", nnz=" << N << std::endl; 
	    }

      set_default_options(&options);
      StatInit(&stat);

      work = 0;
      lwork = 0;
      u = 1.0;
      equil = YES;
      trans = NOTRANS;

	  std::vector<SparseIndexInt> map(N);
	  tripletToCompressedColumn(n, n, N, ridx, cidx, values, Ap, Ai, Az, map); 

	  dCreate_CompCol_Matrix(&A, n, n, N, &Az[0], &Ai[0], &Ap[0], SLU_NC, SLU_D, SLU_GE);
    
      nrhs = 1;
      if ( !(rhsb = doubleMalloc(n * nrhs)) ) ABORT("Malloc fails for rhsb[].");
      if ( !(rhsx = doubleMalloc(n * nrhs)) ) ABORT("Malloc fails for rhsx[].");
      dCreate_Dense_Matrix(&B, n, nrhs, rhsb, n, SLU_DN, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&X, n, nrhs, rhsx, n, SLU_DN, SLU_D, SLU_GE);
      xact = doubleMalloc(n * nrhs);
      ldx = n;
      dGenXtrue(n, nrhs, xact, ldx);
      dFillRHS(trans, nrhs, xact, ldx, &A, &B);

      if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
      if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
      if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
      if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
      if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
      int nrhs = 1;
      if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
      if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

      options.Equil = equil;
      options.DiagPivotThresh = u;
      options.Trans = trans;

      B.ncol = 0;  /* Indicate not to solve the system */
      dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
             &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
             &mem_usage, &stat, &info);
      B.ncol = nrhs;  /* Set the number of right-hand side */

      if ( (info == 0 || info == n+1) && (verbose>0) ) {

	  if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
	  if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
      Lstore = (SCformat *) L.Store;
      Ustore = (NCformat *) U.Store;
	  printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %ld\n", Lstore->nnz + Ustore->nnz - n);
      printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/N);

	  printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	  mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	  fflush(stdout);

      } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %ld bytes\n", info - n);
      } else
      if (info!=0)
        printf("LU factorization: dgssvx() returns info %d\n", info);

      if ( verbose>0 ) StatPrint(&stat);
      StatFree(&stat);
   };

  /// Version of constructor, that destroys input data before factorization: more memory efficient
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  SUPERLUFactorization(SparseIndexInt n_,
                   SparseIndexInt dummy,
                   std::unique_ptr<std::vector<SparseIndexInt> > ridx,
                   std::unique_ptr<std::vector<SparseIndexInt> > cidx,
                   std::unique_ptr<std::vector<Scalar> > values)
	: N(ridx->size()), n(n_)
	{
	  assert(cidx->size()==N && values->size()==N);

	  if (this->getVerbose()>=2)
	    {
	      std::cout << "SuperLU" << " solver, n=" << n << ", nnz=" << N << std::endl; 
	    }

      verbose = Factorization<Scalar>::getVerbose();
      set_default_options(&options);
      StatInit(&stat);

      work = 0;
      lwork = 0;
      u = 1.0;
      equil = YES;
      trans = NOTRANS;

	  std::vector<SparseIndexInt> map(N);
	  tripletToCompressedColumn(n, n, N, ridx, cidx, values, Ap, Ai, Az, map); 

	  dCreate_CompCol_Matrix(&A, n, n, N, &Az[0], &Ap[0], &Ai[0], SLU_NC, SLU_D, SLU_GE);
    
      nrhs = 1;
      if ( !(rhsb = doubleMalloc(n * nrhs)) ) ABORT("Malloc fails for rhsb[].");
      if ( !(rhsx = doubleMalloc(n * nrhs)) ) ABORT("Malloc fails for rhsx[].");
      dCreate_Dense_Matrix(&B, n, nrhs, rhsb, n, SLU_DN, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&X, n, nrhs, rhsx, n, SLU_DN, SLU_D, SLU_GE);
      xact = doubleMalloc(n * nrhs);
      ldx = n;
      dGenXtrue(n, nrhs, xact, ldx);
      dFillRHS(trans, nrhs, xact, ldx, &A, &B);

      if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
      if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
      if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
      if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
      if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
      if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
      if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

      options.Equil = equil;
      options.DiagPivotThresh = u;
      options.Trans = trans;

      SuperMatrix B, X;
      B.ncol = 0;  /* Indicate not to solve the system */
      dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
             &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
             &mem_usage, &stat, &info);
      B.ncol = nrhs;  /* Set the number of right-hand side */

      if ( (info == 0 || info == n+1) && (verbose>0) ) {
	  if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
	  if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
      Lstore = (SCformat *) L.Store;
      Ustore = (NCformat *) U.Store;
	  printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
      printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/N);

	  printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	  mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	  fflush(stdout);

      } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
      }
      if (info!=0)
        printf("LU factorization: dgssvx() returns info %d\n", info);

      if ( verbose>0 ) StatPrint(&stat);
      StatFree(&stat);
	};
                                     

  ~SUPERLUFactorization()
    {
      SUPERLU_FREE (etree);
      SUPERLU_FREE (perm_r);
      SUPERLU_FREE (perm_c);
      SUPERLU_FREE (R);
      SUPERLU_FREE (C);
      SUPERLU_FREE (ferr);
      SUPERLU_FREE (berr);

      SUPERLU_FREE (rhsb);
      SUPERLU_FREE (rhsx);
      SUPERLU_FREE (xact);

//       Destroy_CompCol_Matrix(&A); 
//       Destroy_SuperMatrix_Store(&B); 
//       Destroy_SuperMatrix_Store(&X); 
      Destroy_SuperNode_Matrix(&L); 
      Destroy_CompCol_Matrix(&U);    }

  /**
   * Solves the system for the given right hand side @arg b. @arg x is
   * resized to the number of matrix columns.
   */
  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transpose=false) const
    {
      options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */

    /* Initialize the statistics variables. */
      StatInit(&stat);

      int k;
      for (k=0; k<n; k++) rhsb[k] = b[k];
      dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
             &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
             &mem_usage, &stat, &info);
      for (k=0; k<n; k++) x[k] = rhsx[k];

      if (info!=0)
        printf("Triangular solve: dgssvx() returns info %d\n", info);

      if ( (info == 0 || info == n+1) && (verbose>0) ) {

	    if ( options.IterRefine ) {
          printf("Iterative Refinement:\n");
	      printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	      for (int i = 0; i < nrhs; ++i)
	        printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
	   }
	  fflush(stdout);
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %ld bytes\n", info - n);
      }

      if ( verbose>0 ) StatPrint(&stat);
      StatFree(&stat);
    }

  /**
   * Solves the system for the given right hand side \arg b, which is
   * overwritten with the solution.
   */
  void solve(std::vector<Scalar>& b) const
	{
      options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */

    /* Initialize the statistics variables. */
      StatInit(&stat);

      int k;
      for (k=0; k<n; k++) rhsb[k] = b[k];
      dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
             &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
             &mem_usage, &stat, &info);
      for (k=0; k<n; k++) b[k] = rhsx[k];

      if (info!=0)
        printf("Triangular solve: dgssvx() returns info %d\n", info);

      if ( (info == 0 || info == n+1) && (verbose>0) ) {

	    if ( options.IterRefine ) {
          printf("Iterative Refinement:\n");
	      printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	      for (int i = 0; i < nrhs; ++i)
	        printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i]);
	   }
	  fflush(stdout);
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %ld bytes\n", info - n);
      }

      if ( verbose>0 ) StatPrint(&stat);
      StatFree(&stat);
	}
  
private:
  size_t const N, n; 
  std::vector<SparseIndexInt> Ap, Ai;
  std::vector<Scalar> Az;
  mutable superlu_options_t options;
  mutable int nrhs, ldx;
  mutable SuperMatrix A, L, U;
  mutable NCformat       *Ustore;
  mutable SCformat       *Lstore;
  mutable SuperMatrix    B, X;
  mutable int *perm_c; /* column permutation vector */
  mutable int *perm_r; /* row permutations from partial pivoting */
  mutable int *etree;
  mutable double *R, *C;
  mutable double *ferr, *berr;
  mutable double *rhsb, *rhsx, *xact;
  mutable mem_usage_t mem_usage;
  mutable SuperLUStat_t stat;
  mutable char equed[1];
  mutable void *work;
  mutable int info, lwork;
  mutable double u, rpg, rcond;
  mutable yes_no_t equil;
  mutable trans_t trans;
  int verbose;
};
}  // namespace Kaskade
#endif
