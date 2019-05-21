#ifndef __ACML_TO_STDLAPACK__
#define __ACML_TO_STDLAPACK__

#include <algorithm>
#include <cstdlib>     // malloc, with clang++
#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#if defined(Darwin) || defined(Atlas)
#include <clapack.h>
#endif
#if defined(Cygwin) || defined(Atlas)
#include "lapack/lapacke_utils.h"
#define __CLPK_integer int
int ilaenv_(int*,char*,int,char*,int,int*,int*,int*,int*);
#endif
#ifdef DEBUG
#define VDEBUG 1
#else
#define VDEBUG 0
#endif
#include <cstdio>

void dgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc)
{
  enum CBLAS_TRANSPOSE ntransa=CblasNoTrans, ntransb=CblasNoTrans;
  if ( transa == 't' ) ntransa=CblasTrans;
  if ( transb == 't' ) ntransb=CblasTrans;
  if (VDEBUG) printf("***entering cblas_dgemm ***\n");
  cblas_dgemm(CblasColMajor,ntransa,ntransb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc); 
  if (VDEBUG) printf("***dgemm finished***\n");
};

void dgemv(char transa, int m, int n, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy) 
{
  enum CBLAS_TRANSPOSE ntransa=CblasNoTrans;
  if ( transa == 't' ) ntransa=CblasTrans;
  if (VDEBUG) printf("*** entering cblas_dgemv***\n");
  cblas_dgemv(CblasColMajor,ntransa,m,n,alpha,a,lda,x,incx,beta,y,incy);
  if (VDEBUG) printf("***dgemv finished***\n");
};

void dgels(char trans, int mA, int nA, int nrhsA, double *a, int ldaA, double *b, int ldbA, int *infoA) 
{
  __CLPK_integer m=mA, n=nA, nrhs=nrhsA, lda=ldaA, ldb=ldbA, info=*infoA;
  __CLPK_integer mn = std::min(m,n);
  __CLPK_integer lwork = std::max( 1, mn + std::max( mn, nrhs )*n );
  double *work  = (double*) malloc(lwork*sizeof(double));
  if (VDEBUG) printf("*** entering dgels ***\n");
  dgels_(&trans,&m,&n,&nrhs,a,&lda,b,&ldb,work,&lwork,&info);
  *infoA=info;
  free(work);
  if (VDEBUG) printf("***dgels finished***\n");
};

void dsyevd(char jobz, char uplo, int nA, double *a, int ldaA, double *w, int *infoA) 
{
  __CLPK_integer n=nA, lda=ldaA, info=*infoA;
  __CLPK_integer lwork=1+5*n+5*n*n, liwork=3+5*n;
  double *work = (double*) malloc(lwork*sizeof(double));
  __CLPK_integer *iwork   = (__CLPK_integer*) malloc(liwork*sizeof(int));
  if (VDEBUG) printf("*** entering dsyevd***\n");
  dsyevd_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,iwork,&liwork,&info);
  *infoA=info;
  free(work); free(iwork);
  if (VDEBUG) printf("***dsyevd finished***\n");
};

void dgelsy(int mA, int nA, int nrhsA, double *a, int ldaA, double *b, int ldbA, int *jpvtA, double rcond, int *rankA, int *infoA) 
{
  __CLPK_integer m=mA, n=nA, nrhs=nrhsA, lda=ldaA, ldb=ldbA, *jpvt,
                 rank=*rankA, info=*infoA;
  __CLPK_integer mn = std::min(m,n);
  __CLPK_integer lwork = std::max( mn+2*n+n*(n+1), 2*mn+n*nrhs );
  int i;
  double *work  = (double*) malloc(lwork*sizeof(double));
  jpvt = (__CLPK_integer*)malloc(n*sizeof(__CLPK_integer));
  for (i=0;i<n;i++) jpvt[i]=jpvtA[i];
  if (VDEBUG) printf("*** entering dgelsy ***\n");
  dgelsy_(&m,&n,&nrhs,a,&lda,b,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);
  for (i=0;i<n;i++) jpvtA[i]=jpvt[i];
  *infoA=info;
  *rankA=rank;
  free(work); free(jpvt);
  if (VDEBUG) printf("***dgelsy finished***\n");
};

void dgesv(int nA, int nrhsA, double *a, int ldaA, int *pivA, double *b, int ldbA, int *infoA)
{
  __CLPK_integer n=nA, nrhs=nrhsA, lda=ldaA, *ipiv=pivA, ldb=ldbA;
  __CLPK_integer info;
  info=*infoA;
  if (VDEBUG) printf("*** entering dgesv ***\n");
  dgesv_(&n,&nrhs,a,&lda,ipiv,b,&ldb,&info);
 *infoA=info;
  if (VDEBUG) printf("***dgesv finished***\n");
};

void dgesvd(char jobu, char jobvt, int mA, int nA, double *a, int ldaA, double *sing, double *u, int lduA, double *vt, int ldvtA, int *infoA)
{
  __CLPK_integer m=mA, n=nA, lda=ldaA, ldu=lduA, ldvt=ldvtA,info=*infoA;
  __CLPK_integer lwork = std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n));
  double *work  = (double*) malloc(lwork*sizeof(double));
  if (VDEBUG) printf("*** entering dgesvd ***\n");
  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,sing,u,&ldu,vt,&ldvt,work,&lwork,&info);
  *infoA=info;
  free(work);
  if (VDEBUG) printf("***dgesvd finished***\n");
};

void dstev(char jobz, int nA, double *diag, double *offd, double *z, int ldzA, int *infoA)
{
  __CLPK_integer n=nA, ldz=ldzA, info=*infoA;
  __CLPK_integer lwork = std::max(1,2*n-2);
  double *work  = (double*) malloc(lwork*sizeof(double));
  if (VDEBUG) printf("*** entering dstev ***\n");
  dstev_(&jobz,&n,diag,offd,z,&ldz,work,&info);
  *infoA=info;
  free(work);
  if (VDEBUG) printf("***dstev finished***\n");
};

void dgetrf(int mA, int nA, double *a, int ldaA, int *ipiv, int *infoA)
{
  __CLPK_integer m=mA, n=nA, lda=ldaA,info=*infoA;
  dgetrf_( &m, &n, a, &lda, ipiv, &info );
  *infoA=info;
};

void dgetri(int nA, double *a, int ldaA, int *ipiv, int *infoA)
{
  __CLPK_integer n=nA, lda=ldaA, info=*infoA;
  int minusOne=-1, one=1;
#if defined(Cygwin) || defined(Atlas)
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "dgetri",6, "",0, (int*) &n, &minusOne, &minusOne, &minusOne );
#else
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "dgetri", "", (int*) &n, &minusOne, &minusOne, &minusOne );
#endif
  double *work  = (double*) malloc(lwork*sizeof(double));
  if (VDEBUG) printf("*** entering dgetri ***\n");
  dgetri_( &n, a, &lda, ipiv, work, &lwork, &info );
  *infoA=info;
  free(work);
  if (VDEBUG) printf("***dgetri finished***\n");
};

void sgetrf(int mA, int nA, float *a, int ldaA, int *ipiv, int *infoA)
{
  __CLPK_integer m=mA, n=nA, lda=ldaA,info=*infoA;
  sgetrf_( &m, &n, a, &lda, ipiv, &info );
  *infoA=info;
};

void sgetri(int nA, float *a, int ldaA, int *ipiv, int *infoA)
{
  __CLPK_integer n=nA, lda=ldaA, info=*infoA;
  int minusOne=-1, one=1;
#if defined(Cygwin) || defined(Atlas)
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "sgetri",6, "",0, (int*) &n, &minusOne, &minusOne, &minusOne );
#else
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "sgetri", "", (int*) &n, &minusOne, &minusOne, &minusOne );
#endif
  float *work  = (float*) malloc(lwork*sizeof(float));
  if (VDEBUG) printf("*** entering sgetri ***\n");
  sgetri_( &n, a, &lda, ipiv, work, &lwork, &info );
  *infoA=info;
  free(work);
  if (VDEBUG) printf("***sgetri finished***\n");
};

void sgemv(char transa, int m, int n, float alpha, float *a, int lda, float *x, int incx, float beta, float *y, int incy) 
{
  enum CBLAS_TRANSPOSE ntransa=CblasNoTrans;
  if ( transa == 't' ) ntransa=CblasTrans;
  if (VDEBUG) printf("*** entering cblas_sgemv***\n");
  cblas_sgemv(CblasColMajor,ntransa,m,n,alpha,a,lda,x,incx,beta,y,incy);
  if (VDEBUG) printf("***sgemv finished***\n");
};

// example for a new placeholder routine - use copy/paste, then adapt the name and parameter list,
// and uncomment the code
//
// void sgetrf(int m, int n, float *a, int lda, int *ipiv, int *info)
// {
//   printf("*** sgetrf: Please supply mapping routine to standard LAPACK routine sgetrf_ ***\n");
//   printf("*** by editing this routine in header file linalg/acml_to_stdlapack.hh\n");
//   printf("*** Execution of program aborted! ***\n");
//   exit(10);
// };
// 

#ifdef __cplusplus
}
#endif

#endif
