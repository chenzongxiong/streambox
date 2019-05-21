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
#include <memory>
#include <cstdio>

#include "umfpack.h"

#include "linalg/umfpack_solve.hh"
#include "linalg/pardiso_solve.hh"
#include "linalg/triplet.hh"

#ifdef AIX
#define F77_FUNC(func)  func     
#else
#define F77_FUNC(func)  func ## _
#endif
extern "C"
{
  int F77_FUNC(pardisoinit)
    (void *, int *, int *);
}

extern "C"
{
int F77_FUNC(pardiso)
    (void *, int *, int *, int *, int *, int *, 
     double *, int *, int *, int *, int *, int *, 
     int *, double *, double *, int *);
}

namespace Kaskade
{

template<>
bool PARDISOFactorization<double,int,int>::first(true);
//std::vector<long int> PARDISOFactorization::pt(64);
std::vector<int> iparm0(64);

//---------------------------------------------------------------------


template<>
PARDISOFactorization<double,int,int>::~PARDISOFactorization()
{
    phase = 0;                 /* Release internal memory. */

    for(int i=0; i<Ap.size(); ++i) Ap[i]-=1;
    for(int i=0; i<Ai.size(); ++i) Ai[i]-=1;
     F77_FUNC(pardiso) (&pt[0], &maxfct, &mnum, &mtype, &phase,
 		       &n, &ddum, &Ap[0], &Ai[0], &idum, &nrhs,
 		       &iparm[0], &msglvl, &ddum, &ddum, &error);

    assert(!error);
}


template<>
PARDISOFactorization<double,int,int>::PARDISOFactorization(int n_,
                                           int mtype_,
                                           std::vector<int>const & rid,
                                           std::vector<int>const & cid,
                                           std::vector<double>const & val,
                                           MatrixProperties property)
  : N(rid.size()), n(n_), Ap(n+1), Ai(N), Az(N), maxfct(1), mnum(1), mtype(mtype_), nrhs(1), msglvl(0), error(0)
{
  std::vector<int> ridx(rid), cidx(cid);
  std::vector<double> values(val);
  assert(cidx.size()==N && values.size()==N);
  pt.resize(64);
  iparm.resize(64);
  for(int i=0; i<64; ++i) pt[i]=0;
  for(int i=0; i<64; ++i) iparm[i]=iparm0[i];

//    char *var = getenv("OMP_NUM_THREADS");
//    int num_procs;
//     if(var != NULL)
//         sscanf( var, "%d", &num_procs );
//     else {
//         printf("Set environment OMP_NUM_THREADS to 1");
//         assert(0);
//     }
  iparm[2]  = 1;
  mtype = 11;
  if (property==MatrixProperties::POSITIVEDEFINITE) { mtype = 2; deleteLowerTriangle(ridx,cidx,values); N=ridx.size(); }
  if(first)
  {
    F77_FUNC(pardisoinit) (&pt[0],  &mtype, &iparm0[0]);
  }
  first=false;

// Transform triplet form into compressed row (cidx and ridx swapped)
  std::vector<int> map(N);
  size_t nnz = N;
  tripletToCompressedColumn(n, n, nnz, cidx, ridx, values, Ap, Ai, Az, map); 
  for(int i=0; i<Ap.size(); ++i) Ap[i]+=1;
  for(int i=0; i<Ai.size(); ++i) Ai[i]+=1;
//   for(int i=0; i<N; ++i)
//     { printf("(%d,%d) %e\n", ridx[i]+1, cidx[i]+1, values[i]); }
//   for (int i=0; i<n; i++)
//     {
//       for (int k=Ap[i]; k<Ap[i+1]; k++)
//         printf("(%d,%d) %e\n", i+1, Ai[k-1], Az[k-1]);
//     }

  const char *s = "UNKNOWN";
  if (property==MatrixProperties::GENERAL) { s = "GENERAL"; mtype = 11; }
  else if (property==MatrixProperties::SYMMETRIC) { s = "MatrixProperties::SYMMETRIC"; mtype = 1;}
  else if (property==MatrixProperties::SYMMETRICSTRUCTURE) { s = "SYMMETRICSTRUCTURE"; mtype = 11;}
  else if (property==MatrixProperties::POSITIVEDEFINITE) { s = "MatrixProperties::POSITIVEDEFINITE"; mtype = 2; };
  if (this->getVerbose()>=2)
	{
	  std::cout << "PARDISO" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl; 
	}

  phase = 11; 
  
  F77_FUNC(pardiso) (&pt[0], &maxfct, &mnum, &mtype, &phase,
                     &n, &Az[0], &Ap[0], &Ai[0], &idum, &nrhs,
                     &iparm[0], &msglvl, &ddum, &ddum, &error);
  if (error) printf("phase 11, error=%d\n",error);
  assert(!error);
  phase = 22; 
  
  F77_FUNC(pardiso) (&pt[0], &maxfct, &mnum, &mtype, &phase,
                     &n, &Az[0], &Ap[0], &Ai[0], &idum, &nrhs,
                     &iparm[0], &msglvl, &ddum, &ddum, &error);
  if (error) printf("phase 22, error=%d\n",error);
  assert(!error);
};

template<>
void PARDISOFactorization<double,int,int>::solve(std::vector<double> const& b,
                                                 std::vector<double>& x, int nr, bool transposed) const
{
  assert(b.size()>=n*nr);
  x.resize(n*nr);
  assert(&x != &b);
  phase = 33;


  iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
  if(transposed) iparm[11]=1; else iparm[11]=0;

  F77_FUNC(pardiso) (&pt[0], &maxfct, &mnum, &mtype, &phase,
                     &n, &Az[0], &Ap[0], &Ai[0], &idum, &nr,
                     &iparm[0], &msglvl, const_cast<double*>(&b[0]), &x[0], &error);
  
  if (error) printf("phase 33, error=%d\n",error);
  assert(!error);
  
};


// Creator for DirectType::PARDISO factorizations to be plugged into the factory.

struct PARDISOCreator: public Creator<Factorization<double,int> > 
{
  PARDISOCreator(bool plugin) 
  {
    if (plugin)
      Factory<DirectType,Factorization<double,int> >::plugin(DirectType::PARDISO,this);
  }
  
  virtual std::unique_ptr<Factorization<double,int> >
  create(Creator<Factorization<double,int> >::Argument a) const 
  {
    MatrixAsTriplet<double,int> const& A = std::get<1>(a);
    assert(A.nrows()==A.ncols());
    MatrixProperties property = std::get<0>(a);

    return std::unique_ptr<Factorization<double,int> >
      (new PARDISOFactorization<double,int,int>(A.nrows(),0,A.ridx,A.cidx,A.data,property));
  }
};

PARDISOCreator pardisoCreator(true);

} // End of kaskade namespace ---------------------------------------------------------------------
