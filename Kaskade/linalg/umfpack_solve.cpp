/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>
#include <iostream>

#include "umfpack.h"
#include "linalg/umfpack_solve.hh"
#include "utilities/detailed_exception.hh"

//---------------------------------------------------------------------

namespace Kaskade
{

  template<>
  UMFSymbolic<double,int,int>::UMFSymbolic(std::vector<int> const & Ap, std::vector<int> const & Ai,
    std::vector<double> const & Az, int verbose_) : mem(0), verbose(verbose_)
  {
    if (verbose>=2)
      std::cout << "DirectType::UMFPACK solver, n=" << Ap.size()-1 << ", nnz=" << Ai.size() << std::endl;
    int status = umfpack_di_symbolic(Ap.size()-1,Ap.size()-1,&Ap[0],&Ai[0],&Az[0],&mem,0,0);
    umfpack_di_report_status(0,status);
    assert(mem);
  }
    
template<>
UMFSymbolic<double,int,long>::UMFSymbolic(std::vector<int> const & Ap, std::vector<int> const & Ai,
    std::vector<double> const & Az, int verbose_) : mem(0), verbose(verbose_)
    {
  if (verbose>=2)
    std::cout << "DirectType::UMFPACK3264 solver, n=" << Ap.size()-1 << ", nnz=" << Ai.size() << std::endl;
  lAp.resize(Ap.size());
  lAi.resize(Ai.size());
  int k;
  for (k=0; k<Ap.size(); k++) lAp[k] = Ap[k];
  for (k=0; k<Ai.size(); k++) lAi[k] = Ai[k];
  int status = umfpack_dl_symbolic(lAp.size()-1,lAp.size()-1,&lAp[0],&lAi[0],&Az[0],&mem,0,0);
  umfpack_di_report_status(0,status);
  assert(mem);
    }

template<>
UMFSymbolic<double,int,int>::~UMFSymbolic()
{
  umfpack_di_free_symbolic(&mem);
}

template<>
UMFSymbolic<double,int,long>::~UMFSymbolic()
{
  umfpack_dl_free_symbolic(&mem);
}

template<>
UMFSymbolic<double,long>::UMFSymbolic(std::vector<long> const & Ap, std::vector<long> const & Ai,
    std::vector<double> const & Az, int verbose_) : mem(0), verbose(verbose_)
{
  long status = umfpack_dl_symbolic(Ap.size()-1,Ap.size()-1,&Ap[0],&Ai[0],&Az[0],&mem,0,0);
  umfpack_dl_report_status(0,status);
  assert(mem);
  if (status != UMFPACK_OK)
    throw Kaskade::DirectSolverException("UMFPack failed.\n",__FILE__,__LINE__);
}

template<>
UMFSymbolic<double,long>::~UMFSymbolic()
{
  umfpack_dl_free_symbolic(&mem);
}

template<>
UMFFactorizationSpace<double,int,int>::UMFFactorizationSpace(std::vector<int> const& Ap,
    std::vector<int> const& Ai, std::vector<double> const& Az)
    : symbolic(Ap,Ai,Az), mem(0)
{
  int status = umfpack_di_numeric(&Ap[0],&Ai[0],&Az[0],symbolic.getMem(),&mem,0,0);
  umfpack_di_report_status(0,status);
  assert(mem);
  if (status != UMFPACK_OK)
    throw Kaskade::DirectSolverException("UMFPack failed.\n",__FILE__,__LINE__);
}

template<>
UMFFactorizationSpace<double,int,int>::~UMFFactorizationSpace()
{
  umfpack_di_free_numeric(&mem);
}

template<>
UMFFactorizationSpace<double,int,long>::UMFFactorizationSpace(std::vector<int> const& Ap,
    std::vector<int> const& Ai, std::vector<double> const& Az)
    : symbolic(Ap,Ai,Az), mem(0)
{
  int status = umfpack_dl_numeric(symbolic.getAp(),symbolic.getAi(),&Az[0],symbolic.getMem(),&mem,0,0);
  umfpack_dl_report_status(0,status);
  assert(mem);
  if (status != UMFPACK_OK)
    throw Kaskade::DirectSolverException("UMFPack failed.\n",__FILE__,__LINE__);
}

template<>
UMFFactorizationSpace<double,int,long>::~UMFFactorizationSpace()
{
  umfpack_dl_free_numeric(&mem);
}

template<>
UMFFactorizationSpace<double,long>::UMFFactorizationSpace(std::vector<long> const& Ap,
    std::vector<long> const& Ai, std::vector<double> const& Az)
    : symbolic(Ap,Ai,Az), mem(0)
{
  long status = umfpack_dl_numeric(&Ap[0],&Ai[0],&Az[0],symbolic.getMem(),&mem,0,0);
  umfpack_dl_report_status(0,status);
  assert(mem);
  if (status != UMFPACK_OK)
    throw Kaskade::DirectSolverException("UMFPack failed.\n",__FILE__,__LINE__);
}

template<>
UMFFactorizationSpace<double,long>::~UMFFactorizationSpace()
{
  umfpack_dl_free_numeric(&mem);
}

/// @TODO was ist hiermit??
// template<>
// void UMFFactorization<double,int,int>::solve(std::vector<double> const& b,
//                                              std::vector<double>& x, bool transposed) const
// 	{
// 	  assert(b.size()>=n);
// 	  x.resize(n);
// 	  assert(&x != &b);
// 	  int status = umfpack_di_solve(UMFPACK_A,&Ap[0],&Ai[0],&Az[0],&x[0],&b[0],factorization->getMem(),0,0);
// 	  umfpack_di_report_status(0,status);
// 	};


// template<>
// void UMFFactorization<double,int,long>::solve(std::vector<double> const& b,
//                                               std::vector<double>& x, bool transposed) const
// 	{
// 	  assert(b.size()>=n);
// 	  x.resize(n);
// 	  assert(&x != &b);
// 	  int status = umfpack_dl_solve(UMFPACK_A,factorization->getAp(),factorization->getAi(),
// 	                                &Az[0],&x[0],&b[0],factorization->getMem(),0,0);
// 	  umfpack_dl_report_status(0,status);
// 	};


// template<>
// void UMFFactorization<double,long>::solve(std::vector<double> const& b,
//                                           std::vector<double>& x, bool transposed) const
// 	{
// 	  assert(b.size()>=n);
// 	  x.resize(n);
// 	  assert(&x != &b);
// 	  long status = umfpack_dl_solve(UMFPACK_A,&Ap[0],&Ai[0],&Az[0],&x[0],&b[0],factorization->getMem(),0,0);
// 	  umfpack_dl_report_status(0,status);
// 	};

template<>
void UMFFactorization<double,int,int>::solve_internal(double const* bp, double* xp, bool transposed) const
{
  assert(xp != bp);
  if(transposed)
  {
  int status = umfpack_di_solve(UMFPACK_At,&Ap[0],&Ai[0],&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_di_report_status(0,status);
  } else
  {
  int status = umfpack_di_solve(UMFPACK_A,&Ap[0],&Ai[0],&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_di_report_status(0,status);
  }

}

template<>
void UMFFactorization<double,int,long>::solve_internal(double const* bp, double* xp, bool transposed) const
{
  assert(xp != bp);
  if(transposed)
  {
  int status = umfpack_dl_solve(UMFPACK_At,factorization->getAp(),factorization->getAi(),&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_dl_report_status(0,status);
  } else
  {
  int status = umfpack_dl_solve(UMFPACK_A,factorization->getAp(),factorization->getAi(),&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_dl_report_status(0,status);
  }

}

template<>
void UMFFactorization<double,long>::solve_internal(double const* bp, double* xp, bool transposed) const
{
  assert(xp != bp);
  if(transposed)
  {
  int status = umfpack_dl_solve(UMFPACK_At,&Ap[0],&Ai[0],&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_dl_report_status(0,status);
  } else
  {
  int status = umfpack_dl_solve(UMFPACK_A,&Ap[0],&Ai[0],&Az[0],xp,bp,factorization->getMem(),0,0);
  umfpack_dl_report_status(0,status);
  }
}


namespace {

  // Creator for UMFPack factorizations to be plugged into the factory.

  struct UMFPackCreator64: public Creator<Factorization<double,long> >
  {
    UMFPackCreator64(bool plugin)
    {
      if (plugin)
        Factory<DirectType,Factorization<double,long> >::plugin(DirectType::UMFPACK64,this);
    }

    virtual std::unique_ptr<Factorization<double,long> >
    create(Creator<Factorization<double,long> >::Argument a) const
    {
      MatrixAsTriplet<double,long> const& A = std::get<1>(a);
      assert(A.nrows()==A.ncols());

      return std::unique_ptr<Factorization<double,long> >(new UMFFactorization<double,long>(A.nrows(),0,A.ridx,A.cidx,A.data));
    }
  };

  struct UMFPackCreator: public Creator<Factorization<double,int> >
  {
    UMFPackCreator(bool plugin)
    {
      if (plugin)
        Factory<DirectType,Factorization<double,int> >::plugin(DirectType::UMFPACK,this);
    }

    virtual std::unique_ptr<Factorization<double,int> >
    create(Creator<Factorization<double,int> >::Argument a) const
    {
      MatrixAsTriplet<double,int> const& A = std::get<1>(a);
      assert(A.nrows()==A.ncols());

      return std::unique_ptr<Factorization<double,int> >(new UMFFactorization<double,int>(A.nrows(),0,A.ridx,A.cidx,A.data));
    }
  };

  struct UMFPackCreator3264: public Creator<Factorization<double,int> >
  {
    UMFPackCreator3264(bool plugin)
    {
      if (plugin)
        Factory<DirectType,Factorization<double,int> >::plugin(DirectType::UMFPACK3264,this);
    }

    virtual std::unique_ptr<Factorization<double,int> >
    create(Creator<Factorization<double,int> >::Argument a) const
    {
      MatrixAsTriplet<double,int> const& A = std::get<1>(a);
      assert(A.nrows()==A.ncols());

      return std::unique_ptr<Factorization<double,int> >(new UMFFactorization<double,int,long>(A.nrows(),0,A.ridx,A.cidx,A.data));
    }
  };

  UMFPackCreator64 umfpackCreator64(true);
  UMFPackCreator3264 umfpackCreator3264(true);
  UMFPackCreator umfpackCreator(true);

} // End of anonymous namespace ---------------------------------------------------------------------

void umfpack_solve(std::vector<int> const& ridx,
    std::vector<int> const& cidx,
    std::vector<double> const& values,
    std::vector<double> const& b,
    std::vector<double>& x)
{
  size_t const N = ridx.size();
  size_t const n = b.size();
  assert(cidx.size()==N && values.size()==N);
  assert(x.size()>=n);

  // Transform triplet form into compressed column
  std::vector<int> Ap(n+1), Ai(N), map(N);
  std::vector<double> Az(N);
  umfpack_di_triplet_to_col(n,n,N,
      &ridx[0],&cidx[0],&values[0],
      &Ap[0],&Ai[0],&Az[0],
      &map[0]);

  // Solve system
  UMFFactorizationSpace<double,int> factorization(Ap,Ai,Az);
  umfpack_di_solve(UMFPACK_A,&Ap[0],&Ai[0],&Az[0],&x[0],&b[0],factorization.getMem(),0,0);
}

}  // namespace Kaskade
