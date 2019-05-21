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

#ifndef UMFPACK_SOLVE_HH
#define UMFPACK_SOLVE_HH

#include <vector>
#include <memory>

#include "linalg/factorization.hh"

namespace Kaskade
{

template<class Scalar,class SparseIndexInt, class DIL=int>
class UMFSymbolic
{
public:
  UMFSymbolic(std::vector<SparseIndexInt> const & Ap, std::vector<SparseIndexInt> const & Ai,
    std::vector<Scalar> const & Az, int verbose_ = 0);
  ~UMFSymbolic(); 

  void* getMem() {return mem;}
  DIL* getAp() { return &lAp[0]; }
  DIL* getAi() { return &lAi[0]; }

private:
  void* mem;
  std::vector<DIL> lAp;
  std::vector<DIL> lAi;
  int verbose;
};

/**
 * \cond internals
 */
template<class Scalar,class SparseIndexInt, class DIL=int>
class UMFFactorizationSpace
{
public:
  UMFFactorizationSpace(std::vector<SparseIndexInt> const& Ap, std::vector<SparseIndexInt> const& Ai,
                        std::vector<Scalar> const& Az);
  ~UMFFactorizationSpace();

  void* getMem() {return mem;}
  DIL* getAp() { return symbolic.getAp(); }
  DIL* getAi() { return symbolic.getAi(); }

private:
  UMFSymbolic<Scalar,SparseIndexInt,DIL> symbolic;
  void* mem;
};
/**
 * \endcond
 */

/**
 * \ingroup linalg
 * \brief Factorization of sparse linear systems with UMFPack
 *
 */
template <class Scalar, class SparseIndexInt = int, class DIL = int>
class UMFFactorization: public Factorization<Scalar,SparseIndexInt>
{
  using Factorization<Scalar,SparseIndexInt>::tripletToCompressedColumn;
public:

  ///Version of constructor keeping input data in triplet format (aka coordinate format) constant.
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  UMFFactorization(SparseIndexInt n_,
                   SparseIndexInt dummy,
                   std::vector<SparseIndexInt> const& ridx,
                   std::vector<SparseIndexInt> const& cidx,
                   std::vector<Scalar> const& values,
                   MatrixProperties property_=MatrixProperties::GENERAL) : N(ridx.size()), n(n_), Ap(n+1), Ai(N), Az(N), property(property_)
	{
	  assert(cidx.size()==N && values.size()==N);
	
	  // Transform triplet form into compressed column
	  std::vector<SparseIndexInt> map(N);
	  tripletToCompressedColumn(n, n, N, ridx, cidx, values, Ap, Ai, Az, map); 
	  factorization.reset(new UMFFactorizationSpace<Scalar,SparseIndexInt,DIL>(Ap,Ai,Az));
	}

  /// Version of constructor, that destroys input data before factorization: more memory efficient
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  UMFFactorization(SparseIndexInt n_,
                   std::unique_ptr<std::vector<SparseIndexInt> > const ridx,
                   std::unique_ptr<std::vector<SparseIndexInt> > const cidx,
                   std::unique_ptr<std::vector<Scalar> > const values,
                   MatrixProperties property_=MatrixProperties::GENERAL) : N(ridx->size()), n(n_), Ap(n+1), Ai(N), Az(N), property(property_)
  {
    assert(cidx->size()==N && values->size()==N);
    
    // Transform triplet form into compressed column
    std::vector<SparseIndexInt> map(N);
    tripletToCompressedColumn(n, n, N, ridx, cidx, values, Ap, Ai, Az, map); 
    
    // Clear memory of data that is not needed anymore.
    ridx.reset();
    cidx.reset();
    values.reset();
    { std::vector<SparseIndexInt> dump; std::swap(dump,map); }
    
    factorization.reset(new UMFFactorizationSpace<Scalar,SparseIndexInt,DIL>(Ap,Ai,Az));
  }
                                     

  /**
   * Solves the system for the given right hand side @arg b. @arg x is
   * resized to the number of matrix columns.
   */
  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transposed = false) const
	{
	  assert(b.size()>=n);
	  x.resize(n);
          solve_internal(&b[0],&x[0],transposed);
	}


  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, int nr, bool transposed = false) const
  {
    assert(b.size()>=n*nr);
    x.resize(nr*n);
    for(int i=0; i<nr; ++i)
      solve_internal(&b[i*n],&x[i*n],transposed);
  }

  /**
   * Solves the system for the given right hand side \arg b, which is
   * overwritten with the solution.
   */
  virtual void solve(std::vector<Scalar>& b) const
	{
	  std::vector<Scalar> x(b.size());
	  solve(b,x);
	  std::swap(b,x);
	}
  
private:

  virtual void solve_internal(Scalar const* bp, Scalar* xp, bool transposed) const;

  size_t const N, n; 
  std::vector<SparseIndexInt> Ap, Ai;
  std::vector<Scalar> Az;
  std::unique_ptr<UMFFactorizationSpace<Scalar,SparseIndexInt,DIL>> factorization;
  MatrixProperties property;
};

void umfpack_solve(std::vector<int> const& ridx,
    std::vector<int> const& cidx,
    std::vector<double> const& values,
    std::vector<double> const& b,
    std::vector<double>& x);

}  // namespace Kaskade
#endif
