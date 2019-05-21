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
 * \author Felix Lehmann
 */

#ifndef ICC_0_PRECOND_HH
#define ICC_0_PRECOND_HH

#include <memory>
#include <vector>

#include <boost/fusion/include/at_c.hpp>

#include "dune/grid/common/grid.hh"
#include "dune/istl/matrix.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/common/iteratorfacades.hh"
#include "dune/istl/matrixindexset.hh"

#include "dune/istl/preconditioners.hh"

namespace Kaskade
{

/**
 * \ingroup linalg
 * 
 * Incomplete Cholesky factorization by algorithm from book "Matrix Computations" by Gene Golub & Charles van Loan
 *
 * algorith will create approximate cholesky factorization L by only using the sparsity pattern of A
 *
 * we will only touch the lower triangle of the matrix
 *
 * interface copied from iluprecond.hh
 *
 * iccprecond.hh (using TAUCS Incomplete Cholesky Decomposition) does the same (even faster),
 * if we set its droptolerance parameter sufficiently large
 *
 */

template <class Op>
class ICC_0Preconditioner: public Dune::Preconditioner<typename Op::Range, typename Op::Range>
{
public:
  typedef typename Op::Domain Domain;
  typedef typename Op::Range Range;
  typedef typename Op::Scalar Scalar;

private:
  typedef Dune::FieldMatrix<Scalar,1,1> BlockType;
  typedef Dune::BlockVector<Dune::FieldVector<Scalar,1> > CoeffVector ;
  typedef Dune::BCRSMatrix<BlockType> BCRS_Matrix;
  typedef typename BCRS_Matrix::ConstColIterator ColIter;
  
    BCRS_Matrix mL;
    std::vector<std::vector<int> > mEntries;
    // int const mChoice;
  public:
    static int const category = Dune::SolverCategory::sequential;
    ICC_0Preconditioner( Op& op );
    void pre (Domain&, Range&) {}
    void apply (Domain& x, Range const& y);
    void post (Domain&) {}
};


template <class Op>
ICC_0Preconditioner<Op>::ICC_0Preconditioner(Op& op)
{
  std::unique_ptr<BCRS_Matrix> mA = op.template getPointer<BCRS_Matrix>();
  int N = (*mA).N();
  // create kind of column iterator for BCRS
  // and simultaneously read out the sparsity pattern of mA
  mEntries.resize(N);
  for(int i=0;i<N;i++)
    for (ColIter cI=(*mA)[i].begin(); cI!=(*mA)[i].end() && cI.index()<i; ++cI)
      mEntries[cI.index()].push_back(i);
  
  // backward substitution in ICC_0Preconditioner::apply() with BCRSMatrix mL^T is not
  // efficiently applicable, therefore create two BCRSMatrices mL and mLT := mL^T
  // set up sparsity pattern of lower triangle of mA (upper triangle of mA^T) for mL (mLT)
  Dune::MatrixIndexSet sparsityPattern_mL(N,N);
  for (int i=0; i<N; ++i)
    for (ColIter cI=(*mA)[i].begin(); cI!=(*mA)[i].end() && cI.index()<=i; ++cI)
      sparsityPattern_mL.add(i,cI.index() );
  
  sparsityPattern_mL.exportIdx( mL );
  for(int i=0;i<N;++i)
    for (ColIter cI=mL[i].begin(); cI!=mL[i].end(); ++cI)
      mL[i][cI.index()] = (*mA)[i][cI.index()];
  //find factorization mL of mA
  double tmp;
  for(int k=0;k<N;k++)
  {
    mL[k][k] = sqrt(mL[k][k]); 
    for(std::vector<int>::iterator it = mEntries[k].begin(); it!=mEntries[k].end(); it++)
      mL[*it][k] /= mL[k][k];
    for(std::vector<int>::iterator it = mEntries[k].begin(); it!=mEntries[k].end(); it++)
    {
      tmp = mL[*it][k];
      mL[*it][*it] -= (mL[*it][k])*(mL[*it][k]); // "mEntries" contains not diagonal indices
      for(std::vector<int>::iterator it2 = mEntries[*it].begin(); it2!=mEntries[*it].end(); it2++)
        if( mL.exists(*it2,k) )
	  mL[*it2][*it] -= mL[*it2][k]*tmp;
    }
  }  
}

template <class Op>
void ICC_0Preconditioner<Op>::apply(Domain& x, Range const& b)
{
  int N = x.dim();
  
  CoeffVector &sol = boost::fusion::at_c<0>(x.data);
  CoeffVector const &rhs = boost::fusion::at_c<0>(b.data);
  // for forward substitution
  CoeffVector forward(N); forward = 0;
  
  double tmp;
  // forward substitution
  for(int i=0;i<N;i++)
  {
    tmp = rhs[i];
    for (ColIter cI = mL[i].begin(); cI != mL[i].end(); ++cI)
      if( cI.index()<i )
	tmp -= mL[i][cI.index()]*forward[cI.index()];
    forward[i] = tmp/mL[i][i];
  }
  // backward substitution
  for(int i = N-1;i>=0;i--)
  {
    tmp = forward[i];
    for(std::vector<int>::iterator it = mEntries[i].begin();it!=mEntries[i].end();it++)
      tmp -= mL[*it][i]*sol[*it];
    sol[i] = tmp/mL[i][i]; 
  }
}
}  // namespace Kaskade
#endif
     
