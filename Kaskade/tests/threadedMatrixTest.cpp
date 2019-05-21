/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <numeric>

#include <boost/timer/timer.hpp>

#include <dune/common/fmatrix.hh>
#include "dune/grid/config.h"

#include "linalg/threadedMatrix.hh"


using namespace Kaskade;

// Create a NxN lower triangular matrix with 1x2 entries 
// as a mockup of the divergence discretization in Stokes problems
template <class Index>
void divergence(Index const N)
{
  std::cout << "running divergence with N=" << N << "\n";
  
  // First, sparsity pattern creation. We proceed row by row.
  NumaCRSPatternCreator<Index> creator(N,N,false);
  
  boost::timer::cpu_timer timer;
  std::vector<std::vector<Index>> cols(N), rows(N);
  for (Index k=0; k<cols.size(); ++k)
  {
    cols[k].resize(k+1);
    rows[k].resize(1);
    std::iota(begin(cols[k]),end(cols[k]),0);
    rows[k][0] = k;
  }
  std::cerr << "en bloc filling of indices: " << timer.format() << "\n";
  timer.start();
  creator.addElements(rows,cols);
  std::cerr << "en bloc adding elements:    " << timer.format() << "\n";
  
  // Next, create a sparsity pattern from the creator and a matrix.
  timer.start();
  creator.balance();
  std::cerr << "balance:                    " << timer.format() << "\n";
  timer.start();
  std::shared_ptr<NumaCRSPattern<Index>> pattern(new NumaCRSPattern<Index>(creator));
  std::cerr << "pattern creation:           " << timer.format() << "\n";
  std::cout << "nonzero entries: " << pattern->nonzeroes() << "\n";
  std::cout << "storage size:    " << pattern->storage() << "\n";
  for (int i=0; i<pattern->nodes(); ++i)
  {
    auto const& cp = *pattern->pattern(i);
    std::cout << "chunk " << i << ": [" << cp.first() << "," << cp.last() << "[  nnz=" << cp.nonzeroes() << "  size=" << cp.storage() << "\n";
  }
  timer.start();
  NumaBCRSMatrix<Dune::FieldMatrix<double,1,2>,Index> matrix(pattern);
  std::cerr << "matrix creation:            " << timer.format() << "\n";
  
  // Fill the matrix with ones.
  for (auto r=matrix.begin(); r!=matrix.end(); ++r)
    for (auto c=r->begin(); c!=r->end(); ++c)
      *c = 1.0;
  
  Dune::BlockVector<Dune::FieldVector<double,2>> x(N);
  Dune::BlockVector<Dune::FieldVector<double,1>> y(N);
  for (int i=0; i<N; ++i)
  {
    x[i] = 1;
    y[i] = 0;
  }
  
  timer.start();
  // Compute 1000 matrix-vector products
  for (int i=0; i<1000; ++i) 
    matrix.umv(x,y);
  std::cout << "1000 Ax: " << timer.format() << "\n";
  
  timer.start();
  for (int i=0; i<1000; ++i) 
    matrix.umtv(x,y);
  std::cout << "1000 A^Tx: " << timer.format() << "\n";
  std::cout << "-----------------------------\n\n";
}



// Create a NxN Toeplitz matrix of bandwidth 5
template <class Index>
void toeplitz(Index const N, bool symmetric)
{
  std::cout << "running toeplitz with symmetric=" << symmetric << "\n";
  
  // First, sparsity pattern creation. The pattern is a superposition of 
  // 5x5 blocks on the diagonal. To begin with, we add a couple of those blocks
  // individually.
  NumaCRSPatternCreator<Index> creator(N,N,symmetric,14);
  
  std::vector<int> idx(5);
  for (int k=0; k<100; ++k)
  {
    std::iota(begin(idx),end(idx),k);
    creator.addElements(begin(idx),end(idx),begin(idx),end(idx));
  }
  
  // The rest is filled en bloc.
  boost::timer::cpu_timer timer;
  std::vector<std::vector<Index>> cols(N-104), rows(N-104);
  for (Index k=0; k<cols.size(); ++k)
  {
    cols[k].resize(5);
    rows[k].resize(5);
    std::iota(begin(cols[k]),end(cols[k]),k+100);
    std::iota(begin(rows[k]),end(rows[k]),k+100);
  }
  std::cerr << "en bloc filling of indices: " << timer.format() << "\n";
  timer.start();
  creator.addElements(rows,cols);
  std::cerr << "en bloc adding elements:    " << timer.format() << "\n";
  
  // Next, create a sparsity pattern from the creator and a matrix.
  timer.start();
  creator.balance();
  std::cerr << "balance:                    " << timer.format() << "\n";
  timer.start();
  std::shared_ptr<NumaCRSPattern<Index>> pattern(new NumaCRSPattern<Index>(creator));
  std::cerr << "pattern creation:           " << timer.format() << "\n";
  std::cout << "nonzero entries: " << pattern->nonzeroes() << "\n";
  std::cout << "storage size:    " << pattern->storage() << "\n";
  for (int i=0; i<pattern->nodes(); ++i)
  {
    auto const& cp = *pattern->pattern(i);
    std::cout << "chunk " << i << ": [" << cp.first() << "," << cp.last() << "[  nnz=" << cp.nonzeroes() << "  size=" << cp.storage() << "\n";
  }
  timer.start();
  NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>,Index> matrix(pattern);
  std::cerr << "matrix creation:            " << timer.format() << "\n";
  
  // Fill the matrix with ones.
  for (auto r=matrix.begin(); r!=matrix.end(); ++r)
    for (auto c=r->begin(); c!=r->end(); ++c)
      *c = 1.0;
  
  Dune::BlockVector<Dune::FieldVector<double,1>> x(N), y(N);
  for (int i=0; i<N; ++i)
    x[i] = 1;
  
  timer.start();
  // Compute maximum eigenvalue by power method
  double mx = 0;
  double dp;
  for (int i=0; i<2000; ++i) 
  {
    dp = matrix.smv(.125,x,y);
    if (i%100==0 || i>1995)
    {
      mx = 0;
      for (int i=0; i<N; i+=16)
        mx += std::abs(y[i][0]);
      for (int i=0; i<N; ++i)
        x[i] = y[i]/mx;
    }
  }
  std::cout << "2000 iterations power method: " << timer.format() << "\n";
  std::cout.precision(8);
  std::cout << "max eigenvalue: " << 8*mx << "\n";
  std::cout << "dp: " << dp << "\n";
  
  matrix = 1.0;
  
  if (!symmetric)
  {
    for (int i=0; i<N; ++i)
      x[i] = 1;
    timer.start();
    for (int i=0; i<2000; ++i) 
    {
      matrix.smtv(.125,x,y);
      if (i%100==0 || i>1995)
      {
        mx = 0;
        for (int i=0; i<N; i+=16)
          mx += std::abs(y[i][0]);
        for (int i=0; i<N; ++i)
          x[i] = y[i]/mx;
      }
    }
    std::cout << "2000 iterations power method: " << timer.format() << "\n";
    std::cout.precision(8);
    std::cout << "max eigenvalue: " << 8*mx << "\n";
  }
  std::cout << "-----------------------------\n\n";
}

int main(void) 
{
  int N = 64000;
  toeplitz(N,false);
  toeplitz(static_cast<size_t>(N),true);
  divergence(2000);
  return 0;
}