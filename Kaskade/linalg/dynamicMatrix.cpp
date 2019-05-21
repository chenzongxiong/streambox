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

// Access on acml include file in dependence of the platform HTC-Linux or Mac OS.
// ZIBHTC is to set in the Makefile.Local (copy of Makefile-htc.Local) on HTC platform.
#include <vector>

#ifdef ZIBHTC
#include "acml.h"
#endif

#if defined(Darwin) || defined(Cygwin) || defined(Atlas)
#include "acml_to_stdlapack.hh"
#endif


#include "linalg/dynamicMatrix.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade 
{
  
  namespace DynamicMatrixDetail
  {
    int invertLAPACK(int n, double* A, int lda)
    {
      int info;
      std::vector<int> ipiv(n);
      dgetrf(n,n,A,lda,&ipiv[0],&info);
      assert(info>=0);
      if (info>0)
        throw SingularMatrixException("zero pivot",__FILE__,__LINE__);
      dgetri(n,A,lda,&ipiv[0],&info);
      assert(info>=0);
      return info;
    }
    
    int invertLAPACK(int n, float* A, int lda)
    {
      int info;
      std::vector<int> ipiv(n);
      sgetrf(n,n,A,lda,&ipiv[0],&info);
      assert(info>=0);
      if (info>0)
        throw SingularMatrixException("zero pivot",__FILE__,__LINE__);
      sgetri(n,A,lda,&ipiv[0],&info);
      assert(info>=0);
      return info;
    }
    
    void mv(int n, int m, double const* A, int lda, double const* x, double* y)
    {
      dgemv('N',n,m,1.0,const_cast<double*>(A),lda,const_cast<double*>(x),1,0.0,y,1); 
    }
    
    void mv(int n, int m, float  const* A, int lda, float  const* x, float * y)
    {
      sgemv('N',n,m,1.0,const_cast<float*>(A),lda,const_cast<float*>(x),1,0.0,y,1); 
    }

  }
  
// ----------------------------------------------------------------------------
  
  template <>
  std::tuple<DynamicMatrix<Dune::FieldMatrix<double,1,1>>,std::vector<double>,DynamicMatrix<Dune::FieldMatrix<double,1,1>>>
  svd(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_)
  {
    int info;
    DynamicMatrix<Dune::FieldMatrix<double,1,1>> A(A_), U(A.rows(),A.rows()), VT(A.cols(),A.cols());
    std::vector<double> sing(std::min(A.cols(),A.rows()));
    
    dgesvd('A', 'A', A.rows(), A.cols(), A.data(), A.lda(), &sing[0], U.data(), U.lda(), VT.data(), VT.lda(), &info);
    // todo: check info
    
    // LAPACK computes V^T, but we return V. Transpose.
    for (int i=1; i<VT.rows(); ++i)
      for (int j=0; j<i; ++j)
        std::swap(VT[i][j],VT[j][i]);
      
    return std::make_tuple(U,sing,VT);
  }
  
  // explicit instantiation
  template std::tuple<DynamicMatrix<Dune::FieldMatrix<double,1,1>>,std::vector<double>,DynamicMatrix<Dune::FieldMatrix<double,1,1>>>
  svd(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_);
    
  
// ----------------------------------------------------------------------------
  
  
  template <>
  Dune::DynamicVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                        Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_)
  {
    int info;
    auto A = A_;
    auto b = b_;
    assert(A.rows()==A.cols());
    std::vector<int> ipiv(A.rows());
    dgesv(A.rows(),1,A.data(),A.lda(),&ipiv[0],&b[0][0],1,&info);
    // todo: check info
    return b;
  }
  
  // explicit instantiation
  template Dune::DynamicVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                                 Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_);
  
  
// ----------------------------------------------------------------------------
  
  template <class Scalar>
  void invert(DynamicMatrix<Scalar>& A)
  {
    DynamicMatrixDetail::invertLAPACK(A.N(),A.data(),A.lda());
  }
  // explicit instantiation
  template void invert(DynamicMatrix<double>& A);
  template void invert(DynamicMatrix<float>& A); 
}




// ----------------------------------------------------------------------------

#ifdef UNITTEST

#include <algorithm>
#include <iostream>

int main(void) 
{
  Kaskade::DynamicMatrix<Dune::FieldMatrix<double,1,1>> A(3,5);
  for (int i=0; i<A.rows(); ++i)
    for (int j=0; j<A.cols(); ++j) 
      A[i][j] = i+j;
  
  std::cout << "A = \n" << A << "\n";
  Kaskade::DynamicMatrix<Dune::FieldMatrix<double,1,1>> U, V;
  std::vector<double> sing;
  std::tie(U,sing,V) = svd(A);
  
  std::cout << "U = \n" << U << "\n";
  std::cout << "V = \n" << V << "\n";
  std::cout << "singular values: "; std::copy(begin(sing),end(sing),std::ostream_iterator<double>(std::cout,", ")); std::cout << "\n";
  return 0;
}

#endif