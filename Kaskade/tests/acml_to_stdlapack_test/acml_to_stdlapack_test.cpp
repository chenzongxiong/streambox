/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Access on acml include file in dependence of the platform HTC-Linux or Mac OS.
// ZIBHTC is to set in the Makefile.Local (copy of Makefile-htc.Local) on HTC platform.
#include <iostream>
#include <algorithm>  // std::max
#include <vector>
#include <cassert>

#ifdef ZIBHTC
#include "acml.h"
#endif

#if defined(Darwin) || defined(Cygwin) || defined(Atlas)
#include "linalg/acml_to_stdlapack.hh"
#endif

template <class T> class FortranMatrix
{
  typedef FortranMatrix<T> FMatrix;
  
  public:
  FortranMatrix(int rows, int cols ) { storage.resize(rows*cols); _rows=rows; _cols=cols; }
  FortranMatrix()  {}
  ~FortranMatrix()  {}
  
  int const rows() { return _rows; }
  int const cols() { return _cols; }

  FMatrix& operator = ( std::initializer_list<std::initializer_list<T> > il )
  {
    int i=0;
    for ( auto ita=il.begin();ita!=il.end();ita++ )
    {
      int j=0;
      for ( auto iti=ita->begin();iti!=ita->end();iti++ )
      {
        storage[i+_rows*j] = *iti;
        j++;
      };
      i++;
    };
    return *this;
  }
  
  FMatrix& operator = ( FMatrix& mat )
  {
    if ( _rows*_cols != mat.rows()*mat.cols() )
      storage.resize(mat.rows()*mat.cols());
    _rows=mat.rows();
    _cols=mat.cols();
    for (int i=0;i<_rows*_cols;i++)
      storage[i] = mat[i];
      return *this;
  }
  
  void set( int const i, int const j, T const value )
  {
    storage[i+_rows*j] = value;
  }
  
  T get( int const i, int const j )
  {
    return storage[i+_rows*j];
  }
  
  T& operator [] ( int const row )
  {
    return (T&) storage[row];
  }
    
  void print()
  {
    for (int i=0;i<_rows;i++)
    {
      for (int j=0;j<_cols;j++)
        std::cout << storage[i+_rows*j] << " ";
      std::cout << std::endl;
    }
  };
  
  private:
  int _rows, _cols;
  std::vector<T> storage;
};

void print(std::vector<double> const& a)
{
  for(int i=0;i<a.size();i++)
    std::cout << a[i] << " ";
  std::cout << std::endl;
}

template<class T>
T deviation(int const n, std::vector<T> const& a, std::vector<T> const& b,
                 T* maximum)
{
  T r=0, d;
  assert( (n <= a.size()) &&  (n <= b.size()) );
  for (int j=0;j<n;j++)
  {
    d= a[j]-b[j];
     r += d*d;
  };
  r=sqrt(r);
  *maximum = std::max(*maximum,r);
  return r;
}

template<class T>
T deviation(std::vector<T> const& a, std::vector<T> const& b, T* maximum)
{
  assert( a.size() == b.size() );
  return deviation(a.size(),a,b,maximum);
}

template<class T>
T deviation(FortranMatrix<T>& a, FortranMatrix<T>& b, T* maximum)
{
  T r=0, d;
  assert( a.rows()==b.rows() );
  int const m=a.rows(), n=std::min(a.cols(),b.cols());
  for (int i=0;i<m;i++)
    for (int j=0;i<n;i++)
    {
      int k=i+j*m;
      d= a[k]-b[k];
      r += d*d;
//       if ( d != 0 ) 
//         std::cout << "deviation: " << i << " " << a[k] << " " << b[k] << std::endl;
    };
  r=sqrt(r);
  *maximum = std::max(*maximum,r);
  return r;
}

int main()
{
  int constexpr lda=5, m1=3, m2=5, n=3, k=3, nrhs=1, inc=1;
  double maximum=0, tol=1.0e-6;
  float maximums=0;
  int info;
  // void dgemm(char transa, char transb, int m, int n, int k, double alpha, double *a,
  //           int lda, double *b, int ldb, double beta, double *c, int ldc)
  FortranMatrix<double> A(lda,n),B(lda,n),C(lda,n),CExpected(lda,n);
  A = { {1,2,3},{4,5,6},{7,8,10},{10,11,14},{13,16,25} };
  B = { {9,6,3},{8,5,2},{7,4,1} };
  C = { {-6,-5,-4},{-9,-8,-7},{-3,-2,-1} };
  double alpha = 3, beta=4;
  // C=alpha*A*B'+beta*C
  dgemm('n','t',m2,n,k,alpha,&A[0],lda,&B[0],lda,beta,&C[0],lda);
  CExpected = { {66,52,38},{216,175,134},{411,340,269},{598,497,396},{880,722,564} };
  std::cout << "dgemm test: deviation from expected result is " << 
              deviation(C,CExpected,&maximum) << std::endl;

  // void dgemv(char transa, int m, int n, double alpha, double *a, int lda, double *x,
  //            int incx, double beta, double *y, int incy) 
  std::vector<double> x(m2), x2(n), y(n), yExpected(n);
  x={1,4,9,16,25};
  y={9,4,1};
  // y=alpha*A'*x + beta*y
  dgemv('t', m2,n,alpha,&A[0],lda,&x[0],inc,beta,&y[0],inc);
  yExpected = {1731,2026,2902};
  std::cout << "dgemv test: deviation from expected result is " <<
               deviation(y,yExpected,&maximum) << std::endl;

  // void dgels(char trans, int mA, int nA, int nrhsA, double *a, int ldaA, double *b,
  //            int ldB, int *infoA) 
  FortranMatrix<double> AQR(lda,n);
  AQR=A;
  y=x;
  info=0;
  // y=lsqlin(AQR,y)
  dgels('n',m2,n,nrhs,&AQR[0],lda,&y[0],m2,&info);
  std::cout << "info returned from dgels is " << info << std::endl;
  yExpected = {3.894523326572028,-4.560851926977715,1.897565922920899};
  std::cout << "dgels test: deviation from expected result is " << 
                deviation(n,y,yExpected,&maximum) << std::endl;
  
  // void dsyevd(char jobz, char uplo, int nA, double *a, int ldaA, double *w, int *infoA) 
  AQR={ {1,2,3},{0,3,5},{0,0,7}};
  y.assign(n,0);
  info=0;
  // w=eig(AQR)
  dsyevd('V','U',n,&AQR[0],lda,&y[0],&info);
  yExpected = { -3.888035059e-01,-2.215262886e-01,1.161032979e+01 };
  std::cout << "dsyevd test: deviation from expected result is " << 
               deviation(y,yExpected,&maximum) << std::endl;
  
  // void dgelsy(int mA, int nA, int nrhsA, double *a, int ldaA, double *b, int ldb,
  //             int *jpvtA, double rcond, int *rankA, int *infoA) 
  AQR=A;
  y=x;
  info=0;
  double rcond=1.0e-16;
  int rank;
  std::vector<int> pivot(n,0);
  // y=lsqlin(AQR,y)
  dgelsy(m2,n,nrhs,&AQR[0],lda,&y[0],m2,&pivot[0],rcond,&rank,&info);
  std::cout << "info returned from dgelsy is " << info << std::endl;
  std::cout << "rank returned from dgelsy is " << rank << std::endl;
  yExpected = {3.894523326572028,-4.560851926977715,1.897565922920899};
  std::cout << "dgelsy test: deviation from expected result is " << 
              deviation(n,y,yExpected,&maximum) << std::endl;

  // void dgesv(int nA, int nrhsA, double *a, int ldaA, int *pivA, double *b,
  //            int ldbA, int *infoA)
  FortranMatrix<double> ALU(lda,n);
  ALU=A;
  x2={1,4,9};
  y=x2;
  info=0;
  pivot.assign(n,0);
  // z=linsolve(AQR,x)
  dgesv(n,nrhs,&ALU[0],lda,&pivot[0],&y[0],n,&info);
  std::cout << "info returned from dgesv is " << info << std::endl;
  yExpected = { 3,-4,2 };
  std::cout << "dgesv test: deviation from expected result is " << 
               deviation(y,yExpected,&maximum) << std::endl;

  // void dgesvd(char jobu, char jobvt, int mA, int nA, double *a, int ldaA,
  //             double *sing, double *u, int lduA, double *vt, int ldvtA, int *infoA)
  FortranMatrix<double> SVD(lda,n), U(lda,m2), vt(n,n);
  std::vector<double> sing(n);
  SVD=A;
  info=0;
  // [U,sing,vt] = svd(SVD)
  dgesvd('A','A',m2,n,&SVD[0],lda,&sing[0],&U[0],lda,&vt[0],n,&info);
  std::cout << "info returned from dgesvd is " << info << std::endl;
  FortranMatrix<double> UExpected(lda,m2), vtExpected(n,n);
  std::vector<double> singExpected(n);
  UExpected = {{-0.087677886267887, 0.197967076010598,-0.716565735059678,-0.279964026250644,-0.601061824124043},
               {-0.208290746872814,-0.229306690655220,-0.603301379462155, 0.628079085280756, 0.381546215496717},
               {-0.346501119119734,-0.416983184102440,-0.169047956528042,-0.696230118060224, 0.439031217254653},
               {-0.484711491366654,-0.604659677549660, 0.265205466406070, 0.199898809838595,-0.537724809634173},
               {-0.770661597338707, 0.607238978457326, 0.153785050832248, 0.049405416397172, 0.106069733668949}
              };
  singExpected = {41.992070375666032,2.713341407610269,0.551184153307933};
  vtExpected = {{-0.433702463911394,-0.515602832358961,-0.738955947307148},
                {-0.659936131095202,-0.376619605972657, 0.650109202574534},
                { 0.613503443896024,-0.769617691883716, 0.176924087301021}
               };
  std::cout << "dgesvd test: deviation from expected result for sing is " << 
               deviation(sing,singExpected,&maximum) << std::endl;
  std::cout << "dgesvd test: deviation from expected result for U is " << 
               deviation(U,UExpected,&maximum) << std::endl;
  std::cout << "dgesvd test: deviation from expected result for vt is " << 
               deviation(vt,vtExpected,&maximum) << std::endl;

  // void dstev(char jobz, int nA, double *diag, double *offd, double *z,
  //            int ldzA, int *infoA)
  std::vector<double> diag(n), offd(n-1);
  FortranMatrix<double> evec(n,n);
  // A=[ 1 2 0; 2 3 5; 0 5 7];
  diag = { 1,3,7 };
  offd = { 2,5 };
  info=0;
  dstev('V',n,&diag[0],&offd[0],&evec[0],n,&info);
  std::cout << "info returned from dgesvd is " << info << std::endl;
  std::vector<double> diagExpected(n);
  FortranMatrix<double> evecExpected(n,n);
  evecExpected  = {{ 0.564727113803817,-0.816496580927726, 0.120069231146632},
                   {-0.711781288862587,-0.408248290463863, 0.571577405220368},
                   { 0.417672938745049, 0.408248290463863, 0.811715867513632},
                  };
  diagExpected = {-1.520797289396148,1.999999999999999,10.520797289396150};
  std::cout << "dstev test: deviation from expected result for diag (eigenvalues) is " << 
               deviation(diag,diagExpected,&maximum) << std::endl;
//   std::cout << "eigenvectors computed:\n";
//   evec.print();
//   std::cout << "eigenvectors expected:\n";
//   evecExpected.print();
  for (int j=0;j<n;j++)
    if ( evec[j*n]*evecExpected[j*n]<0 ) // if signs of first component of eigenvector differs ...
      for (int i=j*n;i<(j+1)*n;i++ )
        evecExpected[i]=-evecExpected[i];   // ... multiply the expected eigenvector by -1
  std::cout << "dstev test: deviation from expected result for evec is " << 
               deviation(evec,evecExpected,&maximum) << std::endl;
  
  // void dgetrf(int mA, int nA, double *a, int ldaA, int *ipiv, int *infoA)
  FortranMatrix<double> AInv(lda,n);
  AInv=A;
  info=0;
  pivot.assign(n,0);
// LU Zerlegung von AInv
  dgetrf(m1,n,&AInv[0],lda,&pivot[0],&info);
  std::cout << "info returned from dgetrf is " << info << std::endl;

  // void dgetri(3, double *a, int ldaA, int *ipiv, int *infoA)
  // Ainv=inv(Inv)
  dgetri(n,&AInv[0],lda,&pivot[0],&info);
  std::cout << "info returned from dgetri is " << info << std::endl;
  FortranMatrix<double> AinvExpected(lda,n);
  AinvExpected = { {-0.666666667,-1.333333333,1.0000},
                   {-0.666666667,3.666666667,-2.0000},
                   {1.0000,-2.0000,1.0000} };
  std::cout << "dgetrf/dgetri test: deviation from expected result is " << 
              deviation(AInv,AinvExpected,&maximum) << std::endl;

  // void sgetrf(int mA, int nA, float *a, int ldaA, int *ipiv, int *infoA)
  FortranMatrix<float> AInvS(lda,n);
  AInvS={ {1,2,3},{4,5,6},{7,8,10} };
  info=0;
  pivot.assign(n,0);
// LU Zerlegung von AInvS
  sgetrf(m1,n,&AInvS[0],lda,&pivot[0],&info);
  std::cout << "info returned from sgetrf is " << info << std::endl;

  // void sgetri(int nA, float *a, int ldaA, int *ipiv, int *infoA)
  // AinvS=inv(InvS)
  sgetri(n,&AInvS[0],lda,&pivot[0],&info);
  std::cout << "info returned from sgetri is " << info << std::endl;
  FortranMatrix<float> AinvExpectedS(lda,n);
  AinvExpectedS = { {-0.666667,-1.333333,1.0000},
                    {-0.666667,3.666667,-2.0000},
                    {1.0000,-2.0000,1.0000} };
  std::cout << "sgetrf/sgetri test: deviation from expected result is " <<
               deviation(AInvS,AinvExpectedS,&maximums) << std::endl;

  // void sgemv(char transa, int m, int n, double alpha, double *a, int lda, double *x,
  //            int incx, double beta, double *y, int incy) 
  FortranMatrix<float> As(lda,n);
  As = { {1,2,3},{4,5,6},{7,8,10},{10,11,14},{13,16,25} };
  float alphas = 3, betas=4;
  std::vector<float> xs(m2), x2s(n), ys(n), yExpecteds(n);
  xs={1,4,9,16,25};
  ys={9,4,1};
  // ys=alphas*A'*xs + betas*ys
  sgemv('t', m2,n,alphas,&As[0],lda,&xs[0],inc,betas,&ys[0],inc);
  yExpecteds = {1731,2026,2902};
  std::cout << "sgemv test: deviation from expected result is " <<
               deviation(ys,yExpecteds,&maximums) << std::endl;

// Final message
  maximum = std::max(maximum,(double) maximums);
  std::cout << "Maximum deviation is " << maximum << std::endl;
  if (maximum<tol)
    return 0;
  else
    return 1;
}
