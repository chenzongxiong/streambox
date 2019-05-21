#include <cassert>

// depricated: original include which is not valid on Mac platform
// #include "/htcsoft/amd/sles9/acml/3.6.0/gnu64/include/acml.h"

// Access on acml include file in dependence of the platform HTC-Linux or Mac OS.
// ZIBHTC is to set in the Makefile.Local (copy of Makefile-htc.Local) on HTC platform.
#ifdef ZIBHTC
#include "acml.h"
#endif
#if defined(Darwin) || defined(Cygwin) || defined(Atlas)
#include "acml_to_stdlapack.hh"
#endif

#include "simpleLAPmatrix.hh"

namespace Kaskade
{

void MatrixKT_A_K(SLAPMatrix<double> A, SLAPMatrix<double>& K, SLAPMatrix<double>& out)
{
  assert(K.rows()==A.cols());
  assert(K.cols()==out.cols());
  assert(K.cols()==out.rows());
  SLAPMatrix<double> AK(A.rows(),K.cols());
  dgemm('n','n',AK.rows(),AK.cols(), K.rows(), 1.0, A.ptrToData(), A.LDA(), K.ptrToData(), K.LDA(), 0.0, AK.ptrToData(), AK.LDA());
  dgemm('t','n',out.rows(),out.cols(), AK.rows(), 1.0, K.ptrToData(), K.LDA(), AK.ptrToData(), AK.LDA(), 0.0, out.ptrToData(), out.LDA());
}

void MatrixMultiplication(SLAPMatrix<double> A, SLAPMatrix<double>& B, SLAPMatrix<double>& AB)
{
  assert(B.rows()==A.cols());
  assert(B.cols()==AB.cols());
  assert(A.rows()==AB.rows());
  dgemm('n','n',AB.rows(),AB.cols(), B.rows(), 1.0, A.ptrToData(), A.LDA(), B.ptrToData(), B.LDA(), 0.0, AB.ptrToData(), AB.LDA());
}

void MatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out)
{
  if(in.size() !=A.cols())
  {
    std::cout << "Size Mismatch:" << in.size() << " " << A.cols() << std::endl;
  }
  assert(in.size()==A.cols());
  out.resize(A.rows());
  dgemv('n',A.rows(),A.cols(), 1.0, A.ptrToData(), A.LDA(), &in[0], 1, 0.0, &out[0], 1);
}

void TransposedMatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out)
{
  assert(in.size()==A.rows());
  out.resize(A.cols());
  dgemv('t',A.rows(),A.cols(), 1.0, A.ptrToData(), A.LDA(), &in[0], 1, 0.0, &out[0], 1);
}

void LeastSquares(SLAPMatrix<double> a, std::vector<double> const& b, std::vector<double>& x)
{
  x.reserve(std::max(a.cols(),a.rows()));
  x=b;
  int info;
  dgels('n',a.rows(),a.cols(),1,a.ptrToData(),a.LDA(),&x[0],a.rows(),&info);
  x.resize(a.cols());
}

void SymmetricEigenvalues(SLAPMatrix<double> a, std::vector<double>& eig)
{
  eig.resize(a.rows());
  int info;
  dsyevd('N','U',a.rows(),a.ptrToData(),a.LDA(),&eig[0],&info);
}

//---------------------------------------------------------------------

template <>
void pseudoinverse(SLAPMatrix<double> A, SLAPMatrix<double>& X)
{
  assert(A.rows()==X.cols() && A.cols()==X.rows());

  int n = A.rows();
  SLAPMatrix<double> B(n,n,std::max(A.rows(),A.cols()));
  for (int i=0; i<n; ++i) B(i,i) = 1;

  //int jpvt[A.cols()]; // probably this causes a valgrind error/warning message: use of not inititalized value
  std::vector<int> jpvt(A.cols(),0);
  int rank, info;
  dgelsy(A.rows(),A.cols(),B.cols(),A.ptrToData(),A.LDA(),B.ptrToData(),B.LDA(),
         jpvt.data(),1e-10,&rank,&info);

  if(rank < std::min(A.rows(),A.cols()))
  {
    std::cout << "Warning: rank deficiency detected:" << std::endl;
    std::cout << "A is " << A.rows() << "x" << A.cols() << " but Rank(A) = " << rank << std::endl;
    std::cout << "Info:" << info << std::endl;
    for(int i=0; i<A.rows(); ++i)
    {
      for(int j=0; j<A.cols(); ++j)
        std::cout << A(i,j) << " ";
      std::cout << std::endl;
    }
  }

  for (int j=0; j<X.cols(); ++j)
    for (int i=0; i<X.rows(); ++i)
      X(i,j) = B(i,j);
}

void basisOfKernel(SLAPMatrix<double> A, SLAPMatrix<double>& Basis) 
{
  assert(A.cols()==Basis.rows() && (A.cols()-A.rows())==Basis.cols());

  int info;
  std::vector<double> sing(std::max(A.cols(),A.rows()));
  SLAPMatrix<double> B(A.cols(),A.cols());

  dgesvd('N', 'A', A.rows(), A.cols(), A.ptrToData(), A.LDA(), &sing[0], (double*)0, 1, B.ptrToData(), B.LDA(), &info);
  
  for (int j=0; j<Basis.cols(); ++j)
    for (int i=0; i<Basis.rows(); ++i)
      Basis(i,j) = B(j+A.rows(),i);
}


}  // namespace Kaskade

// explicit instantiation
