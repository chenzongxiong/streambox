#ifndef SIMPLELAPMATRIX_HH
#define SIMPLELAPMATRIX_HH

#include <vector>
#include <dune/grid/config.h>
#include <linalg/linearsystem.hh>

namespace Kaskade
{

/**
 * @file
 * @brief Simple interface to acml (AMD Core Math Library), to be augmented
 * @author Anton Schiela
 *
 * Functions according to the needs of the primary author, to be augmented.
 * Focus was on a quick and simple implementation, usage for auxilliary tasks.
 * Performance could certainly be improved.
 *
 * Functions are not templated to avoid inclusion of the long acml.h - header file into all the code.
 * Anyway, the calls to LAPACK routines for various types would require a total template specialization
 */

/// Very simple vector wrapper for the use of the optimization tool uncmin, which is not capable of using std::vector :(
class UNCMINVector
{
public:
  UNCMINVector(const std::vector<double>& v) : data(v) {};
  UNCMINVector(int n, double v=0.0) : data(n,v) {};
  UNCMINVector() {};
  UNCMINVector(double* start, double* stop) { for(double* s=start; s<stop; ++s) data.push_back(*s); };
  UNCMINVector& operator=(const UNCMINVector &A) { data=A.data; return *this; }
  UNCMINVector& operator=(const double &s) { for(int i=0; i<data.size(); ++i) data[i]=s; return *this; }
  UNCMINVector& operator*=(const double &s) { for(int i=0; i<data.size(); ++i) data[i]*=s; return *this; }
  double& operator()(int i) { return data[i-1]; }
  int size() { return data.size(); }
  void newsize(int n) {data.resize(n); }

  std::vector<double> data;
};

/** 
 * \ingroup linalgbasic
 * \brief Very simple dense matrix class for interfacing with LAPACK routines and the optimization tool uncmin.
 * 
 * Before using this, consider using Kaskade::DynamicMatrix, that implements the Dune dense matrix interface.
 * 
 * Can be constructed from objects that model the Dune Matrix concept, in particular we need
 *  the members
 *
 * - M()     : number of columns
 * - N()     : number of rows
 * - (i,j)   : element of i-th row, j-th column
 *
 * The use of SLAPMatrix is only sensible, if Num is a numeric type, such as float or double
 * 
 * \tparam Num type of matrix entries
 * \tparam offset start index (0 is C-style, 1 is Fortran-style indexing)
 */
template<class Num, int offset=0>
class SLAPMatrix
{
public:
template<class Mat>
  SLAPMatrix(Mat const& mat)
  {
    col = mat.M();
    row = mat.N();
    lda = mat.N();
    data.reserve(col*row);
    // Insert elements in COLUMN-MAJOR ordering
    for(int i=0; i<col; ++i)
      for(int j=0; j<row; ++j)
        data.push_back(mat[j][i]);
  }

  SLAPMatrix(SLAPMatrix<Num,offset> const& mat)
  {
    col=mat.cols();
    row=mat.rows();
    lda=mat.LDA();
    data=mat.data;
  }

  /**
   * \brief Construction by specifying the size of the matrix directly.
   * 
   * The entries are initialized by 0.
   * 
   * \param row_ number of rows
   * \param col_ number of cols
   * \param lda_ leading dimension of the matrix. Values below row_ are implicitly corrected to row_.
   */
  SLAPMatrix(int row_, int col_, int lda_ = -1):
    data(std::max(lda_,row_)*col_,0), col(col_), row(row_), lda(std::max(lda_,row_)) 
  {
    assert(row>=0 && col>=0);
  }
    

  SLAPMatrix(MatrixAsTriplet<Num> const& mat)
  {
    col = mat.ncols();
    row = mat.nrows();
    lda = mat.nrows();
    data.resize(col*row);
    // This should be COLUMN-MAJOR ordering
    for(int i=0; i<mat.data.size(); ++i)
      (*this)(mat.ridx[i],mat.cidx[i]) = mat.data[i];
  }

  void print() const
  {
    std::cout << "[" << std::endl;
    for(int i=0; i<rows(); ++i)
    {
      for(int j=0; j<cols();++j)
        std::cout << (*this)(i,j) << ",";
      std::cout << std::endl;
    }
    std::cout << "]" << std::endl;
  }


  /// Columns
  int cols() const {return col;}
  /// Rows
  int rows() const {return row;}
  
  int size() {return col*row;}
  
  /// Increment for row iteration
  int LDA() const {return lda;}

template<class Mat>
  void toMatrix(Mat& M)
  {
    for(int i=0; i<rows(); ++i)
      for(int k=0; k<cols(); ++k)
        M[i][k] = (*this)(i,k);
  }

template<class RT>
  void toTriplet(MatrixAsTriplet<RT>& M)
  {
    for(int i=0; i<rows(); ++i)
      for(int k=0; k<cols(); ++k)
      {
        M.cidx.push_back(i);
        M.ridx.push_back(k);
        M.data.push_back((*this)(k,i));
      }
  }

  /// Access to the raw data. This is needed by LAPACK routines
  Num* ptrToData() {return &data[0];}

  /**
   * \brief Access to elements.
   */
  Num const& operator()(int r, int c) const { return data[lda*(c-offset)+(r-offset)]; }
  Num&       operator()(int r, int c)       { return data[lda*(c-offset)+(r-offset)]; }
  
  /// Interface to Optimizer
  /// Columns
  int num_columns() const {return col;}
  /// Rows
  int num_rows() const {return row;}

  SLAPMatrix<Num,offset> operator*=(Num alpha)
  {
    for(int i=0; i<rows(); ++i)
      for(int j=0; j<cols(); ++j)
        (*this)(i,j) *= alpha;
    return *this;
  }

  SLAPMatrix<Num,offset> scaleCols(std::vector<Num> scaling)
  {
    for(int i=0; i<rows(); ++i)
      for(int j=0; j<cols(); ++j)
        (*this)(i,j) *= scaling[j];
    return *this;
  }

  SLAPMatrix<Num,offset> scaleRows(std::vector<Num> scaling)
  {
    for(int i=0; i<rows(); ++i)
      for(int j=0; j<cols(); ++j)
        (*this)(i,j) *= scaling[i];
    return *this;
  }

private:
  std::vector<Num> data;
  int col,row,lda;
};

/// Solve linear least squares problem, given by a*x=b
/** Preferrably used like the following:
 *
 * LeastSquares(SLAPMatrix<double>(A),b,x)
 *
 * Here A models the Dune Matrix concept.
 */
void LeastSquares(SLAPMatrix<double> a, std::vector<double> const& b, std::vector<double>& x);

void SymmetricEigenvalues(SLAPMatrix<double> a, std::vector<double>& eig);

/**
 * Computes the pseudoinverse.
 *
 * \arg A The input matrix. Will be overwritten with unspecified data.
 * \arg Ainv the pseudoinverse. Needs to have the shape of A transposed.
 */
template <class Scalar>
void pseudoinverse(SLAPMatrix<Scalar> A, SLAPMatrix<Scalar>& Ainv);

void basisOfKernel(SLAPMatrix<double> A, SLAPMatrix<double>& Basis);


/** 
 * \ingroup linalgbasic
 * \brief Matrix multiplication
 *
 * out = A*in
 * out is resized properly
 */
void MatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out);

void TransposedMatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out);

void MatrixKT_A_K(SLAPMatrix<double> A, SLAPMatrix<double>& K, SLAPMatrix<double>& out);

void MatrixMultiplication(SLAPMatrix<double> A, SLAPMatrix<double>& B, SLAPMatrix<double>& AB);

}  // namespace Kaskade
#endif
