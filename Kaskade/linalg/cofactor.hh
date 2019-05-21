/*
 * cofactor.hh
 *
 *  Created on: June 23, 2014
 *      Author: bzflubko
 */
#ifndef COFACTOR_HH
#define COFACTOR_HH

#include "linalg/wrappedMatrix.hh"
#include "linalg/determinant.hh"

namespace Kaskade
{
  template <int,class> class Determinant;
  template <class,size_t,bool> class WrappedMatrix;

  template<typename Matrix>
  Dune::FieldMatrix<double, Matrix::rows - 1, Matrix::cols - 1> getMinor(Matrix const& A, size_t rowIndex, size_t colIndex)
  {
    static_assert( (Matrix::rows > 1) && (Matrix::cols > 1), "" );
    Dune::FieldMatrix<double, Matrix::rows - 1, Matrix::cols - 1> B(0);
    for(size_t i=0; i<Matrix::rows; ++i)
    {
      if(i<rowIndex)
      {
        for(size_t j=0; j<Matrix::cols; ++j)
        {
          if(j<colIndex) B[i][j] = A[i][j];
          if(j>colIndex) B[i][j-1] = A[i][j];
        }
      }
      if(i>rowIndex)
      {
        for(size_t j=0; j<Matrix::cols; ++j)
        {
          if(j<colIndex) B[i-1][j] = A[i][j];
          if(j>colIndex) B[i-1][j-1] = A[i][j];
        }
      }
    }

    return B;
  }

  template <int dim, class Source = WrappedMatrix<double,dim,false> >
  class Minor
  {
  public:
    typedef Dune::FieldMatrix<double,dim-1,dim-1> ReturnType;
    typedef typename Source::Argument Argument;

    explicit Minor(Source s, size_t row_, size_t col_) : f(s), row(row_), col(col_) {}

    Minor(Minor const&) = default;
    Minor& operator=(Minor const&) = default;

    ReturnType d0() const { return getMinor(f.d0(),row,col); }
    ReturnType d1(Argument const& dF1) const { return getMinor(f.d1(dF1),row,col); }
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return getMinor(f.d2(dF1,dF2),row,col); }
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return getMinor(f.d3(dF1,dF2,dF3),row,col); }

  private:
    Source f;
    size_t row = 0, col = 0;
  };

  template <int dim, class Source = WrappedMatrix<double,dim,false> >
  class Cofactor
  {
  public:
    typedef double ReturnType;
    typedef typename Source::Argument Argument;

    Cofactor(Source f_, size_t row, size_t col) : f(f_), det(Minor<dim,Source>(f,row,col)), sign( (row+col)%2==0 ? 1. : -1. ) {}

    Cofactor(Cofactor const&) = default;
    Cofactor& operator=(Cofactor const&) = default;

    ReturnType d0() const { return sign*det.d0(); }
    ReturnType d1(Argument const& dF1) const { return sign*det.d1(dF1); }
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return sign*det.d2(dF1,dF2); }
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return sign*det.d3(dF1,dF2,dF3); }

  private:
    Source f;
    Determinant<dim-1,Minor<dim,Source> > det;
    double sign = 1;
  };

}

#endif // COFACTOR_HH
