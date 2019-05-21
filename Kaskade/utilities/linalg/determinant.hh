/*
 * determinant.hh
 *
 *  Created on: Sep 21, 2011
 *      Author: bzflubko
 */

#ifndef DETERMINANT_HH_
#define DETERMINANT_HH_

#include <dune/common/fmatrix.hh>
#include <boost/static_assert.hpp>
#include "utilities/get.hh"
#include "fem/fixdune.hh"

/// \cond internals
namespace DeterminantDetail {

  template <typename MatrixType>
  Dune::FieldMatrix<double,MatrixType::rows-1,MatrixType::cols-1> getMinor(MatrixType const& matrix, int rowIndex, int colIndex)
  {
    BOOST_STATIC_ASSERT( (MatrixType::rows > 1) && (MatrixType::cols > 1) );
    Dune::FieldMatrix<double,MatrixType::rows-1,MatrixType::cols-1> subMatrix(0);
    for(int ii=0, iii=0; ii<MatrixType::rows; ++ii)
    {
      if(ii!=rowIndex)
      {
        for(int kk=0, kkk=0; kk<MatrixType::cols; ++kk)
        {
          if(kk!=colIndex)
          {
            subMatrix[iii][kkk]=matrix[ii][kk];
            ++kkk;
          }
        }
        ++iii;
      }
    }

    return subMatrix;
  }

}
/// \endcond

namespace Kaskade {

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* DETERMINANT                                                                                                                     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**
 * \ingroup linalgbasic
 * \brief Class for computing determinants of matrices and (directional) derivatives thereof.
 * 
 * On construction, the Determinant class takes a reference to the matrix (or a functor returning the matrix) and computes
 * determinants and derivatives on demand. That is, if the referenced matrix is modified, the modification is reflected 
 * in the next determinant computation.
 */
template <int dim, class Source=Dune::FieldMatrix<double,dim,dim> >
class Determinant;

template <class Source>
class Determinant<2,Source>
{
public:
  typedef Source source_type;
  typedef Dune::FieldMatrix<double,2,2> MatrixType;

private:
  typedef Get<MatrixType,Source> GetMatrix;

  double composeResult(MatrixType const& m1, MatrixType const& m2) const { return m1[0][0]*m2[1][1] + m1[1][1]*m2[0][0] - (m1[0][1]*m2[1][0] + m1[1][0]*m2[0][1]); }

public:
  /**
   * \brief Constructor.
   * 
   * \param s reference to the source of the matrix of which the determinant is to be computed. This can either 
   *          be a Dune::FieldMatrix<double,n,n> itself or a functor returning such a matrix. Note that the source
   *          is held by reference and needs to exist as long as the Determinant object exists.
   */
  explicit Determinant(Source& s) : matrix(s){}

  double operator()() const { return d0(); }

  double d0() const { return matrix()[0][0]*matrix()[1][1] - matrix()[0][1]*matrix()[1][0]; }

  double d1(MatrixType const& dA) const { return composeResult(matrix(),dA); }

  double d2(MatrixType const& dA, MatrixType const& dB) const { return composeResult(dB,dA); }

  double d3(MatrixType const&, MatrixType const&, MatrixType const&) const { return 0; }

private:
  GetMatrix matrix;
  //MatrixType &matrix;
};

/// Determinant
template <class Source>
class Determinant<3,Source>
{
public:
  typedef Source source_type;
  typedef Dune::FieldMatrix<double,3,3> MatrixType;

private:
  typedef Get<MatrixType,Source> GetMatrix;

  template<int a1, int a2, int b1, int b2, int c1, int c2> double tp(MatrixType const& t) const
  {
    return t[a1][a2]*t[b1][b2]*t[c1][c2];
  }

  template<int a1, int a2, int b1, int b2, int c1, int c2> double circular(MatrixType const& t, MatrixType const& u, MatrixType const& v) const
  {
    return t[a1][a2]*u[b1][b2]*v[c1][c2]+v[a1][a2]*t[b1][b2]*u[c1][c2]+u[a1][a2]*v[b1][b2]*t[c1][c2];
  }

  template<int a1, int a2, int b1, int b2, int c1, int c2> double dtp(MatrixType const& t, MatrixType const& dv) const
  {
    return circular<a1,a2,b1,b2,c1,c2>(t,t,dv);
  }

  template<int a1, int a2, int b1, int b2, int c1, int c2> double ddtp(MatrixType const& t, MatrixType const& dv, MatrixType const& dw) const
  {
    return circular<a1,a2,b1,b2,c1,c2>(t,dv,dw)+circular<a1,a2,b1,b2,c1,c2>(t,dw,dv);
  }

  template<int a1, int a2, int b1, int b2, int c1, int c2> double dddtp(MatrixType const& du, MatrixType const& dv, MatrixType const& dw) const
  {
    return circular<a1,a2,b1,b2,c1,c2>(du,dv,dw)+circular<a1,a2,b1,b2,c1,c2>(du,dw,dv);
  }

public:
  /**
   * \brief Constructor.
   * 
   * \param s reference to the source of the matrix of which the determinant is to be computed. This can either 
   *          be a Dune::FieldMatrix<double,n,n> itself or a functor returning such a matrix. Note that the source
   *          is held by reference and needs to exist as long as the Determinant object exists.
   */
  explicit Determinant(Source &s): matrix(s){}

  double operator()() const { return d0(); }

  double d0() const
  {
    return tp<0,0,1,1,2,2>(matrix())+tp<0,1,1,2,2,0>(matrix())+tp<0,2,1,0,2,1>(matrix())
        - tp<2,0,1,1,0,2>(matrix())-tp<2,1,1,2,0,0>(matrix())-tp<2,2,1,0,0,1>(matrix());
  }

  double d1(MatrixType const& dv) const
  {
    return dtp<0,0,1,1,2,2>(matrix(),dv)+dtp<0,1,1,2,2,0>(matrix(),dv)+dtp<0,2,1,0,2,1>(matrix(),dv)
        -dtp<2,0,1,1,0,2>(matrix(),dv)-dtp<2,1,1,2,0,0>(matrix(),dv)-dtp<2,2,1,0,0,1>(matrix(),dv);
  }

  double d2(MatrixType const& dv, MatrixType const& dw) const
  {
    return ddtp<0,0,1,1,2,2>(matrix(),dv,dw)+ddtp<0,1,1,2,2,0>(matrix(),dv,dw)+ddtp<0,2,1,0,2,1>(matrix(),dv,dw)
        - ddtp<2,0,1,1,0,2>(matrix(),dv,dw)-ddtp<2,1,1,2,0,0>(matrix(),dv,dw)-ddtp<2,2,1,0,0,1>(matrix(),dv,dw);
  }

  double d3(MatrixType const& du, MatrixType const& dv, MatrixType const& dw) const
  {
    return dddtp<0,0,1,1,2,2>(du,dv,dw)+dddtp<0,1,1,2,2,0>(du,dv,dw)+dddtp<0,2,1,0,2,1>(du,dv,dw)
        - dddtp<2,0,1,1,0,2>(du,dv,dw)-dddtp<2,1,1,2,0,0>(du,dv,dw)-dddtp<2,2,1,0,0,1>(du,dv,dw);
  }

private:
  GetMatrix matrix;
//  MatrixType &matrix;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* ADJUGATE                                                                                                                        */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
template <int dim, class Source=Dune::FieldMatrix<double,dim,dim> > class Adjugate;

template <int dim, class Source=Dune::FieldMatrix<double,dim,dim> > class Cofactor;

template <class Source>
class Adjugate<2, Source>{
public:
  typedef Dune::FieldMatrix<double,2,2> MatrixType;
private:
  typedef Get<MatrixType,Source> GetMatrix;

  void composeResult(MatrixType const& m) const
  {
    result[0][0] = m[1][1]; result[0][1] = -m[1][0];
    result[1][0] = -m[0][1]; result[1][1] = m[0][0];
  }

public:
  explicit Adjugate(Source &s) : matrix(s), result(0){}

  MatrixType d0() const
  {
    composeResult(matrix());
    return result;
  }

  MatrixType d1(MatrixType const& dA) const
  {
    composeResult(dA);
    return result;
  }

  MatrixType d2(MatrixType const&, MatrixType const&) const { return 0.0*result; }

  MatrixType d3(MatrixType const&, MatrixType const&, MatrixType const&) { return 0.0*result; }

  //Cofactor<2,Source> transpose() const { return Cofactor<2,Source>(matrix); }

private:
  GetMatrix matrix;
//  MatrixType &matrix;
  mutable MatrixType result;
};

/**
 * Adjugate of a 3x3-matrix
 */
template <class Source>
class Adjugate<3,Source>
{
public:
  enum{ dim=3 };
  typedef Dune::FieldMatrix<double,dim,dim> MatrixType;
  typedef typename GetSubType<double,dim,MatrixType>::type SubMatrixType;
private:
  typedef Get<MatrixType,Source> GetMatrix;

public:
  explicit Adjugate(Source &s) :matrix(s){}

  MatrixType d0() const
  {
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
      {
        SubMatrixType tmp = DeterminantDetail::getMinor(matrix(),i,j);
        int sign = 1;
        if((i+j)%2) sign = -1;
        result[j][i] = sign*Determinant<dim-1,SubMatrixType>(tmp).d0();
      }
    return result;
  }

  MatrixType d1(MatrixType const& dv) const{
    for(int i=0; i<dim; ++i){
      for(int j=0; j<dim; ++j){
        int i1 = (i+1)%dim; int i2 = (i+2)%dim;
        int j1 = (j+1)%dim; int j2 = (j+2)%dim;
        result[j][i] = matrix()[i1][j1]*dv[i2][j2] + matrix()[i2][j2]*dv[i1][j1] -
            matrix()[i1][j2]*dv[i2][j1] - matrix()[i2][j1]*dv[i1][j2];
        //result[j][i] *= ((i+j)%2==0) ? 1. : -1.;
      }
    }
    return result;
  }

  MatrixType d2(MatrixType const& dv, MatrixType const& dw) const {
    result = MatrixType(0);
    for(int i=0; i<dim; ++i){
      for(int j=0; j<dim; ++j){
        int i1 = (i+1)%dim; int i2 = (i+2)%dim;
        int j1 = (j+1)%dim; int j2 = (j+2)%dim;

        result[j][i] += dv[i2][j2]*dw[i1][j1] + dv[i1][j1]*dw[i2][j2];
        result[j][i] -= dv[i2][j1]*dw[i1][j2] + dv[i1][j2]*dw[i2][j1];
        //result[j][i] *= ((i+j)%2==0) ? 1. : -1.;
      }
    }

    return result;
  }

  MatrixType d3(MatrixType const&, MatrixType const&, MatrixType const&) const
  {
    return 0.0*matrix();
  }

  //Cofactor<3,Source> transpose() const { return Cofactor<3,Source>(matrix); }

private:
  GetMatrix matrix;
//  MatrixType &matrix;
  mutable MatrixType result;
};


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Cofactor matrix                                                                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**
 * implemented as adjugate transposed
 */
template <int dimension, class Source>
class Cofactor{	
public:	
  enum{ dim=dimension };
  typedef Dune::FieldMatrix<double,dim,dim> MatrixType;
  typedef Dune::FieldMatrix<double,dim-1,dim-1> SubMatrixType;
private:
  typedef Get<MatrixType,Source> GetMatrix;

public:
  explicit Cofactor(Source& s) : adj(s){}

  MatrixType const d0() const { return transpose( adj.d0() ); }

  MatrixType const d1(MatrixType const& dv) const { return transpose( adj.d1(dv) ); }

  MatrixType const d2(MatrixType const& dv, MatrixType const& dw) const
  { return transpose( adj.d2(dv,dw) ); }
  
  MatrixType const d3(MatrixType const& dv, MatrixType const& dw, MatrixType const& dx) const
  { return transpose( adj.d3(dv,dw,dx) ); }

  //Adjugate<dim,Source> transpose() const { return adj; }

private: 
  Adjugate<dim,Source> const adj;
};

} // end of namespace Kaskade

#endif /* DETERMINANT_HH_ */
