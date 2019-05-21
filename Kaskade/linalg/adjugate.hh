/*
 * adjugate.hh
 *
 *  Created on: June 23, 2014
 *      Author: bzflubko
 */
#ifndef ADJUGATE_HH
#define ADJUGATE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "linalg/cofactor.hh"

namespace Kaskade
{
  template<class,size_t,bool> class WrappedMatrix;

  template <int dim, class Source = WrappedMatrix<double,dim,false> > class Adjugate;

  template<class Source_>
  class Adjugate<2, Source_>
  {
  public:
    static constexpr int dim = 2;
    typedef Source_ Source;
    typedef typename Source::Argument Argument;
    typedef Dune::FieldMatrix<double,dim,dim> ReturnType;

    explicit Adjugate(Source f_) : f(f_) {}

    Adjugate(Adjugate const&) = default;
    Adjugate& operator=(Adjugate const&) = default;

    ReturnType d0() const { return composeResult(f.d0()); }
    ReturnType d1(Argument const& dF1) const { return composeResult(f.d1(dF1)); }
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return composeResult(f.d2(dF1, dF2)); }
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) { return composeResult(f.d3(dF1, dF2, dF3)); }

    Dune::FieldVector<double,dim> operator*(Dune::FieldVector<double,dim> const& v) const { return d0()*v; }

  private:
    ReturnType composeResult(ReturnType const& m) const
    {
      ReturnType result(0);
      result[0][0] = m[1][1];
      result[0][1] = -m[0][1];
      result[1][0] = -m[1][0];
      result[1][1] = m[0][0];
      return result;
    }

    Source f;
  };


  template <class Source_>
  class Adjugate<3,Source_>
  {
    typedef Cofactor<3,Source_> Cof;
  public:
    static constexpr int dim = 3;
    typedef Source_ Source;
    typedef Dune::FieldMatrix<double,dim,dim> ReturnType;
    typedef typename Source::Argument Argument;

    explicit Adjugate(Source const& f_) : f(f_) {}

    Adjugate(Adjugate const&) = default;
    Adjugate& operator=(Adjugate const&) = default;

    ReturnType d0() const { return fillMatrix<0>(); }
    ReturnType d1(Argument const& dF) const { return fillMatrix<1>(dF); }
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return fillMatrix<2>(dF1,dF2); }
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return fillMatrix<3>(dF1,dF2,dF3); }

    Dune::FieldVector<double,dim> operator*(Dune::FieldVector<double,dim> const& v) const { return d0()*v; }

  private:
    template <size_t d, typename... Args>
    ReturnType fillMatrix(const Args&... args) const
    {
      ReturnType A(0);

      for(size_t i=0; i<dim; ++i)
        for(size_t j=0; j<dim; ++j)
          A[i][j] = insert(i,j,std::integral_constant<size_t,d>(),args...);

      return A;
    }

    double insert(size_t i, size_t j, std::integral_constant<size_t,0>) const { return Cof(f,j,i).d0(); }
    double insert(size_t i, size_t j, std::integral_constant<size_t,1>, Argument const& dF) const { return Cof(f,j,i).d1(dF); }
    double insert(size_t i, size_t j, std::integral_constant<size_t,2>, Argument const& dF1, Argument const& dF2) const { return Cof(f,j,i).d2(dF1,dF2); }
    double insert(size_t i, size_t j, std::integral_constant<size_t,3>, Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return Cof(f,j,i).d3(dF1,dF2,dF3); }

    Source f;
  };
}

#endif // ADJUGATE_HH
