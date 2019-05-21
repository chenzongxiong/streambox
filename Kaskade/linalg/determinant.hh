/*
 * determinant.hh
 *
 *  Created on: June 23, 2014
 *      Author: bzflubko
 */
#ifndef DETERMINANT_HH
#define DETERMINANT_HH

#include "linalg/adjugate.hh"

namespace Kaskade
{
  template <class,size_t,bool> class WrappedMatrix;

  template<int dim, class Source = WrappedMatrix<double,dim,false> > class Determinant;


  template<class Source_>
  class Determinant<2, Source_>
  {
  public:
    static constexpr int dim = 2;
    typedef Source_ Source;
    typedef double ReturnType;
    typedef typename Source::Argument Argument;

    explicit Determinant(Source f_) : f(f_) {}

    Determinant(Determinant const&) = default;
    Determinant& operator=(Determinant const&) = default;

    double d0() const
    {
      auto const& m = f.d0();
      return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

    double d1(Argument const& dF1) const { return composeResult(f.d0(), f.d1(dF1)); }

    double d2(Argument const& dF1, Argument const& dF2) const { return composeResult(f.d0(), f.d2(dF1, dF2)) + composeResult(f.d1(dF2), f.d1(dF1)); }

    double d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const
    {
      return composeResult(f.d0(), f.d3(dF1,dF2,dF3)) + composeResult(f.d1(dF3), f.d2(dF1,dF2))
           + composeResult(f.d2(dF2,dF3),f.d1(dF1)) + composeResult(f.d1(dF2), f.d2(dF1,dF3));
    }

  private:
    double composeResult(Dune::FieldMatrix<double,2,2> const& m1, Dune::FieldMatrix<double,2,2> const& m2) const
    {
      return m1[0][0] * m2[1][1] + m1[1][1] * m2[0][0] - (m1[0][1] * m2[1][0] + m1[1][0] * m2[0][1]);
    }

    Source f;
  };

  template <class Source_>
  class Determinant<3,Source_>
  {
  public:
    static constexpr int dim = 3;
    typedef Source_ Source;
    typedef double ReturnType;
    typedef typename Source::Argument Argument;

    Determinant(Source f_) : f(f_), adj(f_) {}

    Determinant(Determinant const&) = default;
    Determinant& operator=(Determinant const&) = default;

    double d0() const { return computeViaLaplaceFormula<0>(); }

    double d1(Argument const& dF1) const { return computeViaLaplaceFormula<1>(dF1); }

    double d2(Argument const& dF1, Argument const& dF2) const { return computeViaLaplaceFormula<2>(dF1,dF2); }

    double d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return computeViaLaplaceFormula<3>(dF1,dF2,dF3); }

  private:
    template <size_t d, typename... Args>
    double computeViaLaplaceFormula(const Args&... args) const
    {
      double det = 0;
      for(size_t i=0; i<dim; ++i) det += computeSummand(i,std::integral_constant<size_t,d>(),args...);
      return det;
    }

    double computeSummand(size_t i, std::integral_constant<size_t,0>) const { return adj.d0()[0][i] * f.d0()[i][0]; }

    double computeSummand(size_t i, std::integral_constant<size_t,1>, Argument const& dF) const
    {
      return adj.d1(dF)[0][i] * f.d0()[i][0] + adj.d0()[0][i]*f.d1(dF)[i][0];
    }

    double computeSummand(size_t i, std::integral_constant<size_t,2>, Argument const& dF1, Argument const& dF2) const
    {
      return adj.d2(dF1,dF2)[0][i] * f.d0()[i][0] + adj.d1(dF1)[0][i] * f.d1(dF2)[i][0]
        + adj.d1(dF2)[0][i] * f.d1(dF1)[i][0] + adj.d0()[0][i] * f.d2(dF1,dF2)[i][0];
    }

    double computeSummand(size_t i, std::integral_constant<size_t,3>, Argument const& dF1, Argument const& dF2, Argument const& dF3) const
    {
      return adj.d3(dF1,dF2,dF3)[0][i] * f.d0()[i][0] + adj.d2(dF1,dF2)[0][i] * f.d1(dF3)[i][0]
        + adj.d2(dF1,dF3)[0][i] * f.d1(dF2)[i][0] + adj.d1(dF1)[0][i] * f.d2(dF2,dF3)[i][0]
        + adj.d2(dF2,dF3)[0][i] * f.d1(dF1)[i][0] + adj.d1(dF2)[0][i] * f.d2(dF1,dF3)[i][0]
        + adj.d1(dF3)[0][i] * f.d2(dF1,dF2)[i][0] + adj.d0()[0][i] * f.d3(dF1,dF2,dF3)[i][0];
    }

    Source f;
    Adjugate<dim,Source> adj;
  };


  template<class Scalar, int dim>
  Scalar determinant(Dune::FieldMatrix<Scalar, dim, dim> const& A)
  {
    return Determinant<dim>(A).d0();
  }
}

#endif // DETERMINANT_HH
