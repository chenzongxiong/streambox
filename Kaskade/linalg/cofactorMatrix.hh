/*
 * cofactorMatrix.hh
 *
 *  Created on: June 23, 2014
 *      Author: bzflubko
 */
#ifndef COFACTORMATRIX_HH
#define COFACTORMATRIX_HH

#include "linalg/adjugate.hh"
#include "fem/fixdune.hh"

namespace Kaskade
{
  template<int dim, class Source = WrappedMatrix<double,dim,false> >
  class CofactorMatrix
  {
  public:
    typedef typename Adjugate<dim,Source>::Argument Argument;
    typedef typename Adjugate<dim,Source>::ReturnType ReturnType;

    explicit CofactorMatrix(Source const& s) : adj(s) {}

    ReturnType const d0() const { return transpose(adj.d0()); }

    ReturnType const d1(Argument const& dv) const { return transpose(adj.d1(dv)); }

    ReturnType const d2(Argument const& dv, Argument const& dw) const { return transpose(adj.d2(dv, dw)); }

    ReturnType const d3(Argument const& dv, Argument const& dw, Argument const& dx) const { return transpose(adj.d3(dv, dw, dx)); }

    Dune::FieldVector<double, dim> operator*(Dune::FieldVector<double, dim> const& v) const
    {
      return d0() * v;
    }

  private:
    Adjugate<dim, Source> const adj;
  };
}

#endif // COFACTORMATRIX_HH
