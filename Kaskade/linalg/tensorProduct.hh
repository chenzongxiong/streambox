/*
 * tensorProduct.hh
 *
 *  Created on: 12.06.2014
 *      Author: bzflubko
 */

#ifndef TENSORPRODUCT_HH_
#define TENSORPRODUCT_HH_

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Kaskade
{
  template <class Scalar, int dim>
  Dune::FieldMatrix<Scalar,dim,dim> tensorProduct(Dune::FieldVector<Scalar,dim> const& v, Dune::FieldVector<Scalar,dim> const& w)
  {
    Dune::FieldMatrix<Scalar,dim,dim> result(0);

    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        result[i][j] = v[i]*w[j];

    return result;
  }
}


#endif /* TENSORPRODUCT_HH_ */
