/*
 * wrappedMatrix.hh
 *
 *  Created on: 06.06.2014
 *      Author: bzflubko
 */

#ifndef WRAPPEDMATRIX_HH_
#define WRAPPEDMATRIX_HH_

#include <type_traits>

namespace Dune { template <class,int,int> class FieldMatrix; }

namespace Kaskade
{
  template <class Scalar, size_t n, bool argByValue=false>
  class WrappedMatrix
  {
  public:
    typedef Dune::FieldMatrix<Scalar,n,n> Argument;
    typedef Dune::FieldMatrix<Scalar,n,n> ReturnType;

    WrappedMatrix(Argument const& F_) : F(F_) {}

    WrappedMatrix(WrappedMatrix const&) = default;
    WrappedMatrix& operator=(WrappedMatrix const&) = default;

    ReturnType const& d0() const { return F; }
    ReturnType const& d1(Argument const& dF) const { return dF; }
    ReturnType d2(Argument const&, Argument const&) const { return Argument(0); }
    ReturnType d3(Argument const&, Argument const&, Argument const&) const { return Argument(0); }

  public:
    typename std::conditional<argByValue,Argument,Argument const&>::type F;
  };
}

#endif /* WRAPPEDMATRIX_HH_ */
