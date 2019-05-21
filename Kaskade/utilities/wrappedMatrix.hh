/*
 * wrappedMatrix.hh
 *
 *  Created on: 06.06.2014
 *      Author: bzflubko
 */

#ifndef WRAPPEDMATRIX_HH_
#define WRAPPEDMATRIX_HH_

namespace Kaskade
{

  template <class Scalar, size_t n, bool byValue=false>
  class WrappedMatrix
  {
  public:
    typedef Dune::FieldMatrix<Scalar,n,n> Argument;
    typedef Dune::FieldMatrix<Scalar,n,n> ReturnType;

    WrappedMatrix(Argument const& F_) : F(F_), zero(0.) {}

    WrappedMatrix(WrappedMatrix const&) = default;
    WrappedMatrix& operator=(WrappedMatrix const&) = default;

    ReturnType const& d0() const { return F; }

    ReturnType const& d1(Argument const& dF) const { return dF; }

    ReturnType d2(Argument const&, Argument const&) const { return zero; }

    ReturnType d3(Argument const&, Argument const&, Argument const&) const { return zero; }

  private:
    typename std::conditional<byValue,Argument,Argument const&>::type F;
    const Argument zero = Argument(0);
  };
}

#endif /* WRAPPEDMATRIX_HH_ */
