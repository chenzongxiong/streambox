/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef STRAINTENSORS_HH
#define STRAINTENSORS_HH

#include <utility>
#include <boost/math/constants/constants.hpp>
#include <dune/common/fmatrix.hh>

#include "fem/fixdune.hh"
#include "fem/diffops/elasto.hh"
#include "linalg/wrappedMatrix.hh"

namespace Kaskade
{
  /// \todo docme
  template <class Scalar, int dim, bool byValue = true>
  class CauchyGreenTensor
  {
  public:
    typedef Dune::FieldMatrix<Scalar,dim,dim> Argument;
    typedef Argument ReturnType;

    explicit CauchyGreenTensor(Argument const& F_) : F(F_) {}

    CauchyGreenTensor(CauchyGreenTensor const&) = default;
    CauchyGreenTensor& operator=(CauchyGreenTensor const&) = default;

    ReturnType d0() const { return transpose(F)*F; }

    ReturnType d1(Argument const& dF1) const
    {
      Argument const& G = transpose(dF1)*F;
      return G + transpose(G);
    }

    ReturnType d2(Argument const& dF1, Argument const& dF2) const
    {
      Argument const& G = transpose(dF1)*dF2;
      return G + transpose(G);
    }

    ReturnType d3(Argument const&, Argument const&, Argument const&) const { return zero; }

  private:
    typename std::conditional<byValue,Argument,Argument const&>::type F;
    Argument zero = Argument(0);
  };

  template <class Scalar,bool byValue>
  class CauchyGreenTensor<Scalar,3,byValue>
  {
    static constexpr int dim = 3;
  public:
    typedef Dune::FieldMatrix<Scalar,dim,dim> Argument;
    typedef Argument ReturnType;

    explicit CauchyGreenTensor(Argument const& F_) : F(F_) {}

    CauchyGreenTensor(CauchyGreenTensor const&) = default;
    CauchyGreenTensor& operator=(CauchyGreenTensor const&) = default;

    __attribute__((always_inline)) ReturnType d0() const
    {
      tmp = 0;
      for(size_t i=0; i<dim; ++i)
        for(size_t j=0; j<dim; ++j)
          for(size_t k=0; k<dim; ++k)
            tmp[i][j] += F[k][i]*F[k][j];

      return tmp;// transpose(F)*F;
    }

    __attribute__((always_inline)) ReturnType d1(Argument const& dF1) const
    {
      tmp = 0;
      for(size_t i=0; i<dim; ++i)
        for(size_t j=0; j<dim; ++j)
          for(size_t k=0; k<dim; ++k)
            tmp[i][j] += F[k][i]*dF1[k][j] + dF1[k][i]*F[k][j];

      return tmp;// G + transpose(G);
    }

    __attribute__((always_inline)) ReturnType d2(Argument const& dF1, Argument const& dF2) const
    {
      tmp = 0;
      for(size_t i=0; i<dim; ++i)
        for(size_t j=0; j<dim; ++j)
          for(size_t k=0; k<dim; ++k)
            tmp[i][j] += dF2[k][i]*dF1[k][j] + dF1[k][i]*dF2[k][j];
      return tmp;
    }

    __attribute__((always_inline)) ReturnType d3(Argument const&, Argument const&, Argument const&) const { return zero; }

  private:
    typename std::conditional<byValue,Argument,Argument const&>::type F;
    Argument zero = Argument(0);
    mutable Argument tmp;
  };

//   /// \todo docme
//   template <class Scalar, int dim>
//   class ElasticLinearizedGreenLagrangeTensor
//   {
//   public:
//     ElasticLinearizedGreenLagrangeTensor() = default;
// 
//     ElasticLinearizedGreenLagrangeTensor(Dune::FieldMatrix<Scalar,dim,dim> const& du, Dune::FieldMatrix<Scalar,dim,dim> const& inelasticStrain_)
//       : elasticStrain(du), inelasticStrain(inelasticStrain_)
//     {}
// 
//     auto d0() const { return elasticStrain.d0() - inelasticStrain; }
// 
//     template <class Arg>
//     auto d1(Arg const& arg) const { return elasticStrain.d1(arg); }
// 
//     template <class Arg>
//     auto d2(Arg const& arg1, Arg const& arg2) const { return elasticStrain.d2(arg1,arg2); }
// 
//     template <class Arg>
//     auto d3(Arg const& arg1, Arg const& arg2, Arg const& arg3) const { return elasticStrain.d3(arg1,arg2,arg3); }
// 
//   private:
//     LinearizedGreenLagrangeTensor<Scalar,dim> elasticStrain;
//     Dune::FieldMatrix<Scalar,dim,dim> inelasticStrain;
//   };
// 

  /// \todo docme
  template <class Scalar, int dim, class Tensor = WrappedMatrix<Scalar,dim> >
  class Deviator
  {
  public:
    typedef typename Tensor::Argument Argument;
    typedef Dune::FieldMatrix<Scalar,dim,dim> ReturnType;

    explicit Deviator(Tensor const& s_) : s(s_), I(unitMatrix<Scalar,dim>()) {}

    ReturnType d0() const
    {
      return s.d0() - boost::math::constants::third<Scalar>() * trace(s.d0()) * I;
    }

    ReturnType d1(Argument const& dF) const
    {
      return s.d1(dF) - boost::math::constants::third<Scalar>() * trace(s.d1(dF)) * I;
    }

    ReturnType d2(Tensor const& dF1, Tensor const& dF2) const
    {
      return s.d2(dF1,dF2) - boost::math::constants::third<Scalar>() * trace(s.d2(dF1,dF2)) * I;
    }

    ReturnType d3(Tensor const& dF1, Tensor const& dF2, Tensor const& dF3) const
    {
      return s.d3(dF1,dF2,dF3) - boost::math::constants::third<Scalar>() * trace(s.d3s(dF1,dF2,dF3)) * I;
    }

  private:
    Tensor s;
    ReturnType zero = ReturnType(0), I;
  };
}

#endif
