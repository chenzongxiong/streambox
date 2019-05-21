/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CONVENIENT_FUNCTIONAL_HH_
#define CONVENIENT_FUNCTIONAL_HH_

#include <type_traits>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include "fem/functionspace.hh"

namespace Kaskade
{
  namespace ConvenientFunctional_Detail
  {
    template <int a, int b>
    struct LessThan
    {
      static constexpr bool value = a<b;
    };

    template <int a, int b>
    struct Equals
    {
      static constexpr bool value = a==b;
    };
  }


  /// Provides overloads of d1 and d2 for (possibly mixed) use of scalar and power spaces.
  /**
   * Implemented using CRTP. Just derive publicly from this functional if you want to use it.
   * Uses VariationalArg as perturbation. If you use Dune::FieldVector/Dune::FieldMatrix use SimpleConvenientCache,
   *
   * \param Functional actual functional (only transports static information)
   * \param Cache Functional::DomainCache or Functional::BoundaryCache, must implement
   * Scalar d1_impl(...);
   * and
   * Scalar d2_impl(...);
   */
  template <class Functional, class Cache>
  class ConvenientCache
  {
  public:
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::AnsatzVars AnsatzVars;
    typedef typename Functional::TestVars TestVars;
    static int constexpr dim = Functional::dim;

  private:
    template <int row> using Vector = Dune::FieldVector<Scalar, TestVars::template Components<row>::m>;
    template <int row, int col> using Matrix = Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>;
    template <int row> using NumberOfComponents = std::integral_constant<int,AnsatzVars::template Components<row>::m>;

  public:
    /// forward to implementation via CRTP
    template <int row>
    Scalar d1(VariationalArg<Scalar,dim,NumberOfComponents<row>::value> const& arg) const
    {
      return static_cast<Cache const*>(this)->template d1_impl<row>(arg);
    }

    /// scalar argument for multi-component variables
    template <int row, class enable = typename std::enable_if<ConvenientFunctional_Detail::LessThan<1,NumberOfComponents<row>::value>::value>::type>
    Vector<row> d1(VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,NumberOfComponents<row>::value> arg;

      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        arg.value = 0; arg.value[i] = scalarArg.value[0];
        arg.gradient = 0; arg.gradient[i] = scalarArg.gradient[0];

        result[i] = d1<row>(arg);
      }

      return result;
    }

    /// forward to implementation via CRTP
    template <int row, int col, int m1, int m2>
    Scalar d2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      return static_cast<Cache const*>(this)->template d2_impl<row,col>(arg1,arg2);
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int m2,
          class enable1 = typename std::enable_if<ConvenientFunctional_Detail::LessThan<1,NumberOfComponents<row>::value>::value>::type,
          class enable2 = typename std::enable_if<ConvenientFunctional_Detail::Equals<m2,NumberOfComponents<row>::value>::value>::type
          >
    Vector<row> d2(VariationalArg<Scalar,dim> const& scalarArg, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,NumberOfComponents<row>::value> arg1;
      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg.value[0];
        arg1.gradient = 0; arg1.gradient[i] = scalarArg.gradient[0];

        result[i] = d2<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int m1,
          class enable1 = typename std::enable_if<ConvenientFunctional_Detail::Equals<m1,NumberOfComponents<row>::value>::value>::type,
          class enable2 = typename std::enable_if<ConvenientFunctional_Detail::LessThan<1,NumberOfComponents<row>::value>::value>::type
          >
    Vector<row> d2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,NumberOfComponents<row>::value> arg2;
      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        arg2.value = 0; arg2.value[i] = scalarArg.value[0];
        arg2.gradient = 0; arg2.gradient[i] = scalarArg.gradient[0];

        result[i] = d2<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col,
          class enable1 = typename std::enable_if<ConvenientFunctional_Detail::LessThan<1,NumberOfComponents<row>::value>::value>::type,
          class enable2 = typename std::enable_if<ConvenientFunctional_Detail::LessThan<1,NumberOfComponents<row>::value>::value>::type
          >
    Matrix<row,col> d2(VariationalArg<Scalar,dim> const& scalarArg1, VariationalArg<Scalar,dim> const& scalarArg2) const
    {
      Matrix<row,col> result(0);
      VariationalArg<Scalar,dim,NumberOfComponents<row>::value> arg1, arg2;
      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg1.value[0];
        arg1.gradient = 0; arg1.gradient[i] = scalarArg1.gradient[0];
        for(int j=0; j<NumberOfComponents<col>::value; ++j)
        {
          arg2.value = 0; arg2.value[j] = scalarArg2.value[0];
          arg2.gradient = 0; arg2.gradient[j] = scalarArg2.gradient[0];

          result[i][j] = d2<row,col>(arg1,arg2);
        }
      }

      return result;
    }
  };

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * *  Simple Convenient Cache  * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  namespace ConvenientFunctional_Detail
  {
    // type of used information of variational argument
    struct Gradient{};
    struct Value{};

    template <class,int,int,class> struct ArgumentTypes;

    template <class Scalar, int dim, int components>
    struct ArgumentTypes<Scalar,dim,components,Value>
    {
      typedef Dune::FieldVector<Scalar,components> type;
      typedef Scalar reduced_type;
    };

    template <class Scalar, int dim, int components>
    struct ArgumentTypes<Scalar,dim,components,Gradient>
    {
      typedef Dune::FieldMatrix<Scalar,components,dim> type;
      typedef Dune::FieldVector<Scalar,dim> reduced_type;
    };
  }

  /// Provides implementations of d1 and d2 for (possibly mixed) use of scalar and power spaces.
  /**
   * In order to use SimpleConvenientCache the cache implementation must provide functions d1_impl(arg), d2_impl(arg1,arg2). These
   * functions implement the first resp. second derivative for given perturbations.
   *
   * Please do NOT provide functions d1, d2 directly in your cache as these would hide the implementations in
   * SimpleConvenientCache.
   *
   * Provides the typedefs
   *  - Scalar
   *  - AnsatzVars
   *  - TestVars
   *  - OriginVars
   *  as well as the static (constexpr) constants
   *  - type
   *  - dim
   *
   * Uses Dune::FieldVector and Dune::FieldMatrix as perturbation. If you use VariationalArg use ConvenientCache.
   * Implemented using CRTP.
   * \param Functional actual functional (only transports type information)
   * \param Cache = Functional::DomainCache or Functional::BoundaryCache, must implement 'Scalar d1(...)' and 'Scalar d2(...)'
   * \param ParameterType = ConvenientFunctional_Detail::Value/ConvenientFunctional_Detail::Gradient type of perturbation
   */
  template <class Functional, class Cache, class ParameterType>
  class SimpleConvenientCache
  {
  public:
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::AnsatzVars AnsatzVars;
    typedef typename Functional::TestVars TestVars;
    typedef typename Functional::OriginVars OriginVars;
    static ProblemType const type = std::is_same<AnsatzVars,TestVars>::value ? VariationalFunctional : WeakFormulation;
    static int constexpr dim = Functional::dim;

  private: // some convenient template aliases (require gcc version >= 4.7)
    template <int row> using Vector = Dune::FieldVector<Scalar, AnsatzVars::template Components<row>::m>;
    template <int row, int col> using Matrix = Dune::FieldMatrix<Scalar, AnsatzVars::template Components<row>::m, AnsatzVars::template Components<col>::m>;
    template <int row> using NumberOfComponents = std::integral_constant<int,AnsatzVars::template Components<row>::m>;

    template <int row> using ReducedType = typename ConvenientFunctional_Detail::ArgumentTypes<Scalar,dim,NumberOfComponents<row>::value,ParameterType>::reduced_type;
    template <int row> using Type = typename ConvenientFunctional_Detail::ArgumentTypes<Scalar,dim,NumberOfComponents<row>::value,ParameterType>::type;

  public:
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * *  D1   * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    template <int row>
    Scalar d1(Type<row> const& du) const { return static_cast<Cache const*>(this)->template d1_impl<row>(du); } // CTRP

    template <int row>
    Vector<row> d1(ReducedType<row> const& dx) const
    {
      Vector<row> result(0);
      Type<row> du;

      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        du = 0; du[i] = dx;
        result[i] = d1<row>(du);
      }

      return result;
    }

    template <int row>
    Vector<row> d1(VariationalArg<Scalar,dim> const& arg) const { return d1_localImplementation<row>(arg,ParameterType()); }

  private:
    template <int row>
    Vector<row> d1_localImplementation(VariationalArg<Scalar,dim> const& arg, ConvenientFunctional_Detail::Gradient /* dummy */) const { return d1<row>(arg.gradient[0]); }

    template <int row>
    Vector<row> d1_localImplementation(VariationalArg<Scalar,dim> const& arg, ConvenientFunctional_Detail::Value /* dummy */) const { return d1<row>(arg.value[0]); }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * *  D2   * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  public:
    template <int row, int col>
    Scalar d2(Type<row> const& du, Type<row> const& dv) const { return static_cast<Cache const*>(this)->template d2_impl<row,col>(du,dv); } // CRTP

    template <int row, int col>
    Vector<row> d2(Type<row> const& du, ReducedType<row> const& dx) const
    {
      Vector<row> result(0);
      Type<row> dv;

      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        dv = 0; dv[i] = dx;
        result[i] = d2<row,col>(du,dv);
      }

      return result;
    }

    template <int row, int col>
    Vector<row> d2(ReducedType<row> const& dx, Type<row> const& dv) const
    {
      Vector<row> result(0);
      Type<row> du;

      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        du = 0; du[i] = dx;
        result[i] = d2<row,col>(du,dv);
      }

      return result;
    }

    template <int row, int col>
    Matrix<row,col> d2(ReducedType<row> const& dx1, ReducedType<row> const& dx2) const
    {
      Matrix<row,col> result(0);
      Type<row> du, dv;

      for(int i=0; i<NumberOfComponents<row>::value; ++i)
      {
        du = 0; du[i] = dx1;
        for(int j=0; j<NumberOfComponents<col>::value; ++j)
        {
          dv = 0; dv[j] = dx2;
          result[i][j] = d2<row,col>(du,dv);
        }
      }

      return result;
    }

    template <int row, int col>
    Matrix<row,col> d2(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2) const { return d2_localImplementation<row,col>(arg1,arg2, ParameterType()); }

  private:
    template <int row, int col>
    Matrix<row,col> d2_localImplementation(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2, ConvenientFunctional_Detail::Gradient /* dummy */) const { return d2<row,col>(arg1.gradient[0],arg2.gradient[0]); }

    template <int row, int col>
    Matrix<row,col> d2_localImplementation(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2, ConvenientFunctional_Detail::Value /* dummy */) const { return d2<row,col>(arg1.value[0],arg2.value[0]); }
  };


  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * Sum Caches  * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  template <class Cache1, class Cache2>
  class SumCache
  {
  public:
    SumCache(Cache1 const& cache1_, Cache2 const& cache2_) : cache1(cache1_), cache2(cache2_)
    {}

    auto d0() const -> decltype(std::declval<Cache1>().d0()+std::declval<Cache2>().d0())
    {
      return cache1.d0() + cache2.d0();
    }

    template <int row, class Arg>
    auto d1(Arg const& arg) const -> decltype(std::declval<Cache1>().template d1<row>(arg)+std::declval<Cache2>().template d1<row>(arg))
    {
      return cache1.template d1<row>(arg) + cache2.template d1<row>(arg);
    }

    template <int row, int col, class Arg1, class Arg2>
    auto d2(Arg1 const& arg1, Arg2 const& arg2) const
    -> decltype(std::declval<Cache1>().template d2<row,col>(arg1,arg2)+std::declval<Cache2>().template d2<row,col>(arg1,arg2))
    {
      return cache1.template d2<row,col>(arg1,arg2) + cache2.template d2<row,col>(arg1,arg2);
    }

  private:
    Cache1 const& cache1;
    Cache2 const& cache2;
  };

} // end of namespace Kaskade

#endif /* CONVENIENT_FUNCTIONAL_HH_ */
