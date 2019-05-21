/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FUNCTIONAL_AUX_HH
#define FUNCTIONAL_AUX_HH

#include <type_traits>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include "fem/functionspace.hh"
/**
 * @file
 * @brief  Utility classes for the definition and use of variational functionals
 * @author Anton Schiela, Lars Lubkoll
 */

namespace Kaskade
{
  /**
   * \ingroup functional
   *
   * \brief Convenience base class for the easy definition of variational
   * functionals and weak formulations.
   *
   * Variational functionals and weak formulations need to provide meta
   * information template classes D1 and D2 that give the assembler
   * static information about the block structure of the problem. It is
   * easy to forget parts of these meta data, in particular if they do
   * not affect your current problem (yet they need to be defined for
   * syntactical reasons).
   *
   * Problem definitions may inherit this convenience base class, which
   * defines safe but maybe inefficient defaults. These can be partially
   * overridden as in the following example.
   *
   * \code
   * class Functional: public FunctionalBase<VariationalFunctional> {
   *   ...
   *   template <int row, int col>
   *   struct D2: public FunctionalBase<VariationalFunctional>::template D2<row,col> {
   *     static bool symmetric = true;
   *   };
   *   ...
   * };
   * \endcode
   * 
   * \seealso VariationalFunctional
   * 
   * \tparam Type the kind of variational problem, either \ref VariationalFunctional or \ref WeakFormulation
   */
  template <ProblemType Type>
  class FunctionalBase
  {
  public:

    static constexpr ProblemType type = Type;

    
    // TODO: docme
    template <class T>
    bool inDomain(T const&) const
    {
      return true;
    }
    
    /**
     * \brief Static a priori default information about the right hand side blocks.
     */
    template <int row>
    struct D1
    {
      /** 
       * \brief Whether the block is to be assembled at all.
       */
      static constexpr bool present = true;
      
      /**
       * \brief The maximum ansatz and test function derivative occuring in this block.
       * 
       * A default value of 1 is provided in the base class, which is appropriate for second order PDEs.
       */
      static constexpr int derivatives = 1;
    };

    /**
     * \brief Static a priori default information about the matrix blocks.
     * 
     * These flags declare structural properties of the individual blocks arising in the 
     * matrix of this problem and can be exploited by the assembler implementation for 
     * improved efficiency.
     */
    template <int row, int col>
    struct D2
    {
      /**
       * \brief Whether the block is to be assembled and stored or not.
       * The default value is the conservative one (but maybe inefficient) that all blocks are structurally nonzero. 
       * Note that for variational problems (which are necessarily symmetric), only the lower triangular part 
       * is accessed and hence in this case the default omits the superdiagonal blocks.
       * 
       * If present is false, the functional's d2 method will never be called for this block.
       */
      static constexpr bool present = type==WeakFormulation || row>=col;
      /**
       * \brief Whether the block is structurally symmetric (or hermitian).
       * The default is true for diagonal blocks of variational problems (which are necessarily symmetric), and
       * the conservative non-symmetric otherwise. Note that actually symmetric blocks are correctly stored if
       * this flag is false, but the storage may be less efficient. The assembler implementation may choose
       * whichever storage it likes, this flag is merely a hint that symmetry-exploiting optimizations are safe.
       */
      static constexpr bool symmetric = type==VariationalFunctional && row==col;
      /**
       * \brief If this flag is true, only the diagonal of this block will be assembled and stored.
       * This is mainly useful for hierarchical error estimators. Note that setting this flag true 
       * for mass matrices does not give a traditional lumped mass matrix (as the integration weights
       * are different).
       */
      static constexpr bool lumped = false;
      /**
       * \brief Whether to use a dense or sparse storage for the block.
       * The default is sparse, which is usually appropriate for PDEs, but in case of constant functions,
       * i.e. scalar parameters, dense blocks do usually appear, in which case a more efficient dense
       * storage may be chosen.
       */
      static constexpr bool dense = false;

      static constexpr bool constant = false;

      /**
       * \brief The maximum ansatz and test function derivative occuring in this block.
       * 
       * A default value of 1 is provided in the base class, which is appropriate for second order PDEs.
       */
      static constexpr int derivatives = 1;
    };

    template <int row, int col>
    struct B2
    {
      static bool const present = true;
      static bool const symmetric = false;
      static bool const lumped = false;
      static bool const constant = false;
    };
  };


  /// \cond internals
  namespace Functional_Aux_Detail
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
  /// \endcond


  /**
   * \ingroup functional
   *
   * \brief Provides implementations of d1 and d2 for (possibly mixed) use of scalar and power spaces.
   *
   * Implemented using CRTP. Just derive publicly from this functional if you want to use it.
   * Uses VariationalArg as perturbation. 
   * 
   * \todo docme better
   * \todo Explain difference between CacheBase and SimpleCacheBase
   * 
   * \tparam Functional actual functional (only transports static information)
   * \tparam Cache Functional::DomainCache or Functional::BoundaryCache, must implement
   * Scalar d1_impl(...);
   * and
   * Scalar d2_impl(...);
   */
  template <class Functional, class Cache>
  class CacheBase
  {
  public:
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::AnsatzVars AnsatzVars;
    typedef typename Functional::TestVars TestVars;

  private:
    template <int row>          using Vector = Dune::FieldVector<Scalar, TestVars::template Components<row>::m>;
    template <int row, int col> using Matrix = Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>;
    template <int row>          using TestComponents =   std::integral_constant<int,TestVars::template Components<row>::m>;
    template <int row>          using AnsatzComponents = std::integral_constant<int,AnsatzVars::template Components<row>::m>;

  public:
    /// default implementation of void moveTo(Cell/FaceIterator) does nothing
    template <class Entity>
    void moveTo(Entity const&){}

    /// default implementation of d0
    Scalar d0() const { return 0.; }

    /// forward to implementation via CRTP
    template <int row, int dimw>
    Scalar d1_local(VariationalArg<Scalar,dimw,TestComponents<row>::value> const& arg) const
    {
      return static_cast<Cache const*>(this)->template d1_impl<row>(arg);
    }

    template <int row, int dimw>
    Dune::FieldVector<Scalar,1> d1(VariationalArg<Scalar,dimw,TestComponents<row>::value> const& arg) const
    {
      return Dune::FieldVector<Scalar,1>(d1_local<row>(arg));
    }  

    /// scalar argument for multi-component variables
    template <int row, int dimw, class enable = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type>
    Vector<row> d1(VariationalArg<Scalar,dimw> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dimw,TestComponents<row>::value> arg;

      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg.value = 0; arg.value[i] = scalarArg.value[0];
        arg.derivative = 0; arg.derivative[i] = scalarArg.derivative[0];

        result[i] = d1_local<row>(arg);
      }

      return result;
    }

    /// forward to implementation via CRTP
    template <int row, int col, int dim>
    Scalar d2_local(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const 
    {
      return static_cast<Cache const*>(this)->template d2_impl<row,col>(arg1,arg2);
    }

    template <int row, int col, int dim,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<row>::value>::value>::type
    >
    Dune::FieldMatrix<Scalar,1,1> d2(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const 
    {
      return d2_local<row,col>(arg1,arg2);
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m2,
              class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
              class enable2 = typename std::enable_if<Functional_Aux_Detail::Equals<m2,AnsatzComponents<col>::value>::value>::type,
              class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> d2(VariationalArg<Scalar,dim> const& scalarArg, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg.derivative[0];

        result[i] = d2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m1,
              class enable1 = typename std::enable_if<Functional_Aux_Detail::Equals<m1,TestComponents<row>::value>::value>::type,
              class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
              class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
             >
    Vector<row> d2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,AnsatzComponents<row>::value> arg2;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg2.value = 0; arg2.value[i] = scalarArg.value[0];
        arg2.derivative = 0; arg2.derivative[i] = scalarArg.derivative[0];

        result[i] = d2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim>
    Matrix<row,col> d2(VariationalArg<Scalar,dim> const& scalarArg1, VariationalArg<Scalar,dim> const& scalarArg2) const
    {
      Matrix<row,col> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      VariationalArg<Scalar,dim,AnsatzComponents<col>::value> arg2;
      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg1.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg1.derivative[0];
        for(int j=0; j<AnsatzComponents<col>::value; ++j)
        {
          arg2.value = 0; arg2.value[j] = scalarArg2.value[0];
          arg2.derivative = 0; arg2.derivative[j] = scalarArg2.derivative[0];

          result[i][j] = d2_local<row,col>(arg1,arg2);
        }
      }

      return result;
    }


    /// forward to implementation via CRTP
    template <int row, int col, int dim>
    Scalar b2_local(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const
    {
      return static_cast<Cache const*>(this)->template b2_impl<row,col>(arg1,arg2);
    }

    template <int row, int col, int dim,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<row>::value>::value>::type
    >
    Dune::FieldMatrix<Scalar,1,1> b2(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const
    {
      return b2_local<row,col>(arg1,arg2);
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m2,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::Equals<m2,AnsatzComponents<col>::value>::value>::type,
    class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> b2(VariationalArg<Scalar,dim> const& scalarArg, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg.derivative[0];

        result[i] = b2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m1,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::Equals<m1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> b2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,AnsatzComponents<row>::value> arg2;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg2.value = 0; arg2.value[i] = scalarArg.value[0];
        arg2.derivative = 0; arg2.derivative[i] = scalarArg.derivative[0];

        result[i] = b2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim>
    Matrix<row,col> b2(VariationalArg<Scalar,dim> const& scalarArg1, VariationalArg<Scalar,dim> const& scalarArg2) const
    {
      Matrix<row,col> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      VariationalArg<Scalar,dim,AnsatzComponents<col>::value> arg2;
      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg1.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg1.derivative[0];
        for(int j=0; j<AnsatzComponents<col>::value; ++j)
        {
          arg2.value = 0; arg2.value[j] = scalarArg2.value[0];
          arg2.derivative = 0; arg2.derivative[j] = scalarArg2.derivative[0];

          result[i][j] = b2_local<row,col>(arg1,arg2);
        }
      }

      return result;
    }
  };


  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * *  Simple Convenient Cache  * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /// \cond internals
  namespace Functional_Aux_Detail
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
  /// \endcond

  /**
   * \ingroup functional
   *
   * \brief Provides implementations of d1 and d2 for (possibly mixed) use of scalar and power spaces.
   *
   * In order to use SimpleCacheBase the cache implementation must provide functions d1_impl(arg), d2_impl(arg1,arg2). These
   * functions implement the first resp. second derivative for given perturbations.
   *
   * Please do NOT provide functions d1, d2 directly in your cache as these would hide the implementations in
   * SimpleCacheBase, thus requiring two additional using directives in order to unhide these implementations.
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
   * \todo Explain difference between CacheBase and SimpleCacheBase
   *
   * Uses Dune::FieldVector and Dune::FieldMatrix as perturbation. If you use VariationalArg use CacheBase.
   * Implemented using CRTP.
   * \param Functional actual functional (only transports type information)
   * \param Cache = Functional::DomainCache or Functional::BoundaryCache, must implement 'Scalar d1_impl(...)' and 'Scalar d2_impl(...)'
   * \param ParameterType = Functional_Aux_Detail::Value/Functional_Aux_Detail::Gradient type of perturbation
   */
  template <class Functional, class Cache, class ParameterType>
  class SimpleCacheBase
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

    template <int row> using ReducedType = typename Functional_Aux_Detail::ArgumentTypes<Scalar,dim,NumberOfComponents<row>::value,ParameterType>::reduced_type;
    template <int row> using Type = typename Functional_Aux_Detail::ArgumentTypes<Scalar,dim,NumberOfComponents<row>::value,ParameterType>::type;

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
    Vector<row> d1_localImplementation(VariationalArg<Scalar,dim> const& arg, Functional_Aux_Detail::Gradient /* dummy */) const { return d1<row>(arg.derivative[0]); }

    template <int row>
    Vector<row> d1_localImplementation(VariationalArg<Scalar,dim> const& arg, Functional_Aux_Detail::Value /* dummy */) const { return d1<row>(arg.value[0]); }

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
    Matrix<row,col> d2_localImplementation(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2, Functional_Aux_Detail::Gradient /* dummy */) const { return d2<row,col>(arg1.derivative[0],arg2.derivative[0]); }

    template <int row, int col>
    Matrix<row,col> d2_localImplementation(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2, Functional_Aux_Detail::Value /* dummy */) const { return d2<row,col>(arg1.value[0],arg2.value[0]); }
  };

  /// \cond internals
  namespace Functional_Aux_Detail
  {
    /// helper function that detects if InnerBoundaryCache is declared in Functional
    template <class Functional>
    constexpr std::true_type hasInnerBoundaryCacheImpl(typename Functional::InnerBoundaryCache*);

    /// helper function that detects if InnerBoundaryCache is not declared in Functional
    template <class Functional>
    constexpr std::false_type hasInnerBoundaryCacheImpl(...);
  }
  /// \endcond

  /// Checks if InnerBoundaryCache is declared in Functional
  template <class Functional>
  constexpr bool hasInnerBoundaryCache()
  {
    typedef decltype(Functional_Aux_Detail::hasInnerBoundaryCacheImpl<Functional>(nullptr)) type;
    return type::value;
  }

  /// \cond internals
  namespace Functional_Aux_Detail
  {
    template <class Functional, class Linearization, bool hasInnerBoundaryCache> struct InnerBoundaryPolicyImpl;

    template <class Functional, class Linearization>
    struct InnerBoundaryPolicyImpl<Functional,Linearization,false>{};

    template <class Functional, class Linearization>
    struct InnerBoundaryPolicyImpl<Functional,Linearization,true>
    {
      typedef typename Functional::InnerBoundaryCache InnerBoundaryCache;
      typedef typename Functional::OriginVars::VariableSet PointOfLinearization;

      /// create inner boundary cache if existent
      InnerBoundaryCache createInnerBoundaryCache(int flags=7) const { return InnerBoundaryCache(static_cast<Linearization const*>(this)->getFunctional(),static_cast<Linearization const*>(this)->getOrigin(),flags); }
    };

    template <class Functional,class Linearization> using InnerBoundaryPolicy = InnerBoundaryPolicyImpl<Functional,Linearization,hasInnerBoundaryCache<Functional>()>;
  }
  /// \endcond

  /** 
   * \ingroup functional
   * \brief Proxy class for the linearization of a nonlinear functional. 
   * 
   * According to the following idea: A
   * VariationalFunctional defines a mapping \f$ F : U \to V \f$, which
   * can be evaluated and linearized in an appropriate ansatz space
   * \f$ U_a\subset U \f$, only if an argument \f$ u \in U_o \f$ in a -
   * possibly different - "origin" subspace \f$ U_o \subset U \f$ is
   * given.
   *
   * Hence, from a functional \f$ f \f$ and a variable \f$ u\in U_o \f$ one has to
   * create a linearization, which is represented by this class.
   */
  template<class Functional_>
  class LinearizationAt : public Functional_Aux_Detail::InnerBoundaryPolicy<Functional_,LinearizationAt<Functional_> >
  {
  public:
    typedef Functional_ Functional;
    static ProblemType const type = Functional::type;

    typedef typename Functional::Scalar Scalar;

    typedef typename Functional::AnsatzVars    AnsatzVars;
    typedef typename Functional::TestVars      TestVars;
    typedef typename Functional::OriginVars    OriginVars;

    typedef typename Functional::DomainCache   DomainCache;
    typedef typename Functional::BoundaryCache BoundaryCache;

    /**
     * \brief Type of variables around which the functional will be linearized.
     *
     * The type of "iterate", i.e. the point around which the functional
     * is linearized in the ansatz and test spaces. Note that this can
     * be, but need not be, the ansatz space itself. In any case, it must
     * refer to the same list of FE spaces.
     * 
     * \TODO: Why that? It must have the
     * same number of variables with the same number of components each,
     * but the underlying FE spaces can be different.
     */
    typedef typename Functional::OriginVars::VariableSet PointOfLinearization;

    /** 
     * \brief Construct a linearization from a given functional
     */
    LinearizationAt(Functional const & f_, PointOfLinearization const& u_) : f(f_), u(u_) {};
    /// create a domain cache
    DomainCache createDomainCache(int flags=7) const { return DomainCache(f,u,flags); }
    /// create a boundary cache
    BoundaryCache createBoundaryCache(int flags=7) const { return BoundaryCache(f,u,flags); }

    /**
     * \brief Rhs block information.
     */
    template <int row>
    struct D1: public Functional::template D1<row> {};

    /// matrix block information
    template <int row, int col>
    struct D2: public Functional::template D2<row,col> {};

    /// integration order
    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return f.template integrationOrder<Cell>(cell, shapeFunctionOrder, boundary);  // clang++
    }

    /// Get point of linearization
    PointOfLinearization const& getOrigin() const {return u;}
    /// Get functional
    Functional const& getFunctional() const {return f;}
    
    /**
     * \brief Gives the highest derivative that arises in the weak formulation.
     */
    int derivatives() const { return f.derivatives(); }

    protected:
    Functional const& f;
    PointOfLinearization const& u;
  };



  /// Convenience routine: construct linearization without having to know the type of the functional
  template<class Functional>
  LinearizationAt<Functional> linearization(Functional const& f, typename Functional::OriginVars::VariableSet const & u)
  {
    return LinearizationAt<Functional>(f,u);
  }

  /** 
   * \ingroup functional
   * \brief Proxy class for evaluating differences of linearizations of a nonlinear functional. 
   * 
   * Acceptance tests in Newton's method for minimization problems usually involves some kind 
   * of comparison \f$ f(u+\delta u) < f(u) \f$, where
   * \f[ f(u) = \sum_{i=1}^N f_i(u). \f]
   * Now taking the difference of both sums (large, because it's assembly over the grid) is 
   * numerically unstable for \f$ \delta u \f$ small. A more stable evaluation is
   * \f[ \sum_{i=1}^N (f_i(u)-f_i(u+\delta u)) > 0. \f]
   * This class supports the evaluation of the sum of differences.
   * 
   * According to the following idea: A VariationalFunctional defines a mapping \f$ F : U \to V \f$, which
   * can be evaluated and linearized in an appropriate ansatz space
   * \f$ U_a\subset U \f$, only if an argument \f$ u \in U_o \f$ in a -
   * possibly different - "origin" subspace \f$ U_o \subset U \f$ is
   * given.
   *
   * Hence, from a functional \f$ f \f$ and two variables \f$ u_1, u_2 \in U_o \f$ one has to
   * create a linearization of \f$ f(u_1)-f(u_2) \f$, which is represented by this class.
   */
  template<class Functional_>
  class LinearizationDifferenceAt : public Functional_Aux_Detail::InnerBoundaryPolicy<Functional_,LinearizationDifferenceAt<Functional_> >
  {
  public:
    typedef Functional_ Functional;
    static ProblemType const type = Functional::type;

    typedef typename Functional::Scalar Scalar;

    typedef typename Functional::AnsatzVars    AnsatzVars;
    typedef typename Functional::TestVars      TestVars;
    typedef typename Functional::OriginVars    OriginVars;


    /**
     * \brief Type of variables around which the functional will be linearized.
     *
     * The type of "iterate", i.e. the point around which the functional
     * is linearized in the ansatz and test spaces. Note that this can
     * be, but need not be, the ansatz space itself. It must have the
     * same number of variables with the same number of components each,
     * but the underlying FE spaces can be different.
     */
    typedef typename Functional::OriginVars::VariableSet PointOfLinearization;

    /** 
     * \brief Construct a difference of linearizations from a given functional
     */
    LinearizationDifferenceAt(Functional const & f_,PointOfLinearization const& u1_, PointOfLinearization const& u2_) : f(f_), u1(u1_), u2(u2_) {};
    
    class DomainCache
    {
    public:
      DomainCache(Functional const& f, PointOfLinearization const& u1, PointOfLinearization const& u2, int flags)
      : dc1(f,u1,flags), dc2(f,u2,flags) {}
      
      template <class Cell>
      void moveTo(Cell const& entity) { dc1.moveTo(entity); dc2.moveTo(entity); }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)  { dc1.evaluateAt(x,evaluators); dc2.evaluateAt(x,evaluators); }

      Scalar d0() const { return dc1.d0()-dc2.d0(); }
    
      template<int row, int dim> 
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
      d1 (VariationalArg<Scalar,dim> const& arg) const { return dc1.d1<row>(arg)-dc2.d1<row>(arg); }

      template<int row, int col, int dim> 
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const { return dc1.d2<row,col>(arg1,arg2)-dc2.d2<row,col>(arg1,arg2); }
      
      
    private:
      typename Functional::DomainCache dc1, dc2;
    };
    
    class BoundaryCache
    {
    public:
      BoundaryCache(Functional const& f, PointOfLinearization const& u1, PointOfLinearization const& u2, int flags)
      : bc1(f,u1,flags), bc2(f,u2,flags) {}
      
      template <class FaceIterator>
      void moveTo(FaceIterator const& entity) { bc1.moveTo(entity); bc2.moveTo(entity); }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators) { bc1.evaluateAt(x,evaluators); bc2.evaluateAt(x,evaluators); }

      Scalar d0() const { return bc1.d0()-bc2.d0(); }
    
      template<int row, int dim> 
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m> 
      d1 (VariationalArg<Scalar,dim> const& arg) const { return bc1.d1<row>(arg)-bc2.d1<row>(arg); }

      template<int row, int col, int dim> 
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const { return bc1.d2<row,col>(arg1,arg2)-bc2.d2<row,col>(arg1,arg2); }
      
    private:
      typename Functional::BoundaryCache bc1, bc2;
    };

    /// create a domain cache
    DomainCache createDomainCache(int flags=7) const {return DomainCache(f,u1,u2,flags);}
    /// create a boundary cache
    BoundaryCache createBoundaryCache(int flags=7) const {return BoundaryCache(f,u1,u2,flags);}

    /**
     * \brief Rhs block information.
     */
    template <int row>
    struct D1: public Functional::template D1<row> {};

    /// matrix block information
    template <int row, int col>
    struct D2: public Functional::template D2<row,col> {};

    /// integration order
    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return f.template integrationOrder<Cell>(cell, shapeFunctionOrder, boundary);
    }

    /// Get functional
    Functional const& getFunctional() const {return f;}

  protected:
    Functional const& f;
    PointOfLinearization const& u1;
    PointOfLinearization const& u2;
  };



  /// Convenience routine: construct linearization without having to know the type of the functional
  template<class Functional>
  LinearizationDifferenceAt<Functional> linearizationDifference(Functional const& f, 
                                                                typename Functional::OriginVars::VariableSet const & u1, 
                                                                typename Functional::OriginVars::VariableSet const & u2)
  {
    return LinearizationDifferenceAt<Functional>(f,u1,u2);
  }

  //---------------------------------------------------------------------

  /**
   * \ingroup functional
   *
   * \brief Proxy class for the linearization of semi-linear time stepping schemes.
   *
   * The sole difference to LinearizationAt is the second linearization
   * point \arg uJ, which is used to compute the Jacobian \f$ J \f$
   * occuring in \f[ B \dot u - J u = f(u) - J u. \f]
   */
  template <class Functional>
  class SemiLinearizationAt: public LinearizationAt<Functional>
  {
  public:
    /// Constructor
    /**
     * \param time dependent functional
     * \param u_ point of linearization for the "abstract" (in the sense of abstract functions, which are not really abstract) time differential operator (i.e. used for Cache::b2)
     * \param uJ_ point of linearization for spatial differential operator (i.e. used for Cache::d2)
     * \param du_ correction of last time step
     */
    SemiLinearizationAt(Functional const& f_,
        typename Functional::AnsatzVars::VariableSet const& u_,
        typename Functional::AnsatzVars::VariableSet const& uJ_,
        typename Functional::AnsatzVars::VariableSet const& du_):
          LinearizationAt<Functional>(f_,u_), uJ(uJ_), du(du_) {}

    typename Functional::DomainCache createDomainCache(int flags) const {
      return typename Functional::DomainCache(this->f,this->u,uJ,du,flags);
    }

    typename Functional::BoundaryCache createBoundaryCache(int flags) const {
      return typename Functional::BoundaryCache(this->f,this->u,uJ,du,flags);

    }
  private:
    typename Functional::AnsatzVars::VariableSet const& uJ;
    typename Functional::AnsatzVars::VariableSet const& du;
  };

  /**
   * \brief Convenience routine for constructing semi-linearizations.
   */
  template<class Functional>
  SemiLinearizationAt<Functional> semiLinearization(Functional const& f,
      typename Functional::AnsatzVars::VariableSet const& u,
      typename Functional::AnsatzVars::VariableSet const& uJ,
      typename Functional::AnsatzVars::VariableSet const& du)
      {
    return SemiLinearizationAt<Functional>(f,u,uJ,du);
      }

  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
  /**
   * Here you may provide easy to use caches that
   */
  //---------------------------------------------------------------------------------------------

  template <class Functional>
  class HomogeneousNeumannBoundary : public CacheBase<Functional,HomogeneousNeumannBoundary<Functional> >
  {
    typedef typename Functional::Scalar Scalar;
    static constexpr int dim = Functional::AnsatzVars::Grid::dimension;
    static constexpr int dimworld = Functional::AnsatzVars::Grid::dimensionworld;

  public:
    HomogeneousNeumannBoundary(Functional const&, typename Functional::OriginVars::VariableSet const&, int){}

    template <class FaceIterator>
    void moveTo(FaceIterator const&) {}

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename Functional::AnsatzVars::Grid::ctype, dim-1> const&, Evaluators const&) {}

    Scalar d0() const { return 0; }

    template <int row>
    Scalar d1_impl(VariationalArg<Scalar,dimworld,Functional::TestVars::template Components<row>::m> const&) const { return 0; }

    template <int row, int col>
    Scalar d2_impl(VariationalArg<Scalar,dimworld,Functional::TestVars::template Components<row>::m> const&,
                   VariationalArg<Scalar,dimworld,Functional::AnsatzVars::template Components<col>::m> const&) const { return 0; }
  };





  namespace Functional_Aux_Detail
  {
    template <bool>
    struct If
    {
      template <class Functor, typename... Args>
      static void evaluate(Functor const&, Args&...)
      {
        return;
      }

      template <class Functor, typename... Args>
      static double evaluateValue(Functor const&, Args&...)
      {
        return 0.0;
      }


      template <class Arg1, class Arg2>
      static auto multiply(Arg1 const&, Arg2 const& arg) -> decltype(arg*arg)
      {
        return 0.0;
      }

      template <class Arg1, class Arg2>
      static Arg2 const& chooseArgument(Arg1 const&, Arg2 const& arg)
      {
        return arg;
      }
    };

    template <>
    struct If<true>
    {
      template <class Functor, typename... Args>
      static void evaluate(Functor const& f, Args&... x)
      {
        f(x...);
      }

      template <class Functor, typename... Args>
      static double evaluateValue(Functor const& f, const Args&... x)
      {
        return f(x...);
      }

      template <class Arg1, class Arg2>
      static auto multiply(Arg1 const& x, Arg2 const& y) -> decltype(x*y)
      {
        return x*y;
      }

      template <class Arg>
      static Arg const& chooseArgument(Arg const& arg, Arg const&)
      {
        return arg;
      }
    };
  }

  /// Optionally evaluate functor/function pointer.
  template <bool enable, class Functor, typename... Args>
  void evaluateIf(Functor const& f, Args&... args)
  {
    Functional_Aux_Detail::If<enable>::evaluate(f,args...);
  }

  /// Optionally evaluate functor/function pointer.
  template <bool enable, class Functor, typename... Args>
  double evaluateValueIf(Functor const& f, Args&... args)
  {
    return Functional_Aux_Detail::If<enable>::evaluateValue(f,args...);
  }


  /// Optionally compute arg1*arg2.
  template <bool enable, class Arg1, class Arg2>
  auto multiplyIf(Arg1 const& arg1, Arg2 const& arg2) -> decltype(Functional_Aux_Detail::If<enable>::multiply(arg1,arg2))
  {
    return Functional_Aux_Detail::If<enable>::multiply(arg1,arg2);
  }

  template <class Arg1, class Arg2>
  Arg2 const& if_(Arg1 const& arg1, Arg2 const& arg2)
  {
    return Functional_Aux_Detail::If<std::is_same<std::remove_reference_t<Arg1>,std::remove_reference_t<Arg2>>::value>::chooseArgument(arg1,arg2);
  }

  template <class Arg1, class Arg2>
  Arg2& if_(Arg1& arg1, Arg2& arg2)
  {
    return Functional_Aux_Detail::If<std::is_same<std::remove_reference_t<Arg1>,std::remove_reference_t<Arg2>>::value>::chooseArgument(arg1,arg2);
  }
  
  template <class T, int n>
  void operator << (Dune::FieldVector<T,n>& target, double source) 
  {
    target = source;
  }
  
  /**
   * \ingroup functional
   * \brief Assignment if types match.
   * 
   * When the return type of a function depends on a template parameter (as often encountered in variational
   * functionals with variables of different components), the return expression has to be
   * legal even if it is statically guaranteed not to be executed. E.g., the following code is illegal:
   * \code
   * if (row==0)
   *   return Dune::FieldVector<double,1>(3.1415);
   * else
   *   return Dune::FieldVector<double,2>(2.71);
   * \endcode 
   * We can circumvent this deficiency of C++ by assigning a value only in case the types match:
   * \code
   * Dune::FieldVector<double,row==0?1:2> result;
   * if (row==0)
   *   result << Dune::FieldVector<double,1>(3.1415);
   * else
   *   result << 2.71;
   * return result;
   * \endcode
   */
  template <class T, int n, int m>
  void operator << (Dune::FieldVector<T,n>& target, Dune::FieldVector<T,m> const& source)
  {
    if (n==m)
      for (int i=0; i<n; ++i)
        target[i] = source[i];
  }
  
  template <class T, int n, int m>
  void operator << (Dune::FieldMatrix<T,n,m>& target, double source) 
  {
    target = source;
  }
  
  template <class T, int n1, int m1, int n2, int m2>
  void operator << (Dune::FieldMatrix<T,n1,m1>& target, Dune::FieldMatrix<T,n2,m2> const& source)
  {
    if (n1==n2 && m1==m2)
      for (int i=0; i<n1; ++i)
        for (int j=0; j<m1; ++j)
          target[i][j] = source[i][j];
  }
  
  
  
} /* end of namespace Kaskade */

#endif
