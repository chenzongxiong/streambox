/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FUNCTIONAL_MANIPULATION_HH
#define FUNCTIONAL_MANIPULATION_HH

namespace Kaskade
{
  template<class Function, class Evaluators, class VarArgs>
  class ComputeTestFunctions
  {
  public:
    ComputeTestFunctions(Function const& u_,
        Evaluators const& evaluators_,
        VarArgs& testfu_):
          u(u_), evaluators(evaluators_), testfu(testfu_)
    {
    }

    template<class Variable>
    void operator()(Variable const& var) const
    {
      int const id = Variable::id;
      int const spc = Variable::spaceIndex;
      using namespace boost::fusion;
      at_c<id>(testfu).value.result_of::value_at_c<VarArgs, id>::type::ValueType::operator=(at_c<id>(u.data).value(at_c<spc>(evaluators)));
      at_c<id>(testfu).gradient.result_of::value_at_c<VarArgs, id>::type::GradientType::operator=(at_c<id>(u.data).gradient(at_c<spc>(evaluators)));
    }
  private:
    Function const& u;
    Evaluators const& evaluators;
    VarArgs& testfu;
  };

  //----------------------------------------------------------------------------------------------------
  /** \ingroup functional
   * \brief Creates a (linear) model to a linearization of a functional.
   *
   * d0_(xs+dx) := d0(xs)+d1(xs)dx+1/2 d2(x)(dx,dx)
   * this is only useful, if d0, d1, and d2 are consistent!
   * in optimization: be careful with Lagrange-Functional, vs. Functional
   *
   * d1_(xs+dx) := d1(xs)+d2(x)dx
   * this is only useful, if d1 and d2 are consistent
   */
  template<class Lin, bool simplified=true>
  class QuadraticModel
  {
  public:

    typedef Lin Linearization;

    static ProblemType const type = Linearization::type;

    typedef typename Linearization::AnsatzVars        AnsatzVars;
    typedef typename Linearization::TestVars          TestVars;
    typedef typename Linearization::OriginVars        OriginVars;
    typedef typename Linearization::RT                RT;
    typedef typename Linearization::RT                Scalar;
    typedef typename AnsatzVars::Grid                 Grid;
    typedef typename AnsatzVars::Spaces               Spaces;
    typedef typename Grid::template Codim<0>::Entity  Entity;
    typedef typename Entity::LeafIntersectionIterator FaceIterator;

    typedef typename Linearization::AnsatzVars::VariableSet DomainElement;
    typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<Spaces, GetEvaluatorTypes>::type>::type Evaluators;

    QuadraticModel(Linearization const& f_, DomainElement const& dx_) : f(f_), ddf(f_), dx(dx_) {}

    QuadraticModel(Linearization const& f_, Linearization const& ddf_, DomainElement const& dx_) : f(f_), ddf(ddf_), dx(dx_) {}

    template<int row> 
    struct RowValueType
    {
      typedef typename boost::fusion::result_of::value_at_c<typename DomainElement::Functions,row>::type::ValueType type;
    };

    // To be used with accumulate.
    template<int row, int dim, class Cache>
    class HessianTimesDirection
    {
    public:
      template <class T> struct result {};

      template<typename T, typename State>
      struct result<HessianTimesDirection<row,dim,Cache>(T,State)>
      {
        typedef  typename RowValueType<row>::type type;
      };

      // arg_  : argument that would be used in a d1 evaluation
      // du_   : direction, in which d2(arg, .) is evaluated
      // cache_: cache_.d2 is evaluated
      // eval: : evaluators to compute shape functions.
      HessianTimesDirection(VariationalArg<RT,dim> const& arg_,  DomainElement const &du_,Cache const& cache_, Evaluators const& eval_):
        arg(arg_), du(du_), cache(cache_), evaluators(eval_)
      {
      }

      template<class Variable>
      typename RowValueType<row>::type operator()(Variable const& var,typename RowValueType<row>::type const& sum) const
      {
        using namespace boost::fusion;
        int const id = Variable::id;
        if(D2<row,id>::present)
        {
          int const csIdx = Variable::spaceIndex;
          int nShapeFcts = at_c<csIdx>(evaluators).size();
          typename RowValueType<row>::type s2(sum);
          for (int i=0; i<nShapeFcts; ++i) 
          {
            Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<id>::m> fm;
            // d2(shapefunction_values, arg)
            fm=cache.template d2<row,id,dim>(at_c<csIdx>(evaluators).evalData[i],arg);
            fm.umv((*at_c<id>(du.data))[at_c<csIdx>(evaluators).globalIndices()[i]],s2);
          }
          return s2;
        }
        else
        {
          return sum;
        }
      }
    private:
      VariationalArg<RT,dim> const& arg;
      DomainElement const& du;
      Cache const& cache;
      Evaluators const& evaluators;
    };

    template<class LinFunctional>
    class LinearFunctionalTimesDirection
    {
    public:
      template <class T> struct result {};

      template<typename T, typename State>
      struct result<LinearFunctionalTimesDirection<LinFunctional>(T,State)>
      {
        typedef RT type;
      };

      LinearFunctionalTimesDirection(DomainElement const &du_, LinFunctional const& fctl_, Evaluators const& eval_):
        du(du_), fctl(fctl_), evaluators(eval_)
      {
      }

      template<class Variable>
      RT operator()(Variable const& var, RT const& sum) const
      {
        using namespace boost::fusion;
        int const id = Variable::id;
        int const csIdx = Variable::spaceIndex;
        int nShapeFcts = at_c<csIdx>(evaluators).size();

        RT s2(sum);
        if(D1<id>::present)
          for (int i=0; i<nShapeFcts; ++i)
          {
            typename RowValueType<id>::type fv;
            fv = fctl.template evaluate<id>(at_c<csIdx>(evaluators).evalData[i]);
            s2 += fv*(*at_c<id>(du.data))[at_c<csIdx>(evaluators).globalIndices()[i]];
          }
        return s2;
      }
    private:
      DomainElement const& du;
      LinFunctional const& fctl;
      Evaluators const& evaluators;
    };



    struct DomainCache
    {
      typedef typename Linearization::DomainCache DomainCache1;
      DomainCache(int flags, Linearization const& f_, Linearization const& ddf_, DomainElement const& du_)
      :  du(du_), 
         fdc(f_.createDomainCache(flags)),
         ffdc(ddf_.createDomainCache(flags))
      {
        if(simplified) ddfdc=&ffdc;
        else
          ddfdc=&fdc;
      }

      template<class EntityOrFace>
      void moveTo(EntityOrFace const& e_)
      {
        fdc.moveTo(e_);
        if(simplified) ddfdc->moveTo(e_);
      }

      template<class Position>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        using namespace boost::fusion;
        fdc.evaluateAt(x,evaluators);
        if(simplified) ddfdc->evaluateAt(x,evaluators);
        eval=&evaluators;
      }

      RT d0() const {
        return boost::fusion::accumulate(vars,fdc.d0(),LinearFunctionalTimesDirection<DomainCache >(du,*this,*eval));
      }

      template <int row, int dim>
      typename RowValueType<row>::type d1(VariationalArg<RT,dim> const& arg) const
      {
        return boost::fusion::accumulate(vars,fdc.d1<row,dim>(arg),HessianTimesDirection<row,dim,DomainCache1>(arg,du,*ddfdc,*eval));
      }

      template <int row, int dim>
      typename RowValueType<row>::type evaluate(VariationalArg<RT,dim> const& arg) const
      {
        return 0.5*boost::fusion::accumulate(vars,2.0*fdc.d1<row,dim>(arg),HessianTimesDirection<row,dim,DomainCache1>(arg,du,*ddfdc,*eval));
      }

      template <int row, int col, int dim>
      Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
      {
        return fdc.d2<row,col,dim>(arg1,arg2);
      }

    private:
      DomainElement const& du;
      DomainCache1 fdc, ffdc;
      DomainCache1* ddfdc;
      typename AnsatzVars::Variables vars;
      Evaluators const* eval;
    };


    struct BoundaryCache
    {
      typedef typename Linearization::BoundaryCache BoundaryCache1;
      static const bool hasInteriorFaces = BoundaryCache1::hasInteriorFaces;

      BoundaryCache(int flags, Linearization const& f_, Linearization const& ddf_, DomainElement const& du_)
      :  du(du_), 
         fdc(f_.createBoundaryCache(flags)),
         ffdc(ddf_.createBoundaryCache(flags))
      {
        if(simplified) ddfdc=&ffdc;
        else
          ddfdc=&fdc;
      }

      template<class EntityOrFace>
      void moveTo(EntityOrFace const& e_)
      {
        fdc.moveTo(e_);
        if(simplified) ddfdc->moveTo(e_);
      }

      template<class Position>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        using namespace boost::fusion;
        fdc.evaluateAt(x,evaluators);
        if(simplified) ddfdc->evaluateAt(x,evaluators);
        eval=&evaluators;
      }

      RT d0() const { assert(!"not implemented"); return fdc.d0(); }

      template <int row, int dim>
      typename RowValueType<row>::type d1(VariationalArg<RT,dim> const& arg) const
      {
        return boost::fusion::accumulate(vars,fdc.d1<row,dim>(arg),HessianTimesDirection<row,dim,BoundaryCache1>(arg,du,*ddfdc,*eval));
      }

      template <int row, int col, int dim>
      Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
      {
        return fdc.d2<row,col,dim>(arg1,arg2);
      }

    private:
      DomainElement const& du;
      BoundaryCache1 fdc, ffdc;
      BoundaryCache1* ddfdc;
      typename AnsatzVars::Variables vars;
      Evaluators const* eval;
    };



    DomainCache createDomainCache(int flags=7) const {return DomainCache(flags,f,ddf,dx);}

    BoundaryCache createBoundaryCache(int flags=7) const {return BoundaryCache(flags,f,ddf,dx);}

    template <class Cell>
    int integrationOrder(Cell const& cell , int shapeFunctionOrder, bool boundary) const
    {
      return f.integrationOrder(cell,shapeFunctionOrder,boundary);
    }

    DomainElement const& getOrigin() const {return f.getOrigin(); }

    template <int row, int col>
    struct D2 
    {
      static bool const present = Linearization::template D2<row,col>::present;
      static bool const lumped = Linearization::template D2<row,col>::lumped;
      static bool const symmetric = true;
    };

    template <int row>
    struct D1 
    {
      static bool const present = Linearization::template D1<row>::present;
    };

    private:
    Linearization const& f;
    Linearization const& ddf;
    DomainElement const&  dx;
  };

  //----------------------------------------------------------------------------------------------------
  /** \ingroup functional
   * \brief Creates a (linear) model to a linearization of a functional.
   *
   * f := f1+f2
   */

  namespace ifff{
    template <bool c>
    struct if_
    {
      template <class A, class B>
      static A value(A const& a, B const& b) { return a; }
    };

    template <>
    struct if_<false>
    {
      template <class A, class B>
      static B value(A const& a, B const& b) { return b; }
    };

  };

  template<class Fu1, class Fu2>
  class SumFunctional
  {
  public:

    typedef Fu1 Functional;
    typedef Fu2 Functional2;

    static ProblemType const type = Functional::type;

    typedef typename Functional::AnsatzVars    AnsatzVars;
    typedef typename Functional::TestVars      TestVars;
    typedef typename Functional::OriginVars    OriginVars;
    typedef typename Functional::RT RT;
    typedef typename Functional::RT Scalar;
    typedef typename AnsatzVars::Grid Grid;
    typedef typename AnsatzVars::Spaces    Spaces;


    typedef typename Grid::template Codim<0>::Entity Entity;
    typedef typename Entity::LeafIntersectionIterator FaceIterator;
    typedef typename Functional::AnsatzVars::VariableSet Function;
    typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<Spaces, GetEvaluatorTypes>::type>::type Evaluators;


    SumFunctional(Functional const& f_, Functional2 const& f2_) : f(f_), f2(f2_) {}

    template<int row>
    struct RowValueType
    {
      typedef typename boost::fusion::result_of::value_at_c<typename Function::Functions,row>::type::ValueType type;
    };

    struct DomainCache
    {
      typedef typename Functional::DomainCache DomainCache1;
      typedef typename Functional2::DomainCache DomainCache2;
      DomainCache(int flags, Functional const& f_, Functional2 const& f2_)
      :  fdc(f_.createDomainCache(flags)),
         f2dc(f2_.createDomainCache(flags))
      {
      }

      template<class EntityOrFace>
      void moveTo(EntityOrFace const& e_)
      {
        fdc.moveTo(e_);
        f2dc.moveTo(e_);
      }

      template<class Position>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        using namespace boost::fusion;
        fdc.evaluateAt(x,evaluators);
        f2dc.evaluateAt(x,evaluators);
      }

      RT d0() const { return fdc.d0()+f2dc.d0(); }

      template <int row, int dim>
      typename RowValueType<row>::type d1(VariationalArg<RT,dim> const& arg) const
      {
        return fdc.d1<row,dim>(arg)+f2dc.d1<row,dim>(arg);
      }

      template <int row, int col, int dim>
      Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
      {

        Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);
        if(Functional::template D2<row,col>::present) result +=
            ifff::if_<Functional::template D2<row,col>::present>::value(fdc.d2<row,col,dim>(arg1,arg2),result);
        if(Functional2::template D2<row,col>::present) result +=
            ifff::if_<Functional2::template D2<row,col>::present>::value(f2dc.d2<row,col,dim>(arg1,arg2),result);


        //       if(Functional::template D2<row,col>::present) result += fdc.d2<row,col,dim>(arg1,arg2);
        //       if(Functional2::template D2<row,col>::present) result += f2dc.d2<row,col,dim>(arg1,arg2);

        return result;
      }

    private:
      DomainCache1 fdc;
      DomainCache2 f2dc;
    };


    struct BoundaryCache
    {
      typedef typename Functional::BoundaryCache BoundaryCache1;
      typedef typename Functional2::BoundaryCache BoundaryCache2;
      static const bool hasInteriorFaces = BoundaryCache1::hasInteriorFaces || BoundaryCache2::hasInteriorFaces;

      BoundaryCache(int flags, Functional const& f_, Functional2 const& f2_)
      :  fdc(f_.createBoundaryCache(flags)),
         f2dc(f2_.createBoundaryCache(flags))
      {
      }

      template<class EntityOrFace>
      void moveTo(EntityOrFace const& e_)
      {
        fdc.moveTo(e_);
        f2dc.moveTo(e_);
      }

      template<class Position>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        using namespace boost::fusion;
        fdc.evaluateAt(x,evaluators);
        f2dc.evaluateAt(x,evaluators);
      }

      RT d0() const { return fdc.d0()+f2dc.d0(); }

      template <int row, int dim>
      typename RowValueType<row>::type d1(VariationalArg<RT,dim> const& arg) const
      {
        typename RowValueType<row>::type result=fdc.d1<row,dim>(arg);

        result +=f2dc.d1<row,dim>(arg); ;
        return result;
      }

      template <int row, int col, int dim>
      Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
      {

        Dune::FieldMatrix<RT, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);
        if(Functional::template D2<row,col>::present) result +=
            ifff::if_<Functional::template D2<row,col>::present>::value(fdc.d2<row,col,dim>(arg1,arg2),result);
        if(Functional2::template D2<row,col>::present) result +=
            ifff::if_<Functional2::template D2<row,col>::present>::value(f2dc.d2<row,col,dim>(arg1,arg2),result);

        return result;
        //       if(!Functional::template D2<row,col>::present) return f2dc.d2<row,col,dim>(arg1,arg2);
        //       if(!Functional2::template D2<row,col>::present) return fdc.d2<row,col,dim>(arg1,arg2);
        //       return fdc.d2<row,col,dim>(arg1,arg2)+f2dc.d2<row,col,dim>(arg1,arg2);
      }

    private:
      BoundaryCache1 fdc;
      BoundaryCache2 f2dc;
    };



    DomainCache createDomainCache(int flags=7) const {return DomainCache(flags,f,f2);}

    BoundaryCache createBoundaryCache(int flags=7) const {return BoundaryCache(flags,f,f2);}

    template <class Cell>
    int integrationOrder(Cell const& cell , int shapeFunctionOrder, bool boundary) const
    {
      return std::max(f.integrationOrder(cell,shapeFunctionOrder,boundary),f2.integrationOrder(cell,shapeFunctionOrder,boundary));
    }

    Function const& getOrigin() const {return f.getOrigin(); }

    template <int row>
    struct D1 
    {
      static bool const present = Functional::template D1<row>::present || Functional2::template D1<row>::present;
    };

    template <int row, int col>
    struct D2 
    {
      static bool const present = Functional::template D2<row,col>::present || Functional2::template D2<row,col>::present;
      static bool const lumped = Functional::template D2<row,col>::lumped;
      static bool const symmetric = Functional::template D2<row,col>::symmetric && Functional2::template D2<row,col>::symmetric;
    };

    private:
    Functional const& f;
    Functional2 const& f2;
  };

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * Sum Caches  * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /**
   * Sum caches. As on the level of domain and boundary cache no information on symmetry and
   * existence of blocks is available it is assumed that this information coincides for both caches.
   */
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
} /* end of namespace Kaskade */

#endif
