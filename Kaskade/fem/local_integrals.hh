/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LOCAL_INTEGRALS_HH
#define LOCAL_INTEGRALS_HH

#include "fem/assemble.hh"

/** \ingroup functional
 * \brief Helps to evaluate integral functionals over subregions
 */
namespace Kaskade
{
  template<class TestVariables, class Functional, class NewFktVariables, class Function>
  class LocalValueFunctional
  {
  public:
    typedef typename TestVariables::Grid Grid;
    typedef typename Functional::RT RT;
    typedef typename Grid::template Codim<0>::Entity Entity;

    typedef TestVariables AnsatzVars;
    typedef TestVariables TestVars;
    static ProblemType const type = WeakFormulation;
    typedef typename Functional::AnsatzVars::VariableSet ArgOfFunctional;

    typedef typename AnsatzVars::Spaces    Spaces;

    typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<Spaces, GetEvaluatorTypes>::type>::type Evaluators;

    typedef typename boost::fusion::result_of::as_vector
    <typename  boost::fusion::result_of::transform<NewFktVariables,AssemblyDetail::RhsLocalData<RT> >::type>::type LocalRhs;


    template<class Arg>
    bool inDomain(Arg const& a) const { return true; }

    LocalValueFunctional(Functional const& f_, Function const& dx_) : f(f_), dx(dx_)
    {
    }

    struct DomainCache : public EvalCacheBase
    {
      DomainCache(Functional const& f_, Function const& u_) : f(f_), u(u_), fdc(f.createDomainCache())
      {
      }

      void moveTo(Entity const& e_)
      {
        fdc.moveTo(e_);
      }

      template<class Position>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        using namespace boost::fusion;
        using namespace AssemblyDetail;
        for_each(insertFuVariables,ClearLocalRhs<LocalRhs,Evaluators>(localRhs,evaluators));
        fdc.evaluateAt(x,evaluators);
        for_each(insertFuVariables,UpdateLocalRhs<LocalRhs, Evaluators, RT, typename Functional::DomainCache>(localRhs,evaluators,1.0,fdc));
        eval=&evaluators;
      }

      RT d0() const { return 0; }

      template <int row, int dim>
      typename boost::fusion::result_of::value_at_c<typename Function::Functions,row>::type::ValueType d1(VariationalArg<RT,dim> const& arg) const
      {
        using namespace boost::fusion;
        int const sIdx = result_of::value_at_c<NewFktVariables, row>::type::spaceIndex;
        int const m = result_of::value_at_c<NewFktVariables, row>::type::m;
        typename result_of::value_at_c<Evaluators,sIdx>::type const& reval = at_c<sIdx>(*eval);
        typename result_of::value_at_c<typename Function::Functions,row>::type::ValueType sum(0.0);
        for (size_t i=0; i<reval.size(); i++)
        {
          for(int j=0; j<m; ++j)
            sum[j] += (*at_c<row>(u.vars))[reval.globalIndices()[i]][j]*at_c<row>(localRhs)[i][j];
        }
        sum*=arg.value[0];
        return sum;
        // Cannot work for multi-component shapefunctions
      }

      template <int row, int col, int dim>
      Dune::FieldMatrix<RT,1,1> d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
      {
        return 0;
      }

    private:
      LocalRhs localRhs;
      Functional const& f;
      Function const& u;
      typename Functional::DomainCache fdc;
      NewFktVariables insertFuVariables;
      Evaluators const* eval;
    };

    struct BoundaryCache : public HomNeumannBoundaryCache<RT> {};

    DomainCache createDomainCache(int flags=7) const {return DomainCache(f,dx);}
    BoundaryCache createBoundaryCache(int flags=7) const {return BoundaryCache();}

    Functional const& f;
    Function const&  dx;

    public:
    template <int row, int col>
    struct D2
    {
      static bool const present = Functional::template D2<row,col>::present;
      static bool const lumped = Functional::template D2<row,col>::lumped;
      static bool const symmetric = true;
    };

    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
    {
      int matrixOrder =  boundary? 2*shapeFunctionOrder: 2*(shapeFunctionOrder);
      int lastIterateOrder = 4;

      return matrixOrder+lastIterateOrder;
    }
  };


  template <template <class TVrs, class Fctl, class NFkt, class Fkt> class Operator, class Functional, class FktVariableSet, class Quadrature, class TestSpace>
  class LocalValues
  {
    typedef typename FktVariableSet::VariableSet Function;
    typedef typename Functional::AnsatzVars::Grid Grid;
  public:

  private:
    static int const spaceSize = boost::fusion::result_of::size<typename Functional::AnsatzVars::Spaces>::type::value;

    typedef typename  boost::fusion::result_of::join<typename Functional::AnsatzVars::Spaces, typename FktVariableSet::Spaces>::type Spaces2Seq;
    typedef typename boost::fusion::result_of::push_back<Spaces2Seq, TestSpace*>::type Spaces3Seq;
    typedef typename boost::fusion::result_of::as_vector<Spaces3Seq>::type Spaces3;

    template<int sidx>
    struct SpaceIdSet {
      template <class T> struct result {};

      template<class VD>
      struct result<SpaceIdSet<sidx>(VD)> {
        typedef VariableDescription<sidx,VD::m,VD::id> type;
      };
    };

    template<int sidx>
    struct SpaceIdAdd {
      template <class T> struct result {};

      template<class VD>
      struct result<SpaceIdAdd<sidx>(VD)> {
        typedef VariableDescription<sidx+VD::spaceIndex,VD::m,VD::id> type;
      };
    };

    static int const spaceSize3 = boost::fusion::result_of::size<Spaces3>::type::value;
    typedef typename boost::fusion::result_of::transform<typename Functional::TestVars::Variables, SpaceIdSet<spaceSize3-1> >::type TestVarView;
    typedef typename boost::fusion::result_of::as_vector<TestVarView>::type NewTestVariables;

    typedef typename boost::fusion::result_of::transform<typename FktVariableSet::Variables, SpaceIdAdd<spaceSize> >::type FktVarView;
    typedef typename boost::fusion::result_of::as_vector<FktVarView>::type NewFktVariables;

  public:

    typedef VariableSetDescription<Spaces3,NewTestVariables> TestVariables;
    typedef VariableSetDescription<Spaces3,typename Functional::AnsatzVars::Variables> AnsatzVariables;

    LocalValues(Functional const& f_, Function const& dx_, int order_): f(f_),dx(dx_), order(order_)
    {
    }

    void getFunction(std::vector<double>& rhs)
    {
      using namespace boost::fusion;

      typedef Operator<TestVariables,Functional,NewFktVariables,Function> FKT;

      FKT LVF(f,dx);

      typename Functional::AnsatzVars::Spaces const& s1(f.getOrigin().descriptions.spaces);
      typename FktVariableSet::Spaces const& s2(dx.descriptions.spaces);
      TestSpace dspace(at_c<0>(dx.vars).space().gridManager(), at_c<0>(dx.vars).space().indexSet(),order);
      Spaces3 spaces3(as_vector(push_back(join(s1,s2),&dspace)));

      typename Functional::AnsatzVars::Variables v1;

      NewTestVariables v2;

      AnsatzVariables ansatzSet(spaces3);
      TestVariables testSet(spaces3);
      // construct Galerkin representation
      typedef VariationalFunctionalAssembler<FKT,
      CachingBoundaryDetector<typename FKT::AnsatzVars::Grid,
      typename FKT::AnsatzVars::IndexSet>,
      Quadrature> Assembler;

      Assembler gop(spaces3);

      size_t size = testSet.degreesOfFreedom(0,TestVariables::noOfVariables);

      rhs.resize(size);

      gop.assemble(LVF,2);
      gop.toSequence(0,TestVariables::noOfVariables,rhs.begin());
    }


  private:

    Functional const& f;
    Function const&  dx;
    int const order;

  };
} // end of namespace Kaskade

#endif
