#ifndef IPFUNCTIONAL_HH
#define IPFUNCTIONAL_HH

namespace Kaskade
{

/**
 * @file 
 * @brief Functional for Barrier/Interior Point methods
 * @author Anton Schiela
 *
 */
#include "fem/functional_aux.hh"

//---------------------------------------------------------------------------
// Barrier functions of various orders


/// Barrier<q,p> = sum_{i=p}^{q} Barrier<i,i> 
template<int q, int p>
struct Barrier
{
  static const int order = q;
  static const int loworder = p;

  static double b(double mu, double mudx)
  {
    return Barrier<q,q>::b(mu,mudx)+Barrier<q-1,p>::b(mu,mudx);
  }

  static double db(double mu, double mudx)
  {
    return Barrier<q,q>::db(mu,mudx)+Barrier<q-1,p>::db(mu,mudx);
  }

  static double ddb(double mu, double mudx)
  {
    return Barrier<q,q>::ddb(mu,mudx)+Barrier<q-1,p>::ddb(mu,mudx);
  }

  static double bmu(double mu, double mudx)
  {
    return Barrier<q,q>::bmu(mu,mudx)+Barrier<q-1,p>::bmu(mu,mudx);
  }

  static double dbmu(double mu, double mudx)
  {
    return Barrier<q,q>::dbmu(mu,mudx)+Barrier<q-1,p>::dbmu(mu,mudx);
  }
};

/// Rational Barrier functions of order q/2+1 (gradient has rational order q/2+1)
template<int q>
struct Barrier<q,q>
{
  static const int order = q;
  static const int loworder = q;

  static double b(double mu, double mudx)
  {
    return mu*std::pow(mudx,q/2.0)/(q/2.0);
  }

  static double db(double mu, double mudx)
  {
    return -std::pow(mudx,q/2.0+1.0);
  }

  static double ddb(double mu, double mudx)
  {
    return (q/2.0+1.0)*std::pow(mudx,q/2.0+2.0)/mu;
  }

  static double bmu(double mu, double mudx)
  {
    return std::pow(mudx,q/2.0)*(2.0/q+1);
  }

  static double dbmu(double mu, double mudx)
  {
    return -(q/2.0+1.0)*std::pow(mudx,q/2.0+1.0)/mu;
  }
};


/// Logarithmic Barrier functions
template<> struct Barrier<0,0>
{
  static const int order = 0;
  static const int loworder = 0;

  static double b(double mu, double mudx)
  {
    return -mu*std::log(mu/mudx);
  }

  static double db(double mu, double mudx)
  {
    return -mudx;
  }

  static double ddb(double mu, double mudx)
  {
    return mudx*mudx/mu;
  }

  static double bmu(double mu, double mudx)
  {
    return -std::log(mu/mudx);
  }

  static double dbmu(double mu, double mudx)
  {
    return -mudx/mu;
  }
};

/// Functional that adds barrier terms to VarFu
template <class VarFu, class BarrierFu, int paralin=0>
class IPFunctional
{
public:
  static const int parameterLin = paralin;
  typedef IPFunctional<VarFu,BarrierFu,2> ParameterLinearization;
  typedef IPFunctional<typename VarFu::FunctionalDiagonal,BarrierFu,1> ParameterValueLinearization;
  typedef double Parameter;
private:
  template<class VarF, int parmeterln> 
  struct FunctionalTraits
  {
    typedef typename VarF::Functional OptimalityFunctional;
  };

  template<class VarF>  struct FunctionalTraits<VarF,1>
  {
    typedef typename VarF::FunctionalDiagonal OptimalityFunctional;
  };

public:

  typedef typename VarFu::ConstraintsCache ConstraintsCache;

  typedef typename VarFu::Scalar Scalar;

  typedef typename FunctionalTraits<VarFu, paralin>::OptimalityFunctional OptimalityFunctional;

  typedef typename VarFu::OriginVars OriginVars;
  typedef typename VarFu::AnsatzVars AnsatzVars;
  typedef typename VarFu::TestVars TestVars;
  static int const nrows = TestVars::noOfVariables;
  static ProblemType const type = VariationalFunctional;

  typedef BarrierFu BarrierFunction;

  static const int order = BarrierFu::order;
  static const int loworder = BarrierFu::loworder;

  typedef IPFunctional<OptimalityFunctional,BarrierFu,paralin> Functional;
  typedef typename VarFu::Entity Entity;


  static int const yIdx = VarFu::yIdx;
  static int const ySIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,yIdx>::type::spaceIndex;

  IPFunctional(Parameter mu_, OptimalityFunctional& f) : mu(mu_), unconstrainedFunctional(f) 
  { 
  }

  template<class DomainVector>
  void prepareConstraintsCache(DomainVector const& x) const
  {
    unconstrainedFunctional.prepareConstraintsCache(x);
  }

  template<class Arg>
  bool inDomain(Arg const& a) const 
  { 
    return true;
  }

private:

  template<class Evaluators>
  class ComputeBoundDifference
  {
  public:
    ComputeBoundDifference(typename AnsatzVars::VariableSet const& u_, 
                           ConstraintsCache const& dc_,
                           Evaluators const& evaluators_, 
                           Dune::FieldVector<Scalar,AnsatzVars::noOfVariables>& dylower_, 
                           Dune::FieldVector<Scalar,AnsatzVars::noOfVariables>& dyupper_):
      u(u_), dc(dc_), evaluators(evaluators_), dylower(dylower_), dyupper(dyupper_)
    {
    }
    
    template<class Variable>
    void operator()(Variable const& var) const
    {
      int const row = Variable::id;
      int const spc = Variable::spaceIndex;
      using namespace boost::fusion;   
      if(ConstraintsCache::template bounds<row>::upper)
      {
        dyupper[row]=dc.upperbound()-at_c<row>(u.data).value(at_c<spc>(evaluators))[0];
      }
      if(ConstraintsCache::template bounds<row>::lower)
      {
        dylower[row]=at_c<row>(u.data).value(at_c<spc>(evaluators))[0]-dc.lowerbound();
      }
    }
  private:
    typename AnsatzVars::VariableSet const& u;
    ConstraintsCache const& dc;
    Evaluators const& evaluators;
    Dune::FieldVector<Scalar,AnsatzVars::noOfVariables> &dylower, &dyupper;
  };

  template<int present, int component, int paraln>
  struct boundTraits
  {
    static Scalar b0(Parameter mu, Scalar yp) {return 0;}
    static Scalar b1(Parameter mu, Scalar yp) {return 0;}
    static Scalar b2(Parameter mu, Scalar yp) {return 0;}
  };

  template<int component>
  struct boundTraits<1, component, 0>
  {      
    static Scalar b0(Parameter mu, Scalar yp) { return BarrierFu::b(mu,mu/yp); }
    static Scalar b1(Parameter mu, Scalar yp) { return BarrierFu::db(mu,mu/yp);}
    static Scalar b2(Parameter mu, Scalar yp) { return BarrierFu::ddb(mu,mu/yp);}
  };

  template<int component>
  struct boundTraits<1, component, 1>
  {      
    static Scalar b0(Parameter mu, Scalar yp) { return BarrierFu::bmu(mu,mu/yp); }
    static Scalar b1(Parameter mu, Scalar yp) { return BarrierFu::db(mu,mu/yp);}
    static Scalar b2(Parameter mu, Scalar yp) { return BarrierFu::ddb(mu,mu/yp);}
  };

  template<int component>
  struct boundTraits<1, component, 2>
  {      
    static Scalar b0(Parameter mu, Scalar yp) { return BarrierFu::b(mu,mu/yp); }
    static Scalar b1(Parameter mu, Scalar yp) { return BarrierFu::dbmu(mu,mu/yp);}
    static Scalar b2(Parameter mu, Scalar yp) { return BarrierFu::ddb(mu,mu/yp);}
  };


public:
  
  struct DomainCache : public EvalCacheBase
  {    
    
    DomainCache(Functional const& f_, typename AnsatzVars::VariableSet const& x_, int flags_=7) : 
      x(x_), 
      f(f_), 
      fcc(f.getFunctional().createConstraintsCache(x_)),
      pmu(f.mu)
    {
    }

    void moveTo(Entity const& e_)
    {
      e=&e_;
      fcc.moveTo(e_);
    }

    template<class Position, class Evaluators>
    void evaluateAt(Position const& y, Evaluators const& evaluators) 
    {
      fcc.evaluateAt(y, evaluators);
      boost::fusion::for_each(typename AnsatzVars::Variables(),ComputeBoundDifference<Evaluators>(x,fcc,evaluators,dylower,dyupper));

      b1[0]=boundTraits<ConstraintsCache::template bounds<0>::lower, 0, paralin>::b1(pmu,dylower[0])
        -boundTraits<ConstraintsCache::template bounds<0>::upper, 0, paralin>::b1(pmu,dyupper[0]);
      b1[1]=boundTraits<ConstraintsCache::template bounds<1>::lower, 1, paralin>::b1(pmu,dylower[1])
        -boundTraits<ConstraintsCache::template bounds<1>::upper, 1, paralin>::b1(pmu,dyupper[1]);

      b2[0]=boundTraits<ConstraintsCache::template bounds<0>::lower, 0, paralin>::b2(pmu,dylower[0])
        +boundTraits<ConstraintsCache::template bounds<0>::upper, 0, paralin>::b2(pmu,dyupper[0]);
      b2[1]=boundTraits<ConstraintsCache::template bounds<1>::lower, 1, paralin>::b2(pmu,dylower[1])
        +boundTraits<ConstraintsCache::template bounds<1>::upper, 1, paralin>::b2(pmu,dyupper[1]);
    }

    Scalar d0() const { 
      Scalar sum=0.0;
      sum += boundTraits<ConstraintsCache::template bounds<1>::upper, 1, paralin>::b0(pmu,dyupper[1]);
      sum += boundTraits<ConstraintsCache::template bounds<2>::upper, 2, paralin>::b0(pmu,dyupper[2]);
      return sum;
    }

    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const
    {
      Dune::FieldVector<Scalar,TestVars::template Components<row>::m> result(0);
      if(TestVars::template Components<row>::m==1) result[0]=b1[row]*arg.value[0];
      return result;
    }

    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2) const
    {
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);
      if(row==col && TestVars::template Components<row>::m==1 &&TestVars::template Components<col>::m==1)
        result[0][0]=b2[row]*arg1.value[0]*arg2.value[0];
      return result;
    }

//     template <int row, int dim>
//     Dune::FieldVector<Scalar,1> d1(VariationalArg<Scalar,dim> const& arg) const
//     {
//       return (boundTraits<ConstraintsCache::template bounds<row>::lower, row, paralin>::b1(pmu,dylower[row])-
//               boundTraits<ConstraintsCache::template bounds<row>::upper, row, paralin>::b1(pmu,dyupper[row]))
//         *arg.value[0];
//     }

//     template <int row, int col, int dim>
//     Dune::FieldMatrix<Scalar,1,1> d2(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2) const
//     {
//       return (boundTraits<ConstraintsCache::template bounds<col>::upper, col, paralin>::b2(pmu,dyupper[col])+ 
//               boundTraits<ConstraintsCache::template bounds<col>::lower, col, paralin>::b2(pmu,dylower[col]))
//         *arg1.value[0]*arg2.value[0];
//     }
   
  private:
    Dune::FieldVector<Scalar,2> b1,b2;
    typename AnsatzVars::VariableSet const& x;
    Functional const& f;
    ConstraintsCache fcc;
    Entity const* e;
    Parameter pmu;
    Dune::FieldVector<Scalar,AnsatzVars::noOfVariables> dyupper,dylower;
  };


  struct BoundaryCache : public EvalCacheBase
  {
    static const bool hasInteriorFaces = false;

    BoundaryCache(Functional const& , typename AnsatzVars::VariableSet const&, int flags=7) {};
    
    template<class Position, class Evaluators>
    void evaluateAt(Position const& /* x */, Evaluators const& /* evaluators */) 
    {
    }    

    Scalar d0() const 
    {
        return 0;
    }

    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const
    {
      Dune::FieldVector<Scalar,TestVars::template Components<row>::m> result(0);
        return result;
    }

    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2) const
    {
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);
      return result;
    }
  };
  

  template <int row>
  struct D1
  {
    static bool const present = (ConstraintsCache::template bounds<row>::upper || ConstraintsCache::template bounds<row>::lower);
  };
     
  template <int row, int col>
  struct D2 
  {
    static int const present = (row == col) && (ConstraintsCache::template bounds<row>::upper || ConstraintsCache::template bounds<row>::lower);
    static int const symmetric = true;
    static bool const lumped = VarFu::template D2<row,col>::lumped;

  };

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool  boundary ) const 
  {
    if(boundary) return 0;
    // matrix and rhs do depend on gradient and value, i.e. the polynomial
    // order of grad(arg1)*grad(arg2) + arg1*arg1 is
    int matrixOrder = shapeFunctionOrder;
    int lastIterateOrder = 0;
    return matrixOrder+lastIterateOrder;
  }

  OptimalityFunctional const& getFunctional() const {return unconstrainedFunctional;}

  Parameter mu;
  OptimalityFunctional& unconstrainedFunctional;
};

}  // namespace Kaskade
#endif
