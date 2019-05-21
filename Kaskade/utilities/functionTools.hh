#ifndef FUNCTION_TOOLS_HH
#define FUNCTION_TOOLS_HH

#include <functional> // std::declval
#include "utilities/straintensors.hh"

namespace Dune { template <class,int,int> class FieldMatrix; }

namespace Kaskade
{
  template <class Scalar, int dim>
  struct StrainBase
  {
    StrainBase() : F(0), S(F) {}

    StrainBase(StrainBase const&) = default;
    StrainBase& operator=(StrainBase const&) = default;

    void updateDisplacementGradient(Dune::FieldMatrix<Scalar,dim,dim> const& G) 
    { 
      F = G + unitMatrix<Scalar,dim>();
      S = CauchyGreenTensor<Scalar,dim>(F);
    }
    
    void updateDeformationGradient(Dune::FieldMatrix<Scalar,dim,dim> const& F_)
    { 
      F = F_;
      S = CauchyGreenTensor<Scalar,dim>(F);
    }
    
  public:
    Dune::FieldMatrix<Scalar,dim,dim> F;
    CauchyGreenTensor<Scalar,dim> S;
  };


  template <class F, class G>
  struct Product
  {
    typedef typename F::Argument Argument;

    Product(F const& f_, G const& g_) : f(f_), g(g_) {}

    Product(Product const&) = default;
    Product& operator=(Product const&) = default;

    auto d0() const { return f.d0()*g.d0(); }

    auto d1(Argument const& dF1) const { return f.d1(dF1)*g.d0() + f.d0()*g.d1(dF1); }

    auto d2(Argument const& dF1, Argument const& dF2) const
    {
      return f.d2(dF1,dF2)*g.d0() + f.d1(dF1)*g.d1(dF2) + f.d1(dF2)*g.d1(dF1) + f.d0()*g.d2(dF1,dF2);
    }

    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const
    {
      return f.d3(dF1,dF2,dF3)*g.d0() + f.d2(dF1,dF2)*g.d1(dF3) + f.d2(dF1,dF3)*g.d1(dF2) + f.d1(dF1)*g.d2(dF2,dF3)
        + f.d2(dF2,dF3)*g.d1(dF1) + f.d1(dF2)*g.d2(dF1,dF3) + f.d1(dF3)*g.d2(dF1,dF2) + f.d0()*g.d3(dF1,dF2,dF3);
    }
    
  private:
    F f;
    G g;
  };
 
  template <class F>
  struct Squared
  {
    typedef typename F::Argument Argument;

    Squared(F const& f_) : f(f_){}

    Squared(Squared const&) = default;
    Squared& operator=(Squared const&) = default;

    auto d0() const { return f.d0()*f.d0(); }

    auto d1(Argument const& dF1) const { return 2.*f.d1(dF1)*f.d0(); }

    auto d2(Argument const& dF1, Argument const& dF2) const
    {
      return 2.*(f.d2(dF1,dF2)*f.d0() + f.d1(dF1)*f.d1(dF2));
    }

    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const
    {
      return 2.*( f.d3(dF1,dF2,dF3)*f.d0() + f.d2(dF1,dF2)*f.d1(dF3) + f.d2(dF1,dF3)*f.d1(dF2) + f.d1(dF1)*f.d2(dF2,dF3) );
    }

  private:
    F f;
  };

  template <class F, class Scaling=double>
  struct Scaled
  {
  public:
    typedef typename F::Argument Argument;

    Scaled(Scaling c_, typename F::Source const& s) : c(c_), f(s) {}
    Scaled(Scaling c_, F const& f_) : c(c_), f(f_) {}
    
    Scaled(Scaled const&) = default;
    Scaled& operator=(Scaled const&) = default;
    
    auto d0() const { return c()*f.d0(); }
    
    auto d1(Argument const& dF1) const { return c()*f.d1(dF1); }
    
    auto d2(Argument const& dF1, Argument const& dF2) const { return c()*f.d2(dF1,dF2); }
    
    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return c()*f.d3(dF1,dF2,dF3); }
  
  private:
    Scaling c;
    F f;
  };
  
  template <class F>
  struct Scaled<F,double>
  {
  public:
    typedef typename F::Argument Argument;
    
//    Scaled(double c_, typename F::Source const& s) : c(c_), f(s) {}
    Scaled(double c_, F const& f_) : c(c_), f(f_) {}
    
    Scaled(Scaled const&) = default;
    Scaled& operator=(Scaled const&) = default;
    
    auto d0() const { return c*f.d0(); }
    
    auto d1(Argument const& dF1) const { return c*f.d1(dF1); }
    
    auto d2(Argument const& dF1, Argument const& dF2) const { return c*f.d2(dF1,dF2); }
    
    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return c*f.d3(dF1,dF2,dF3); }
  
  private:
    double c;
    F f;
  };
  
  
  template <typename...> class Sum;
  
  template <class F, typename... Args>
  struct Sum<F,Args...> : public Sum<F,Sum<Args...> >
  {
    Sum(F const& f, const Args&... args) : Sum<F,Sum<Args...> >(f,Sum<Args...>(args...)) {}
    
    Sum(Sum const&) = default;
    Sum& operator=(Sum const&) = default;
  };

  template <class F, class G>
  class Sum<F,G>
  {
  public:
    typedef typename F::Argument Argument;

    Sum(F const& f_, G const& g_) : f(f_), g(g_) {}

    Sum(Sum const&) = default;
    Sum& operator=(Sum const&) = default;

    auto d0() const { return f.d0() + g.d0(); }

    auto d1(Argument const& dF1) const { return f.d1(dF1) + g.d1(dF1); }

    auto d2(Argument const& dF1, Argument const& dF2) const { return f.d2(dF1,dF2) + g.d2(dF1,dF2); }

    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return f.d3(dF1,dF2,dF3) + g.d3(dF1,dF2,dF3); }

  private:
    F f;
    G g;
  };

  
  template <class Argument_=double>
  class ScalarConstant
  {
  public:    
    typedef Argument_ Argument;
    explicit ScalarConstant(double v_) : v(v_) {}
    
    ScalarConstant(ScalarConstant const&) = default;
    ScalarConstant& operator=(ScalarConstant const&) = default;
    
    double d0() const { return v; }
    double d1(Argument const&) const { return 0; }
    double d2(Argument const&,Argument const&) const { return 0; }
    double d3(Argument const&,Argument const&,Argument const&) const { return 0; }
    
  private:
    double v;
  };
  
  template <class Argument_, class Source_=double> 
  class Constant
  {
  public:
    typedef Source_ Source;
    typedef Argument_ Argument;
    
    explicit Constant(Source v_) : v(v_) {}
    
    Constant(Constant const&) = default;
    Constant& operator=(Constant const&) = default;
    
    auto d0() const { return v.d0(); }
    auto d1(Argument const&) const { return 0; }
    auto d2(Argument const&,Argument const&) const { return 0; }
    auto d3(Argument const&,Argument const&,Argument const&) const { return 0; }
    
    auto operator()() const { return d0(); }
    
  private:
    Source v;
  };
  
  template <class Argument_>
  class Constant<Argument_,double>
  {
  public:
    typedef double Source;
    typedef Argument_ Argument;
    
    explicit Constant(auto const& v_) : v(v_) {}
    
    Constant(Constant const&) = default;
    Constant& operator=(Constant const&) = default;
    
    auto d0() const { return v; }
    auto d1(Argument const&) const { return 0; }
    auto d2(Argument const&,Argument const&) const { return 0; }
    auto d3(Argument const&,Argument const&,Argument const&) const { return 0; }
    
    auto operator()() const { return d0(); }

  private:
    double v;
  };


  template <class F, class Source>
  class Sum<F,Constant<typename F::Argument,Source> >
  {
    typedef Constant<typename F::Argument,Source> G;
  public:
    typedef typename F::Argument Argument;

    Sum(F const& f_, G const& g_) : f(f_), g(g_) {}

    Sum(Sum const&) = default;
    Sum& operator=(Sum const&) = default;

    auto d0() const { return f.d0() + g.d0(); }

    auto d1(Argument const& dF1) const { return f.d1(dF1); }
      //                                             if(std::fabs(val) > 1e6) std::cout << "val: " << f.d1(dF1) << " + " << g.d1(dF1) << std::endl;
      //                                     return val;}

    auto d2(Argument const& dF1, Argument const& dF2) const { return f.d2(dF1,dF2); }

    auto d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return f.d3(dF1,dF2,dF3); }

  private:
    F f;
    G g;
  };

  template <class ScalarFunction, class MatrixFunction>
  class MatrixToScalarFunction
  {
  public:
//    typedef typename MatrixFunction::Source Source;
    typedef typename MatrixFunction::Argument Argument;

    explicit MatrixToScalarFunction(MatrixFunction const& m_) : m(m_), f() {}
        
    MatrixToScalarFunction(MatrixFunction const& m_, ScalarFunction const& f_) : m(m_), f(f_) {}
    MatrixToScalarFunction(ScalarFunction const& f_, MatrixFunction const& m_) : m(m_), f(f_) {}

    MatrixToScalarFunction(MatrixToScalarFunction const&) = default;
    MatrixToScalarFunction& operator=(MatrixToScalarFunction const&) = default;

    double d0() const { return f.d0(m.d0());  }

    double d1(Argument const& dF) const { return f.d1(m.d0())*m.d1(dF); }

    double d2(Argument const& dF1, Argument const& dF2) const { return f.d2(m.d0())*m.d1(dF1)*m.d1(dF2) + f.d1(m.d0())*m.d2(dF1,dF2); }

    double d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const
    {
      return f.d3(m.d0())*m.d1(dF1)*m.d1(dF2)*m.d1(dF3) + f.d2(m.d0())*m.d2(dF1,dF3)*m.d1(dF2) + f.d2(m.d0())*m.d1(dF1)*m.d2(dF2,dF3)
	+ f.d2(m.d0())*m.d2(dF1,dF2)*m.d1(dF3) + f.d1(m.d0())*m.d3(dF1,dF2,dF3);
    }

  private:
    MatrixFunction m;
    ScalarFunction f;
  };

  template <class Invariant>
  struct InvariantIndicator
  {
    typedef typename Invariant::Argument Argument;
    InvariantIndicator(Invariant i_, double cutoff_=1) : i(i_), cutoff(cutoff_)
    {}

    InvariantIndicator(InvariantIndicator const&) = default;
    InvariantIndicator& operator=(InvariantIndicator const&) = default;

    size_t d0() const { return ( i.d0() < cutoff ) ? 0 : 1; }

    size_t d1(Argument const&) const { return 0; }
    size_t d2(Argument const&, Argument const&) const { return 0; }
    size_t d3(Argument const&, Argument const&, Argument const&) const { return 0; }

  private:
    Invariant i;
    double cutoff;
  };

  template <class Function, class Indicator>
  class MultiplyWithIndicator
  {
  public:
    typedef typename Function::Argument Argument;

    MultiplyWithIndicator(MultiplyWithIndicator const&) = default;
    MultiplyWithIndicator& operator=(MultiplyWithIndicator const&) = default;

    MultiplyWithIndicator(Function f_, Indicator xi_) : f(f_), xi(xi_)
    {}

    double d0() const { return xi.d0() * f.d0(); }
    double d1(Argument const& dF) const { return xi.d0() * f.d1(dF); }
    double d2(Argument const& dF1, Argument const& dF2) const { return xi.d0() * f.d2(dF1,dF2); }
    double d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return xi.d0() *  f.d3(dF1,dF2,dF3); }

  private:
    Function f;
    Indicator xi;
  };

  template <size_t dim>
  struct SimpleFiberDistribution
  {
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    explicit SimpleFiberDistribution(double kappa) : principalCoefficient(kappa), mixedCoefficient(1-3*kappa) { assert ( kappa >= 0 && kappa <= 1./3. ); }

    SimpleFiberDistribution(SimpleFiberDistribution const&) = default;
    SimpleFiberDistribution& operator=(SimpleFiberDistribution const&) = default;

    Constant<Argument> getPrincipalCoefficient() const { return principalCoefficient; }

    Constant<Argument> getMixedCoefficient() const { return mixedCoefficient; }

    void setKappa(double kappa)
    {
      principalCoefficient = Constant<Argument>(kappa);
      mixedCoefficient = Constant<Argument>(1-3*kappa);
    }

  private:
    Constant<Argument> principalCoefficient, mixedCoefficient;
  };
}

#endif
