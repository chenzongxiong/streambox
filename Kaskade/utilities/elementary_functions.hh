#ifndef ELEMENTARY_FUNCTIONS_HH
#define ELEMENTARY_FUNCTIONS_HH

#include <cmath>
#include <limits>

namespace Kaskade
{
    /// Sets a function to std::numeric_limits<double>::max() for negative arguments.
    template <class Base>
    class PenaltyFunction : public Base
    {
    public:
      template <typename... Args>
      PenaltyFunction(const Args&... args) : Base(args...){}
      
      PenaltyFunction(PenaltyFunction const&) = default;
      PenaltyFunction& operator=(PenaltyFunction const&) = default;

      double d0(double const& x) const
      {
	if(x < 1e-16) std::cout << "negative value x in penalty function d0: " << x << std::endl;
        if(x < 1e-16) return std::numeric_limits<double>::max();
        else return Base::d0(x);
      }
      double d1(double const& x) const
      {
	if(x < 1e-16) std::cout << "negative value x in penalty function d1: " << x << std::endl;
        if(x < 1e-16) return std::numeric_limits<double>::max();
        else return Base::d1(x);
      }
      double d2(double const& x) const
      {
	if(x < 1e-16) std::cout << "negative value x in penalty function d2: " << x << std::endl;
        if(x < 1e-16) return std::numeric_limits<double>::max();
        else return Base::d2(x);
      }
      double d3(double const& x) const
      {
	if(x < 1e-16) std::cout << "negative value x in penalty function d3: " << x << std::endl;
        if(x < 1e-16) return std::numeric_limits<double>::max();
        else return Base::d3(x);
      }
    };

  /// Monomial a*x^k of order k
  template <int k>
  class Monomial{
  public:
    enum{ K=k };
    /// Constructor
    /**
     * \param a_ scaling factor
     */
    explicit Monomial(double const a_=1.0) : a(a_)
    {
#ifdef TESTOUTPUT
      std::cout << "Created monomial of order " << k << std::endl;
#endif
    }

    /// Function value
    /**
     * \param x evaluation point
     */
    double d0(double const x) const
    {
      return a*pow(x,k);
    }

    /// First derivative
    /**
     * \param x evaluation point
     */
    double d1(double const x) const
    {
      return a*k*pow(x,k-1);
    }

    /// Second derivative
    /**
     * \param x evaluation point
     */
    double d2(double const x) const
    {
      return a*k*(k-1)*pow(x,k-2);
    }

    /// Third derivative
    /**
     * \param x evaluation point
     */
    double d3(double const x) const
    {
      return a*k*(k-1)*(k-2)*pow(x,k-3);
    }

    /// Set a
    void setA(double const a_) { a=a_; }

  private:
    double a;
  };

  template <int k> using PenaltyMonomial = PenaltyFunction<Monomial<k> >;

  /// Natural logarithm of the form a*\ln(x)
  class LN{
  public:
    /// Constructor
    /**
     * \param a_ scaling factor
     */
    explicit LN(double const a_=1.0) :
    a(a_)
    {
#ifdef TESTOUTPUT
      std::cout << "Created natural logarithm" << std::endl;
#endif
    }

    /// Function value
    /**
     * \param x evaluation point
     */
    double d0(double const x) const
    {
      if(x<0) std::cout << "Warning: negative argument for log" << std::endl;
      return a*log(x);
    }

    /// First derivative
    /**
     * \param x evaluation point
     */
    double d1(double const x) const
    {
      if(x<0) std::cout << "Warning: negative argument for a/x" << std::endl;
      return a/x;
    }

    /// Second derivative
    /**
     * \param x evaluation point
     */
    double d2(double const x) const
    {
      return -a/(x*x);
    }

    /// Third derivative
    /**
     * \param x evaluation point
     */
    double d3(double const x) const
    {
      return 2.0*a/(x*x*x);
    }

    /// Set a
    void setA(double const a_) { a=a_; }

  private:
    double a;
  };

  typedef PenaltyFunction<LN> PenaltyLN;

  /// Exponential function of the form a*exp(b*x)
  class Exp
  {
  public:
    /// Constructor
    /**
     * \param a_ scaling factor
     * \param b_ scaling factor for x
     */
    explicit Exp(double a_=1.0) : a(a_){}

    /// Function value
    /**
     * \param x evaluation point
     */
    double d0(double x) const
    {
      return a*exp(x);
    }

    /// First derivative
    /**
     * \param x evaluation point
     */
    double d1(double x) const
    {
      return a*exp(x);
    }

    /// Second derivative
    /**
     * \param x evaluation point
     */
    double d2(double x) const
    {
      return a*exp(x);
    }

    /// Third derivative
    /**
     * \param x evaluation point
     */
    double d3(double x) const
    {
      return a*exp(x);
    }

    /// Set a
    void setA(double a_) { a=a_; }

  private:
    double a;
  };

  class Sqrt
  {
  public:

    double d0(double x) const noexcept
    {
      return sqrt(x);
    }

    double d1(double x) const noexcept
    {
      return 1./(2*sqrt(x));
    }

    double d2(double x) const noexcept
    {
      return -1./(4*pow(x,1.5));
    }

    double d3(double x) const noexcept
    {
      return 3./(8*pow(x,2.5));
    }

  private:
  };
}

#endif // ELEMENTARY_FUNCTIONS_HH
