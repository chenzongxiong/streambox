#ifndef ELEMENTARY_FUNCTIONS_HH
#define ELEMENTARY_FUNCTIONS_HH

#include <cmath>

/// Monomial a*x^k of order k
template <int k>
class Monomial{
public:
  enum{ K=k };
  /// Constructor
  /**
    * \param a_ scaling factor
    */
  explicit Monomial(double const a_=1.0) : a(a_){}

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

/// Natural logarithm of the form a*\ln(x)
class LN{
public:
  /// Constructor
  /**
    * \param a_ scaling factor
    */
  explicit LN(double const a_=1.0) :
              a(a_)
  {}

  /// Function value
  /**
    * \param x evaluation point
    */
  double d0(double const x) const
  {
    return a*log(x);
  }

  /// First derivative
  /**
    * \param x evaluation point
    */
  double d1(double const x) const
  {
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

/// Exponential function of the form a*exp(b*x)
class Exponential{
public:
  /// Constructor
  /**
    * \param a_ scaling factor
    * \param b_ scaling factor for x
    */
  explicit Exponential(double const a_=1.0) : a(a_){}

  /// Function value
  /**
    * \param x evaluation point
    */
  double d0(double const x) const
  {
    return a*exp(x);
  }

  /// First derivative
  /**
    * \param x evaluation point
    */
  double d1(double const x) const
  {
    return a*exp(x);
  }

  /// Second derivative
  /**
    * \param x evaluation point
    */
  double d2(double const x) const
  {
    return a*exp(x);
  }

  /// Third derivative
  /**
    * \param x evaluation point
    */
  double d3(double const x) const
  {
    return a*exp(x);
  }

  /// Set a
  void setA(double const a_) { a=a_; }

private:
  double a;
};


#endif // ELEMENTARY_FUNCTIONS_HH
