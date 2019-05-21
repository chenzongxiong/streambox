#ifndef PENALTY_FUNCTIONS_HH
#define PENALTY_FUNCTIONS_HH

#include <cmath>

#include "utilities/elementary_functions.hh"

namespace Kaskade
{
  class QuadraticAndLogBase
  {
  public:
    QuadraticAndLogBase(double c_, double d_)
    {
      init(c_,d_);
    }
    
    QuadraticAndLogBase() : QuadraticAndLogBase(0.,1.) {}
    
    QuadraticAndLogBase(QuadraticAndLogBase const&) = default;
    QuadraticAndLogBase& operator=(QuadraticAndLogBase const&) = default;
    
    double d0(double t) const { return c*t*t + d*log(t);}// + e; }
    double d1(double t) const { return 2*c*t + d/t; }
    double d2(double t) const { return 2*c - d/(t*t); }
    double d3(double t) const { return 2*d/(t*t*t); }
    
  private:
void init(double c_, double d_)
    {
      //      assert(c_>0);
      c = c_;
      d = d_;
      e = -d0(1);
    }

    double c = 1, d = 1, e = 0;
  };
  


  template <size_t k1, int k2>
  class Polynomial
  {
//    static_assert(k2<0, "power not admissible");
  public:
    explicit Polynomial(double c_, double d_) : c(c_), d(d_)
    { 
      assert(c>0); 
      e = -d0(1);
    }

    Polynomial() : Polynomial(1.,1.) {}
    Polynomial(Polynomial const&) = default;
    Polynomial& operator=(Polynomial const&) = default;
    
    double d0(double t) const { return c*pow(t,k1) + d*pow(t,k2) + e; }
    double d1(double t) const { return c*k1*pow(t,k1-1) + d*k2*pow(t,k2-1); }
    double d2(double t) const { return c*k1*(k1-1)*pow(t,k1-2) + d*k2*(k2-1)*pow(t,k2-2); }
    double d3(double t) const { return c*k1*(k1-1)*(k1-2)*pow(t,k1-3) + d*k2*(k2-1)*(k2-2)*pow(t,k2-3); }
    
  private:
    double c, d, e = 0;
  };
  
  typedef PenaltyFunction<Polynomial<5,-5> > Hartmann03;
  typedef PenaltyFunction<Polynomial<2,-2> > QuadQuad;
  typedef PenaltyFunction<QuadraticAndLogBase> QuadraticAndLog;

  template <int k> using QuadraticAndPolynomial = PenaltyFunction<Polynomial<2,k> >;
}

#endif
