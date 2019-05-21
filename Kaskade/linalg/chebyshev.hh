/*
 * chebyshev.hh
 *
 * \brief Chebyshev semi-iteration.
 *
 *  Created on: 14.12.2013
 *      Author: Lars Lubkoll
 */

#ifndef CHEBYSHEV_HH_
#define CHEBYSHEV_HH_

#include <cstdlib>
#include <iostream>

#include "linalg/jacobiPreconditioner.hh"

namespace Kaskade
{
  /// \cond internals
  namespace ChebyshevDetail
  {
    double chebyshev(double x, size_t m)
    {
      if(m==0) return 1;
      if(m==1) return x;

      double n0 = 1, n1 = x;

      for(size_t i=2; i<m+1; ++i)
      {
        double n = 2*x*n1 - n0;
        n0 = n1;
        n1 = n;
      }
      return n1;
    }
  }
  /// \endcond

  /**
   * \brief Preconditioned chebyshev semi iteration.
   *
   * Standard implementation based on a three-term recurrence and explicit computation of the residuals to avoid
   * accumulation of round of errors in the computation of the residuals. In contrast to Krylov methods this does not
   * slow down convergence (see Gutknecht,Röllin: The Chebyshev iteration revisited).
   *
   * The Chebyshev semi-iteration requires bounds on the spectrum! Currently we only provide these bounds for
   * tetrahedral elements (3D!) with piecewise linear ansatz functions. For more bounds see a.o. "Wathen: Realistic eigenvalue bounds for the Galerkin mass matrix the spectrum of the preconditioned mass matrix".
   * For reason of consistency with the published bounds the implementation is restricted to a one-step Jacobi preconditioner.
   */
  template <class Operator, bool isPreconditioner=false>
  class ChebyshevSemiIteration : public Dune::Preconditioner<typename Operator::Domain,typename Operator::Range>
  {
  public:
    typedef typename Operator::Domain Domain;
    typedef typename Operator::Range Range;

    ChebyshevSemiIteration(Operator const& A_, double tol, size_t steps, size_t verbose_=0)
    : A(A_), jac(A), a(0), b(0), eta(0), tolerance(tol), maxSteps(steps), verbose(verbose_), isInitialized(false)
    {}

    virtual ~ChebyshevSemiIteration(){}

    virtual void pre(Domain&, Range&){}
    virtual void post(Domain&){}

    virtual void apply(Domain& x, Range const& y)
    {
      if(!isInitialized)
      {
        std::cerr << "In file " << __FILE__ << ": line " << __LINE__;
        std::cerr << ": No spectral bounds have been provided for chebyshev semi-iteration. Aborting." << std::endl;
        abort();
      }

      double beta = 0;
      double gammainv = 0;
      Domain x0(x), x1(x); x1 *= 0;
      Range r(y), jr(y), tmp(y), res0(y), res1(y); res1 *= 0;
      A.applyscaleadd(-1.,x,tmp);
      jac.apply(res0,tmp);
      double initialResidual = res0*tmp, currentResidual = 0, desiredResidual = tolerance*tolerance*initialResidual;
      if(verbose > 1) std::cout << "initial residual: " << initialResidual << std::endl;
      for(size_t i=0; i<maxSteps; ++i)
      {
        if(i==0)
        {
          beta = 0;
          gammainv = -1./a;
        }
        else
        {
          if(i==1) beta = -0.5*b*b/a;
          else beta = 0.25*b*b*gammainv;
          gammainv = -1./(a+beta);
        }


        x*=0;
        x.axpy(-1.,jr);
        x.axpy(-a,x0);
        x.axpy(-beta,x1);
        x *= gammainv;

        x1 = x0;
        x0 = x;

        r = y;
        A.applyscaleadd(-1.,x,r);
        jac.apply(jr,r);

        currentResidual = r*jr;
        if((currentResidual < desiredResidual) && (isPreconditioner==false))
        {
          if(verbose>0)
          {
            std::cout << "Terminating chebyshev semi iteration after " << i << " iterations with final residual: " << currentResidual << "." << std::endl;
            std::cout << "Residual reduction: " << currentResidual/initialResidual << " (initialResidual=" << initialResidual << ")" << std::endl;
          }
          return;
        }
      }
    }

    /**
     * \brief Init spectral bounds for the mass matrix arising from tetrahedral discretization of the domain and linear elements.
     *
     * According to "Wathen: Realistic eigenvalue bounds for the Galerkin mass matrix" the spectrum of the preconditioned mass matrix is contained in \f$[\frac{1}{2},\frac{5}{2}]\f$.
     * For chebyshev semi-iteration the bounds on the spectrum are in general given in the form \f$ [a-b,a+b] \f$ yielding the coefficients \f$a=1.5,\ b=1\f$.
     * See Wathen,Rees: "Chebyshev semi-iteration in preconditioning for problems including the mass matrix".
     */
    void initForMassMatrix_TetrahedralQ1Elements()
    {
      // parameters for describing the spectral bounds
      // i.e. \f$\lambda\in[\lambda_{min},\lambda_{max}]\f$ with \f$\lambda_{min}=a-b\f$ and \f$\lambda_{max}=a+b\f$
      a = 1.5; // center of spectral bounds
      b = 1;  // half diameter of the spectral bounds
      eta = -a/b;
      isInitialized = true;
    }

    void setSteps(size_t steps) { maxSteps = steps; }

  private:
    Operator const& A;
    JacobiPreconditioner<Operator> jac;
    double a, b, eta, tolerance;
    size_t maxSteps, verbose;
    bool isInitialized;
  };

  /**
   * \brief Preconditioner based on the Chebyshev semi-iteration. In order to provide a linear preconditioner the termination
   * criterion based on relative descent in the energy-norm induced by the preconditioner is disabled here. Instead the preconditioner
   * performs a fixed number of iterations thus yielding a linear preconditioner.
   */
  template <class Operator>
  class ChebyshevPreconditioner : public ChebyshevSemiIteration<Operator,true>
  {
  public:
    ChebyshevPreconditioner(Operator const& A, size_t steps, size_t verbose=0) : ChebyshevSemiIteration<Operator,true>(A,0,steps,verbose) {}
  };
}

#endif /* CHEBYSHEV_HH_ */
