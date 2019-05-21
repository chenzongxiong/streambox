/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TCG_HH
#define TCG_HH

#include <numeric>

#include <boost/circular_buffer.hpp>

#include "dune/common/timer.hh"
#include "dune/istl/istlexception.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"
#include "dune/istl/solvers.hh"

#include "utilities/geometric_sequence.hh"
#include "utilities/scalar.hh"
#include "utilities/duneInterface.hh"

#include "linalg/apcg.hh"

namespace Kaskade
{
  struct NoRegularization
  {
    template <class X, class Xstar>
    void operator()(const X&, Xstar&){}
  };

  template <class Linearization, class VariableSet>
  struct DoRegularization
  {
    DoRegularization(Linearization const& lin_, typename VariableSet::Descriptions const& description) : lin(lin_), v(description)
    {}

    template <class X, class Xstar>
    void operator()(const X& x, Xstar& r)
    {
      v = x; //
      Bridge::Vector<VariableSet> y(v), z(v);
      y *= 0;
      lin.ddxpy(y,z,0,2,2,3);
      y *= -1;
      v = r; //

      boost::fusion::at_c<0>(v.data) += boost::fusion::at_c<0>(y.get().data);
      boost::fusion::at_c<1>(v.data) += boost::fusion::at_c<1>(y.get().data);

      r = v; //
    }

    private:
      Linearization const& lin;
      VariableSet v;
  };

  /**
   * \ingroup linalgsolution
   * \brief preconditioned conjugate gradient method
   *
   * This implements a preconditioned IterateType::CG iteration for an operator \f$ A: X\to x^* \f$, preconditioned by a
   * preconditioner \f$ B^{-1}: X^* \to X \f$. The termination is based on an estimate of the absolute energy error.
   *
   * The implementation follows Deuflhard/Weiser, Section 5.3.3.
   *
   */
  template<class X, class Xstar, class Regularization=NoRegularization>
  class TCG : public Dune::InverseOperator<X,Xstar> {
  public:
    enum class Result { Converged, Failed, EncounteredNonConvexity };
    /**
     * \brief the real field type corresponding to X::field_type
     */
    typedef typename ScalarTraits<typename GetScalar<X>::type >::Real Real;

    /** 
     * \brief Set up conjugate gradient solver with absolute energy error termination criterion.
     * 
     * \param verbose
     */
    template <class Int=int, class enable = typename std::enable_if<!std::is_same<Regularization,NoRegularization>::value && std::is_same<Int,int>::value>::type>
    TCG(Dune::LinearOperator<X,Xstar>& op_, Dune::Preconditioner<X,Xstar>& prec_, DualPairing<X,Xstar> const& dp_,
        Regularization& regularization_, double relTol_=1e-3, size_t maxSteps_=100, Int verbose_=0) :
          op(op_), prec(prec_), dp(dp_), regularization(regularization_), relTol(relTol_), absTol(1e-9), maxSteps(maxSteps_), verbose(verbose_)
    {
      // Do we need this in Kaskade7?
      //         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
      //                             "L and P must have the same category!");
      //         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
      //                             "L must be sequential!");
    }

    /**
     * \brief Set up conjugate gradient solver with absolute energy error termination criterion.
     *
     * \param verbose
     */
    template <class Int=int, class enable = typename std::enable_if<std::is_same<Regularization,NoRegularization>::value && std::is_same<Int,int>::value>::type>
    TCG(Dune::LinearOperator<X,Xstar>& op_, Dune::Preconditioner<X,Xstar>& prec_, DualPairing<X,Xstar> const& dp_,
        double relTol_=1e-3, size_t maxSteps_=100, Int verbose_=0) :
          op(op_), prec(prec_), dp(dp_), relTol(relTol_), absTol(1e-9), maxSteps(maxSteps_), verbose(verbose_)
    {
      // Do we need this in Kaskade7?
      //         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
      //                             "L and P must have the same category!");
      //         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
      //                             "L must be sequential!");
    }

    /**
     * \brief Apply inverse operator.
     * 
     * \param u the initial value (starting iterate)
     * \param b the right hand side
     */
    virtual int apply (X& u, X& q, Xstar& b, Dune::InverseOperatorResult& res) 
    {
      result = Result::Failed;
      int numberOfResults = 1;
      res.clear();                // clear solver statistics
      //terminate.clear();          // clear termination criterion
      Dune::Timer watch;          // start a timer

      prec.pre(u,b);

      Xstar r(b); 
      op.applyscaleadd(-1.0,u,r);  // r = b-Au
      regularization(u,r);

      X rq(u), du(u); rq = 0;

      prec.apply(rq,r); // rq = B^{-1} r

      q = rq;

      X Aq(b);


      // some local variables
      Real alpha,beta,sigma,gamma;
      sigma = dp(rq,r); // preconditioned residual norm squared

      //terminate.residual(sigma);
      double const sigma0 = std::abs(sigma);

      if (verbose>0) {            // printing
        std::cout << "=== Kaskade7 TCG" << std::endl;
        if (verbose>1) {
          this->printHeader(std::cout);
          this->printOutput(std::cout,0,sigma0);
        }
      }

      // the loop
      int i=0;
      while(true)
      {
        ++i;
        // minimize in given search direction p
        op.apply(q,Aq);             // h = Aq
        Real qAq = dp(q,Aq);

        if(qAq <=0)
        {
          result = Result::EncounteredNonConvexity;
          if(verbose > 0) std::cout << "Nonconvexity at iteration " << i << ", qAq=" << qAq << ", ||q||=" << sqrt(q*q) << std::endl;
          if(i==1)
          {
            u.axpy(1.0,q);
            break;
          }
          numberOfResults = 2;
          break;      // Abbruch wg. fehlender Konvexität
        }

        alpha = sigma/qAq;
        u.axpy(alpha,q);
        du = q;
        du *= 1./alpha;
        if(dp(du,du) < absTol)
        {
          result = Result::Converged;
          break;
        }

        gamma = sigma*alpha;
        //terminate.step(ScalarTraits<Real>::real(gamma));
        resDebug.push_back(std::sqrt(sigma));

        // convergence test
        //if (terminate)
          //break;

        r.axpy(-alpha,Aq); // r = r - alpha*A*q
        regularization(u,r);
        rq = 0;
        prec.apply(rq,r); // rq = B^{-1}r

        // determine new search direction
        double sigmaNew = dp(rq,r); 
        //terminate.residual(ScalarTraits<Real>::real(sigmaNew));
        if (verbose>1)             // print
          this->printOutput(std::cout,i,std::abs(sigmaNew),std::abs(sigma));

        // convergence check to prevent division by zero if we obtained an exact solution
        // (which may happen for low dimensional systems)
        if (std::abs(sigmaNew) < relTol*sigma0)
        {
          result = Result::Converged;
          break;
        }
        beta = sigmaNew/sigma;
        if(verbose > 0) std::cout << "step reduction: " << (sigmaNew/sigma) << ", overall reduction: " << (sigmaNew/sigma0) << std::endl;
        sigma = sigmaNew;
        q *= beta;                  // scale old search direction
        q += rq;                     // orthogonalization with correction
        
        if(i > maxSteps) break;
      }

      res.iterations = i;               // fill statistics
      res.reduction = std::sqrt(std::abs(sigma)/sigma0);
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();

      if (verbose>0)                 // final print 
      {
        std::cout << "=== rate=" << res.conv_rate
            << ", time=" << res.elapsed
            << ", iterations=" << i << std::endl;
      } 

      prec.post(u);

      return numberOfResults;
    }

    /**
     * \brief Apply tcg and possibly forget second descent direction
     */
    virtual void apply (X& u, Xstar& b, Dune::InverseOperatorResult& res) {
      X q(u);
      apply(u,q,b,res);
    }

    void apply (X& u, Xstar& b)
    {
      Dune::InverseOperatorResult res;
      apply(u,b,res);
    }

    /** 
     * \brief Apply inverse operator with given absolute tolerance.
     */
    virtual void apply (X& x, X& b, double relTol, Dune::InverseOperatorResult& res) {
      //terminate.relTol(relTol);
      (*this).apply(x,b,res);
    }

    void setRelativeAccuracy(double relTol_) { relTol = relTol_; }

    void setMaxSteps(size_t maxSteps_) { maxSteps = maxSteps_; }
    
    bool localConvergenceLikely() const { return result == Result::Converged; }

    bool encounteredNonConvexity() const { return result == Result::EncounteredNonConvexity; }

  private:
    Dune::LinearOperator<X,Xstar>& op;
    Dune::Preconditioner<X,Xstar>& prec;
    DualPairing<X,Xstar> const& dp;
    typename std::conditional<std::is_same<Regularization,NoRegularization>::value,NoRegularization,Regularization&>::type regularization;
    double relTol, absTol;
    size_t maxSteps;
    int verbose;
    std::vector<double> resDebug;
    Result result;
  };

  template <class X, class Xstar, class Regularization = NoRegularization> using TPCG = TCG<X,Xstar,Regularization>;
} // namespace Kaskade
#endif
