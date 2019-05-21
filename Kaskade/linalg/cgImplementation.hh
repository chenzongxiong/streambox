/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CG_IMPLEMENTATION_HH
#define CG_IMPLEMENTATION_HH

#include <iostream>
#include <memory>

#include "algorithm/opt_aux/src/include/Fmin.h"

#include "boost/timer/timer.hpp"
#include "utilities/scalar.hh"
#include "linalg/cgTerminationCriteria.hh"

namespace Kaskade
{
  enum class CGImplementationType { STANDARD, TRUNCATED, REGULARIZED, HYBRID };

  namespace CG_Detail
  {
    template <class X, class Xstar, class Derived>
    struct Regularization
    {
      typedef typename ScalarTraits<typename GetScalar<X>::type>::Real Real;

      Regularization(double eps_, int verbose_) : eps(eps_), verbose(verbose_)
      {}

      void initializeRegularization()
      {
        theta = 0;
      }

      void regularize(double& qAq, double qPq)
      {
        qAq += theta*qPq;
      }

      void updateRegularization(double qAq, double qPq)
      {
        double oldTheta = theta > 0 ? theta : eps;
        theta += (1-qAq)/std::fabs(qPq);
        if( verbose > 0 ) std::cout << "Regularization: Computed theta: " << theta << std::endl;
        theta = std::min(std::max(2*oldTheta,theta),1e2*oldTheta);
        if( verbose > 0 ) std::cout << "Regularization: Updating theta from " << oldTheta << " to " << theta << std::endl;
      }

      void adjustResidual(Xstar& r, double alpha, Xstar const& Pq)
      {
        r.axpy(-alpha*theta,Pq);
      }

      Real theta = 0, eps = 1e-16;
      int verbose = 0;
    };

    template <class X, class Xstar, class>
    struct NoRegularization
    {
      NoRegularization(double,int) {}

      void initializeRegularization() {}
      void regularize(double&, double) {}
      void updateRegularization(double,double) {}
      void adjustResidual(Xstar&,double,Xstar const&) {}

    };

    template <class X, class Xstar, class Derived, CGImplementationType impl>
    struct ChooseRegularization
    {
      typedef typename std::conditional< (impl == CGImplementationType::STANDARD || impl == CGImplementationType::TRUNCATED ),
      NoRegularization<X,Xstar,Derived>,
      Regularization<X,Xstar,Derived>
      >::type type;
    };
  }


  struct DoNotMeasureTime
  {
    void resumeMVTimer() {}
    void resumeDPTimer() {}
    void resumeAPTimer() {}
    void resumePRTimer() {}

    void stopMVTimer() {}
    void stopDPTimer() {}
    void stopAPTimer() {}
    void stopPRTimer() {}

    void printTimers() {}
  };

  struct MeasureTime
  {
    MeasureTime()
    {
      mvTimer.stop();
      dpTimer.stop();
      apTimer.stop();
      prTimer.stop();
    }

    void resumeMVTimer() { mvTimer.resume(); }
    void resumeDPTimer() { dpTimer.resume(); }
    void resumeAPTimer() { apTimer.resume(); }
    void resumePRTimer() { prTimer.resume(); }

    void stopMVTimer() { mvTimer.stop(); }
    void stopDPTimer() { dpTimer.stop(); }
    void stopAPTimer() { apTimer.stop(); }
    void stopPRTimer() { prTimer.stop(); }

    void printTimers()
    {

      std::cout << "Matrix-Vector time: " << mvTimer.format()
                << "dual pairings     : " << dpTimer.format()
                << "vector ops        : " << apTimer.format()
                << "preconditioner    : " << prTimer.format() << '\n';
    }

  private:
    boost::timer::cpu_timer mvTimer;
    boost::timer::cpu_timer dpTimer;
    boost::timer::cpu_timer apTimer;
    boost::timer::cpu_timer prTimer;
  };

  struct DoNotAlwaysRegularize
  {
    template <typename... Args>
    void init(const Args&...){}

    double operator()() { return 1; }
  };

  struct CGFminModel
  {
    CGFminModel() = default;
    CGFminModel(CGFminModel const&) = default;
    CGFminModel& operator=(CGFminModel const&) = default;

    CGFminModel(double quadraticPart_, double linearPart_, double omegaL_,
                double dndn_, double dnd_, double dd_)
      : linearPart(linearPart_), quadraticPart(quadraticPart_), omegaL(omegaL_),
        dndn(dndn_), dnd(dnd_), dd(dd_)
    {}

    double operator()(double t)
    {
      double normSquared = dndn + t*dnd + t*t*dd;
      return t*linearPart + 0.5*t*t*quadraticPart + pow(normSquared,1.5);
    }

    double linearPart = 0, quadraticPart = 0, omegaL = 0, dndn = 0, dnd = 0, dd = 0;
  };

  /**
   * In test phase. Unclear if it makes sense. This class was motivated by the idea to always regularize possibly non-convex problems in such a way that
   * the regularization automatically goes to zero close to the solution. Not available from "cg.hh". If you want to play around with it use CGBase.
   */
  template <class POperator, class Vector>
  struct AlwaysRegularizeWithThetaGuess
  {
    AlwaysRegularizeWithThetaGuess(POperator const& P_, Vector const& dn_, double nu_, double omegaL_) : P(P_), dn(dn_), omegaL(omegaL_), nu(nu_)
    {}

    template <class Operator, class Rhs, class DualPairing>
    void init(Operator const& A, Rhs const& rhs, Vector const& d, DualPairing const& dp)
    {
      Rhs Ad(rhs), Adn(rhs);
      A.apply(d,Ad);
      A.apply(dn,Adn);
      double dAd = dp(d,Ad);
      double lin = dp(d,rhs) + nu*dp(d,Adn);

      Rhs Pd(rhs), Pdn(rhs);
      P.apply(d,Pd);
      P.apply(dn,Pdn);
      double dndn = nu*nu*dp(dn,dn),
             dnd = nu*dp(d,dn),
             dd = dp(d,d);


      model = CGFminModel(dAd,lin,omegaL,dndn,dnd,dd);
    }

    double operator()()
    {
      return Fmin(0.,1.,model,1e-3);
    }

    POperator const& P;
    Vector const& dn;
    double omegaL, nu;
    CGFminModel model;
  };

  /**
   * \ingroup linalgsolution
   * \brief regularized preconditioned conjugate gradient method
   *
   * This implements a preconditioned IterateType::CG iteration for an operator \f$ A: X\to x^* \f$, preconditioned by a
   * preconditioner \f$ B^{-1}: X^* \to X \f$. The termination is based on an estimate of the absolute energy error.
   *
   * The implementation follows Deuflhard/Weiser, Section 5.3.3.
   *
   */
  template<class X, class Xstar, CGImplementationType impl = CGImplementationType::STANDARD, class TimerPolicy = DoNotMeasureTime, class Functor = DoNotAlwaysRegularize>
  class CGBase : public Dune::InverseOperator<X,Xstar>, public TimerPolicy, public CG_Detail::ChooseRegularization<X,Xstar,CGBase<X,Xstar,impl,TimerPolicy>,impl>::type
  {
    enum class Result { Converged, Failed, EncounteredNonConvexity, TruncatedAtNonConvexity };
  public:
    /**
     * \brief the real field type corresponding to X::field_type or X::Scalar
     */
    typedef typename ScalarTraits<typename GetScalar<X>::type >::Real Real;

    /**
     * \brief Set up conjugate gradient solver with termination criterion.
     *
     * \param verbose
     */
    CGBase(Dune::LinearOperator<X,Xstar>& op_, Dune::Preconditioner<X,Xstar>& prec_, DualPairing<X,Xstar> const& dp_, PCGTerminationCriterion<Real>& terminate_, int verbose_ = 0, double eps_ = 1e-15)
      : CG_Detail::ChooseRegularization<X,Xstar,CGBase<X,Xstar,impl,TimerPolicy>,impl>::type(eps_,verbose_),
        op(op_), prec(prec_), dp(dp_), terminate(terminate_), verbose(verbose_), eps(eps_)
    {}

    CGBase(Dune::LinearOperator<X,Xstar>& op_, Dune::Preconditioner<X,Xstar>& prec_, DualPairing<X,Xstar> const& dp_, PCGTerminationCriterion<Real>& terminate_, Functor f_, int verbose_ = 0, double eps_ = 1e-15)
      : CG_Detail::ChooseRegularization<X,Xstar,CGBase<X,Xstar,impl,TimerPolicy>,impl>::type(eps_,verbose_),
        op(op_), prec(prec_), dp(dp_), terminate(terminate_), verbose(verbose_), eps(eps_), f(f_)
    {}


    void apply (X& u, Xstar& b)
    {
      Dune::InverseOperatorResult res;
      apply(u,b,res);
    }

    /**
     * @param u initial guess
     * @param b right hand side
     * @param tolerance tolerance of the termination criterion
     * @param res contains, number of iterations, residual reduction, ...
     */
    virtual void apply(X& u, Xstar&b, Real tolerance, Dune::InverseOperatorResult& res)
    {
      terminate.setTolerance(tolerance);
      apply(u,b,res);
    }

    virtual void apply(X& u, Xstar& b, Dune::InverseOperatorResult& res)
    {
      int step = 0;
      if( impl == CGImplementationType::STANDARD || impl == CGImplementationType::TRUNCATED ) step = cgLoop(u,b,res);
      else
      {
        X u0(u);
        Xstar b0(b);
        result = Result::Failed;
        while( result != Result::Converged && result != Result::TruncatedAtNonConvexity )
        {
          u = u0;
          b = b0;
          step = cgLoop(u,b,res);
        }
      }

      if (verbose>0)                 // final print
      {
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/step
                  << ", IT=" << step << std::endl;

        this->printTimers();
      }
    }

    int cgLoop (X& u, Xstar& b, Dune::InverseOperatorResult& res)
    {
      Dune::Timer watch;
      terminate.clear();          // clear termination criterion

      result = Result::Failed;
      res.converged = false;
      prec.pre(u,b);

      int step = 0;
      const int maxNoiterations=terminate.getMaxIterationSteps();

      Xstar r(b);
      op.applyscaleadd(-1,u,r);  // r = b-Au
      std::unique_ptr<X> rq(new X(u)); *rq = 0;

      this->resumePRTimer();
      prec.apply(*rq,r); // rq = B^{-1} r
      this->stopPRTimer();
      std::unique_ptr<X> q(new X(*rq));
      X Aq(b);

      Xstar Pq(r);

      // some local variables
      Real alpha = 0, beta = 0, energyNorm2 = 0;
      Real sigma = std::abs(dp(*rq,r)); // preconditioned residual norm squared
      double const sigma0 = sigma;

      if (verbose>0) {            // printing
        std::cout << "=== Kaskade7 IterateType::CG" << std::endl;
        if (verbose>0) {
          this->printHeader(std::cout);
          this->printOutput(std::cout,(double)0,sigma0);
        }
      }

      // monitor proposal for initial theta
      if(!std::is_same<Functor,DoNotAlwaysRegularize>::value)
      {
        f.init(op,b,*q,dp);
        double tau = f();
        Xstar Hd(b);
        op.apply(*q,Hd);
        double dHd = dp(*q,Hd);

        double initialThetaGuess = 1./tau - dHd/sigma;
        std::cout << "CG: TAU = " << tau << std::endl;
        std::cout << "CG: INITIAL GUESS FOR THETA = " << initialThetaGuess << std::endl;
      }
         // (sigma/tau - dHd)/dPd;
      
      // the loop
      for ( step = 1; step<=maxNoiterations; step++ )
      {
        // minimize in given search direction p
        this->resumeMVTimer();
        op.apply(*q,Aq);
        this->stopMVTimer();

        this->resumeDPTimer();
        Real qAq = dp(*q,Aq);
        Real qPq = dp(*q,Pq);
        this->stopDPTimer();
        this->regularize(qAq,qPq);

        alpha = sigma/qAq;
        energyNorm2 += alpha*qAq;

        //  don't trust small numbers
        if( std::fabs(qAq) < eps*eps )
        {
          if( verbose > 0 ) std::cout << pre << "Terminating due to numerically almost vanishing step." << std::endl;
          result = Result::Converged;
          res.converged = true;
          break;
        }

        if( qAq < 0 )
        {
          if( impl == CGImplementationType::STANDARD )
          {
            std::cout << pre << "Direction of negative curvature encountered in standard IterateType::CG implementation!" << std::endl;
            std::cout << pre << "Either something is wrong with your operator or you should use TCG, RCG or HCG. Terminating IterateType::CG!" << std::endl;
//            result = Result::EncounteredNonConvexity;
//            break;
          }

          if( impl == CGImplementationType::TRUNCATED || ( impl == CGImplementationType::HYBRID && terminate.minimalDecreaseAchieved() ) )
          {
            // At least do something to retain a little chance to get out of the nonconvexity. In fact if a nonconvexity is encountered in the first step something probably went wrong
            // elsewhere. Chances that a way out of the nonconvexity can be found are small in this case.
            if( step == 1 ) u.axpy(1.0,*q);
            if( verbose > 0 ) std::cout << pre << "Truncating at nonconvexity in iteration " << step << ": " << qAq << std::endl;
            result = Result::TruncatedAtNonConvexity;
            break;
          }

          if( impl == CGImplementationType::HYBRID || impl == CGImplementationType::REGULARIZED )
          {
            this->updateRegularization(qAq,qPq);
            if( verbose > 0 ) std::cout << pre << "Regularizing at nonconvexity in iteration " << step << "." << std::endl;
            result = Result::EncounteredNonConvexity;
            break;
          }
        }

        terminate.addStepQuantities(alpha,qAq,qPq,sigma);

        this->resumeAPTimer();
        u.axpy(alpha,*q);
        this->stopAPTimer();

        // convergence test
        if (terminate)
        {
          result = Result::Converged;
          res.converged = true;
          break;
        }

        this->resumeAPTimer();
        r.axpy(-alpha,Aq); // r = r - alpha*A*q
        this->adjustResidual(r,alpha,Pq);
        this->stopAPTimer();

        this->resumePRTimer();
        prec.apply(*rq,r);
        this->stopPRTimer();

        // determine new search direction
        this->resumeDPTimer();
        double sigmaNew = std::abs(dp(*rq,r));
        this->stopDPTimer();
        if (verbose>1) this->printOutput(std::cout,(double)step,std::abs(sigmaNew),std::abs(sigma));

        beta = sigmaNew/sigma;
        sigma = sigmaNew;
        this->resumeAPTimer();
        rq->axpy(beta,*q); // compute q = rq + beta*q
        this->stopAPTimer();
        Pq *= beta;
        Pq += r;
        std::swap(rq,q);
      }

      prec.post(u);

      if (verbose>0) this->printOutput(std::cout,(double)step,std::abs(sigma));
      res.iterations = step;               // fill statistics
      res.reduction = std::sqrt(std::abs(sigma)/sigma0);
      res.conv_rate  = pow(res.reduction,1.0/step);
      res.elapsed = watch.elapsed();

      return step;
    }

    bool encounteredNonConvexity() const { return result==Result::EncounteredNonConvexity || result==Result::TruncatedAtNonConvexity; }

    double getSolutionEnergyNormEstimate() { return sqrt(energyNorm2); }

  private:
    Dune::LinearOperator<X,Xstar>& op;
    Dune::Preconditioner<X,Xstar>& prec;
    DualPairing<X,Xstar> const& dp;
    PCGTerminationCriterion<double>& terminate;
    int verbose;
    Result result;
    double eps, energyNorm2;
    std::string pre = std::string("KASKADE IterateType::PCG: ");
    Functor f;
  };
} // namespace Kaskade
#endif
