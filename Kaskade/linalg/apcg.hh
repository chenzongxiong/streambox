/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef NMIII_APCG_HH
#define NMIII_APCG_HH

#include <numeric>
#include <string>

#include <boost/circular_buffer.hpp>
#include <boost/timer/timer.hpp>

#include "dune/common/timer.hh"
#include "dune/istl/istlexception.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"
#include "dune/istl/solvers.hh"

#include "linalg/symmetricOperators.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/geometric_sequence.hh"
#include "utilities/scalar.hh"

namespace Kaskade
{
  /**
   * \ingroup iterative
   * \brief Whether a preconditioner needs zero initialized result vector or not.
   * 
   * For a couple of known preconditioners, this function returns the actual information.
   * For unknown preconditioners, it errs on the conservative side, i.e. reports true.
   * 
   * The function is moderately expensive, hence don't use it inside of loops.
   */
  template <class Domain, class Range>
  bool requiresInitializedInput(Dune::Preconditioner<Domain,Range> const* p)
  {
    if (dynamic_cast<Dune::Richardson<Domain,Range> const*>(p))
      return false;
    if (auto sp = dynamic_cast<SymmetricPreconditioner<Domain,Range> const*>(p))
      return sp->requiresInitializedInput();
    return true;
  }
  
  /**
   * \ingroup iterative
   * \brief Interface for IterateType::PCG termination criterion policy classes
   * \tparam R a floating point type for real numbers
   */
  template <class R>
  class PCGTerminationCriterion 
  {
  public:
    /**
     * \brief real field type
     */
    typedef R Real;
    
    /**
     * \brief re-initializes the termination criterion for a new IterateType::CG run
     */
    virtual void clear() = 0;
    
    /**
     * \brief set requested tolerance
     * 
     * \param tol the requested tolerance (nonnegative)
     */
    virtual PCGTerminationCriterion<R>& tolerance(Real tol) = 0;
    
    /**
     * \brief supplies energy of step to the termination criterion
     * \param gamma2 the energy of the step
     */
    virtual void step(Real gamma2) = 0;
    
    /**
     * \brief supplies the preconditioned residual to the termination criterion
     * \param sigma the preconditioned residual norm
     */
    virtual void residual(Real sigma) = 0;
    
    /**
     * \brief termination decision
     * \return true if the iteration has reached the required accuracy
     */
    virtual operator bool() const = 0;
  };
  
  /**
   * \ingroup iterative
   * \brief TerminationCriterion based on an absolute energy error estimate
   * 
   * This termination criterion assumes the error reduction occurs in two phases: First a smoothing phase where
   * \f$ \epsilon_k^2 \approx a^2/(1+\sqrt{k})^2 + b^2 \f$, and a linear convergence phase with \f$ \epsilon_k \approx q^k \epsilon_0 \f$.
   * The termination criterion tries to fit the model parameters \f$ a, b, q \f$ and to decide where the crossover of both 
   * phases occurs. Based on this decision, the energy error \f$ \epsilon_k \f$ is estimated as \f$ [\epsilon_k] \f$.
   * 
   * For greater robustnes, the actual error estimate used is \f[ [\epsilon_{k-l}]^2 := [\epsilon_k]^2 + \sum_{i=k-l}^{k-1} \gamma_i^2. \f] with
   * \f$ \gamma_i = \|u_{i+1}-u_i\|_A \f$.The lookahead value \f$ l \f$ can be set using \ref lookahead.
   * 
   * \tparam R the real field type. Currently instantiated for float and double in file apcg.cpp.
   */
  template <class R>
  class PCGEnergyErrorTerminationCriterion: public PCGTerminationCriterion<R> {
  public:
    typedef R Real;
    
    /**
     * \brief constructor
     * 
     * The pcg iteration is terminated as soon as either \f$ [\epsilon] \le \mathrm{atol} \f$ or 
     * the number of iterations exceeds the limit maxit.
     * 
     * \param atol the absolute error tolerance for termination
     * \param maxit the maximum number of iterations
     */
    PCGEnergyErrorTerminationCriterion(Real atol, int maxit);
    
    /**
     * \brief re-initializes the termination criterion for a new IterateType::CG run
     */
    virtual void clear();
    
    /**
     * \brief set requested absolute tolerance
     * 
     * \param atol the requested tolerance (nonnegative)
     */
    virtual PCGEnergyErrorTerminationCriterion<Real>& tolerance(Real atol);
    
    /**
     * \brief set requested lookahead value
     * 
     * \param lah the requested lookahead (nonnegative)
     * 
     * The default value is 6.
     */
     PCGEnergyErrorTerminationCriterion<Real>& lookahead(int lah);
    
    /**
     * \brief supplies energy of step to the termination criterion
     * \param gamma2 the energy \f$ \gamma^2 = \| u_{k+1}-u_k \|_A^2 \f$ of the step
     */
    virtual void step(Real gamma2);
    
    /**
     * \brief supplies the preconditioned residual to the termination criterion
     * \param sigma the preconditioned residual norm
     */
    virtual void residual(Real sigma) {}
    
    /**
     * \brief termination decision
     * \return true if the iteration has reached the required accuracy
     */
    virtual operator bool() const;
    
    /**
     * \brief returns the estimated absolute energy error
     */
    Real error() const;
    
    std::vector<Real> const& gamma2() const {return gammas2; }
    
  private:
    Real tol;
    int maxit;
    // squared gammas
    std::vector<Real> gammas2;
    int lookah;
  };
  

 /**
  * \ingroup iterative
  * \brief preconditioned conjugate gradient method
  * 
  * This implements a preconditioned IterateType::CG iteration for an operator \f$ A: X\to X^* \f$, preconditioned by a
  * preconditioner \f$ B^{-1}: X^* \to X \f$. The termination is policy-based.
  * 
  * The implementation follows Deuflhard/Weiser, Section 5.3.3.
  * 
  * \tparam X the type of vectors from the primal space
  * \tparam Xstar the type of vectors from the dual space
  * 
  * \see PCGEnergyErrorTerminationCriterion
  */
  template<class X, class Xstar>
  class Pcg : public Dune::InverseOperator<X,Xstar> {
  public:
    /// \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    /// \brief The range type of the operator to be inverted.
    typedef Xstar range_type;
    /// \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;
    
    /**
     * \brief the real field type corresponding to field_type
     */
    typedef typename ScalarTraits<field_type>::Real Real;

    /** 
     * \brief Set up conjugate gradient solver.
     * 
     * \param op the operator
     * \param prec the preconditioner
     * \param terminate the termination criterion. The object has to exist during the lifetime of the pcg object as it is referenced.
     * \param verbose controls the verbosity of logging to std::cout. 0 means no output at all (default). Values in 0,1,2 are valid.
     */
    Pcg(SymmetricLinearOperator<X,Xstar>& op_, SymmetricPreconditioner<X,Xstar>& prec_, 
        PCGTerminationCriterion<Real>& terminate_, int verbose_=0) : 
      proxyOp(op_,zeroDp), op(op_), 
      proxyPrec(prec_,zeroDp), prec(prec_), 
      terminate(terminate_), verbose(verbose_), requiresInit(requiresInitializedInput(&prec))
    {
// Do we need this in Kaskade7? 
//         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
//                             "L and P must have the same category!");
//         dune_static_assert( static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
//                             "L must be sequential!");
    }

    /** 
     * \brief Set up conjugate gradient solver.
     * 
     * \param op the operator
     * \param prec the preconditioner
     * \param dp the dual pairing with respect to which the operator is symmetric
     * \param terminate the termination criterion. The object has to exist during the lifetime of the pcg object as it is referenced.
     * \param verbose controls the verbosity of logging to std::cout. 0 means no output at all (default). Values in 0,1,2 are valid.
     */
    Pcg(Dune::LinearOperator<X,Xstar>& op_, Dune::Preconditioner<X,Xstar>& prec_, DualPairing<X,Xstar> const& dp_,
        PCGTerminationCriterion<Real>& terminate_, int verbose_=0) : 
      proxyOp(op_,dp_), op(proxyOp), 
      proxyPrec(prec_,dp_), prec(proxyPrec), 
      terminate(terminate_), verbose(verbose_), requiresInit(requiresInitializedInput(&prec))
    {    }
    
    /**
     * \brief Apply inverse operator by performing a number of IterateType::PCG iterations.
     * 
     * \param u the initial value (starting iterate)
     * \param b the right hand side (which will not be modified)
     */
    virtual void apply (X& u, Xstar& b, Dune::InverseOperatorResult& res) {
      res.clear();                // clear solver statistics
      terminate.clear();          // clear termination criterion
      Dune::Timer watch;          // start a timer
      
      boost::timer::cpu_timer mvTimer;
      mvTimer.stop();
      boost::timer::cpu_timer apTimer;
      apTimer.stop();
      boost::timer::cpu_timer prTimer;
      prTimer.stop();
    
      Xstar r(b); 
      op.applyscaleadd(-1,u,r);  // r = b-Au

      // Keeping rq and q in unique_ptr allows efficient computation of q = rq + beta*q below by swapping the pointers.
      // More elegant would be a swap method directly on X, but this is hard to implement as there is no swap either 
      // for boost::fusion vectors nor for Dune::BlockVectors. Even more elegant would be a loop fusion expression
      // framework for BLAS level 1 operations.
      std::unique_ptr<X> rq(new X(u)); *rq = 0;
      prec.apply(*rq,r); // rq = B^{-1} r
      
      std::unique_ptr<X> q(new X(*rq));

      X Aq(b);
    

      // some local variables
      field_type alpha,beta,sigma,gamma;
      sigma = op.dp(*rq,r); // preconditioned residual norm squard
      terminate.residual(ScalarTraits<field_type>::real(sigma));
      double const sigma0 = std::abs(sigma);
      
      // Check for trivial case
      if (sigma0 == 0)
        return;

      if (verbose>0) {            // printing
        std::cout << "=== Kaskade7 IterateType::PCG" << std::endl;
        if (verbose>1) {
          this->printHeader(std::cout);
          this->printOutput(std::cout,0.0,sigma0);
        }
      }

      // the loop
      int i=1; 
      for ( ; true; i++ ) {
        // minimize in given search direction p
        mvTimer.resume();
        field_type qAq = op.applyDp(*q,Aq);             // compute Aq = A*q and qAq = <q,A*q>
if (std::isnan(qAq)) std:: cerr << "qAq is nan\n";        
        mvTimer.stop();
        	
	if (std::abs(qAq) < 1e-28*sigma) 
	  break; // oops.
	  
        if (qAq <= 0)
          throw NonpositiveMatrixException("encountered nonpositive energy product " + std::to_string(qAq) + " vs sigma " + std::to_string(sigma) + " in IterateType::PCG.",__FILE__,__LINE__);
          
        alpha = sigma/qAq;
        apTimer.resume();
        u.axpy(alpha,*q);
        apTimer.stop();
        gamma = sigma*alpha;
	terminate.step(ScalarTraits<field_type>::real(gamma));
	resDebug.push_back(std::sqrt(sigma));

        // convergence test
        if (terminate)
	  break;

        apTimer.resume();
        r.axpy(-alpha,Aq); // r = r - alpha*A*q
        apTimer.stop();
        
        prTimer.resume();
        if (requiresInit)
          *rq = 0;
        double sigmaNew = prec.applyDp(*rq,r); // compute rq = B^{-1}r and <rq,r>
        prTimer.stop();

        // determine new search direction
	terminate.residual(ScalarTraits<field_type>::real(sigmaNew));
        if (verbose>1)             // print
          this->printOutput(std::cout,static_cast<double>(i),std::abs(sigmaNew),std::abs(sigma));

	// convergence check to prevent division by zero if we obtained an exact solution
	// (which may happen for low dimensional systems)
	if (std::abs(sigmaNew) < 1e-28*sigma0)
	  break;
	
        beta = sigmaNew/sigma;
        sigma = sigmaNew;
        apTimer.resume();
        rq->axpy(beta,*q); // compute q = rq + beta*q
        std::swap(rq,q);
        apTimer.stop();
      }

      if (verbose>0)                // printing for non verbose
        this->printOutput(std::cout,static_cast<double>(i),std::abs(sigma));

      res.iterations = i;               // fill statistics
      res.reduction = std::sqrt(std::abs(sigma)/sigma0);
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();

      if (verbose>0)                 // final print 
      {
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/i
                  << ", IT=" << i << std::endl;
        
        std::cout << "Matrix-Vector time: " << mvTimer.format() 
                  << "vector ops        : " << apTimer.format() 
                  << "preconditioner    : " << prTimer.format() << '\n';
      } 
  }

    /** 
     * \brief Apply inverse operator with given tolerance.
     * 
     * This method is equivalent to setting first the tolerance in the termination criterion, then calling apply().
     */
  virtual void apply (X& x, X& b, double tol, Dune::InverseOperatorResult& res) {
    terminate.tolerance(tol);
    (*this).apply(x,b,res);
  }

  private:
    // The solver may be constructed with either symmetric operator/preconditioner interface provided,
    // or with Dune interfaces provided. The implementation relies on the symmetric interface
    // in order to exploit performance gains with simultaneous evaluation of matrix-vector products
    // and dual pairings.
    // Therefore we hold "working" references to symmetric operator/preconditioners. Those reference
    // either the provided ones, or wrapper objects relaying the evaluation to provided Dune interfaces.
    // For the second case, we hold the (small) wrapper objects ourselves. In the first case, those
    // wrapper objects have to be initialized correctly (even though they are never accessed). For this 
    // we need a dual pairing, hence the dummy zeroDp.
    ZeroDualPairing<X,Xstar>                zeroDp;  // just a dummy
    
    SymmetricLinearOperatorWrapper<X,Xstar> proxyOp; // proxy in case a Dune::LinearOperator is provided
    SymmetricLinearOperator<X,Xstar>& op;  // working reference to operator

    SymmetricPreconditionerWrapper<X,Xstar>  proxyPrec; // proxy in case a Dune::Preconditioner is provided
    SymmetricPreconditioner<X,Xstar>&        prec;      // working reference to preconditioner
    
    PCGTerminationCriterion<Real>& terminate;
    int verbose;
    bool requiresInit;
    
    std::vector<double> resDebug;
  };
} // namespace Kaskade
#endif
