/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copied from dune/istl/solvers.hh and adapted to NM III notation          */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef NMIII_PCG_HH
#define NMIII_PCG_HH

#include "dune/common/static_assert.hh"
#include "dune/istl/solvers.hh"
#include "dune/istl/istlexception.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"

namespace Kaskade {

 //! \brief conjugate gradient method
  template<class X>
  class NMIIIPCGSolver : public Dune::InverseOperator<X,X> {
  public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;

    /*! 
      \brief Set up conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&,double,int,int)
    */
    template<class L, class P>
    NMIIIPCGSolver (L& op, P& prec, double reduction, int maxit, int verbose) : 
      ssp(), _op(op), _prec(prec), _sp(ssp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
        dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
                            "L and P must have the same category!");
        dune_static_assert( static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                            "L must be sequential!");
    }

    /*! 
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& u, X& b, Dune::InverseOperatorResult& res)
    {
      res.clear();                  // clear solver statistics
      Dune::Timer watch;                // start a timer
    
      X r(u); 
      X rq(u);
      X q(u);
      X h(u);
    
      r = b;
      _op.applyscaleadd(-1,u,r);  // r = r-Au
      rq = 0;
      _prec.apply(rq,r);
      q = rq;

      double def0 = _sp.norm(rq);// compute norm
      if (def0<1E-30)    // convergence check 
      {
        res.converged  = true;
        res.iterations = 0;               // fill statistics
        res.reduction = 0;
        res.conv_rate  = 0;
        res.elapsed=0;
        if (_verbose>0)                 // final print 
                        std::cout << "=== rate=" << res.conv_rate
                                  << ", T=" << res.elapsed << ", TIT=" << res.elapsed
                                  << ", IT=0" << std::endl;
        return;
      }

      if (_verbose>0)             // printing
      {
        std::cout << "=== NMIIIPCGSolver" << std::endl;
        if (_verbose>1) {
          this->printHeader(std::cout);
          this->printOutput(std::cout,(double) 0,def0);
        }
      }

      // some local variables
      double def=def0, defnew=def0;   // loop variables
      field_type alpha,beta,sigma;
      sigma = _sp.dot(r,rq); 

      // the loop
      int i=1; 
      for ( ; i<=_maxit; i++ )
      {
        // minimize in given search direction p
        _op.apply(q,h);             // h=Aq
        alpha = sigma/_sp.dot(q,h);       // scalar product
        u.axpy(alpha,q);           // update solution

        // convergence test

        if (_verbose>1)             // print
          this->printOutput(std::cout,(double) i,defnew,def);

        if (defnew<sqrt(_reduction) || def<1E-30)    // convergence check
        {
          res.converged  = true;
          break;
        }

        r.axpy(-alpha,h);
        rq = 0;
        _prec.apply(rq,r);

        defnew=_sp.norm(rq);        // comp defect norm

        // determine new search direction
        double sigmaNew = _sp.dot(r,rq); 
        beta = sigmaNew/sigma;
        sigma = sigmaNew;
        q *= beta;                  // scale old search direction
        q += rq;                     // orthogonalization with correction
        def = defnew;               // update norm
      }

      if (_verbose>0)                // printing for non verbose
        this->printOutput(std::cout,(double) i,def);

      res.iterations = i;               // fill statistics
      res.reduction = def/def0;
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();

      if (_verbose>0)                 // final print 
      {
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/i
                  << ", IT=" << i << std::endl;
      }
  }

    /*! 
      \brief Apply inverse operator with given reduction factor.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
  virtual void apply (X& x, X& b, double reduction, Dune::InverseOperatorResult& res)
  {
    _reduction = reduction;
    (*this).apply(x,b,res);
  }

  private:
  Dune::SeqScalarProduct<X> ssp;
  Dune::LinearOperator<X,X>& _op;
  Dune::Preconditioner<X,X>& _prec;
  Dune::ScalarProduct<X>& _sp;
  double _reduction;
  int _maxit;
  int _verbose;
  };
}
#endif
