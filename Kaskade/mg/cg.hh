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

#ifndef NMII_CG_HH
#define NMII_CG_HH

#include "dune/istl/istlexception.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"

namespace Kaskade {

 //! \brief conjugate gradient method
  template<class X>
  class NMIIICGSolver : public InverseOperator<X,X> {
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
    NMIIICGSolver (L& op, P& prec, double reduction, int maxit, int verbose) : 
      ssp(), _op(op), _sp(ssp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
        dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
                            "L and P must have the same category!");
        dune_static_assert( static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
                            "L must be sequential!");
    }

    /*! 
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& u, X& b, InverseOperatorResult& res)
    {
      res.clear();                  // clear solver statistics
      Timer watch;                // start a timer
    
      X r(u);              // the search direction
      X p(u);              // the search direction
      X q(u);              // a temporary vector
    
      r = b;
      _op.applyscaleadd(-1,u,r);  // r = r-Au
      p = r;

      double def0 = _sp.norm(r);// compute norm
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
        std::cout << "=== NMIIICGSolver" << std::endl;
        if (_verbose>1) {
          this->printHeader(std::cout);
          this->printOutput(std::cout,0,def0);
        }
      }

      // some local variables
      double def=def0, defnew=def0;   // loop variables
      field_type alpha,beta;

      // the loop
      int i=1; 
      for ( ; i<=_maxit; i++ )
      {
        // minimize in given search direction p
        _op.apply(p,q);             // q=Ap
        alpha = def/_sp.dot(p,q);       // scalar product
        u.axpy(alpha,p);           // update solution

        // convergence test

        if (_verbose>1)             // print
          this->printOutput(std::cout,i,defnew,def);

        if (defnew<def0*_reduction || def<1E-30)    // convergence check
        {
          res.converged  = true;
          break;
        }

        r.axpy(-alpha,q);          // update defect
        defnew=_sp.norm(r);        // comp defect norm

        // determine new search direction
        beta = defnew/def;         // scaling factor
        p *= beta;                  // scale old search direction
        p += r;                     // orthogonalization with correction
        def = defnew;               // update norm
      }

      if (_verbose>0)                // printing for non verbose
        this->printOutput(std::cout,i,def);

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
  virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res)
  {
    _reduction = reduction;
    (*this).apply(x,b,res);
  }

  private:
  SeqScalarProduct<X> ssp;
  LinearOperator<X,X>& _op;
  ScalarProduct<X>& _sp;
  double _reduction;
  int _maxit;
  int _verbose;
  };
}
#endif
