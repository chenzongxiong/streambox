/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef UZAWA_H
#define UZAWA_H

#include <boost/fusion/sequence.hpp>

#include "dune/common/static_assert.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvers.hh"

#include "fem/linearspace.hh"

namespace Kaskade
{
  /**
   * Convenience function for simple application of
   * preconditioners. Calls pre(), apply(), and post() sequentially.
   */
  template <class X, class Y>
  void applyPreconditioner(Dune::Preconditioner<X,Y>& p, X& x, Y& y) {
    p.pre(x,y);
    p.apply(x,y);
    p.post(x);
  }

  /**
   * An inexact Uzawa solver for solving symmetric saddle point systems
   *
   * \f[ \begin{array}{cc} A & B^T \\ B & \end{array}\begin{matrix}{c} u \\ \lambda \end{matrix} =
   * \begin{matrix}{c} f \\ g \end{matrix} \f]
   *
   * It implements an inexact preconditioned conjugate gradient method for the
   * Schur complement \f$ BA^{-1}B^T\f$.
   *
   * Template parameters:
   * - X: the linear space of the primal variable \f$ u \f$
   * - Y: the linear space of the dual variable \f$ \lambda \f$
   */
  template<class X, class Y>
  class UzawaSolver : public Dune::InverseOperator<LinearProductSpace<typename X::field_type,boost::fusion::vector<X,Y> >,
  LinearProductSpace<typename X::field_type,boost::fusion::vector<X,Y> > > {
  public:
    typedef LinearProductSpace<typename X::field_type,boost::fusion::vector<X,Y> > domain_type;
    typedef domain_type                                                            range_type;
    typedef typename X::field_type                                                 field_type;
    typedef field_type Scalar;
    typedef domain_type Domain;
    typedef range_type Range;

    /**
     * Constructor.
     * \param opA_ the linear operator \f$ A \f$
     * \param solA_ the approximate inverse operator \f$ A^{-1} \f$, needs to be a subclass of Dune::InverseOperator<X,X>
     * \param opB_ the linear operaotr \f$ B \f$
     * \param opBt_ the linear operator \f$ B^T \f$
     * \param pb_ the Schur complement preconditioner
     * \param reduction_ the default required reduction in preconditioned Schur complement's residual
     * \param maxit_ the maximum number of iterations
     * \param verbose_ printing level to stdout (0=silent, 1=result, 2=every iteration)
     */
    template<class A, class SolA, class B, class Bt, class PB>
    UzawaSolver (A& opA_, SolA& solA_, B& opB_, Bt& opBt_, PB& pb_, double reduction_, int maxit_, int verbose_) :
    opA(opA_), solA(solA_), opB(opB_), opBt(opBt_), precB(pb_), reduction(reduction_), maxit(maxit_), verbose(verbose_) {
      using namespace Dune;
      dune_static_assert( static_cast<int>(PB::category) == static_cast<int>(SolverCategory::sequential) , "Preconditioner PB should be sequential.");
      dune_static_assert( static_cast<int>(A::category) == static_cast<int>(SolverCategory::sequential)  , "Operator A should be sequential");
      dune_static_assert( static_cast<int>(B::category) == static_cast<int>(SolverCategory::sequential)  , "Operator B should be sequential");
      dune_static_assert( static_cast<int>(Bt::category) == static_cast<int>(SolverCategory::sequential) , "Operator Bt should be sequential");
    }


    virtual void apply (Domain& x, Range& b, Dune::InverseOperatorResult& res) {
      using namespace boost::fusion;
      using namespace Dune;


      Timer watch;                // start a timer

      InverseOperatorResult resA;
      if (verbose>0) printf("=== UzawaSolver\n");

      // Declare variables
      X& u(at_c<0>(x.data));
      X& f(at_c<0>(b.data));
      X r(u), h(u), rSave(u);
      Y& lambda(at_c<1>(x.data));
      Y g(at_c<1>(b.data)), d(g), p(g);
      Y tmp(g); // WARNING: needed only for applying P twice!

      // initialize
      r = f;
        opBt.applyscaleadd(-1,lambda,r);  // r = f-B'*u
      solA.apply(u,r,resA);                    // Au = r
      opB.applyscaleadd(-1,u,g);               // g = g-Bu
      applyPreconditioner(precB,tmp,g);
      applyPreconditioner(precB,d,tmp); // d = Pg // WARNING: assumes precB^2 is preco
      double sigma = spY.dot(g,d);            // sigma = g'd
      double sigma0 = sigma;

      //std::cout << "sigma = " << sigma << "\n g=" << at_c<0>(g.data) << "d=" << at_c<0>(d.data);

      // Perform stationary iteration.
      int i=0;
      for ( ; i<maxit && sigma>reduction*reduction*sigma0; i++) {

        // std::cout << "d=P(g-Bu)=" << at_c<0>(d.data);

        opBt.apply(d,r);                       // r = B'd
        rSave = r;
        solA.apply(h,r,resA);                  // Ah = r
        r = rSave;

        double alpha = sigma / spX.dot(h,r);   // alpha = sigma / h'r
        lambda.axpy(-alpha,d);                  // lambda = lambda-alpha d
        u.axpy(alpha,h);                      // u = u+alpha h
        opB.applyscaleadd(-alpha,h,g);          // g = g-alpha Bh

        applyPreconditioner(precB,tmp,g);        // p = Pg
        applyPreconditioner(precB,p,tmp); // WARNING

        double sigmaOld = sigma;
        sigma = spY.dot(g,p);
        double beta = sigma / sigmaOld;


        d *= beta; d.axpy(1,p);               // d = beta d + p


        if (verbose>1) {
          if (i%30==0)
            std::printf("%5s %14s %14s\n","Iter","Preco. Res.","Rate");
          std::printf("%5d %14.4e %14.4e\n",i,std::sqrt(sigma),std::sqrt(beta));
        }
      }

      // Fill statistics
      res.clear();
      res.iterations = i;
      res.reduction = std::sqrt(sigma/sigma0);
      res.converged = res.reduction<=reduction;
      res.conv_rate = std::pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();


      if (verbose>0)                 // final print
        printf("=== rate=%g, time=%g, time/it=%g, iter=%d\n",res.conv_rate,res.elapsed,res.elapsed/i,i);
    }

    virtual void apply (Domain& x, Range& b, double reduction_, Dune::InverseOperatorResult& res)
    {
      reduction = reduction_;
      this->apply(x,b,res);
    }

    void apply(Domain& x, Range& b)
    {
      Dune::InverseOperatorResult tmpResult;
      apply(x,b,tmpResult);
    }

  private:
    Dune::SeqScalarProduct<X> spX;
    Dune::SeqScalarProduct<Y> spY;

    Dune::LinearOperator<X,X>& opA;
    Dune::InverseOperator<X,X>& solA;

    Dune::LinearOperator<X,Y>& opB;
    Dune::LinearOperator<Y,X>& opBt;
    Dune::Preconditioner<Y,Y>& precB;

    double reduction;
    int    maxit;
    int    verbose;
  };
}

#endif
