#ifndef FGMRES_HH
#define FGMRES_HH

#include<cmath>
#include<complex>
#include<iostream>
#include<iomanip>
#include<string>

#include "dune/istl/solvers.hh"
#include "dune/istl/istlexception.hh"
#include "dune/istl/operators.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/scalarproducts.hh"
#include <dune/common/timer.hh>
#include <dune/common/static_assert.hh>


  /**
     \brief implements the Generalized Minimal Residual (GMRes) method

     GMRes solves the unsymmetric linear system Ax = b using the
     Generalized Minimal Residual method as described the SIAM Templates
     book (http://www.netlib.org/templates/templates.pdf).
     
     \todo F durch rebind erzeugen und nur den field_type für F übergeben
     
  */

  template<class X>
  class RestartedFGMResSolver : public Dune::InverseOperator<X,X>
  {
  public:
    typedef X Y;
    typedef X F;

    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef Y range_type;
    //! \brief The field type of the operator to be inverted
    typedef typename X::field_type field_type;
    //! \brief The field type of the basis vectors
    typedef F basis_type;
    
    /*! 
      \brief Set up solver.
      
      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
      \param restart number of GMRes cycles before restart
      \param recalc_defect recalculate the defect after everey restart or not [default=false]
    */
    template<class L, class P>
    RestartedFGMResSolver (L& op, P& prec, double reduction, int restart, int maxit, int verbose, bool recalc_defect = false) : 
      _A_(op), _M(prec),
      ssp(), _sp(ssp), _restart(restart),
      _reduction(reduction), _maxit(maxit), _verbose(verbose),
      _recalc_defect(recalc_defect)
    {
      dune_static_assert(static_cast<int>(P::category) == static_cast<int>(SolverCategory::sequential),
        "P must be sequential!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
        "L must be sequential!");
    }

    /*!
      \brief Set up solver.
      
      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
      \param restart number of GMRes cycles before restart
      \param recalc_defect recalculate the defect after everey restart or not [default=false]
    */
    template<class L, class S, class P>
    RestartedFGMResSolver (L& op, S& sp, P& prec, double reduction, int restart, int maxit, int verbose, bool recalc_defect = false) :  
      _A_(op), _M(prec),
      _sp(sp), _restart(restart),
      _reduction(reduction), _maxit(maxit), _verbose(verbose),
      _recalc_defect(recalc_defect)
    {
      dune_static_assert(static_cast<int>(P::category) == static_cast<int>(SolverCategory::sequential),
        "P must be sequential!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(SolverCategory::sequential),
        "L must be sequential!");
    }

    //! \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    virtual void apply (X& x, X& b, Dune::InverseOperatorResult& res)
    {
      apply(x,b,_reduction,res);
    };
    
    /*! 
      \brief Apply inverse operator.
      
      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, Y& b, double reduction, Dune::InverseOperatorResult& res)
    {
      int m = _restart;
      field_type norm;
      field_type norm_old = 0.0;
      field_type norm_0;
      field_type beta;
      int i, j = 1, k;
      std::vector<field_type> s(m+1), cs(m), sn(m);
      std::vector<field_type> sw(m+1), csw(m), snw(m);
      // helper vector
      X tmp(b.size());
      std::vector< std::vector<field_type> > H(m+1,s), Hw(m+1,s);
      std::vector<F> v(1,b);
      std::vector<X> w(1,b);

      // start timer
      Timer watch;                // start a timer

      // clear solver statistics
      res.clear();

      _M.pre(x,b);
        
        norm_0 = _sp.norm(b);
        if(x.size()!=0)
          _A_.applyscaleadd(-1,x,b);  // overwrite b with defect
        else
        {
          x.resize(b.size());
          x=0.0;
        }
        v[0]=b;
        beta = _sp.norm(v[0]);


      // avoid division by zero
      if (norm_0 == 0.0)
        norm_0 = 1.0;

      norm = norm_old = _sp.norm(v[0]);

      // print header
      if (_verbose > 0)
      {
        std::cout << "=== RestartedGMResSolver" << std::endl;
        if (_verbose > 1) 
        {
          this->printHeader(std::cout);
          this->printOutput(std::cout,0,norm);
        }
      }

      // check convergence
      if (norm <= reduction * norm_0) {
        _M.post(x);                  // postprocess preconditioner
        res.converged  = true;
        if (_verbose > 0)                 // final print
          print_result(res);
        return;
      }

      while (j <= _maxit && res.converged != true) {
        v[0] *= (1.0 / beta);
        for (i=1; i<=m; i++) s[i] = 0.0;
        s[0] = beta;
    
        for (i = 0; i < m && j <= _maxit && res.converged != true; i++, j++) {
          if(w.size()<i+1) w.push_back(b);
          w[i] = 0.0;
          if(v.size()<i+2) v.push_back(b);
          v[i+1] = 0.0; 
          _M.apply(w[i], v[i]);
          _A_.apply(w[i], /* => */ v[i+1]);
          for (k = 0; k <= i; k++) {
            H[k][i] = _sp.dot(v[i+1], v[k]);
            // w -= H[k][i] * v[k];
            v[i+1].axpy(-H[k][i], v[k]);
          }
          H[i+1][i] = _sp.norm(v[i+1]);
          if (H[i+1][i] == 0.0)
            DUNE_THROW(ISTLError,"breakdown in GMRes - |w| "
              << w[i] << " == 0.0 after " << j << " iterations");
          // v[i+1] = w * (1.0 / H[i+1][i]);
          v[i+1] *= (1.0 / H[i+1][i]);

          for (k = 0; k < i; k++)
            applyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);
      
          generatePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
          applyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
          applyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

          norm = std::abs(s[i+1]);
                
          if (_verbose > 1)             // print
          {
            this->printOutput(std::cout,j,norm,norm_old);
          }

          norm_old = norm;
                
          if (norm < reduction * norm_0) {
            res.converged = true;
          }
        }

        tmp = 0.0;
        // calc update vector
        update(tmp, i - 1, H, s, w);
        if (_verbose > 1)
        {
          std::cout << "|dx|:" << _sp.norm(tmp) << " " << std::flush;
        };
        
        // r = (b - A * x);
        // update defect
        _A_.applyscaleadd(-1,tmp, /* => */ b);

        x += tmp;

        beta = _sp.norm(b);
        norm = beta;
        v[0]=b;

        res.converged = false;
                        
        norm_old = norm;
                
        if (norm < reduction * norm_0) {
          // fill statistics
          res.converged = true;
        }
            
        if (res.converged != true && _verbose > 0)
          std::cout << "=== GMRes::restart\n";
      }
  
      _M.post(x);                  // postprocess preconditioner
        
      res.iterations = j;
      res.reduction = norm / norm_0;
      res.conv_rate  = pow(res.reduction,1.0/j);
      res.elapsed = watch.elapsed();

      if (_verbose>0)
        print_result(res);
    }
  private:

    void
    print_result (const Dune::InverseOperatorResult & res) const
    {
      int j = res.iterations>0?res.iterations:1;
      std::cout << "=== rate=" << res.conv_rate
                << ", T=" << res.elapsed
                << ", TIT=" << res.elapsed/j
                << ", IT=" << res.iterations
                << std::endl;
    }
    
    static void 
    update(X &x, int k,
      std::vector< std::vector<field_type> > & h,
      std::vector<field_type> & s, std::vector<F> v)
    {
      std::vector<field_type> y(s);
        
      // Backsolve:  
      for (int i = k; i >= 0; i--) {
        y[i] /= h[i][i];
        for (int j = i - 1; j >= 0; j--)
          y[j] -= h[j][i] * y[i];
      }
        
      for (int j = 0; j <= k; j++)
        // x += v[j] * y[j];
        x.axpy(y[j],v[j]);
    }
    
    void
    generatePlaneRotation(field_type &dx, field_type &dy, field_type &cs, field_type &sn)
    {
      if (dy == 0.0) {
        cs = 1.0;
        sn = 0.0;
      } else if (std::abs(dy) > std::abs(dx)) {
        field_type temp = dx / dy;
        sn = 1.0 / std::sqrt( 1.0 + temp*temp );
        cs = temp * sn;
      } else {
        field_type temp = dy / dx;
        cs = 1.0 / std::sqrt( 1.0 + temp*temp );
        sn = temp * cs;
      }
    }


    void
    applyPlaneRotation(field_type &dx, field_type &dy, field_type &cs, field_type &sn)
    {
      field_type temp  =  cs * dx + sn * dy;
      dy = -sn * dx + cs * dy;
      dx = temp;
    }
    
    Dune::LinearOperator<X,X>& _A_;
    Dune::Preconditioner<X,X>& _M;
    Dune::SeqScalarProduct<X> ssp;
    Dune::ScalarProduct<X>& _sp;
    int _restart;
    double _reduction;
    int _maxit;
    int _verbose;
    bool _recalc_defect;
  };

  /** @} end documentation */

#endif
