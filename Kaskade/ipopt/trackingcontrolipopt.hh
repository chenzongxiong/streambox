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

#ifndef TRACKINGCONTROLIPOPT_HH
#define TRACKINGCONTROLIPOPT_HH

#include <memory>
#include <algorithm>
#include <cstdlib>

#include <boost/timer.hpp>

#include "IpTNLP.hpp"

#include "assemble.hh"
#include "functional_aux.hh"


/**
 * Interface to Ipopt.
 *
 * The VariationalFunctional provided has to implement Ipopt
 * semantics, i.e., in evaluating the right hand side and the
 * functional the value of lambda is implicitly assumed to be zero,
 * independent of the provided value. This means that the gradient and
 * value of the objective have to be computed instead of gradient and
 * value of the Lagrange functional.
 */
template <class VariationalFunctional>
class TrackingControlInterface: public Ipopt::TNLP 
{
  typedef VariationalFunctionalAssembler<LinearizationAt<VariationalFunctional> >  Gop;
  typedef typename Gop::RT RT;
  typedef typename Gop::AnsatzVariableSet VarSet;
  typedef typename VarSet::Grid Grid;
  typedef VariationalFunctional Functional;
  
  typedef Ipopt::Index Index;
  typedef Ipopt::TNLP::IndexStyleEnum IndexStyleEnum;
  typedef Ipopt::SolverReturn SolverReturn;
  typedef Ipopt::Number Number;
  

  enum { U = 0, Y = 1, LAMBDA = 2, END = 3 };
  typedef int AssemblyParts;
  
  
  
public:
  TrackingControlInterface(GridSignals& signals,
                           VarSet const& varSet_,
                           Functional& f_):
    xx(varSet_),
    varSet(varSet_),
    gop(signals,varSet.spaces),
    f(f_),
    nu(varSet.degreesOfFreedom(U,U+1)),
    ny(varSet.degreesOfFreedom(Y,Y+1)),
    nl(varSet.degreesOfFreedom(LAMBDA,LAMBDA+1)),
    n(nu+ny)
  {
    f.initVars(xx);
  }
  
  
  /** overload this method to return the number of variables
   *  and constraints, and the number of non-zeros in the jacobian and
   *  the hessian. The index_style parameter lets you specify C or Fortran
   *  style indexing for the sparse matrix iRow and jCol parameters.
   *  C_STYLE is 0-based, and FORTRAN_STYLE is 1-based.
   */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) 
  {

    n = this->n;
    m = nl;

    nnz_jac_g = gop.nnz(LAMBDA,END,U,Y+1,false);
    
    int nnz = gop.nnz(U,LAMBDA,U,LAMBDA,true);
    std::vector<int> rows(nnz), cols(nnz);
    std::vector<RT> data(nnz);
    typename std::vector<int>::iterator ir = rows.begin(), ic = cols.begin();
    typename std::vector<RT>::iterator id = data.begin();
    gop.toTriplet(U,LAMBDA,U,LAMBDA,ir,ic,id,true);
    nnz_h_lag = 0;
    for (int i=0; i<nnz; ++i)
      if (rows[i]>=cols[i])
        ++nnz_h_lag;

    index_style = C_STYLE;

    std::cout << "nu=" << nu << " ny=" << ny << " nl=" << nl << " nnz=" << nnz << " nnz_jac_g=" << nnz_jac_g << '\n';
    std::cout.flush();
    
    return true;
  }
  
  /** overload this method to return the information about the bound
   *  on the variables and constraints. The value that indicates
   *  that a bound does not exist is specified in the parameters
   *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
   *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
   *  1e19. (see TNLPAdapter) */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u) 
  {
    assert(n == ny+nu);
    assert(m == nl);

    f.boxConstraints(x_l+nu,x_u+nu,x_l,x_u);
    
    std::fill_n(g_l,nl,0);
    std::fill_n(g_u,nl,0);
    
    return true;
  }
  
  /** overload this method to return the starting point. The bools
   *  init_x and init_lambda are both inputs and outputs. As inputs,
   *  they indicate whether or not the algorithm wants you to
   *  initialize x and lambda respectively. If, for some reason, the
   *  algorithm wants you to initialize these and you cannot, set
   *  the respective bool to false.
   */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda) 
  {
    if (init_x) {
      std::vector<RT> data(n+m);
      xx.write(data.begin());
      // std::cout << "Initial u:\n";
      // std::copy(data.begin(),data.begin()+nu,std::ostream_iterator<RT>(std::cout," "));
      // std::cout << "\n";
      // std::cout.flush();
      
      std::copy(data.begin(),data.begin()+n,x);
    }
    if (init_z) {
      std::fill_n(z_L,n,1);
      std::fill_n(z_U,n,1);
    }
    if (init_lambda) {
      std::fill_n(lambda,m,0);
    }
    

    return true;
  }
  

  /** overload this method to return the value of the objective function */
  virtual bool eval_f(Index n, const Number* x, bool new_x,
                      Number& obj_value) 
  {
    update(new_x,x,false,0);
    obj_value = gop.functional();
    
    return true;
  }
  

  /** overload this method to return the vector of the gradient of
   *  the objective w.r.t. x */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                           Number* grad_f) 
  {
    update(new_x,x,false,0);
    gop.toSequence(U,LAMBDA,grad_f);

    return true;
  }
  

  /** overload this method to return the vector of constraint values */
  virtual bool eval_g(Index n, const Number* x, bool new_x,
                      Index m, Number* g) 
  {
    update(new_x,x,false,0);
    gop.toSequence(LAMBDA,END,g);

    return true;
  }
  
    
  /** overload this method to return the jacobian of the
   *  constraints. The vectors iRow and jCol only need to be set
   *  once. The first call is used to set the structure only (iRow
   *  and jCol will be non-NULL, and values will be NULL) For
   *  subsequent calls, iRow and jCol will be NULL. */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow,
                          Index *jCol, Number* values) 
  {
    if (iRow != 0) {
      // provide sparsity structure
      std::vector<RT> data(nele_jac);
      typename std::vector<RT>::iterator v = data.begin();
      gop.toTriplet(LAMBDA,END,U,LAMBDA,iRow,jCol,v,false);
    } else {
      assert(values);
      
      update(new_x,x,false,0);

      std::vector<Index> idx(nele_jac);
      typename std::vector<Index>::iterator i = idx.begin(), j = idx.begin();
      gop.toTriplet(LAMBDA,END,U,LAMBDA,i,j,values,false);
    }
    

    return true;
  }
  
    

  /** overload this method to return the hessian of the
   *  lagrangian. The vectors iRow and jCol only need to be set once
   *  (during the first call). The first call is used to set the
   *  structure only (iRow and jCol will be non-NULL, and values
   *  will be NULL) For subsequent calls, iRow and jCol will be
   *  NULL. This matrix is symmetric - specify the lower diagonal
   *  only.  A default implementation is provided, in case the user
   *  wants to se quasi-Newton approximations to estimate the second
   *  derivatives and doesn't not neet to implement this method. */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess,
                      Index* iRow, Index* jCol, Number* values) 
  {
    int nnz = gop.nnz(U,LAMBDA,U,LAMBDA,true);
    assert(nnz==nele_hess);

    if (values) {
      if (obj_factor != f.ipoptSigma) {
        f.ipoptSigma = obj_factor;
        assert(lambda);
        new_lambda = true;
        std::cout << "new sigma: " << obj_factor << "\n";
      }
      update(new_x,x,new_lambda,lambda);
      
      std::vector<int> unused(nnz);
      std::vector<int>::iterator ri = unused.begin(), ci = unused.begin();
      gop.toTriplet(U,LAMBDA,U,LAMBDA,ri,ci,values,true);
    } else {
      std::vector<RT> unused(nnz);
      gop.toTriplet(U,LAMBDA,U,LAMBDA,iRow,jCol,unused.begin(),true);
    }

    return true;
  }
  
    

  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value) 
  {
    int nnz = gop.nnz(U,END,U,END,false);
    std::vector<size_t> ri(nnz), ci(nnz);
    std::vector<double> di(nnz);
    gop.toTriplet(U,END,U,END,ri.begin(),ci.begin(),di.begin(),false);
    
    std::ofstream out("result.m");
    out << "function [A] = result()\n"
        << "ri = [\n";
    std::copy(ri.begin(),ri.end(),std::ostream_iterator<size_t>(out,"\n"));
    out << "];\nci = [\n";
    std::copy(ci.begin(),ci.end(),std::ostream_iterator<size_t>(out,"\n"));
    out << "];\ndi = [\n";
    std::copy(di.begin(),di.end(),std::ostream_iterator<double>(out,"\n"));
    out << "];\n"
        << "A = spconvert([ri+1 ci+1 di]);\n";
    


    std::vector<RT> data(n+nl,0);
    std::copy(x,x+n,data.begin());
    if (lambda)
      std::copy(lambda,lambda+nl,data.begin()+n);
    xx.read(data.begin());
    std::cout << "\nControl: ";
    std::copy(x,x+24,std::ostream_iterator<double>(std::cout," "));
    std::cout << "\n";
    
    writeVTKFile(varSet,xx,"solution");
  }
  
    

  // /** Intermediate Callback method for the user.  Providing dummy
  //  *  default implementation.  For details see IntermediateCallBack
  //  *  in IpNLP.hpp. */
  // virtual bool intermediate_callback(AlgorithmMode mode,
  //                                    Index iter, Number obj_value,
  //                                    Number inf_pr, Number inf_du,
  //                                    Number mu, Number d_norm,
  //                                    Number regularization_size,
  //                                    Number alpha_du, Number alpha_pr,
  //                                    Index ls_trials,
  //                                    const IpoptData* ip_data,
  //                                    IpoptCalculatedQuantities* ip_cq);
  

  /**
   * The current iterate.
   */
  typename VarSet::VariableSet xx;
  
private:
  void update(bool new_x, Number const* x, bool new_lambda, Number const* lambda) 
  {
    if (new_x || new_lambda) {
      std::vector<RT> data(n+nl);
      xx.write(data.begin());
      
      if (new_x) {
        assert(x);
        std::copy(x,x+n,data.begin());
      }
      
      if (new_lambda) {
        assert(lambda);
        std::copy(lambda,lambda+nl,data.begin()+n);
      }
      
      xx.read(data.begin());
      
      // it's more efficient to simultaneously assemble functional/rhs/matrix
      boost::timer timer;
      gop.assemble(linearization(f,xx));
    }
  }
  

    
  VarSet const& varSet;
  Gop gop;  
  Functional& f;
  int nu, ny, nl, n;
};



#endif
