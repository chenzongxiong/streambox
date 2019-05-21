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
#include "variationalproblemutil.hh"

struct FirstLess 
{
  template <class T>
  bool operator()(std::pair<T,T> const& p) const { return p.first < p.second; }
};



/**
 * Ipopt interface implementation for optimal control problems. This
 * adapts Galerkin operator representations from variational
 * functionals to Ipopt.  Since Ipopt requires not only value and
 * derivatives of the Lagrange functional, the primal variables are
 * split into two quasi-identical groups: those that enter the cost
 * functional (u,y) and those that enter the constraints (v,z).
 */
template <class GOP>
class TrackingControlInterface: public Ipopt::TNLP 
{
  typedef GOP Gop;
  typedef typename Gop::RT RT;
  typedef typename Gop::VariableSet VarSet;
  typedef typename VarSet::Grid Grid;
  typedef typename Gop::Functional Functional;
  
  typedef Ipopt::Index Index;
  typedef Ipopt::TNLP::IndexStyleEnum IndexStyleEnum;
  typedef Ipopt::SolverReturn SolverReturn;
  typedef Ipopt::Number Number;
  

  enum { U=0, Y, V, Z, LAMBDA, END};
  typedef int AssemblyParts;
  
  
  
public:
  TrackingControlInterface(VarSet const& varSet_,
                           Functional& f_):
    varSet(varSet_),
    gop(varSet),
    f(f_),
    nu(varSet.dimension(U,U+1)),
    ny(varSet.dimension(Y,Y+1)),
    nl(varSet.dimension(LAMBDA,LAMBDA+1)),
    n(nu+ny)
  {}
  
  
  /** overload this method to return the number of variables
   *  and constraints, and the number of non-zeros in the jacobian and
   *  the hessian. The index_style parameter lets you specify C or Fortran
   *  style indexing for the sparse matrix iRow and jCol parameters.
   *  C_STYLE is 0-based, and FORTRAN_STYLE is 1-based.
   */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) 
  {
    std::cout << "nu=" << nu << " ny=" << ny << " nl=" << nl << '\n';

    n = this->n;
    m = nl;

    nnz_jac_g = gop.nnz(LAMBDA,LAMBDA+1,V,Z+1,false);
    
    int nnz1 = gop.nnz(U,Y+1,U,Y+1,true);
    int nnz2 = gop.nnz(V,Z+1,V,Z+1,true);
    int nnzMax = std::max(nnz1,nnz2);
    
    std::vector<int> rows(nnzMax), cols(nnzMax);
    std::vector<RT> data(nnzMax);
    std::vector<std::pair<int,int> > entries(nnz1+nnz2);
    
    gop.toTriplet(U,Y+1,U,Y+1,rows.begin(),cols.begin(),data.begin(),true);
    for (int i=0; i<nnz1; ++i)
      entries[i] = std::make_pair(rows[i],cols[i]);
    gop.toTriplet(V,Z+1,V,Z+1,rows.begin(),cols.begin(),data.begin(),true);
    for (int i=0; i<nnz2; ++i)
      entries[i+nnz1] = std::make_pair(rows[i],cols[i]);

    // remove upper diagonal entries
    entries.erase(std::remove_if(entries.begin(),entries.end(),FirstLess()),entries.end());
    // remove duplicate entries
    std::sort(entries.begin(),entries.end());
    entries.erase(std::unique(entries.begin(),entries.end()),entries.end());
    

    nnz_h_lag = entries.size();

    index_style = C_STYLE;

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
      std::fill_n(x,n,0);
      for (size_t i=0; i<nu; ++i)
        x[i] = static_cast<Number>(std::rand())/RAND_MAX;
    }
    if (init_z) {
      std::fill_n(z_L,n,1);
      std::fill_n(z_U,n,1);
    }
    if (init_lambda) {
      std::fill_n(lambda,m,1);
    }
    

    return true;
  }
  

  /** overload this method to return the value of the objective function */
  virtual bool eval_f(Index n, const Number* x, bool new_x,
                      Number& obj_value) 
  {
    std::vector<Number> lambda(nl,0);
    assemble(new_x,x,true,&lambda[0],GOP::VALUE);

    obj_value = gop.functional();
    
    return true;
  }
  

  /** overload this method to return the vector of the gradient of
   *  the objective w.r.t. x */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                           Number* grad_f) 
  {
    std::vector<Number> lambda(nl,0);
    assemble(new_x,x,true,&lambda[0],GOP::RHS);

    ToSequence<U,LAMBDA>::call(gop,grad_f);

    return true;
  }
  

  /** overload this method to return the vector of constraint values */
  virtual bool eval_g(Index n, const Number* x, bool new_x,
                      Index m, Number* g) 
  {
    assemble(new_x,x,false,0,GOP::RHS);

    ToSequence<LAMBDA,END>::call(gop,g);

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
      assert(x);
      assert(values);
      
      assemble(new_x,x,false,0,GOP::MATRIX);

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
      if (x || lambda) {
        // if (lambda) {
        //   std::cerr << "set lambda: ";
        //   std::copy(lambda,lambda+m,std::ostream_iterator<Number>(std::cerr," "));
        //   std::cerr << '\n';
        // }
        
        assemble(new_x,x,new_lambda,lambda,GOP::MATRIX,obj_factor);
      } else
        std::cerr << "eval_h called without x or lambda\n";
      
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
    f.vars.read(data.begin());
    writeVTKFile(varSet,f.vars,"solution");
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

private:
  void assemble(bool new_x, Number const* x, bool new_lambda, Number const* lambda, AssemblyParts ap, double sigma = 1.0) 
  {
    if (new_x || new_lambda || (sigma!=f.ipoptSigma)) {
      
      std::vector<RT> data(n+nl);
      f.vars.write(data.begin());
      
      if (new_x) {
        assert(x);
        std::copy(x,x+n,data.begin());
        // std::cout << "u = " << x[0] << '\n';
        
      }
      
      if (new_lambda) {
        assert(lambda);
        std::copy(lambda,lambda+nl,data.begin()+n);
      }
      
      f.vars.read(data.begin());
      f.ipoptSigma = sigma;
      
      // it's more efficient to simultaneously assemble functional/rhs/matrix
      boost::timer timer;
      gop.assemble(f);

      // std::cout << "\nAssembled (sigma=" << sigma << "):\nu=";
      // std::copy(data.begin(),data.begin()+nu,std::ostream_iterator<RT>(std::cout," ")); 
      // std::cout << "\ny="; std::copy(data.begin()+nu,data.begin()+n,std::ostream_iterator<RT>(std::cout," ")); 
      // std::cout << "\nl="; std::copy(data.begin()+n,data.begin()+n+nl,std::ostream_iterator<RT>(std::cout," "));
      // toSequence<0,3>(gop,data.begin());
      // std::cout << "\nLagrange gradient:\ndu="; std::copy(data.begin(),data.begin()+nu,std::ostream_iterator<RT>(std::cout," ")); 
      // std::cout << "\ndy="; std::copy(data.begin()+nu,data.begin()+n,std::ostream_iterator<RT>(std::cout," ")); 
      // std::cout << "\ndl="; std::copy(data.begin()+n,data.begin()+n+nl,std::ostream_iterator<RT>(std::cout," "));
      // std::cout << "\n";      

      //std::cout << "assembly time: " << timer.elapsed() << "s\n";
      // std::cout << "As. u=[";
      // std::cout.precision(2);
      // std::copy(x+ny,x+ny+nu,std::ostream_iterator<Number>(std::cout,","));
      // std::cout << "]\n";
    }
  }
  

    
  VarSet const& varSet;
  Gop gop;
  Functional& f;
  int nu, ny, nl, n;
};



#endif
