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

#ifndef TRIVIALPRECONDITIONER_HH
#define TRIVIALPRECONDITIONER_HH

#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvercategory.hh"

namespace Kaskade
{

/**
 * The trivial preconditioner: this is simply the identity that does
 * exactly nothing. Use this as dummy preconditioner in case you don't
 * want preconditioning for some reason.
 */
template <class Operator>
class TrivialPreconditioner: public Dune::Preconditioner<typename Operator::domain_type,typename Operator::range_type>
{
public:
  typedef typename Operator::domain_type domain_type;
  typedef typename Operator::range_type range_type;

  static int const category = Dune::SolverCategory::sequential;
  
  virtual void pre (domain_type&, range_type&) {}
  virtual void apply (domain_type& x, range_type const& y) { x = y; }
  virtual void post (domain_type&) {}
};

} // namespace Kaskade
#endif
