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

#include <cassert>
#include <memory>

#include "linalg/superlu_solve.hh"

namespace Kaskade
{

// Creator for DirectType::MUMPS factorizations to be plugged into the factory.

struct SUPERLUCreator: public Creator<Factorization<double,int> > 
{
  SUPERLUCreator(bool plugin) 
  {
    if (plugin)
      Factory<DirectType,Factorization<double,int> >::plugin(DirectType::SUPERLU,this);
  }
  
  virtual std::unique_ptr<Factorization<double,int> >
  create(Creator<Factorization<double,int> >::Argument a) const 
  {
    MatrixAsTriplet<double,int> const& A = std::get<1>(a);
    assert(A.nrows()==A.ncols());

    return std::unique_ptr<Factorization<double,int> >(new SUPERLUFactorization<double,int>(A.nrows(),0,A.ridx,A.cidx,A.data));
  }
};

SUPERLUCreator superluCreator(true);

} // End of Kaskade namespace ---------------------------------------------------------------------

