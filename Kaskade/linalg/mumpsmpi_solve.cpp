/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>
#include <memory>

#include "mpi.h"
#include "linalg/mumpsmpi_solve.hh"

namespace Kaskade
{

// Creator for DirectType::MUMPS factorizations to be plugged into the factory.

struct MUMPSCreator: public Creator<Factorization<double,int> > 
{
  MUMPSCreator(bool plugin) 
  {
    if (plugin)
      Factory<DirectType,Factorization<double,int> >::plugin(DirectType::MUMPS,this);
  }
  
  virtual std::unique_ptr<Factorization<double,int> >
  create(Creator<Factorization<double,int> >::Argument a) const 
  {
    MatrixAsTriplet<double,int> const& A = std::get<1>(a);
    assert(A.nrows()==A.ncols());
    MatrixProperties property = std::get<0>(a);

    return std::unique_ptr<Factorization<double,int> >
      (new MUMPSFactorization<double,int>(const_cast<MatrixAsTriplet<double,int>&>(A).nrows(),0,
                                          const_cast<MatrixAsTriplet<double,int>&>(A).ridx,
                                          const_cast<MatrixAsTriplet<double,int>&>(A).cidx,
                                          const_cast<MatrixAsTriplet<double,int>&>(A).data,
                                          property));
  }
};

MUMPSCreator mumpsCreator(true);

#undef ICNTL
#undef JOB_INIT
#undef JOB_END
#undef USE_COMM_WORLD
} // namespace Kaskade 

