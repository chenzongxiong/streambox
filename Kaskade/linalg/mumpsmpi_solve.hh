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

#ifndef MUMPS_SOLVE_HH
#define MUMPS_SOLVE_HH

#include <vector>
#include <iostream>
#include <memory>

#include "dmumps_c.h"
#define ICNTL(I) icntl[(I)-1]
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#include "linalg/factorization.hh"

namespace Kaskade
{

/**
 * \ingroup linalg
 * \brief Factorization of sparse linear systems with mumps
 *
 */
template <class Scalar,class SparseIndexInt=int, class DIL=int>
class MUMPSFactorization: public Factorization<Scalar,SparseIndexInt>
{
public:

  ///Version of constructor keeping input data in triplet format (aka coordinate format) constant.
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  MUMPSFactorization(SparseIndexInt n_,
                     SparseIndexInt dummy,
                     std::vector<SparseIndexInt> & ridx,
                     std::vector<SparseIndexInt> & cidx,
                     std::vector<Scalar> & values,
                     MatrixProperties property = MatrixProperties::GENERAL,
                     int verb = 0)
  : N(ridx.size()), n(n_)
  {
    assert(cidx.size()==N && values.size()==N);

    id = MPI::COMM_WORLD.Get_rank ( );
    this->setVerbose(verb);
     if ( (this->getVerbose()>=2) && (id == 0) )
      {
        const char *s = "UNKNOWN";
        if (property==MatrixProperties::GENERAL) s = "GENERAL";
        else if (property==MatrixProperties::SYMMETRIC) s = "SYMMETRIC";
        else if (property==MatrixProperties::SYMMETRICSTRUCTURE) s = "SYMMETRICSTRUCTURE";
        else if (property==MatrixProperties::POSITIVEDEFINITE) s = "POSITIVEDEFINITE";
        std::cout << "MUMPS" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl; 
      }

    init(property);
    analyze(&ridx[0], &cidx[0]);
    factorize(&values[0]);
   };

  /// Version of constructor, that destroys input data before factorization: more memory efficient
  /**
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  MUMPSFactorization(SparseIndexInt n_,
                     SparseIndexInt dummy,
                     std::unique_ptr<std::vector<SparseIndexInt> > ridx,
                     std::unique_ptr<std::vector<SparseIndexInt> > cidx,
                     std::unique_ptr<std::vector<Scalar> > values,
                     MatrixProperties property = MatrixProperties::GENERAL,
                     int verb = 0)
  : N(ridx->size()), n(n_)
  {
    assert(cidx->size()==N && values->size()==N);

    id = MPI::COMM_WORLD.Get_rank ( );
    this->setVerbose(verb);
    if ( (this->getVerbose()>=2) && (id == 0) )
      {
        char *s = "UNKNOWN";
        if (property==MatrixProperties::GENERAL) s = "GENERAL";
        else if (property==MatrixProperties::SYMMETRIC) s = "SYMMETRIC";
        else if (property==MatrixProperties::SYMMETRICSTRUCTURE) s = "SYMMETRICSTRUCTURE";
        else if (property==MatrixProperties::POSITIVEDEFINITE) s = "POSITIVEDEFINITE";
        std::cout << "MUMPS" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl; 
      }

    init(property);
    analyze(&ridx[0], &cidx[0]);
    factorize(&values[0]);
  };
                                     

  ~MUMPSFactorization()
  {
    mumps_par.job = JOB_END;
    mumps_par.ICNTL(3) = 0;
    dmumps_c(&mumps_par);
  }

  void init(MatrixProperties property)
    {
    mumps_par.job = JOB_INIT;
    mumps_par.par = 1;
    switch (property)
      {
        case MatrixProperties::GENERAL:
          mumps_par.sym = 0;
          break;
        case MatrixProperties::SYMMETRICSTRUCTURE:
          mumps_par.sym = 0;
          break;
        case MatrixProperties::SYMMETRIC:
          mumps_par.sym = 2;
          break;
        case MatrixProperties::POSITIVEDEFINITE:
          mumps_par.sym = 1;
          break;
      }
    mumps_par.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&mumps_par);
    if (Factorization<Scalar>::getVerbose()<1)
      mumps_par.ICNTL(3) = 0;
    mumps_par.ICNTL(14) += 5; // 5% additonal workspace
    }

  void analyze(SparseIndexInt *irn, SparseIndexInt *jcn)
    {
    mumps_par.job = 1;
    mumps_par.n = n;
    mumps_par.nz = N;
    if ( id == 0 )
    {
      for (int k=0; k<N; k++)
      {
        irn[k]++;
        jcn[k]++;
      };
    };
    MPI_Barrier(MPI_COMM_WORLD);
    mumps_par.irn = irn;
    mumps_par.jcn = jcn;
    dmumps_c(&mumps_par);
    }

  void factorize(Scalar *a)
    {
    mumps_par.job = 2;
    mumps_par.a = a;
    dmumps_c(&mumps_par);
    }

  /**
   * Solves the system for the given right hand side @arg b. @arg x is
   * resized to the number of matrix columns.
   */
  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transpose=false) const
    {
    assert(b.size()>=n);
    x.resize(n);
    assert(&x != &b);
    x = b;
    mumps_par.job = 3;
    mumps_par.rhs = &x[0];
    dmumps_c(&mumps_par);
    }

  /**
   * Solves the system for the given right hand side \arg b, which is
   * overwritten with the solution.
   */
  void solve(std::vector<Scalar>& b) const
  {
    mumps_par.job = 3;
      mumps_par.rhs = &b[0];
      dmumps_c(&mumps_par);
  }
  
private:
  int id;
  size_t const N; 
  size_t const n;
  mutable DMUMPS_STRUC_C mumps_par;
};

}  // namespace Kaskade
#endif
