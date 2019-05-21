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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "dmumps_c.h"
#define ICNTL(I) icntl[(I)-1]
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#include "linalg/factorization.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{

/**
 * \ingroup direct
 * \brief Factorization of sparse linear systems with mumps
 *
 */
template <class Scalar,class SparseIndexInt=int, class DIL=int>
class MUMPSFactorization: public Factorization<Scalar,SparseIndexInt>
{
public:

  /**
   * \brief Constructor keeping input data in triplet format (aka coordinate format) constant.
   *
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  MUMPSFactorization(SparseIndexInt n_,
                     SparseIndexInt dummy,
                     std::vector<SparseIndexInt>& ridx,
                     std::vector<SparseIndexInt>& cidx,
                     std::vector<Scalar>& values,
                     MatrixProperties property = MatrixProperties::GENERAL,
                     int verb = 0)
  : N(ridx.size()), n(n_)
  {
    assert(cidx.size()==N && values.size()==N);

    this->setVerbose(verb);
    if (this->getVerbose()>=2)
    {
      const char *s = "UNKNOWN";
      if (property==MatrixProperties::GENERAL) s = "GENERAL";
      else if (property==MatrixProperties::SYMMETRIC) s = "MatrixProperties::SYMMETRIC";
      else if (property==MatrixProperties::SYMMETRICSTRUCTURE) s = "SYMMETRICSTRUCTURE";
      else if (property==MatrixProperties::POSITIVEDEFINITE) s = "MatrixProperties::POSITIVEDEFINITE";
      std::cout << "MUMPS" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl;
    }

    init(property);
    analyze(&ridx[0], &cidx[0]);
    factorize(&values[0]);
  }

  /**
   * \brief Constructor destroying input data before factorization: more memory efficient
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param dummy unused, present for compatibility with DirectType::PARDISO interface. Specifies problem type there
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  MUMPSFactorization(SparseIndexInt n_,
                     SparseIndexInt dummy,
                     std::unique_ptr<std::vector<SparseIndexInt>> ridx,
                     std::unique_ptr<std::vector<SparseIndexInt>> cidx,
                     std::unique_ptr<std::vector<Scalar>> values,
                     MatrixProperties property = MatrixProperties::GENERAL,
                     int verb = 0)
  : N(ridx->size()), n(n_)
  {
    assert(cidx->size()==N && values->size()==N);

    this->setVerbose(verb);
    if (this->getVerbose()>=2)
    {
      char *s = "UNKNOWN";
      if (property==MatrixProperties::GENERAL) s = "GENERAL";
      else if (property==MatrixProperties::SYMMETRIC) s = "MatrixProperties::SYMMETRIC";
      else if (property==MatrixProperties::SYMMETRICSTRUCTURE) s = "MatrixProperties::SYMMETRICSTRUCTURE";
      else if (property==MatrixProperties::POSITIVEDEFINITE) s = "MatrixProperties::POSITIVEDEFINITE";
      std::cout << "MUMPS" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl;
    }

    init(property);
    analyze(&ridx[0], &cidx[0]);
    factorize(&values[0]);
  }


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
    if (Factorization<Scalar>::getVerbose() < 1)        // If verbose reporting is asked for
      mumps_par.ICNTL(3) = 0;                           // set the output stream to stdout
    mumps_par.ICNTL(14) += 5; // 5% additonal workspace
  }

  void analyze(SparseIndexInt *irn, SparseIndexInt *jcn)
  {
    mumps_par.job = 1;
    mumps_par.n = n;
    mumps_par.nz = N;
    for (int k=0; k<N; k++)
    {
      irn[k]++;
      jcn[k]++;
    }
    mumps_par.irn = irn;
    mumps_par.jcn = jcn;
    dmumps_c(&mumps_par);
    if (mumps_par.info[0] < 0)
      throw DirectSolverException("MUMPS analyze failed with error " + std::to_string(mumps_par.info[0]),__FILE__,__LINE__);
  }

  void factorize(Scalar *a)
  {
    mumps_par.job = 2;
    mumps_par.a = a;
    dmumps_c(&mumps_par);
    if (mumps_par.info[0] < 0)
    {
      std::string errmsg;
      switch(mumps_par.info[0])
      {
        case -13:
        {
          size_t entries = mumps_par.info[1]<0? -1000000*mumps_par.info[1]: mumps_par.info[1];
          errmsg = "\nAllocation of " + std::to_string(entries) + " scalar entries (" + std::to_string(entries*sizeof(Scalar)/1024/1024) + "MB) failed.";
          break;
        }
        default:
          break;
      }
      throw DirectSolverException("MUMPS factorization failed with error " + std::to_string(mumps_par.info[0])+errmsg,__FILE__,__LINE__);
    }
  }

  /**
   * \brief Solves the system for the given right hand side @arg b.
   *
   * @arg x is resized to the number of matrix columns.
   */
  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transpose=false) const
  {
    assert(b.size()>=n);
    x.resize(n);
    assert(&x != &b);
    x = b;
    solve(x);
  }

  /**
   * \brief Solves the system for the given right hand side \arg b, which is
   * overwritten with the solution.
   */
  void solve(std::vector<Scalar>& b) const
  {
    mumps_par.job = 3;
    mumps_par.rhs = &b[0];
    dmumps_c(&mumps_par);
    if (mumps_par.info[0] < 0)
      throw DirectSolverException("MUMPS solve failed with error " + std::to_string(mumps_par.info[0]),__FILE__,__LINE__);
  }

private:
  size_t const N;
  size_t const n;
  mutable DMUMPS_STRUC_C mumps_par;
};

}  // namespace Kaskade
#endif
