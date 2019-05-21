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

#ifndef PARDISO_SOLVE_HH
#define PARDISO_SOLVE_HH


#include <vector>
#include <memory>
/** \ingroup linalg
 * \brief Factorization with DirectType::PARDISO
 *
 * Important: can only be used with an appropriate license file!
 * Set PARDISO_LIC_PATH to the directory where it is located.
 */
class PARDISOFactorization
{
public:

  /**
   * Matrix structures accepted by the solver. Note that for symmetric
   * matrices (SYMMETRIC_POSITIVE_DEFINITE, SYMMETRIC_INDEFINITE) the
   * upper triangular part has to be given.
   */
  enum MatrixStructure {
    STRUCTURALLY_SYMMETRIC = 1,
    SYMMETRIC_POSITIVE_DEFINITE = 2,
    SYMMETRIC_INDEFINITE = -2,
    NONSYMMETRIC = 11 
  };
  

  ///Version of constructor leaving input data in triplet format (aka coordinate format) unchanged.
  /** 
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param mtype_ Matrix Type: several options available: 
   * - 1 real structurally symmetric
   * - 2 real symmetric positive definite
   * - -2 real symmetric indefinite
   * - 11 real nonsymmetric
   * - 3 complex structurally symmetric
   * - 4 complex symmetric positive definite
   * - -4 complex symmetric indefinite
   * - 13 complex nonsymmetric
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  PARDISOFactorization(int n,
                       int mtype_,
                       std::vector<int> const& ridx,
                       std::vector<int> const& cidx,
                       std::vector<double> const& values);
                                     
  ~PARDISOFactorization();

  /**
   * Solves the system for the given right hand side @arg b. @arg x is
   * resized to the number of matrix columns.
   */
  void solve(std::vector<double> const& b,
             std::vector<double>& x) const;
private:
  mutable int N; 
  mutable int n; 
  mutable std::vector<int> Ap, Ai;
  mutable std::vector<double> Az;
  mutable int phase;

  mutable int maxfct, mnum, mtype, nrhs, msglvl, error,  idum;
  mutable double ddum;
  static bool first;
  mutable std::vector<int> iparm;
  mutable std::vector<long int> pt;
};


#endif
