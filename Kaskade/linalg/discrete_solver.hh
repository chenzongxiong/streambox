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

#ifndef DISCRETE_SOLVER_HH
#define DISCRETE_SOLVER_HH

/**
 * @file
 * @brief  Routines for the solution of (sparse) linear systems
 * @author Anton Schiela
 */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "linalg/factorization.hh"
#include "linalg/linearsystem.hh"
#include "linalg/umfpack_solve.hh"
//#include "linalg/direct_solver_factory.hh"

namespace Kaskade
{
  /** \ingroup linalg
   * \brief Direct solver.Implements Direct Solver Interface of Algorithms
   *
   *
   *  Factorizations:
   * - UMF
   * - DirectType::MUMPS
   * - DirectType::SUPERLU
   */
  //template<class Factorization>
  template <class Scalar, class SparseInt=int>
  class DirectLinearSolver
  {
  public:
    /// needs a matrix
    static const bool needMatrix = true;
    /// default constructor
    explicit DirectLinearSolver(std::string const& solverName_, bool verbose = false) : report(false), solverName(solverName_)/*, factory(verbose) */{}

    ~DirectLinearSolver() { flushFactorization(); }

    /// solve a system, keep factorization
    void solve(std::vector<double>& sol, SparseLinearSystem& lin)
    {
      flushFactorization();
      boost::timer::cpu_timer timer;

      MatrixAsTriplet<double> mat;
      lin.getMatrix(mat);

      //   mat.print();
      factorization.reset(new UMFFactorization<Scalar,int>(mat.nrows(), 1, mat.ridx, mat.cidx, mat.data)/*factory.create(solverName, mat.nrows(), mat.ridx, mat.cidx, mat.data) */);
      //    factorization.reset(new Factorization(mat.nrows(),
      //                                          matrixtype,
      //                                          mat.ridx,
      //                                          mat.cidx,
      //                                          mat.data));
      if(report && (double)(timer.elapsed().user)/1e9 > 2.0) {
        std::cout << "Factorization: " << (double)(timer.elapsed().user)/1e9 << " sec." << std::endl;
      }
      resolve(sol,lin);
    }
    /// solve a system again
    void resolve(std::vector<double>& sol,
        SparseLinearSystem const& lin) const
    {
      boost::timer::cpu_timer timer;
      if(report >= 2) std::cout << "Starting Substitution...";
      std::vector<double> rhs;
      lin.getRHS(rhs);
      factorization->solve(rhs,sol);
      if(report >=2)  std::cout << "Finished: " << (double)(timer.elapsed().user)/1e9 << " sec." << std::endl;

    }

    /// Solves always exactly
    void setRelativeAccuracy(double) {}

    /// Always exact solution
    double getRelativeAccuracy() {return 0.0;}

    /// Always exact solution
    double getAbsoluteAccuracy() {return 0.0;}

    bool improvementPossible() { return false; }


    /// Do sth if linearization changed
    void onChangedLinearization() {flushFactorization(); }

    /// flush factorization
    void flushFactorization()
    {
      if(factorization.get()!=nullptr) factorization.reset();
    }
    /// Report of progress (0=no report, 1=brief, 2=verbose)
    int report;

  private:
    std::string const& solverName;
    //DirectSolverFactory<Scalar,SparseInt> factory;
    std::unique_ptr<Factorization<Scalar,SparseInt> > factorization;
  };
}
#endif
