/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SCHURBLOCKLU_SOLVER_HH
#define SCHURBLOCKLU_SOLVER_HH

/**
 * @file
 * @brief  Routines for the solution of (sparse) linear systems
 * @author Anton Schiela
 */

#include <iostream>
#include <memory> // std::unique_ptr

#include "dune/common/fmatrix.hh"
#include "dune/istl/matrix.hh"

#include "linalg/umfpack_solve.hh"

#include "linalg/linearsystem.hh"
#include "linalg/simpleLAPmatrix.hh"

#include <boost/timer/timer.hpp>

namespace Kaskade
{

/** \ingroup linalg 
 *\ brief Adapter class for DUNE::IterativeSolver
 */


void printVec(std::vector<double> const&v, int vend=1000000);


template<class Factorization=UMFFactorization<double> >
class BlockLUFactorization
{
public:
/// needs a matrix
  static const bool needMatrix = true;
template<class Sys>
BlockLUFactorization(Sys const& lin, int start2, int end2, int start3, int end3)
  {
    MatrixAsTriplet<double> matL;
    lin.getMatrixBlocks(matL,start3,end3,start2,end2);
    matA.resize(0);
    lin.getMatrixBlocks(matA,start2,end2,start2,end2);
    std::cout << "Inner, " << std::flush;
    boost::timer::cpu_timer timer;
    factoredL.reset(new Factorization(matL.nrows(),
                                      2,
                                      matL.ridx, 
                                      matL.cidx, 
                                      matL.data,
                                      MatrixProperties::GENERAL));
    std::cout << "Finished: " << (double)(timer.elapsed().user)/1e9 << " sec." << std::endl;
  }

  void resetBlock22(MatrixAsTriplet<double>const & matA_) { matA = matA_; };
    
  void solve(std::vector<double>const& rhs, std::vector<double>& sol, int nr=1);
  MatrixAsTriplet<double> matA;
  std::unique_ptr<Factorization> factoredL;
};

struct BlockSchurParameters
{
  BlockSchurParameters(bool reg_)
    : refactorizeInner(true), refactorizeOuter(true)
  {
    if(reg_) regularizationMethod=AddId; else regularizationMethod=None;
  }

  BlockSchurParameters(bool reg_, bool innerf_, bool outerf_)
    : refactorizeInner(innerf_), refactorizeOuter(outerf_)
  {
    if(reg_) regularizationMethod=AddId; else regularizationMethod=None;
  }

  bool refactorizeInner, refactorizeOuter;
  typedef enum { None=0, AddId=1, CG = 2} RegularizationMethod;

  RegularizationMethod regularizationMethod;
};


/// Solver, which is especially designed for the hyperthermia planning problem. 
template<class Factorization>
class DirectBlockSchurSolver
{
public:
  int start1, end1, start2, end2, start3, end3;

/// needs a matrix
  static const bool needMatrix = true;

  DirectBlockSchurSolver(bool doregularize = false) : 
    start1(0), end1(1), 
    start2(1), end2(2), 
    start3(2), end3(3),
    report(false), paras(doregularize)
  {}

  void ax(std::vector<double>& sol, std::vector<double>const &r) const;

  void resolve(std::vector<double>& sol, SparseLinearSystem const& lin) const;

  void resolveAdjAndNormal(std::vector<double>& sol, SparseLinearSystem const& lin) const;

  void resolveNormal(std::vector<double>& sol, SparseLinearSystem const& lin,std::vector<double>const *addrhs=0);

  void solve(std::vector<double>& sol,
             SparseLinearSystem const& lin);

  void solveAdjAndNormal(std::vector<double>& sol,
             SparseLinearSystem const& lin);

  void solveTCG(std::vector<double>& sol1, std::vector<double>& sol2,
                SparseLinearSystem const& linT, SparseLinearSystem const& linN, std::vector<double>const & normalStep, double nu0);

  /// Solves always exactly
  void setRelativeAccuracy(double) {}

  /// Always exact solution
  double getRelativeAccuracy() {return 0.0;}

  void onChangedLinearization() {flushFactorization(); }

  void flushFactorization() 
  { 
    if(factorization.get() && paras.refactorizeInner) factorization.reset(); 
    mC.setSize(0,0);
    B.resize(0);
    AinvB.resize(0);
    matANormal.flush();
  }

  bool report;

  /// Always exact solution
  double getAbsoluteAccuracy() {return 0.0;}

  bool improvementPossible() { return false; }

  void resetParameters(BlockSchurParameters const& p_) { paras=p_; }

private:
  void resolveN(std::vector<double>& sol, std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const;

  void tsolve(std::vector<double>& sol1,std::vector<double>& sol2, 
              std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const;

  void fwd(std::vector<double>& sol, std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const;
  void bwd(std::vector<double>& sol, std::vector<double>const &x2,std::vector<double>const &s,std::vector<double>const &t) const;


  void buildNewSchurComplement(SparseLinearSystem const& lin,int task);

  std::unique_ptr<BlockLUFactorization<Factorization> > factorization;
  MatrixAsTriplet<double> matANormal;
  Dune::Matrix<Dune::FieldMatrix<double,1,1> > mC,mCNormal;
  std::vector<double> B,AinvB;

  int rowsB, colsBC, rowsC;
  int rows1, rows2, rows3;

  BlockSchurParameters paras;

};


class ModifiedSparseSystem : public SparseLinearSystem
{
public:
  ModifiedSparseSystem(SparseLinearSystem const &lin_, MatrixAsTriplet<double> const& mat_, MatrixAsTriplet<double> const& mat2_
                       , std::vector<double> const& scaling_) :
    scaling(scaling_), lin(&lin_), mat(mat_), mat2(mat2_)
  {
  }

  virtual int rows(int rbegin, int rend) const { return lin->rows(rbegin,rend);}
  virtual int cols(int colbegin, int colend) const { return lin->cols(colbegin,colend);}

  /// Return matrix blocks of the linear system in triplet format 
  virtual void getMatrixBlocks(MatrixAsTriplet<double>& m, int rbegin, int rend, int colbegin, int colend) const
  {
    lin->getMatrixBlocks(m,rbegin,rend,colbegin,colend);
    MatrixAsTriplet<double> mm(mat);
    mm *= -1.0;
    if(rbegin==0 && colbegin == 0) m+=mm;
  }

  void resetLin(SparseLinearSystem const& lin_)
  {
    lin=&lin_;
  }

  /// value of function
  virtual double getValue() const { return lin->getValue();}

  /// Return components of the right hand side of the linear system
  virtual void getRHSBlocks(std::vector<double>& rhs, int rbegin, int rend) const 
  {
    rhs.resize(0);
    lin->getRHSBlocks(rhs,rbegin,rend);
    if(rbegin==0)
    {
      std::vector<double> t(lin->rows(3,4),0.0);
      std::vector<double> Mt(lin->rows(0,1),0.0);
      lin->getRHSBlocks(t,3,4);
      for(int i=0; i<t.size();++i)
      {
        t[i] *=scaling[i];
      }
      mat2.ax(Mt,t);
      for(int i=0; i<Mt.size();++i)
      {
//        std::cout << Mt[i] << std::endl;
      }
      for(int i=0; i<lin->rows(0,1); ++i)
        rhs[i] -= Mt[i];
    }
  }


  /// number of column blocks
  virtual int nColBlocks() const { return lin->nColBlocks();};

  /// number of row blocks
  virtual int nRowBlocks() const { return lin->nRowBlocks();};
private:
  std::vector<double> scaling;
  SparseLinearSystem const* lin;
  MatrixAsTriplet<double> mat;
  MatrixAsTriplet<double> mat2;
};


/// Solver, which is especially designed for the hyperthermia planning problem with amplitude ratio
template<class Factorization>
class ARDirectBlockSchurSolver
{
public:
/// needs a matrix
  static const bool needMatrix = true;

  ARDirectBlockSchurSolver(bool doregularize = false) : report(false), DBSSolver(doregularize), justsolved(false) 
{
}

  void resolve(std::vector<double>& sol, SparseLinearSystem const& lin) const;

  void solve(std::vector<double>& sol,
             SparseLinearSystem const& lin);

  /// Solves always exactly
  void setRelativeAccuracy(double) {}

  /// Always exact solution
  double getRelativeAccuracy() {return 0.0;}

  void onChangedLinearization() {flushFactorization(); }

  void flushFactorization() 
  { 
    DBSSolver.flushFactorization();
    F.setSize(0,0);
    FT.setSize(0,0);
    FTVinv.setSize(0,0);
    FTVinvF.setSize(0,0);
    Vinv.setSize(0,0);
  }

  bool report;

  /// Always exact solution
  double getAbsoluteAccuracy() {return 0.0;}

  bool improvementPossible() { return false; }


private:
  std::vector<double> scaling;
  DirectBlockSchurSolver<Factorization> DBSSolver;
  Dune::Matrix<Dune::FieldMatrix<double,1,1> > Vinv,F,FT, FTVinv, FTVinvF;
  std::unique_ptr<MatrixAsTriplet<double> > L,FTVi;
  std::unique_ptr<ModifiedSparseSystem> linMod;
  bool justsolved;
};

}  // namespace Kaskade
#endif
