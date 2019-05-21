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

#ifndef SCHUR_SOLVER_HH
#define SCHUR_SOLVER_HH

/**
 * @file
 * @brief  Routines for the solution of (sparse) linear systems
 * @author Anton Schiela
 */

#include <iostream>
#include <memory> // std::unique_ptr

#include <boost/timer/timer.hpp>


#include "dune/common/fmatrix.hh"
#include "dune/istl/matrix.hh"

#include "linalg/umfpack_solve.hh"
#include "linalg/linearsystem.hh"
#include "linalg/simpleLAPmatrix.hh"
#include "algorithm/opt_interface.hh"
#include "algorithm/newton_bridge.hh"

namespace Kaskade
{

/** 
 * \ingroup linalg 
 * \brief Adapter class for DUNE::IterativeSolver
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
    boost::timer::cpu_timer timer;
    factoredL.reset(new Factorization(matL.nrows(),
                                      2,
                                      matL.ridx, 
                                      matL.cidx, 
                                      matL.data,
                                      MatrixProperties::GENERAL));
  }

  void resetBlock22(MatrixAsTriplet<double>const & matA_) { matA = matA_; };
    
  void solve(std::vector<double>const& rhs, std::vector<double>& sol, int nr=1);
  MatrixAsTriplet<double> matA;
  std::unique_ptr<Factorization> factoredL;
};

void printVec(std::vector<double> const&v, int vend)
{
  int endv=v.size();
  if(vend < v.size()) endv=vend;
  for(int i=0; i< endv; ++i)
  {
    if(i%5 == 0)  std::cout << ".  " << v[i] << std::endl;
    else std::cout << "   " << v[i] << std::endl;

  }
}

template<class Factorization>
void BlockLUFactorization<Factorization>::solve(std::vector<double>const& rhs, std::vector<double>& sol, int nr)
{
  int r=rhs.size()/2/nr;
  std::vector<double> rhs1(r*nr);
  std::vector<double> rhs2(r*nr);
  std::vector<double> sol1(r*nr);
  std::vector<double> sol2(r*nr);
  for(int i=0; i<nr; ++i)
  for(int j=0; j<r; ++j)
  {
    rhs1[i*r+j]=rhs[i*2*r+j];
    rhs2[i*r+j]=-rhs[i*2*r+r+j];
  }
  factoredL->solve(rhs2,sol2,false);  // solve system
  //  factoredL->solve(rhs2,sol2,nr,false);  // solve system
  matA.axpy(rhs1,sol2,1.0,nr);
  factoredL->solve(rhs1,sol1,true);  //solve transposed system
  // factoredL->solve(rhs1,sol1,nr,true);  //solve transposed system
  sol.resize(2*r*nr);
  for(int i=0; i<nr; ++i)
  for(int j=0; j<r; ++j)
  {
    sol[i*2*r+j]=-sol2[i*r+j];
    sol[i*2*r+r+j]=sol1[i*r+j];
  }
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
  typedef enum { None=0, AddId=1, IterateType::CG = 2} RegularizationMethod;

  RegularizationMethod regularizationMethod;
};



template<class Factorization, class VariableSet>
class DirectBlockSchurSolver : public AbstractNormalDirection, public Dune::Preconditioner<typename VariableSet::Descriptions::template CoefficientVectorRepresentation<>::type,
 typename VariableSet::Descriptions::template CoefficientVectorRepresentation<>::type> 
{
public:
  typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<>::type Domain;
  typedef Domain Range;

  int start1, end1, start2, end2, start3, end3;

/// needs a matrix
  static const bool needMatrix = true;

  DirectBlockSchurSolver(bool doregularize = false) : 
    start1(0), end1(1), 
    start2(1), end2(2), 
    start3(2), end3(3),
    report(false), paras(doregularize)
  {}

  virtual void pre(Domain &x, Range &b) {}

  virtual void apply (Domain &v, const Range &d) {
    
    std::vector<double> r(rows1),s(rows2),t(rows3), sol1(rows1+rows2+rows3);

    d.write(sol1.begin());

    for(int i=0; i<rows1; ++i) r[i]=sol1[i];
    for(int i=0; i<rows2; ++i) s[i]=sol1[i+rows1];
    for(int i=0; i<rows3; ++i) t[i]=sol1[i+rows1+rows2];
    
    resolveN(sol1,r,s,t);

    v.read(sol1.begin());
  }
  virtual void post (Domain &x) {}

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

  void resetParameters(BlockSchurParameters const& p_) { paras=p_; }

  virtual void computeCorrectionAndAdjointCorrection(AbstractVector& correction, AbstractVector& adjointCorrection, AbstractLinearization& linearization)
  {
    SparseLinearSystem& lins=dynamic_cast<SparseLinearSystem &>(linearization);
    flushFactorization();
    buildNewSchurComplement(lins);

    std::vector<double> r,s,t, sol1, sol2;
    
    lins.getRHSBlocks(r,start1,end1);
    lins.getRHSBlocks(s,start2,end2);
    lins.getRHSBlocks(t,start3,end3);

    std::vector<double> r0(r.size(),0.0),s0(s.size(),0.0),t0(t.size(),0.0);
    resolveN(sol1,r0,s0,t);
    resolveN(sol2,r,s,t0);

    dynamic_cast<Bridge::Vector<VariableSet>& >(correction).read(sol1);
    dynamic_cast<Bridge::Vector<VariableSet>& >(adjointCorrection).read(sol2);
    correction *= -1.0;
    adjointCorrection *= -1.0;
  }

  virtual void computeSimplifiedCorrection(AbstractVector& correction, AbstractLinearization const& lin) const
  {
    SparseLinearSystem const& lins=dynamic_cast<SparseLinearSystem const&>(lin);

    std::vector<double> t, sol1;
    lins.getRHSBlocks(t,start3,end3);
    std::vector<double> r0(rows1,0.0),s0(rows2,0.0);
    resolveN(sol1,r0,s0,t);
    dynamic_cast<Bridge::Vector<VariableSet>& >(correction).read(sol1);
    correction *= -1.0;
  }

private:

  void resolveN(std::vector<double>& sol, std::vector<double>const &rhs0,std::vector<double>const &rhs1, std::vector<double>const &rhs2) const;

  void fwd(std::vector<double>& sol, std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const;
  void bwd(std::vector<double>& sol, std::vector<double>const &x2,std::vector<double>const &s,std::vector<double>const &t) const;


  void buildNewSchurComplement(SparseLinearSystem const& lin,int task);

  std::unique_ptr<BlockLUFactorization<Factorization> > factorization; //
  MatrixAsTriplet<double> matANormal;
  Dune::Matrix<Dune::FieldMatrix<double,1,1> > mC,mCNormal;
  std::vector<double> B,AinvB;

  int rowsB, colsBC, rowsC;
  int rows1, rows2, rows3;

  BlockSchurParameters paras;

};

//  UU UY BU* ...  C B* B*
//  YU YY AY* ...  B A  A
//  BU AY 00  ...  B A  A
//  .. .. ..  ...
template<class Factorization, class VariableSet>
void DirectBlockSchurSolver<Factorization, VariableSet>::buildNewSchurComplement(SparseLinearSystem const& lin,int task=0)
{
    boost::timer::cpu_timer timer;
    if(report) std::cout << "Schur Complement: " << std::flush;

    MatrixAsTriplet<double> matB, matC;
    rowsB=lin.rows(start2,end3);
    colsBC=lin.cols(start1,end1);
    rowsC=lin.rows(start1,end1);
    if(task==0)
    {
      rows1=lin.rows(start1,end1);
      rows2=lin.rows(start2,end2);
      rows3=lin.rows(start3,end3);
    }

    if(paras.refactorizeOuter || !factorization.get() || B.size()==0)
    {
      if(report) std::cout << "Outer, " << std::flush;
      if(task==0)
      {
        flushFactorization();
      }

      if((task==0 && paras.refactorizeInner) || !factorization.get())
      {
        factorization.reset(new BlockLUFactorization<Factorization>(lin,start2,end2,start3,end3));
        lin.getMatrixBlocks(matANormal,start2,end2, start2, end2);
      }
      else
      {
        MatrixAsTriplet<double> matA;
        lin.getMatrixBlocks(matA,start2,end2, start2, end2);
        factorization->resetBlock22(matA);
      }

  
// Matrix = [C B^T; B A]; rhs=[r2, r1]; sol=[x_2,x_1]
// Compute sol via the Schur complement
// B is stored in column-first formal in a vector B(i,j)=B[i+j*rowsB]
      if(task==0 || B.size()==0)
      {
        matB.resize(0);
        lin.getMatrixBlocks(matB,start2,end3, start1, end1);
        B.resize(0);
        matB.toVector(B);
      }

// Compute A^{-1}B
      AinvB.resize(0);
      factorization->solve(B,AinvB,colsBC);

    }


// Compute C as a full matrix
    matC.resize(0);
    lin.getMatrixBlocks(matC,start1,end1, start1, end1);
    mC.setSize(rowsC, colsBC);
    mC=0.0;
    matC.addToMatrix(mC);

// Compute S := B^T A^{-1}B-C
    for(int i=0; i<rowsC; ++i)
      for(int j=0; j<colsBC; ++j)
      {
        mC[i][j] *= -1.0;
        for(int k=0; k<rowsB; ++k) mC[i][j] += B[i*rowsB+k]*AinvB[j*rowsB+k];
      }
    if(task==0)
    {
      mCNormal.setSize(rowsC, colsBC);
      mCNormal = mC;
    }
    if(report) std::cout << "Finished: " << (double)(timer.elapsed().user)/1e9 << " sec." << std::endl;
}

template<class Factorization, class VariableSet>
void DirectBlockSchurSolver<Factorization,VariableSet>::resolveN(std::vector<double>& sol, std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const
{
  std::vector<double> xx2;

  fwd(xx2,r,s,t);

  std::vector<double> x2(xx2.size(),0.0);

  // Compute x2 := S^{-1}xx_2

  LeastSquares(SLAPMatrix<double>(mCNormal), xx2, x2);  

  bwd(sol,x2,s,t);
}

template<class Factorization, class VariableSet>
void DirectBlockSchurSolver<Factorization,VariableSet>::fwd(std::vector<double>& sol, std::vector<double>const &r,std::vector<double>const &s,std::vector<double>const &t) const
{
    std::vector<double> xx1, r1;

    r1.reserve(s.size()+t.size());

    for(int i=0; i<s.size();++i) r1.push_back(s[i]);
    for(int i=0; i<t.size();++i) r1.push_back(t[i]);
// Compute xx1=A^{-1}r1

    xx1.resize(r1.size());
    factorization->solve(r1,xx1);

// Compute sol := -r2+B^T xx1

    sol.resize(r.size());
    for(int i=0; i<r.size(); ++i)
    {
      sol[i]=-r[i];
      for(int k=0; k<xx1.size(); ++k)
        sol[i]+=B[i*rowsB+k]*xx1[k];
    }
}

template<class Factorization, class VariableSet>
void DirectBlockSchurSolver<Factorization,VariableSet>::bwd
(std::vector<double>& sol, std::vector<double>const &x2,std::vector<double>const &s,std::vector<double>const &t) const
{
  std::vector<double> xx1,x1;

   xx1.reserve(s.size()+t.size());

    for(int i=0; i<s.size();++i) xx1.push_back(s[i]);
    for(int i=0; i<t.size();++i) xx1.push_back(t[i]);

// Compute xx1 := r1-B x2
    for(int i=0; i<xx1.size(); ++i)
      for(int k=0; k<x2.size(); ++k)
        xx1[i]-=B[k*rowsB+i]*x2[k];

// Compute x1 = A^{-1}xx1
    factorization->solve(xx1,x1);    

// sol=[x2;x1]
    sol.reserve(x2.size()+x1.size());//x2.size()+x1.size());
    sol.resize(x2.size()+x1.size(),0.0);
    int k=0;
    for(int i=0; i< x2.size(); ++i)
    {
      sol[k]=x2[i];
      ++k;
    }
    for(int i=0; i< x1.size(); ++i)
    {
      sol[k]=x1[i];
      ++k;
    }
}


}  // namespace Kaskade
#endif
