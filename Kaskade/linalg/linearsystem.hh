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

#ifndef LINEARSYSTEM_HH
#define LINEARSYSTEM_HH

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "linalg/triplet.hh"

namespace Kaskade 
{
/** \ingroup linalg
 *
 *  \brief Abstract base class for a sparse linear system
 */
class SparseLinearSystem
{
public:
  /// Return the matrix of the linear system in triplet format 
  virtual void getMatrix(MatrixAsTriplet<double>& mat) const { getMatrixBlocks(mat,0,nRowBlocks(),0,nColBlocks());}

  /// Return the right hand side of the linear system
  virtual void  getRHS(std::vector<double>& rhs) const { getRHSBlocks(rhs,0,nRowBlocks()); }

  /// Return the number of variables of the linear system
  virtual int size() const {return cols(0,nColBlocks());};

  virtual int rows(int rbegin, int rend) const = 0;
  virtual int cols(int colbegin, int colend) const = 0;

  /// Return matrix blocks of the linear system in triplet format 
  virtual void getMatrixBlocks(MatrixAsTriplet<double>& mat, int rbegin, int rend, int colbegin, int colend) const = 0;

  /// value of function
  virtual double getValue() const = 0;

  /// Return components of the right hand side of the linear system
  virtual void getRHSBlocks(std::vector<double>& rhs, int rbegin, int rend) const = 0;

  /// number of column blocks
  virtual int nColBlocks() const = 0;

  /// number of row blocks
  virtual int nRowBlocks() const = 0;

  virtual ~SparseLinearSystem() {};
};

/** \ingroup linalg
 *
 *  \brief Simple Implementation for SparseLinearSystem
 */
class SimpleSparseLinearSystem : public SparseLinearSystem
{
public:
/// Construction
  SimpleSparseLinearSystem(std::vector<int>const& ridx, std::vector<int>const& cidx, std::vector<double>const& data, std::vector<double>const& rhs_)
  {
    mat.ridx = ridx;
    mat.cidx = cidx;
    mat.data = data;
    rhs = rhs_;
  }

  virtual void getMatrixBlocks(MatrixAsTriplet<double>& mat_, int rbegin, int rend, int colbegin, int colend) 
    const {mat_=mat;}
  virtual void getRHSBlocks(std::vector<double>& rhs_, int rbegin, int rend) const { rhs_=rhs;}
  virtual int  cols(int colbegin, int colend) const {return mat.ncols(); }
  virtual int  rows(int rowbegin, int rowend) const {return mat.nrows(); }
  
  virtual double getValue() const {return 0;}

  virtual int nColBlocks() const {return 1;}
  virtual int nRowBlocks() const {return 1;}

private:
  MatrixAsTriplet<double> mat;
  std::vector<double> rhs;
};

namespace Bridge {

  class UnknownLinearization;

/// Traits class to choose the right linearization class
  template<class Vector, class Functional>
  class LinearizationTraits
  {
  public:
    typedef UnknownLinearization Linearization;
  };
}
} // end of namespace Kaskade

#endif
