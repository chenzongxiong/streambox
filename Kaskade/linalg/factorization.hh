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

#ifndef FACTORIZATION_HH
#define FACTORIZATION_HH

#include <sstream>
#include <tuple>
#include <vector>

#include <dune/grid/config.h>

#include "utilities/enums.hh"
#include "dune/istl/bcrsmatrix.hh"

#include "linalg/factory.hh"
#include "linalg/triplet.hh"

namespace Kaskade
{

//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Abstract base class for matrix factorizations.
 *
 * \tparam Scalar the underlying field type of the matrix elements.
 * \tparam SparseIndexInt the integral type used for indices
 */
template <class Scalar,class SparseIndexInt=int>
class Factorization 
{
public:
  /**
   * \brief The type of matrix elements (a field type).
   */
  typedef Scalar field_type;

  /**
   * \brief Solves the system \f$ Ax=b \f$ for the given right hand
   * side \f$ b \f$.
   * \arg b is overwritten with the solution \f$ x \f$.
   */
  virtual void solve(std::vector<field_type>& b) const = 0;
  virtual void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transposed=false) const = 0;
  
  Factorization() : verbose(0) {}
  virtual ~Factorization() {}
  void setVerbose(int verbose_) { verbose= verbose_; }
  int getVerbose() { return verbose; }

protected:
  /**
   * \brief Converts a triplet to a compressed column format using umfpack.
   */
  void tripletToCompressedColumn(SparseIndexInt nRows, SparseIndexInt nCols,
    size_t nNonZeros, std::vector<SparseIndexInt> const &ridx,
    std::vector<SparseIndexInt> const &cidx, std::vector<Scalar> const &values,
    std::vector<SparseIndexInt> &Ap, std::vector<SparseIndexInt> &Ai,
    std::vector<Scalar> &Az, std::vector<SparseIndexInt> &map);
  int verbose;
};

//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Abstract base class for factorization creators to be plugged into a factorization factory.
 *
 * Derive from this class if you need to plug a direct linear solver
 * into the factory.
 */
template <class Scalar,class SparseIndexInt>
struct Creator<Factorization<Scalar,SparseIndexInt> > 
{
  typedef std::tuple<MatrixProperties,MatrixAsTriplet<Scalar,SparseIndexInt> const&> Argument;

  virtual std::unique_ptr<Factorization<Scalar,SparseIndexInt> > create(Argument a) const = 0;

  virtual ~Creator() {}
  
  static std::string lookupFailureHint(DirectType direct)
  {
    std::ostringstream hint;
    hint << "Solver key " << (int)direct << " not registered.\n"
         << "Has the corresponding object file been linked in?\n"
         << "Did you specify the right index type (int/long) in your matrix?\n";
    return hint.str();
  }
};


//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Creates a factorization of the given triplet matrix.
 */
template <class Scalar,class SparseIndexInt>
std::unique_ptr<Factorization<Scalar,SparseIndexInt> > getFactorization(DirectType directType, MatrixProperties properties,
                                                                        MatrixAsTriplet<Scalar,SparseIndexInt> const& A) 
{
  assert(A.nrows()==A.ncols());
  if (A.nnz()==0)
    throw SingularMatrixException("No nonzero entries.",__FILE__,__LINE__);
    std::tuple<MatrixProperties,MatrixAsTriplet<Scalar,SparseIndexInt> > args = std::make_tuple(properties,A);
  if (directType == DirectType::ANY)
    return Factory<DirectType,Factorization<Scalar,SparseIndexInt> >::createAny(args);
  else
    return Factory<DirectType,Factorization<Scalar,SparseIndexInt> >::create(directType,args);
}

/**
 * \ingroup direct
 * \brief Creates a factorization of the given BCRSmatrix.
 */
template <class Scalar,class SparseIndexInt=int /* typename Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::size_type*/ >  // TODO: what is the best default?
std::unique_ptr<Factorization<Scalar,SparseIndexInt> >  
getFactorization(DirectType directType, MatrixProperties properties, Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& A) 
{
  return getFactorization(directType,properties,MatrixAsTriplet<Scalar,SparseIndexInt>(A));
}

/**
 * \ingroup direct
 * \brief Creates a factorization of the given NumaBCRS matrix.
 */
template <class Scalar, int n, int m, class Index, class SparseIndexInt=int>  // TODO: what is the best default?
std::unique_ptr<Factorization<Scalar,SparseIndexInt> >  
getFactorization(DirectType directType, MatrixProperties properties, NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Index> const& A) 
{
  return getFactorization(directType,properties,MatrixAsTriplet<Scalar,SparseIndexInt>(A));
}

}  // namespace Kaskade

#endif 
