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

#ifndef PARTIAL_DIRECT_PRECONDITIONER_HH
#define PARTIAL_DIRECT_PRECONDITIONER_HH

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvercategory.hh>

#include "fem/istlinterface.hh"
#include "utilities/duneInterface.hh"

namespace Kaskade
{

  /**
   * \ingroup linalg
   * \brief DEPRECATED A partial direct preconditioner applicable to assembled operators.
   *
   * THIS PRECONDITIONER IS DEPRECATED. PREFER A P-MULTIGRID PRECONDITIONER, WHICH IS MUCH MORE
   * EFFECTIVE FOR HIGHER POLYNOMIAL DEGREES AND WORKS WITH SPACES OTHER THAN HIERARCHICAL LAGRANGE.
   *
   * The specified direct solver is used for preconditioning
   * the submatrix with the given range, and a Jacobi preconditioner is
   * used for the remaining variables.
   *
   * This approach is effective for elliptic problems discretized with hierarchical FE spaces
   * (e.g., using \ref ContinuousHierarchicMapper) where the first order ansatz functions
   * span the usual linear FE space. This subspace forms a rather good coarse space with
   * relatively small dimension (in 2D with quadratic elements 1/4, for 3D and/or higher
   * polynomial degree less). The small dimension of the linear FE subspace reduces the
   * direct solver fill-in and hence its memory consumption significantly.
   *
   * The quality of the preconditioner deteriorates with increasing polynomial degree of the
   * ansatz functions. Thus, the preconditioner is not asymptotically optimal for increasing degree.
   * Nevertheless it is very effective in practice for moderate ansatz order (2-5). A (small) benefit
   * is that its relative cost decreases with increasing order.
   *
   * \tparam Op the assembled operator type, a AssembledGalerkinOperator<...> type.
   */
  template <class Op>
  class PartialDirectPreconditioner: public Dune::Preconditioner<typename Op::Domain, typename Op::Range>
  {
    typedef typename Op::Domain Domain;
    typedef typename Op::Domain Range;
    typedef typename GetScalar<Domain>::type Scalar;

  public:
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \arg op the assembled operator
     * \arg first,last a half-open index range denoting the range where the direct solver will be applied.
     * \arg directType the type of direct solver to be used
     * \arg matprop the properties of the matrix
     */
    PartialDirectPreconditioner(Op const& op, size_t first, size_t last,
                                DirectType directType=DirectType::UMFPACK, MatrixProperties matprop=MatrixProperties::SYMMETRIC)
    {
      MatrixAsTriplet<Scalar> A = op.template get<MatrixAsTriplet<Scalar> >();

      // Delete all the off-diagonal entries whose indices are not
      // contained in [first,last)x[first,last).
      size_t out = 0;
      for (size_t in=0; in<A.nnz(); ++in)
        if (A.ridx[in]==A.cidx[in] || (A.ridx[in]>=first && A.ridx[in]<last && A.cidx[in]>=first && A.cidx[in]<last)) {
          A.ridx[out] = A.ridx[in];
          A.cidx[out] = A.cidx[in];
          A.data[out] = A.data[in];
          ++out;
        }
      A.ridx.erase(A.ridx.begin()+out,A.ridx.end());
      A.cidx.erase(A.cidx.begin()+out,A.cidx.end());
      A.data.erase(A.data.begin()+out,A.data.end());

      solver.reset(getFactorization(directType,matprop,A).release());
    }

    virtual void pre(Domain&, Range&) {}
    virtual void post (Domain&) {}

    virtual void apply (Domain& x, Range const& y) {
      std::vector<Scalar> b(y.dim());
      y.write(b.begin());
      solver->solve(b);
      x.read(b.begin());
    }

  private:
    std::unique_ptr<Factorization<Scalar> > solver;
  };

} // namespace Kaskade
#endif
