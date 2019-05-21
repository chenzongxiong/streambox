/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MG_MULTIGRID_HH
#define MG_MULTIGRID_HH

#include <memory>
#include <type_traits>

#include "fem/spaces.hh"
#include "mg/additiveMultigrid.hh"
#include "mg/multiplicativeMultigrid.hh"
#include "utilities/memory.hh"

namespace Kaskade
{

  /**
   * \ingroup multigrid
   * \brief Convenience function for creating a (multiplicative) multigrid preconditioner for elliptic equations.
   *
   * This function chooses sensible default values for the structure of the preconditioner depending in particular on the
   * finite element ansatz order and the grid refinement depth. As the resulting preconditioners have different types
   * depending on runtime information, we return a unique pointer to SymmetricPreconditioner. Use this if you don't have
   * more specific ideas on how you would like the multigrid to work.
   *
   * The function assumes the space consists of globally continuous finite element functions on simplicial grids
   * with scalar shape functions.
   */
  template <class Entry, class Index, class Space>
  auto makeMultigrid(NumaBCRSMatrix<Entry,Index> const& A, Space const& space, bool onlyLowerTriangle=false)
  {
    using Matrix = NumaBCRSMatrix<Entry,Index>;
    using Preconditioner = SymmetricPreconditioner<typename MatrixTraits<Matrix>::NaturalDomain, typename MatrixTraits<Matrix>::NaturalRange>;

    static_assert(Space::continuity >= 0,"makeMultigrid assumes global continuity of finite element functions");

    // TODO: Consider using cheaper additive multigrid for elastomechanics problems (the difference in contraction rates
    //       appears to be smaller for elastomechanics, and hence the cheaper application pays off).

    // Choose number of pre- and post-smoothings depending on the number of components. Scalar (Laplace) problems tend to
    // have a better contraction, and hence more smoother iterations pay off compared to, say, elastomechanics. In any case
    // we use the same number of pre- and post-smoothings.
    int const nSmooth = Entry::rows > 1? 2: 3;

    int const order = space.mapper().maxOrder();
    if (order == 1)  // these are (assumed to be) P1 elements
      return std::unique_ptr<Preconditioner>(moveUnique(makeJacobiMultiGrid(A,space.gridManager(),nSmooth,nSmooth,onlyLowerTriangle)));

    if (order == 2)  // P2 elements. For those, a Jacobi smoother is cheaper than and almost as efficient as a patch smoother.
    {
      auto mg = makeJacobiPMultiGrid(A,space,nSmooth,nSmooth,onlyLowerTriangle);
      mg.setSmoothings(3,3);
      return std::unique_ptr<Preconditioner>(moveUnique(std::move(mg)));
    }

    // Higher order elements. Here, a patch smoother is much more efficient than a Jacobi smoother
    // (at least from order 4 on, order 3 is probably break even).
    // The patch smoother is quite effective, but rather expensive. Hence we use a smaller number of smoothing iterations than
    // for the Jacobi smoother in the h-multigrid part below.
    auto mg = makeBlockJacobiPMultiGrid(A,space,nSmooth,nSmooth,onlyLowerTriangle);
    mg.setSmoothings(3,3);
    return std::unique_ptr<Preconditioner>(moveUnique(std::move(mg)));
  }

  // ---------------------------------------------------------------------------------------------------------


}

#endif
