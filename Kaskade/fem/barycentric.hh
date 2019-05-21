/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BARYCENTRIC_HH
#define BARYCENTRIC_HH

#include <cassert>
#include <numeric>   // std::accumulate, with clang++

#include "dune/common/fvector.hh"

namespace Kaskade
{
  /**
   * \brief Computes the barycentric coordinates of a point in the unit simplex.
   * 
   * The barycentric coordinates of \f$ x \f$ are just the cartesian ones with \f$ 1-\|x\|_{\infty}\f$ appended.
   */
  template <class CoordType, int dim>
  Dune::FieldVector<CoordType,dim+1> barycentric(Dune::FieldVector<CoordType,dim> const& x)
  {
    Dune::FieldVector<CoordType,dim+1> zeta;
    zeta[dim] = 1;
    for (int i=0; i<dim; ++i) {
      zeta[i] = x[i];
      zeta[dim] -= x[i];
    }

    // Shape functions vanishing on the boundary will probably rely on
    // one barycentric coordinate being exactly zero, in order not to be
    // affected by large penalty values from Dirichlet boundary
    // conditions. Hence we here enforce exactly zero values for very
    // small barcyentric coordinates.
    for (int i=0; i<=dim; ++i)
      if (std::abs(zeta[i]) < 1e-12) // well, that's probably *exactly* on the boundary
        zeta[i] = 0;

    return zeta;
  }
  
  /**
   * \ingroup fem
   * \brief Converts integer coordinate index to barycentric indices.
   * \tparam dim the spatial dimension
   * \param x the coordinate index
   * \param bsum the barycentric sum
   * 
   * The entries of x need to be between 0 and bsum, and their sum 
   * shall not exceed bsum. Then the returned index contains in the 
   * first dim entries just x, and in the last one bsum-sum.
   */
  template <size_t dim>
  std::array<int,dim+1> barycentric(std::array<int,dim> const& x, int bsum)
  {
    std::array<int,dim+1> zeta;
    std::copy(begin(x),end(x),begin(zeta));
    zeta[dim] = bsum - std::accumulate(begin(x),end(x),0);
    return zeta;
  }
  
} /* end of namespace Kaskade */

#endif
