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

#ifndef COMBINATORICS_HH
#define COMBINATORICS_HH

#include <algorithm>
#include <cassert>
#include <functional>

namespace Kaskade
{
  /**
   * \brief Computes the binomial coefficient \f$ \binom{n}{k} \f$.
   */
  inline int binomial(int n, int k)
  {
    assert(0<=k && k<=n);

    if (2*k>n)
      k = n-k;

    if (k==0) return 1;

    // Use the product form for efficiency and overflow avoidance.
    int r = n;
    for (int i=2; i<=k; ++i) {
      r *= n+1-i;
      r /= i;
    }

    return r;
  }

  /**
   * \brief Computes the next nonnegative multiindex of order m. 
   * 
   * Using this function one can cycle through all integer tuples that sum up to at
   * most m.
   */
  template <class It>
  void next_multiindex(It first, It last, int const m)
  {
    // If the sum limit is reached, we need to "make room" for an increment.
    if (std::accumulate(first,last,0)==m) {
      // Find the first (i.e. least significant) nonzero entry and set it to zero.
      first = std::find_if(first,last,std::bind2nd(std::not_equal_to<int>(),0));
      assert(first!=last);
      *first = 0;
      // Move to the next higher (i.e. more significant) entry.
      ++first;
    }

    // Increment the relevant index.
    if (first!=last)
      ++(*first);
  }
} /* end of namespace Kaskade */


#endif
