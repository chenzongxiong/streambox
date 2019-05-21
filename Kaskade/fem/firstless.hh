/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FIRSTLESS_HH
#define FIRSTLESS_HH

namespace Kaskade
{
  /**
   * \brief A comparator functor that supports sorting std::pair by their first component.
   */
  struct FirstLess
  {
    template <class Pair>
    bool operator()(Pair const& p1, Pair const& p2) const
    {
      return p1.first < p2.first;
    }
  };
  
  /**
   * \brief A comparator functor that supports sorting std::pair by their first component.
   */
  struct FirstGreater
  {
    template <class Pair>
    bool operator()(Pair const& p1, Pair const& p2) const 
    {
      return p1.first > p2.first;
    }
  };
} /* end of namespace Kaskade */

#endif
 