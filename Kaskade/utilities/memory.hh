/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef UTILITIES_MEMORY
#define UTILITIES_MEMORY

namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Creates a dynamically allocated object from an r-value reference.
   *
   * Use this if you get an object of a runtime-polymorphic class hierarchy
   * only as an r-value, e.g. from a makeSomething() function, and you don't
   * care for its actual type:
   *
   * \code
   * std::unique_ptr<Base> = moveUnique(makeComplexDerivedClass(arguments));
   * \endcode
   *
   * \tparam A needs to be move (preferably) or copy constructible.
   */
  template <class A>
  auto moveUnique(A&& a)
  {
    return std::unique_ptr<A>(new A(std::move(a)));
  }

  /**
   * \ingroup utilities
   * \brief Creates a copy of a given object.
   *
   * Use this if you need a temporary copy that can be moved from.
   *
   * \tparam A needs to be copy constructible and (preferably) moveable.
   */
  template <class A>
  A duplicate(A const& a)
  {
    return A(a);
  }
}


#endif
