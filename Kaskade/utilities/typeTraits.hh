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

#ifndef UTILITIES_TYPE_TRAITS
#define UTILITIES_TYPE_TRAITS

namespace Kaskade
{
  /**
   * \ingroup utilities 
   * \brief A convenience template for removing the const qualifier from references and pointers.
   * 
   * Use this instead of 
   * \code
   * const_cast<T&>(t)
   * \endcode
   * in polymorphic lambda functions where the type T is not named (just inferred from auto) 
   * and hence more involved to obtain, e.g. using decltype.
   */
  template <class T>
  T& removeConst(T const& t)
  {
    return const_cast<T&>(t);
  }
  
  template <class T>
  T* removeConst(T const* t)
  {
    return const_cast<T*>(t);
  }
  
}

#endif