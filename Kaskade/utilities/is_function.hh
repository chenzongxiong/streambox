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
/*
 * is_function.hh
 *
 *  Created on: 20.04.2012
 *      Author: Lars Lubkoll
 */

#ifndef IS_FUNCTION_HH_
#define IS_FUNCTION_HH_

namespace Kaskade
{
  /* * * * * * * * * * * * * * * * * no function * * * * * * * * * * * * * * * */
  /// Checks if type is a function or member function
  /**
   * static member variables:
   * \return isMemberFunction
   * \return value true if T is (member) function, else false
   */
  template <typename T>
  struct IsFunction
  {
    static bool const isMemberFunction = false;
    static bool const value = false;
  };

  /* * * * * * * * * * * * * * * * member function * * * * * * * * * * * * * * */
  /**
   * \cond internals
   */
  template <typename ReturnType, typename ClassName, typename... Arguments>
  struct IsFunction<ReturnType (ClassName::*)(Arguments...)>
  {
    static bool const isMemberFunction = true;
    static bool const value = true;
  };

  template <typename ReturnType, typename ClassName, typename... Arguments>
  struct IsFunction<ReturnType (ClassName::*)(Arguments...) const>
  {
    static bool const isMemberFunction = true;
    static bool const value = true;
  };

  template <typename ReturnType, typename ClassName, typename... Arguments>
  struct IsFunction<ReturnType (ClassName::*)(Arguments...) volatile>
  {
    static bool const isMemberFunction = true;
    static bool const value = true;
  };

  template <typename ReturnType, typename ClassName, typename... Arguments>
  struct IsFunction<ReturnType (ClassName::*)(Arguments...) const volatile>
  {
    static bool const isMemberFunction = true;
    static bool const value = true;
  };
  /* * * * * * * * * * * * * * * *  free function  * * * * * * * * * * * * * * */
  template <typename ReturnType, typename... Arguments>
  struct IsFunction<ReturnType (*)(Arguments...)>{
    static bool const isMemberFunction = false;
    static bool const value = true;
  };
  /**
   * \endcond
   */
}

#endif /* IS_FUNCTION_HH_ */
