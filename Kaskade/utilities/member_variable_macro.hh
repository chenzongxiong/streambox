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
 * member_variable_macro.hh
 *
 *  Created on: 18.04.2012
 *      Author: Lars Lubkoll
 */

#ifndef MEMBER_VARIABLE_MACRO_HH_
#define MEMBER_VARIABLE_MACRO_HH_

#include <boost/type_traits/is_same.hpp>
#include <boost/preprocessor/cat.hpp>

/*
 * Creates template function that checks if its template parameter has a member variable called
 * VARIABLENAME with type VARIABLETYPE. The struct is named NAME.
 */
#define KASKADE_CREATE_MEMBER_VARIABLE_CHECK(VARIABLE_TYPE, VARIABLE_NAME, NAME) template<class T, bool variableExists> \
  struct BOOST_PP_CAT(Kaskade_,BOOST_PP_CAT(NAME,_VariableTypeCheck{)) \
    static bool const value = false; \
  };\
  \
  template<class T> \
  struct BOOST_PP_CAT(Kaskade_,BOOST_PP_CAT(NAME,_VariableTypeCheck))<T,true>{ \
    static bool const value = boost::is_same<VARIABLE_TYPE,decltype(T::VARIABLE_NAME)>::value; \
  };\
  \
template<typename T> \
struct NAME { \
    struct Fallback { \
        VARIABLE_TYPE VARIABLE_NAME; \
    }; \
    \
    struct Derived : T, Fallback { }; \
    template<typename C, C> \
    struct ChT; \
    \
    template<typename C> \
    static std::false_type f(ChT< VARIABLE_TYPE Fallback::*, &C::VARIABLE_NAME >*); \
    template<typename C> \
    static std::true_type f(...); \
    \
    typedef decltype(f<Derived>(0)) type; \
    static bool const value = BOOST_PP_CAT(Kaskade_,BOOST_PP_CAT(NAME,_VariableTypeCheck))<T,type::value>::value; \
};

/*
 * Creates template function that checks if its template parameter has a member variable or function
 * called VARIABLENAME. The signature/type will not be checked.
 * The struct is named NAME.
 */
#define KASKADE_CREATE_MEMBER_NAME_CHECK(VARIABLE_NAME, NAME) template<typename T> \
struct NAME { \
    struct Fallback { \
        bool VARIABLE_NAME; \
    }; \
    \
    struct Derived : T, Fallback { }; \
    template<typename C, C> \
    struct ChT; \
    \
    template<typename C> \
    static std::false_type f(ChT< bool Fallback::*, &C::VARIABLE_NAME >*); \
    template<typename C> \
    static std::true_type f(...); \
    \
    typedef decltype(f<Derived>(0)) type; \
    static bool const value = type::value; \
};

#endif /* MEMBER_VARIABLE_MACRO_HH_ */
