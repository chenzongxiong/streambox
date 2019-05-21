/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef KASKADE_DUNE_INTERFACE_HH_
#define KASKADE_DUNE_INTERFACE_HH_

#include <type_traits>

#include <boost/mpl/plus.hpp>

#include "is_function.hh"
#include "member_variable_macro.hh"

namespace Kaskade
{
  /// Get dimension as nested integral static constexpr.
  template <class Type>
  struct Dim
  {  
    struct Dummy
    { 
      static constexpr int dim = 0;
      static constexpr int dimension = 0;
    };
  
    template <class T, bool enable>
    struct EvaluateDim{ static constexpr int value = 0; };
    
    template <class T>
    struct EvaluateDim<T,true>{ static constexpr int value = std::conditional<!IsFunction<decltype(&(T::dim))>::value,T,Dummy>::type::dim; };
    
    template <class T, bool enable>
    struct EvaluateDimension{ static constexpr int value = 0; };
    
    template <class T>
    struct EvaluateDimension<T,true>{ static constexpr int value = std::conditional<!IsFunction<decltype(&T::dimension)>::value,T,Dummy>::type::dimension; };
    
    
    KASKADE_CREATE_MEMBER_NAME_CHECK(dimension, HasDimension)
    KASKADE_CREATE_MEMBER_NAME_CHECK(dim, HasDim)

    static_assert(!HasDimension<Type>::value && !HasDim<Type>::value, "Type must provide one of the integral static constexpr members \'dimension\' (DUNE) or \'dim\' (Kaskade7)!");
 
    static constexpr int value = std::conditional<
                                                  HasDimension<Type>::value && HasDim<Type>::value, 
                                                  EvaluateDimension<Type,true>,
                                                  boost::mpl::plus<
                                                    EvaluateDimension<Type,HasDimension<Type>::value>,
                                                    EvaluateDim<Type,HasDim<Type>::value> 
                                                  > 
                                 >::type::value;
  };

  /// Get scalar type from Dune (-> field_type) or Kaskade (-> Scalar) member
  template <class Type>
  struct GetScalar
  {
    struct TypeNotFound{};

    template <class LocalType>
    static typename LocalType::field_type hasFieldType(typename LocalType::field_type);

    template <class LocalType>
    static TypeNotFound hasFieldType(...);

    template <class LocalType>
    static typename LocalType::Scalar hasScalar(typename LocalType::Scalar);
    
    template <class LocalType>
    static TypeNotFound hasScalar(...);

    template <class LocalType, bool hasScalar_, bool hasFieldType_> struct ExtractScalar;

    template <class LocalType, bool hasFieldType_>
    struct ExtractScalar<LocalType,true,hasFieldType_>
    {
      typedef typename LocalType::Scalar type;
    };
    
    template <class LocalType>
    struct ExtractScalar<LocalType,false,true>
    {
      typedef typename LocalType::field_type type;
    };

    typedef typename ExtractScalar<Type,
                          !std::is_same<decltype(hasScalar<Type>(0)),TypeNotFound>::value,
                          !std::is_same<decltype(hasFieldType<Type>(0)),TypeNotFound>::value
                          >::type type;
    
    static_assert(!std::is_same<type,TypeNotFound>::value, "No scalar type found. Type must provide one of the nested types \'field_type\' (DUNE) or \'Scalar\' (Kaskade7)!");
  };
  
  /**
   * \brief Extracts the scalar field type from linear algebra data types.
   * 
   * Different conventions for specifying the scalar field type as member type alias are in use in Dune (field_type) and Kaskade (Scalar).
   * This type alias extracts the scalar type from any of those.
   */
  template <class Type>
  using ScalarType = typename GetScalar<Type>::type;
  

}

#endif // KASKADE_DUNE_INTERFACE_HH_
