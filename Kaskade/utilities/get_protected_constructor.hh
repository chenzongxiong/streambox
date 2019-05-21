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
/*
 * get_protected_constructor.hh
 *
 *  Created on: 06.01.2012
 *      Author: Lars Lubkoll
 */

#ifndef GET_PROTECTED_CONSTRUCTOR_HH_
#define GET_PROTECTED_CONSTRUCTOR_HH_

/// A wrapper class for access to protected constructors
template <class Type>
class GetProtectedConstructor  : public Type{
public:
  GetProtectedConstructor(Type const& type) : Type(type){}
};

#endif /* GET_PROTECTED_CONSTRUCTOR_HH_ */
