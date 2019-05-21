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

#ifndef SIGN_HH
#define SIGN_HH

/**
 * Computes the sign of a scalar value. T must be comparable (support
 * < and >). There must be an implicit conversion from int to T.
 */
template <class T>
T sign(T x) 
{
  if (x>0) return 1;
  if (x<0) return -1;
  return 0;
}


#endif
