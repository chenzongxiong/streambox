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

#include <cstdint>
#include "check_endianness.hh"

Kaskade::Endian Kaskade::endianness()
{
  uint8_t buffer[4] = {0x00, 0x01, 0x02, 0x03};

  switch (*((uint32_t *)buffer)) 
  {
  case 0x00010203: return Endian::Big;
  case 0x03020100: return Endian::Little;
  case 0x02030001: return Endian::BigWord;
  case 0x01000302: return Endian::LittleWord;
  default:         return Endian::Unknown;
  }

  return Endian::Unknown;
}
