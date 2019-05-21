/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <ios>
#include <sstream>

#include "io/iobase.hh"

namespace Kaskade 
{
  std::string paddedString(int n, int places)
  {
    std::ostringstream s;
    s.width(places);
    s.fill('0');
    s.setf(std::ios_base::right,std::ios_base::adjustfield);
    s << n;
    s.flush();
    return s.str();
  }
  
  IoOptions ioOptions_default;
}