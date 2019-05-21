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

#ifndef DATE_AND_TIME_HH_
#define DATE_AND_TIME_HH_

#include <string>

namespace Kaskade
{
//  std::vector<const char[4]> months = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
//  //std::vector<const char[3]> numericalMonths = { "01", "02" };
//
//
//  std::string numericalMonth(std::string const& month);
//
//  std::string replaceMonth(std::string const& str);

  std::string getDateAndTime();

  std::string getDate();

//  std::string getNumDate();

  /// append current time and date to string str
  std::string& appendDateAndTime(std::string& str);

  /// append current date to string str
  std::string& appendDate(std::string& str);
}
#endif 
