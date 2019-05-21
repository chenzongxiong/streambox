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

#include <iostream>
#include "kaskopt.hh"
#include <boost/property_tree/info_parser.hpp>

int main(int argc, char *argv[])
  {
    int verbosity = 1;
    std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosity, true);

    std::cout << "Start property tree test program (r." << pt->get<std::string>("run.version") << ")" << std::endl;

	write_info(std::cout, *pt);

    std::cout << "End property tree test program" << std::endl;

  }
