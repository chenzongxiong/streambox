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

#ifndef READPOLY_HH
#define READPOLY_HH

#include <memory>
#include <string>

#include <dune/grid/config.h>
#undef HAVE_UG
#define HAVE_UG 1
#include "dune/grid/uggrid.hh"

/**
 * \ingroup gridInput
 * \brief Construct an unstructured UG grid from .node and .ele files (Triangle output).
 *
 * \param name the base name of the files (without the .node/.ele extension).
 * \param heapSize UG heap size (default 500)
 * \param envSize UG environment size (default 10)
 */
std::unique_ptr<Dune::UGGrid<2>> readPolyData(std::string const& name, int heapSize=500, int envSize=10);

#endif
