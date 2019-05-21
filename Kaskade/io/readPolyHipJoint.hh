/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef READPOLY_HH
#define READPOLY_HH

#include <fstream>
#include <memory>
#include <string>
#include <array>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/config.h>
#undef HAVE_UG
#define HAVE_UG 1
#include "dune/grid/uggrid.hh"

/**
 * \brief Construct an unstructured UG grid from .node and .ele files (Triangle output) and offers
 * the functionality of rotating a grid before handing it over to the grid manager. Obstacle or mortar
 * information (depending on mortarCheck) can be written into obsField (needed for FuFEM contact solver).
 *
 * \param name the base name of the files (without the .node/.ele extension)
 * \param heapSize UG heap size (default 500)
 * \param alpha Rotation angle
 * \param deviation Stores the deviation of the upper boundary of the object after rotation
 * \param barycenter To be determined (and also used later) in order to turn around
 * \param mortarCheck Tells the programm whether we look at an obstacle or mortar object (needed for FuFEM)
 * \param obsField stores whether a vertex belongs to an obstical field or not (needed for FuFEM, boundary contact problem)
 * \param pelvis stores barycenter and radius information (possibly to be deleted later on)
 */
std::unique_ptr<Dune::UGGrid<2> > readPolyData(std::string const& name, int heapSize=500, double alpha = 0, std::array<double,2>* barycenter = nullptr,
                                               double* deviation = nullptr, int mortarCheck = 0,
                                               Dune::BitSetVector<1>* obsField = nullptr, std::array<double,3>* pelvis=nullptr);

#endif
