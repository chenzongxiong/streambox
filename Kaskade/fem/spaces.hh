/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SPACES_HH
#define SPACES_HH

/**
 * \file
 * \brief  Convenience typedefs and creation functions for common function spaces
 * 
 * This is pure convenience. If compile time is the dominant concern, consider whether using 
 * FEFunctionSpace and the required mappers directly, and include only one of fem/hierarchicspace.hh
 * and fem/lagrangespace.hh as required.
 */

#include "fem/functionspace.hh"
#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"

namespace Kaskade
{
  /** 
   * \ingroup fem
   * \brief An \f$ H^1 \f$ conforming finite element space on the leaf grid view with Lagrange basis.
   */
  template <class Grid, class Scalar=double>
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<Scalar,typename Grid::LeafGridView>>;
  
  /**
   * \ingroup fem
   * \brief An \f$ H^1 \f$ conforming finite element space on the leaf grid view with hierarchic basis.
   */
  template <class Grid, class Scalar=double>
  using H1HierarchicSpace = FEFunctionSpace<ContinuousHierarchicMapper<Scalar,typename Grid::LeafGridView>>;
  
  /**
   * \ingroup fem
   * \brief An \f$ H^1 \f$ conforming finite element extension space on the leaf grid view with hierarchic basis.
   */
  template <class Grid, class Scalar=double>
  using H1HierarchicExtensionSpace = FEFunctionSpace<ContinuousHierarchicExtensionMapper<Scalar,typename Grid::LeafGridView>>;
  
  /**
   * \ingroup fem
   * \brief An \f$ L^2 \f$ conforming finite element spaceon the leaf grid view  with Lagrange basis.
   */
  template <class Grid, class Scalar=double>
  using L2Space = FEFunctionSpace<DiscontinuousLagrangeMapper<Scalar,typename Grid::LeafGridView>>;
  
  /**
   * \ingroup fem
   * \brief An \f$ L^2 \f$ conforming finite element space on the leaf grid view with hierarchic basis.
   */
  template <class Grid, class Scalar=double>
  using L2HierarchicSpace = FEFunctionSpace<DiscontinuousHierarchicMapper<Scalar,typename Grid::LeafGridView>>;

  /**
   * \ingroup fem
   * \brief An \f$ L^2 \f$ conforming finite element extension space on the leaf grid view with hierarchic basis.
   */
  template <class Grid, class Scalar=double>
  using L2HierarchicExtensionSpace = FEFunctionSpace<DiscontinuousHierarchicExtensionMapper<Scalar,typename Grid::LeafGridView>>;

}

#endif
