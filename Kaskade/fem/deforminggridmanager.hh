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
/*
 * deforminggridmanager.hh
 *
 *  Created on: 30.07.2012
 *      Author: Lars Lubkoll
 */

#ifndef DEFORMINGGRIDMANAGER_HH_
#define DEFORMINGGRIDMANAGER_HH_

#include <dune/grid/uggrid.hh>

#include "fem/fetransfer.hh"
#include "utilities/interpolation/boundarynormalcollector.hh"
#include "utilities/interpolation/hermite/hermite.hh"
#include "utilities/interpolation/polynomialcollection.hh"
#include "utilities/interpolation/tools.hh"

namespace Kaskade
{
  /**
   * After grid refinement this grid manager deforms the refined mesh using cell-wise specified third-order
   * hermite interpolation. The required slopes are computed using given or estimated surface normals defined
   * at vertices. This will result in a smoother domain boundary and inner boundaries. Moreover the
   * discretization error may be reduced.
   *
   * There exist three policies that can be used to adjust the DeformingGridManager to specific needs:
   *
   * \param OuterBoundaryPolicy specifies if (possibly parts of) the outer boundary is (are) considered
   * \param InnerBoundaryPolicy specifies if (possibly parts of) inner boundaries are considered
   * \param ThresholdPolicy specify which slopes are admissible (equivalently which given/estimated surface
   * normals are ignored). Use this policy for feature preservation.
   */
  template <class Scalar, int dimension,
  class OuterBoundaryPolicy=Policy::ConsiderOuterBoundary,
  class InnerBoundaryPolicy=Policy::IgnoreInnerBoundary,
  template <class> class ThresholdPolicy=Policy::NoGradientThreshold>
  class DeformingGridManager : public GridManagerBase<Dune::UGGrid<dimension>>
  {
    typedef typename Dune::UGGrid<dimension>::LeafGridView GridView;
    typedef HermiteInterpolation<Scalar,GridView,OuterBoundaryPolicy,ThresholdPolicy> InterpolationPolynomial;
    typedef InterpolationPolynomialCollection<InterpolationPolynomial> Interpolation;

  public:
    /// Constructor.
    /**
     * Use this constructor if you need the grid manager before being able to provide the necessary
     * data for the computation of the deformation (i.e. if you have to construct a function space element for storing phase ids)
     */
    explicit DeformingGridManager(std::unique_ptr<Dune::UGGrid<dimension> >&& grid, bool verbose = false)
    : GridManagerBase<Dune::UGGrid<dimension> >(std::move(grid), verbose), deformation(nullptr)
    {}

    /// Constructor
    /**
     * Use this constructor if no phase information is provided.
     */
    explicit DeformingGridManager(std::unique_ptr<Dune::UGGrid<dimension>>&& grid,
        OuterBoundaryPolicy const& outerBoundaryPolicy,
        InnerBoundaryPolicy const& innerBoundaryPolicy,
        ThresholdPolicy<Scalar> const& thresholdPolicy,
        bool verbose = false)
    : GridManagerBase<Dune::UGGrid<dimension>>(std::move(grid), verbose),
      deformation(new Interpolation(this->gridptr->leafView(), BoundaryNormalCollector<GridView,OuterBoundaryPolicy,InnerBoundaryPolicy>(this->gridptr->leafView(), outerBoundaryPolicy, innerBoundaryPolicy).meanNormals(), outerBoundaryPolicy, thresholdPolicy) )
      {}

    /// Constructor
    /**
     * Use this constructor if for some reasons you have a function space element storing phase ids.
     * Keep in mind that the grid the function space element is using should be consistent with the given grid
     */
    template <class Space>
    DeformingGridManager(std::unique_ptr<Dune::UGGrid<dimension>>&& grid, FunctionSpaceElement<Space,1> const& phase,
        OuterBoundaryPolicy const& outerBoundaryPolicy, InnerBoundaryPolicy const& innerBoundaryPolicy,
        ThresholdPolicy<Scalar> const& thresholdPolicy, bool verbose = false)
    : GridManagerBase<Dune::UGGrid<dimension>>(std::move(grid),verbose),
      deformation(new Interpolation(this->gridptr->leafView(), phase, BoundaryNormalCollector<GridView, OuterBoundaryPolicy, InnerBoundaryPolicy, FunctionSpaceElement<Space,1>, Policy::PhaseAsFSElement>(this->gridPtr->leafView(), phase, outerBoundaryPolicy, innerBoundaryPolicy).meanNormals(), outerBoundaryPolicy, thresholdPolicy) )
      {}

    void initializeDeformation(InterpolationTools::NormalContainer<Scalar,dimension> const& normalContainer,
                                 OuterBoundaryPolicy const& outerBoundaryPolicy, ThresholdPolicy<Scalar> const& thresholdPolicy)
    {
      deformation.reset( new Interpolation( this->gridptr->leafView(), normalContainer, outerBoundaryPolicy, thresholdPolicy ) );
    }

    template <class Space>
    void initializeDeformation(InterpolationTools::NormalContainer<Scalar,dimension> const& normalContainer, FunctionSpaceElement<Space,1> const& phase,
                                 OuterBoundaryPolicy const& outerBoundaryPolicy, ThresholdPolicy<Scalar> const& thresholdPolicy)
    {
      deformation.reset( new Interpolation( this->gridptr->leafView(), normalContainer, phase, outerBoundaryPolicy, thresholdPolicy ) );
    }

    void initializeDeformation(OuterBoundaryPolicy const& outerBoundaryPolicy, InnerBoundaryPolicy const& innerBoundaryPolicy,
                                 ThresholdPolicy<Scalar> const& thresholdPolicy)
    {
      deformation.reset( new Interpolation( this->gridptr->leafView(), BoundaryNormalCollector<GridView,OuterBoundaryPolicy,InnerBoundaryPolicy>(this->gridptr->leafView(), outerBoundaryPolicy, innerBoundaryPolicy).meanNormals(), outerBoundaryPolicy, thresholdPolicy) );
    }

    template <class Space>
    void initializeDeformation(FunctionSpaceElement<Space,1> const& phase, OuterBoundaryPolicy const& outerBoundaryPolicy,
                                 InnerBoundaryPolicy const& innerBoundaryPolicy,  ThresholdPolicy<Scalar> const& thresholdPolicy)
    {
      BoundaryNormalCollector<GridView,OuterBoundaryPolicy,InnerBoundaryPolicy,FunctionSpaceElement<Space,1>,Policy::PhaseAsFSElement> collector(this->gridptr->leafView(), phase, outerBoundaryPolicy, innerBoundaryPolicy);
      deformation.reset( new Interpolation( this->gridptr->leafView(), collector.meanNormals(), phase, outerBoundaryPolicy, thresholdPolicy ) );
    }

    void deleteDeformation() { deformation = nullptr; }

    virtual void update()
    {
      //this->gridptr->update();
    }

    virtual void refineGrid(int refCount)
    {
      if(deformation==nullptr)
      {
        std::cout << "SKIPPING MESH REFINEMENT!\n Initialize deformation before refinement!!" << std::endl;
        return;
      }

      if (this->verbose) std::cout << std::endl << "mesh is refined: " << this->gridptr->size(0) << " -> " << std::flush;

      deformGrid(this->grid_non_const(), *deformation, true, true);
      this->gridptr->globalRefine(refCount);
      deformGrid(this->grid_non_const(), *deformation, true, false);

      if (this->verbose) std::cout << this->gridptr->size(0) << " cells" << std::endl;
    }

    virtual bool adaptGrid()
    {
      if(deformation==nullptr)
      {
        std::cout << "SKIPPING MESH REFINEMENT!\n Initialize deformation before refinement!!" << std::endl;
        return false;
      }
      if (this->verbose) std::cout << std::endl << "mesh is refined: "  << this->gridptr->size(0) << " -> " << std::flush;
      deformGrid(this->grid_non_const(), *deformation, true, true);
      bool res = this->gridptr->adapt();
      deformGrid(this->grid_non_const(), *deformation, true, false);
      if (this->verbose) std::cout << this->gridptr->size(0) << " cells" << std::endl;
      return res;
    }

  private:
    std::unique_ptr<Interpolation> deformation;
  };
}

#endif /* DEFORMINGGRIDMANAGER_HH_ */
