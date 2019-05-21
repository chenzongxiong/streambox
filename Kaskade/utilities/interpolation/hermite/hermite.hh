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

#ifndef HERMITE_INTERPOLATION_HH
#define HERMITE_INTERPOLATION_HH

#include <memory>
#include <utility>

#include "utilities/interpolation/tools.hh"
#include "utilities/interpolation/hermite/hermite2D.hh"
#include "utilities/interpolation/hermite/hermite3D.hh"

namespace Kaskade
{
  /**
   * @file hermiteinterpolation.hh
   * @brief  Wrapper choosing the correct template specialization for hermite interpolation.
   * @author Lars Lubkoll
   */

  /// Hermite interpolation class.
  /*
   * Hermite interpolation on a simplicial grid. It is assumed that the function vanishes
   * ot the vertices. Moreover the tangent planes at each vertex and in 3D also
   * at each edge center are given through its outer normals.
   *
   * \param InterpolationPolicy policy type for the specification of interpolation details
   * \param Container data container holding vectors associated to vertices and in 3D also edges
   */
  template<class Scalar_, class GridView_, class OuterBoundaryPolicy=Policy::IgnoreOuterBoundary, template <typename> class ThresholdPolicy=Policy::NoGradientThreshold>
  class HermiteInterpolation
  {
  public:
    typedef Scalar_ Scalar;
    static int const dim = GridView_::dimension;
    static bool const needsShapeFunctionSet = true;
    typedef GridView_ GridView;
    typedef InterpolationTools::NormalContainer<Scalar, dim> Container;

    typedef Dune::FieldVector<Scalar,dim> range_type;
    typedef range_type domain_type;
    typedef range_type Range;
    typedef domain_type Domain;
    typedef Scalar field_type;

    /// Default constructor.
    HermiteInterpolation() : implementation(0)
    {}

    /// Copy constructor.
    HermiteInterpolation(HermiteInterpolation const& other) : implementation(other.implementation)
    {}

    /// Move constructor.
    HermiteInterpolation(HermiteInterpolation&& other) : implementation(0)
    {
      implementation.swap(other.implementation);
    }

    /// Copy assignment
    HermiteInterpolation& operator=(HermiteInterpolation const& other)
    {
      if(this!=&other) implementation(other.implementation);
      return *this;
    }

    /// Move assignment
    HermiteInterpolation& operator=(HermiteInterpolation&& other)
    {
      if(this != &other)
      {
        implementation.swap(other.implementation);
        other.implementation.reset();
      }
      return *this;
    }

    /// Constructor
    /*
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals and, in 3D, edge normals in container.vertexNormals/container.edgeNormals
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     */
    template <class Entity, class ShapeFunctionSet>
    HermiteInterpolation(Entity const &entity, GridView const &gridView,
        ShapeFunctionSet const& shapeFunctions, InterpolationTools::NormalContainer<Scalar,GridView::dimension> const& container,
        OuterBoundaryPolicy const& outerBoundaryPolicy=OuterBoundaryPolicy(),
        ThresholdPolicy<Scalar> const& thresholdPolicy=ThresholdPolicy<Scalar>())
        : implementation(new HIPImpl<Scalar,GridView,Policy::MergedPolicy<OuterBoundaryPolicy,ThresholdPolicy<Scalar>,Policy::NoPhaseInfo<Scalar> >,dim>
          (
            entity, gridView, shapeFunctions, container,
            Policy::MergedPolicy<OuterBoundaryPolicy,ThresholdPolicy<Scalar>,Policy::NoPhaseInfo<Scalar> >( outerBoundaryPolicy, thresholdPolicy, Policy::NoPhaseInfo<Scalar>() )
          )
        )
        {}

    /// Constructor
    /*
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals and, in 3D, edge normals in container.vertexNormals/container.edgeNormals
     * \param phaseElement (discontinuous) function space element containing information on different phases via associated ids
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     */
    template <class Entity, class ShapeFunctionSet, class PhaseElement>
    HermiteInterpolation(Entity const &entity, GridView const &gridView,
        ShapeFunctionSet const& shapeFunctions, InterpolationTools::NormalContainer<Scalar,GridView::dimension> const& container, PhaseElement const& phaseElement,
        OuterBoundaryPolicy const& outerBoundaryPolicy=OuterBoundaryPolicy(),
        ThresholdPolicy<Scalar> const& thresholdPolicy=ThresholdPolicy<Scalar>())
        : implementation(new HIPImpl<Scalar,GridView,Policy::MergedPolicy<OuterBoundaryPolicy,ThresholdPolicy<Scalar>,Policy::PhaseAsFSElement<PhaseElement> >,dim>
          (
            entity, gridView, shapeFunctions, container,
            Policy::MergedPolicy<OuterBoundaryPolicy,ThresholdPolicy<Scalar>,Policy::PhaseAsFSElement<PhaseElement> >(outerBoundaryPolicy, thresholdPolicy, Policy::PhaseAsFSElement<PhaseElement>(phaseElement))
          )
        )
        {}

    /// Evaluate interpolation polynomial at position x in local coordinates.
    template <class Vector, class ShapeFunctionSet>
    Range evaluate(Vector const &x, ShapeFunctionSet const& shapeFunctionSet) const
    {
      return implementation->evaluate(x,shapeFunctionSet);
    }

  private:
    mutable std::shared_ptr<HIPBase<Scalar,dim> > implementation;
  };
}
#endif
