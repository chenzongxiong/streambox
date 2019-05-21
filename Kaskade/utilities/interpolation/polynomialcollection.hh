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

#ifndef INTERPOLATION_POLYNOMIAL_COLLECTION_HH
#define INTERPOLATION_POLYNOMIAL_COLLECTION_HH

#include <vector>

#include "utilities/interpolation/tools.hh"

namespace Kaskade
{
  /**
   * @file interpolationpolynomialcollection.hh
   * @brief  Container class holding a polynomial for each cell. Can be evaluated like WeakFunctionViews.
   * @author Lars Lubkoll
   */
  template<class InterpolationPolynomial, template <typename> class PolynomialPolicy = Policy::ShapeFunctionSetPolicy>
  class InterpolationPolynomialCollection
  : private PolynomialPolicy<InterpolationPolynomial>
  {
    typedef typename InterpolationPolynomial::GridView GridView;
    typedef typename GridView::template Codim<0>::Entity Entity;
    typedef typename GridView::template Codim<0>::Iterator CellIterator;
    typedef PolynomialPolicy<InterpolationPolynomial> Policy;
  public:
    typedef typename InterpolationPolynomial::range_type ValueType;

    // no copying
    InterpolationPolynomialCollection(InterpolationPolynomialCollection const&) = delete;
    InterpolationPolynomialCollection& operator=(InterpolationPolynomialCollection const&) = delete;

    /// Constructor
    /**
     * Constructs an interpolation polynomial for each entity. The constructors signature must be:
     * a) InterpolationPolynomial(GridView const&, Entity const&, Parameters const&...) (default)
     * b) InterpolationPolynomial(GridView const&, Entity const&, ShapeFunctionSet const&, Parameters const&...)
     *    (if PolynomialPolicy=UseShapeFunctionSet)
     *
     * For polynomials that can not be constructed this way you must implement your own policy.
     * It is assumed that policies used here take the type 'InterpolationPolynomial' as template parameter.
     * For examples see Policy::NoShapeFunctionSet, Policy::UseShapeFunctionSet.
     *
     * \param gridView_
     * \param parameters
     */
    template <typename... Parameters>
    InterpolationPolynomialCollection(GridView const& gridView_, Parameters const&... params) :
    gridView(gridView_),
    interpolationPolynomials(gridView.size(0))
    {
      CellIterator end = gridView.template end<0>();
      for(CellIterator ci = gridView.template begin<0>(); ci!=end; ++ci) interpolationPolynomials[gridView.indexSet().index(*ci)] = Policy::init(*ci, gridView, params...);
    }

    /// Evaluation.
    /**
     * This method is called by interpolateGloballyWeak in fetransfer.hh
     * \param Cell codim 0 entity type
     * \param cell cell containing the evaluation point
     * \param localCoordinate coordinates of the evaluation point in the cell
     * \return interpolated value at localCoordinate
     */
    template <class Cell>
    ValueType value(Cell const& cell, Dune::FieldVector<typename Cell::ctype, Cell::dimension> const& localCoordinate) const
    {
      return Policy::evaluate( interpolationPolynomials[ gridView.indexSet().index(cell) ], localCoordinate );
    }

  private:
    GridView const& gridView;
    std::vector<InterpolationPolynomial> interpolationPolynomials;
    using PolynomialPolicy<InterpolationPolynomial>::init;
  };
} /* end of namespace Kaskade */
#endif
