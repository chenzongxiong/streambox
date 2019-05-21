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

#ifndef FEM_INTEGRATION_HH
#define FEM_INTEGRATION_HH

/**
 * @file
 * @brief  Functionalities for integration of FunctionSpaceElement s or FunctionViews.
 * @author Anton Schiela, Jakob Schneck
 */

#include <dune/geometry/quadraturerules.hh>

#include "fem/functionspace.hh"
#include "fem/forEach.hh"

namespace Kaskade
{
/// \cond internals
namespace Integration_Detail {
  // C++ 14 has polymorphic lambdas..., but without method "order".
  template <class Function, class Functor>
  struct G {
    G(Function const& f_, Functor& g_): f(f_), g(g_) {}
    template <class Cell, class QuadPoint, class Evaluator>
    void operator()(Cell const& cell, QuadPoint const& qp, Evaluator const& eval) const {
      g(cell,qp,f.value(eval));
    }
    template <class Evaluator>
    int order(Evaluator const& eval) const { return f.order(eval); }
    Function const& f;
    Functor& g;
  };

}
/// \endcond

/**
 * \brief Loops over the whole domain occupied by the grid and calls the
 * functor g for each integration point with the value of the function
 * f. This version of forEach requires the function f to be
 * evaluatable by an Evaluator of the associated space of f.
 *
 * g has to provide a (possibly non-const) operator()(CellPointer,
 * QuadPoint const&, Value) with Cell=Grid::Codim<0>::Entity,
 * QuadPoint=QuadraturePoint<Grid::ctype,Grid::dimension>, and
 * Value=Function::ValueType.
 *
 */
template <class Function, class Functor>
void forEachQuadPoint(Function const& f, Functor& g)
{
  Integration_Detail::G<Function,Functor> gWrapper(f,g);
// this lambda would be great, but it doesn't provide method order
//     auto gWrapper = [&] (auto const& cellPointer, auto const& quadPoint, auto w, auto const& eval)
//     {
//       g(cellPointer,quadPoint,f.value(eval));
//     };
  forEachQP(gWrapper,f.space());
}

/**
 * \ingroup fem
 * \brief Loops over the whole domain occupied by the grid and calls the
 * functor g for each integration point with an evaluator of the
 * associated space.
 *
 * g has to provide a (possibly non-const) operator()(CellPointer,
 * QuadPoint const&, Grid::ctype, Evaluator const&) with Cell=Grid::Codim<0>::Entity,
 * QuadPoint=QuadraturePoint<Grid::ctype,Grid::dimension>
 * and a method int order(Evaluator const&).
 *
 * WARNING: This function is deprecated in favor of an alternative interface
 * that does not require g to have a method "order" but takes a second functor
 * providing this order.
 */
template <class Functor, class Space>
void forEachQP(Functor& g, Space const& space)
{
  typedef typename Space::Grid Grid;
  typedef typename Grid::LeafGridView GridView ;
  typedef typename Space::Evaluator Evaluator;

  typedef ShapeFunctionCache<Grid,typename Space::Scalar> SfCache;
  SfCache sfCache;

  GridView const& gridView = space.grid().leafGridView() ;

  Evaluator evaluator(space,&sfCache);
  for(auto const& cell : elements(gridView)) {

    evaluator.moveTo(cell);
    Dune::GeometryType gt = cell.type();
    auto const& qr = Dune::QuadratureRules<typename Grid::ctype,Grid::dimension>::rule(gt,g.order(evaluator));
    evaluator.useQuadratureRule(qr,0);

    size_t nQuadPos = qr.size();
    for (size_t i=0; i<nQuadPos; ++i) {
      evaluator.evaluateAt(qr[i].position(),qr,i,0);
      // @todo: think whether passing the iterator ci to the functor g
      // is the best interface design. This has been introduced by
      // Anton for storing EntityPointers inside the functor during
      // error estimation. But maybe an IndexSet based version is
      // cleaner?
      //g(typename GridView::template Codim<0>::EntityPointer(ci),qr[i],qr[i].weight()*ci->geometry().integrationElement(qr[i].position()),sfs);

      //hopefully better:
      g(cell,qr[i],evaluator);
    }
  }
}

/**
 * \ingroup fem
 * \brief Loops over the whole domain occupied by the grid and calls the
 * functor g for each integration point with the value of the function f.
 *
 * This version of forEach requires the function f to be
 * evaluatable by cell and local coordinate.
 *
 * g has to provide a (possibly non-const) operator()(Cell const&,
 * QuadPoint const&, Value) with Cell=Grid::Codim<0>::Entity,
 * QuadPoint=QuadraturePoint<Grid::ctype,Grid::dimension>, and
 * Value=Function::ValueType.
 */
template <class Functor, class Function, class Space>
void forEachQuadPoint(Function const& f, Functor& g, Space& space)
{
  typedef typename Space::Grid Grid;
  typedef typename Grid::LeafGridView GridView;
  typedef typename Space::Evaluator Evaluator;

  typedef ShapeFunctionCache<Grid,typename Space::Scalar> SfCache;
  SfCache sfCache;

  //   typedef typename IndexSet::template Codim<0>::template Partition<Dune::All_Partition> Entity;
  //   typedef typename Entity::Iterator CellIterator;
  //   IndexSet const& indexSet = space.grid().leafIndexSet();
  typedef typename GridView::template Codim<0>::Iterator CellIterator ;
  GridView const& gridView = f.space().grid().leafGridView() ;

  Evaluator sfs(space, &sfCache);

  for (CellIterator ci=gridView.template begin<0,Dune::All_Partition>();
      ci!=gridView.template end<0,Dune::All_Partition>(); ++ci) {
    sfs.moveTo(*ci);

    Dune::GeometryType gt = ci->type();
    Dune::QuadratureRule<typename Grid::ctype,Grid::dimension> const& qr
    = Dune::QuadratureRules<typename Grid::ctype,Grid::dimension>::rule(gt,f.order(sfs));

    size_t nQuadPos = qr.size();
    for (size_t i=0; i<nQuadPos; ++i) {
      g(typename GridView::template Codim<0>::EntityPointer(ci),qr[i],f.value(*ci,qr[i].position()));
    }
  }
}

/**
 * \ingroup fem
 * @brief integrateOverIntersection computes the integral of an FE function over a grid intersection (face).
 *
 * @param function is the FE function to be integrated.
 * @param intersection codim one domain of integration.
 * @param evaluator is suitable evaluator for the function.
 * @tparam FEFunction is the type of the integrand.
 * @return value of integral
 */
template <class FEFunction>
typename FEFunction::ValueType integrateOverIntersection(FEFunction const& function, typename FEFunction::Space::GridView::Intersection const& intersection, typename FEFunction::Space::Evaluator & evaluator) {
  typename FEFunction::Space::Grid::template Codim<0>::Entity cell = intersection.inside(); //it is important to assign the temporary object to a local variable since the evaluator stores the address of the cell (as pointer) and dereferences it later
  evaluator.moveTo(cell);
  typename FEFunction::ValueType integral(0.0);
  auto const& refElem = evaluator.shapeFunctions().referenceElement();
  auto const& qr = Dune::QuadratureRules<typename FEFunction::Space::ctype,FEFunction::Space::dim-1>::rule(intersection.type(),function.order(evaluator));
  int nQuadPos = qr.size();
  for (int i=0; i<nQuadPos; ++i) {
    auto localPositionInCell = refElem.template geometry<1>(intersection.indexInInside()).global(qr[i].position());
    evaluator.evaluateAt(localPositionInCell,qr,i,intersection.indexInInside());
    integral += qr[i].weight()*intersection.geometry().integrationElement(qr[i].position())*function.value(evaluator);
  }
  return integral;
}
/**
 * \ingroup fem
 * @brief integrateOverBoundary computes the integral of an FE function over the whole boundary of the underlying grid.
 *
 * @param function is the FE function to be integrated.
 * @tparam FEFunction is the type of the integrand.
 * @return value of integral
 *
 */
template <class FEFunction>
typename FEFunction::ValueType integrateOverBoundary(FEFunction const& function) {
  ShapeFunctionCache<typename FEFunction::Space::Grid, typename FEFunction::Space::Scalar> sfCache;
  typename FEFunction::Space::Evaluator evaluator(function.space(),&sfCache);
  typename FEFunction::Space::GridView const& gridView = function.space().gridView();
  typename FEFunction::ValueType integral(0.0);
  forEachBoundaryFace(gridView,[&](typename FEFunction::Space::GridView::Intersection const& intersection) {
    integral += integrateOverIntersection(function,intersection,evaluator);
  });
  return integral;
}
/**
 * \ingroup fem
 * @brief integrateOverBoundary computes the integral of an FE function which is restricted to (parts of) the boundary over the boundary of the underlying grid.
 *
 * This specialized overload does integration only on those faces where the boundary FE function is nonzero.
 *
 * @param function is the FE function to be integrated.
 * @return value of integral
 *
 */
template <template <class, class> class DomainMapper, typename Scalar, typename GridView, int m>
typename FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper, Scalar, GridView>>,m>::ValueType
                integrateOverBoundary(FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper, Scalar, GridView>>,m> const& function) {
  using FEFunction = FunctionSpaceElement<FEFunctionSpace<BoundaryMapper<DomainMapper, Scalar, GridView>>,m>;
  ShapeFunctionCache<typename GridView::Grid, typename FEFunction::Space::Scalar> sfCache;
  typename FEFunction::Space::Evaluator evaluator(function.space(),&sfCache);
  GridView const& gridView = function.space().gridView();
  BoundaryMapper<DomainMapper, Scalar, GridView> const& mapper = function.space().mapper();
  typename FEFunction::ValueType integral(0.0);
  forEachBoundaryFace(gridView,[&](typename GridView::Intersection const& intersection) {
    if(mapper.inDomain(intersection)) integral += integrateOverIntersection(function,intersection,evaluator);
  });
  return integral;
}

//---------------------------------------------------------------------
/// Evaluation class for integrals
template <typename Space>
class Integral
{
class Integrator
{
  typedef typename Space::Grid Grid;
  typedef typename Space::GridView GridView;

  typedef typename Grid::template Codim<0>::Entity Cell;
  typedef Dune::QuadraturePoint<typename Grid::ctype,Grid::dimension> QuadPoint;
  typedef typename Space::template Element<1>::type::ValueType ValueType;

public:
  Integrator() : sum(0.0) {}

//  void operator()(std::conditional_t<std::is_same<Cell,typename Cell::EntityPointer>::value,Cell const*,typename Cell::EntityPointer> const& c, QuadPoint const& q, ValueType v)
//  {
//    v *= q.weight()*c->geometry().integrationElement(q.position());
//    sum += v;
//  }

//  void operator()(typename GridView::template Codim<0>::Iterator const& ci, QuadPoint const& q, ValueType v)
//  {
//    v *= q.weight()*ci->geometry().integrationElement(q.position());
//    sum += v;
//  }
  template <class CellPointer>
  void operator()(CellPointer const& c, QuadPoint const& q, ValueType v)
  {
    v *= q.weight()*c->geometry().integrationElement(q.position());
    sum += v;
  }

  void operator()(Cell const& cell, QuadPoint const& q, ValueType v)
  {
    v *= q.weight()*cell.geometry().integrationElement(q.position());
    sum += v;
  }

  ValueType getIntegral() { return sum; }
private:
  ValueType sum;
};

public:
  template<typename Function>
  typename Function::ValueType operator()(Function const& f)
  {
    Integrator integrator;
    forEachQuadPoint(f,integrator);
    return integrator.getIntegral();
  }

  template<typename Function>
  typename Function::ValueType operator()(Function const& f, typename Function::Space const& space)
  {
    Integrator integrator;
    forEachQuadPoint(f,integrator,space);
    return integrator.getIntegral();
  }
};


/// Integral
struct IntegralF
{
/// Evaluation of integral
  template<typename Function>
  double id(Function const& f) const
  {
    Integral<typename Function::Space> integral;
    return integral(f);
  }

/// Evaluation of norm (??? that's not true)
  template<typename Function>
  double operator()(Function const& f) const
  {
    return id(f);
  }
};


} // end of namespace Kaskade
#endif //FEM_INTEGRATION_HH
