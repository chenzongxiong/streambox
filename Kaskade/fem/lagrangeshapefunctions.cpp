/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dune/grid/config.h"

#include "fem/lagrangeshapefunctions.hh"


namespace Kaskade
{
  
  template <class ctype, int dimension, class Scalar, int O>
  LagrangeCubeShapeFunction<ctype,dimension,Scalar,O>::LagrangeCubeShapeFunction(Dune::FieldVector<int,dimension> const& xi_)
  : xi(xi_)
  {
    for (int i=0; i<dim; ++i)
      assert(xi[i]>=0 && xi[i]<=order);

    // construct interpolation nodes (here: equidistant nodes). ?: in
    // denominator is just to prevent division by 0 for order 0.
    for (int i=0; i<=order; ++i)
      if (order==0)
        node[i] = 0.5;
      else
        node[i] = static_cast<CoordType>(i)/(order? order: 1);

    // compute normalization factor such that the shape functions have
    // value 1 at their node
    Dune::FieldVector<CoordType,dim> x;
    for (int i=0; i<dim; ++i)
      x[i] = node[xi[i]];
    normalization = 1.0 / evaluateNonNormalized(x);

    // compute index local in element
    local_ = 0;
    for (int i=0; i<dim; ++i)
      local_ = (order+1)*local_ + xi[i];

    // compute codimension of subentity in which the node resides. If
    // an index is 0 or order, the node lies on a face in that
    // coordinate direction, hence the codimension has to be increased
    // by 1.
    codim_ = 0;
    for (int i=0; i<dim; ++i)
      if (xi[i]==0 || xi[i]==order)
        ++codim_;

    // compute subentity index. See the Dune reference element page for subentity numbering.
    entity_ = 0;
    if (codim_ == 0) // the cell itself
      entity_ = 0;
    else if (codim_ == dim) {
      // vertices are numbered according to the binary number of their
      // reversed coordinates. All of xi's entries are either 0 or
      // order.
      entity_ = 0;
      for (int i=dim-1; i>=0; --i)
        entity_ = 2*entity_ + (xi[i]==0? 0: 1);
    } else if (codim_ == 1) {
      // faces are numbered according to which coordinate direction they delimit.
      entity_ = 0;
      for (int i=0; i<dim; ++i)
        if (xi[i]==0 || xi[i]==order) // only one i passes this test!
          entity_ = 2*i + (xi[i]==0? 0: 1);
    } else if (codim_ == dim-1) {
      entity_ = 0;
      int off = 0;
      for (int i=0; i<dim; ++i) {
        if (xi[i]>0 && xi[i]<order) // exactly one i passes this test!
          off = (1<<(dim-1))*i;
        else
          entity_ = 2*entity_ + (xi[i]==0? 0: 1);
      }
      entity_ += off;
    } else
      assert("Unknown codimension!"==0);

    // Precompute completely local subentity indices. Computation of
    // those subentity indices which depend on global orientation of
    // the subentity are postponed to actual call.
    if (codim()==0) {
      // We are in the interior of the cell. Interpete the interior
      // nodes as indices to the reduced order order-1
      // by subtracting 1 from all coordinates.
      entityIndex_ = 0;
      for (int i=0; i<dim; ++i)
        if (xi[i]>0 && xi[i]<order)
          entityIndex_ = (order-1)*entityIndex_ + xi[i]-1;
    } else if (order<=2) {
      // Ok, since there is only one degree of freedom on this
      // subentity, a globally unique numbering is trivial...
      entityIndex_ = 0;
    } else
      // Oops, this is the hard case. Here the numbering depends on
      // the actual element. Schedule the computation to be performed
      // lateron.
      entityIndex_ = -1;
  }
  // explicit instantiation
  template LagrangeCubeShapeFunction<double,2,double,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,double,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,double,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  
  template LagrangeCubeShapeFunction<double,3,double,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,double,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,double,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  
  
  template <class ctype, int dimension, class Scalar, int O>
  Scalar LagrangeCubeShapeFunction<ctype,dimension,Scalar,O>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<ctype,dimension> const& x) const 
  {
    ResultType independentFactor = 1;
    for (int i=0; i<dim; ++i)
      if (i!=dir)
        for (int k=0; k<=order; ++k)
          if (k!=xi[i])
            independentFactor *= x[i]-node[k];

    ResultType result = 0;
    for (int j=0; j<=order; ++j) {
      if (j!=xi[dir]) {
        ResultType dphi = 1;
        for (int  k=0; k<=order; ++k)
          if (k!=xi[dir] && k!=j)
            dphi *= x[dir]-node[k];
        result += dphi;
      }
    }

    return independentFactor * result;
  }
  // explicit instantiation
  template double LagrangeCubeShapeFunction<double,2,double,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template double LagrangeCubeShapeFunction<double,2,double,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template double LagrangeCubeShapeFunction<double,2,double,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  
  template double LagrangeCubeShapeFunction<double,3,double,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template double LagrangeCubeShapeFunction<double,3,double,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template double LagrangeCubeShapeFunction<double,3,double,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  
  
  template <class ctype, int dimension, class Scalar, int Ord>
  LagrangeCubeShapeFunctionSet<ctype,dimension,Scalar,Ord>::LagrangeCubeShapeFunctionSet()
  : ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::cube,dimension))
  {
    Dune::FieldVector<int,dimension> xi(0);
    for (int i=0; i<power(Ord+1,dimension); ++i) 
    {
      sf.push_back(value_type(xi));
      this->iNodes.push_back(sf.back().position());
      
      // shift to next node index
      for (int j=0; j<dimension; ++j)
        if (xi[j]<Ord) {
          ++xi[j];
          break;
        } else
          xi[j] = 0;
    }
    
    this->size_ = sf.size();
    this->order_ = dimension * Ord;
  }
  // explicit instantiation
  template LagrangeCubeShapeFunctionSet<double,1,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,double,2>::LagrangeCubeShapeFunctionSet();
  
  template LagrangeCubeShapeFunctionSet<double,2,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,double,2>::LagrangeCubeShapeFunctionSet();
  
  template LagrangeCubeShapeFunctionSet<double,3,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,double,2>::LagrangeCubeShapeFunctionSet();
  
}
