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

#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

double const PI = 3.141592653589793238462;

/// Example simple poisson equation
///

template <class RType, class VarSet>
class PoissonFunctional: public FunctionalBase<VariationalFunctional>
{
public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  //static ProblemType const type = VariationalFunctional;
  //static ProblemType const type = WeakFormulation;

/// \class DomainCache
///

  class DomainCache 
  {
  public:
    DomainCache(PoissonFunctional<RType,AnsatzVars> const&,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      data(vars_) 
    {}

    template <class Entity>
    void moveTo(Entity const &entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      int const uIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;

      u  = at_c<0>(data.data).value(at_c<uIdx>(evaluators));
      du = at_c<0>(data.data).gradient(at_c<uIdx>(evaluators))[0];
      
      xglob = e->geometry().global(x);
      double x_1, x_2;
      x_1 = xglob[0];
      x_2 = xglob[1];
      //right hand side f = -laplace(u) , where u is known solution
      f = PI*PI*(9*std::cos(3*PI*x_1)+16*std::sin(4*PI*x_2));
    }

    Scalar
    d0() const 
    {
      return du*du/2 - f*u;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return -f*arg.value; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1,
                               VariationalArg<Scalar,dim> const &arg2) const 
    {
      return arg1.gradient[0]*arg2.gradient[0];
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Scalar u, f;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> du;
    
    Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> xglob;
  };

/// \class BoundaryCache
///


  class BoundaryCache 
  {
  public:
    static const bool hasInteriorFaces = false;
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(PoissonFunctional<RType,AnsatzVars> const&,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      data(vars_), e(0)
    {}

    void moveTo(FaceIterator const& entity)
    {
      e = &entity;
      penalty = 1.0e14;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                                      AnsatzVars::Grid::dimension-1>
                   const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      
      int const uIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      xglob = (*e)->geometry().global(x);

      u = at_c<0>(data.data).value(at_c<uIdx>(evaluators));
      
      double x_1 = xglob[0];
      double x_2 = xglob[1];
      // term for dirichlet boundary condition
      u0 = std::cos(3*PI*x_1) + std::sin(4*PI*x_2);
    }

    // inhomogeneous Dirichlet condition on one part of the boundary and homogeneous Neumann condition on the rest
    Scalar
    d0() const 
    {
      if ( (xglob[1]<=1e-12) || (xglob[1]>=(1-1e-12)) )
	return penalty*(u-u0)*(u-u0)/2;
      else return 0.0;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      if ( (xglob[1]<=1e-12) || (xglob[1]>=(1-1e-12)) )
        return penalty*(u-u0)*arg.value[0];
      else return 0.0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, 
                       VariationalArg<Scalar,dim> const &arg2) const 
    {
      if ( (xglob[1]<=1e-12) || (xglob[1]>=(1-1e-12)) )
        return penalty*arg1.value*arg2.value;
      else return 0.0;
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    FaceIterator const* e;
    
    Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> xglob;
    
    Scalar penalty, u, u0;
  };

/// \struct D2
///

  template <int row>
  struct D1: public FunctionalBase<VariationalFunctional>::D1<row> 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };
  
public:
  template <int row, int col>
  struct D2: public FunctionalBase<VariationalFunctional>::D2<row,col> 
  {
    static bool const present = true;
    static bool const symmetric = true;
    static bool const lumped = false;
  };

/// \fn integrationOrder
///

  template <class Cell>
  int integrationOrder(Cell const& /* cell */,
                       int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder-1;
  }
};

/// \example poisson.cpp
/// show the usage of PoissonFunctional describing a stationary heat transfer problem,
/// no adaptive grid refinement.
///
#endif
