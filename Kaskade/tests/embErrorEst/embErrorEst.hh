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

#ifndef PEAKSOURCE_HH
#define PEAKSOURCE_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

/// Example simple stationary heat transfer equation
/// using embedded error estimation for adaptive mesh refinement

template <class RType, class VarSet>
class EmbErrorEstFunctional: public FunctionalBase<VariationalFunctional>
{
public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;

/// \class DomainCache
///

  class DomainCache 
  {
  public:
    DomainCache(EmbErrorEstFunctional<RType,AnsatzVars> const&,
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
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>
                          xglob = e->geometry().global(x);

      T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
      dT = at_c<0>(data.data).gradient(at_c<TIdx>(evaluators))[0];
      
      double g,h,gX,hX,gXX,hXX,gY,hY,gYY,hYY, p;
      p=100.0;
      g   = xglob[0]*(xglob[0]-1)*xglob[1]*(xglob[1]-1);
      gX  = (2*xglob[0]-1)*xglob[1]*(xglob[1]-1);
      gXX = 2*xglob[1]*(xglob[1]-1);
      gY  = (2*xglob[1]-1)*xglob[0]*(xglob[0]-1);
      gYY = 2*xglob[0]*(xglob[0]-1);
      h   = exp( -p*( (xglob[0]-0.5)*(xglob[0]-0.5) + (xglob[1]-0.5)*(xglob[1]-0.5) ) );
      hX  = -2*p*(xglob[0]-0.5)*h;
      hXX = -2*p*h - 2*p*(xglob[0]-0.5)*hX;
      hY  = -2*p*(xglob[1]-0.5)*h;
      hYY = -2*p*h - 2*p*(xglob[1]-0.5)*hY;  
      f = gXX*h + 2*gX*hX + g*hXX  + gYY*h + 2*gY*hY + g*hYY;
      f=-f;
    }

    Scalar
    d0() const 
    {
      return dT*dT/2-f*T;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return dT*arg.gradient[0] - f*arg.value; 
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
    Scalar T, f;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> dT;
  };

/// \class BoundaryCache
///

  class BoundaryCache 
  {
  public:
    static const bool hasInteriorFaces = false;
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(EmbErrorEstFunctional<RType,AnsatzVars> const&,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      data(vars_), e(0)
    {}

    void moveTo(FaceIterator const& entity)
    {
      e = &entity;
      penalty = 1.0e15;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                                      AnsatzVars::Grid::dimension-1>
                   const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;

      T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
      T0 = 0;   
    }

    Scalar
    d0() const 
    {
      return penalty*(T-T0)*(T-T0)/2;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return penalty*(T-T0)*arg.value[0];
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, 
                       VariationalArg<Scalar,dim> const &arg2) const 
    {
      return penalty*arg1.value*arg2.value;
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    FaceIterator const* e;
    Scalar penalty, T, T0;
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

/// \example embErrorEst.cpp
/// show the usage of EmbErrorEstFunctional describing a stationary heat transfer problem with peak formed source, 
/// and adapted grids by embedded error estimator
///
#endif
