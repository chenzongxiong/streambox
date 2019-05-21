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

#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

/// Example simple heat transfer equation
///

template <class RType, class VarSet>
class ATPFunctional : public FunctionalBase<WeakFormulation>
{
public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;

/// \class DomainCache

  class DomainCache : public CacheBase<ATPFunctional,DomainCache> 
  {
  public:
    DomainCache(ATPFunctional<RType,VarSet> const& F_,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      F(F_),
      data(vars_) 
    {}

    template <class Entity>
    void moveTo(Entity const &entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& z, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>
                          xglob = e->geometry().global(z);

      U = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
      dU = at_c<0>(data.data).gradient(at_c<TIdx>(evaluators))[0];
      x = xglob[0];
      q = x*x;
      Uh = 0.9*exp(-q)+0.1*U;
      Uxxh = Uh*(4.0*x*x-2.0);
      dUxxh = 0.4*x*x-0.2;
      
      f = -(Uxxh+F.fsign*(exp(U)-exp(exp(-q))));
      df = -(dUxxh+F.fsign*exp(U));
      
    }

    Scalar
    d0() const 
    {
      return dU*dU/2-f*U;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return dU*arg.gradient[0] - f*arg.value; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1,
                               VariationalArg<Scalar,dim> const &arg2) const 
    {
      return arg1.gradient[0]*arg2.gradient[0]-df*arg1.value*arg2.value;
    }

  private:
    ATPFunctional<RType,VarSet> const& F;
    typename AnsatzVars::VariableSet const& data;
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Scalar U, f, df;
    Scalar q, Uh, Uxxh, dUxxh, x, fsign;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> dU;
  };
/// \class BoundaryCache

  class BoundaryCache : public CacheBase<ATPFunctional,BoundaryCache> 
  {
  public:
    static const bool hasInteriorFaces = false;

    BoundaryCache(ATPFunctional<RType,VarSet> const& F_,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      F(F_),
      data(vars_)
    {}

    template <class Entity>
    void moveTo(Entity const &entity)
    {
      penalty = 1.0e9; 
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                                      AnsatzVars::Grid::dimension-1>
                   const& z, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>

      U = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
      U0 = 0.0;   
    }

    Scalar
    d0() const 
    {
      return penalty*(U-U0)*(U-U0)/2;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return penalty*(U-U0)*arg.value[0];
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
    ATPFunctional<RType,VarSet> const& F;
    typename AnsatzVars::VariableSet const& data;
    Scalar penalty, U, U0;
  };
/// \struct D2

  template <int row>
  struct D1: public FunctionalBase<WeakFormulation>::D1<row> 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };
  
public:
  template <int row, int col>
  struct D2: public FunctionalBase<WeakFormulation>::D2<row,col>
  {
    static bool const present = true;
    static bool const symmetric = false;
    static bool const lumped = false;
  };

/// \fn integrationOrder

  template <class Cell>
  int integrationOrder(Cell const& /* cell */,
                       int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder-1;
  }

ATPFunctional(Scalar fsign_): fsign(fsign_) {};

private:
  Scalar fsign;
};

/// \example atp.cpp
/// show the usage of ATPFunctional describing an artificial nonlinear test problem,
/// no adaptive grid refinement.
///
#endif
