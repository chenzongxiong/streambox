/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MOVINGSOURCE_HH
#define MOVINGSOURCE_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

/// Example simple moving source equation
///

int printCount = 10;

template <class RType, class VarSet>
class MovingSourceEquation: public FunctionalBase<WeakFormulation>
{
  using Self = MovingSourceEquation<RType,VarSet>;
public:
  using Scalar = RType;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using OriginVars = VarSet;


  using Cell = typename AnsatzVars::Grid::template Codim<0>::Entity;
  using Position = Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension>;

/// \class DomainCache
///
  class DomainCache 
  {
  public:
    DomainCache(Self const& F_,
                typename AnsatzVars::VariableSet const& vars_,int flags=7):
      F(F_), vars(vars_) 
    {}

    void moveTo(Cell const &entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,0>::type::spaceIndex;
      auto xglob = e->geometry().global(x);

      u = at_c<0>(vars.data).value(at_c<TIdx>(evaluators));
      du = at_c<0>(vars.data).gradient(at_c<TIdx>(evaluators))[0];
      t = F.time();

      double vx = xglob[0]-0.25*(2.0 + sin(3.14159265358979323846*t));
      double vy = xglob[1]-0.25*(2.0 + cos(3.14159265358979323846*t));
      double h = 0.8*exp(-80.0*(vx*vx + vy*vy));
      double vt = vx * cos(3.14159265358979323846*t) - vy * sin(3.14159265358979323846*t);
      double g  = -25600.0 * ( vx*vx + vy*vy) + 320.0 + 3.14159265358979323846*40.0*vt;
      f = h*g;
    }

    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return -du*arg.gradient[0] + f*arg.value;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return -arg1.gradient[0]*arg2.gradient[0] + 0.0*arg1.value*arg2.value;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    b2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return arg1.value[0]*arg2.value[0];
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;

    Cell const* e;

    Scalar u, t, f;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> du;
  };

/// \class BoundaryCache
///
  class BoundaryCache 
  {
  public:
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(Self const& F_,
                  typename AnsatzVars::VariableSet const& vars_,int flags=7):
      F(F_), vars(vars_), e(0)
    {}

    void moveTo(FaceIterator const& entity)
    {
      e = &entity;
      penalty = 1.0e9;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> xglob = (*e)->geometry().global(x);
      u = at_c<0>(vars.data).value(at_c<TIdx>(evaluators));
      t = F.time();
	  
	  double vx = xglob[0]-0.25*(2.0 + sin(3.14159265358979323846*t));
	  double vy = xglob[1]-0.25*(2.0 + cos(3.14159265358979323846*t));

      u0 = 0.8*exp(-80.0*(vx*vx + vy*vy));
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return penalty*(u-u0)*arg.value[0];
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return penalty*arg1.value*arg2.value;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    b2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return 0;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* e;
    Scalar penalty, u, u0, t;
  };

  class Scaling 
  {
  public:
    Scaling(Self const& self_): self(&self_) {}

    template <class Vals>
    void scale(Cell const&, Position const&, Vals& vals) const 
    {
    }

  private:
    Self const* self;
  };

/// \struct D2
///
  
public:

  void newTime(Scalar dt)
    {
      t += dt;
    }
  Scalar time() const { return t; }
  void time(Scalar tnew)  { t=tnew; }
  Scaling scaling() const { return Scaling(*this); }
  void temporalEvaluationRange(double t0, double t1) { tau = t1-t0; }

  template <int row>
  struct D1: public FunctionalBase<WeakFormulation>::D1<row> 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };

  template <int row, int col>
  struct D2: public FunctionalBase<WeakFormulation>::D2<row,col> 
  {
    static bool const present = true;
    static bool const symmetric = true;
    static bool const lumped = false;
  };

  template <int row, int col>
  struct B2 
  {
    static bool const present   = true;
    static bool const symmetric = true;
    static bool const constant  = false;
    static bool const lumped = false;
  };

  /**
   * Given an initial value, this is transfered to a properly scaled
   * finite element iterate.
   */
  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& u0, typename AnsatzVars::VariableSet& u) const 
  {
    interpolateGloballyWeak<Volume>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(true,u0,*this));
  }

/// \fn integrationOrder
///

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder+1;
  }

private:
  Scalar t, tau;  

  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static int const components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& u0_,
                   MovingSourceEquation<RType,AnsatzVars> const& f_): doScaling(doScaling_), u0(u0_), f(f_) {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
    {
      ValueType u = u0.value(cell,localCoordinate);
      return u;
    }

  private:
    bool doScaling;
    WeakFunctionView const& u0;
    MovingSourceEquation<RType,AnsatzVars> const& f;
  };

};
/// \example movingsource.cpp
/// show the usage of MovingSourceEquation describing an instationary heat transfer problem,
/// time discretization is done by LIMEX, spatial discretization by FEM using an embedded
/// error estimation.
///
#endif
