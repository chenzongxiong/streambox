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
//   * \todo:       This file currently depends on having a UUGrid !!! 

#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include <cmath>

//#include "fem/fixdune.hh"
#include "fem/variables.hh"
#include "fem/celldata.hh"

/// Example simple heat transfer equation
///

template <class RType, class VarSet>
class SSTFunctional: public FunctionalBase<WeakFormulation>
{
  using Self = SSTFunctional<RType,VarSet>;
  
public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;

  static constexpr int dim = AnsatzVars::Grid::dimension;
  using Cell = typename AnsatzVars::Grid::template Codim<0>::Entity;
  using Position = Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension>;

private:
  static int const u0Idx = 0;
  static int const u1Idx = 1;
  static int const u2Idx = 2;
  static int const u3Idx = 3;
  static int const u0SIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,u0Idx>::type::spaceIndex;
  static int const u1SIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,u1Idx>::type::spaceIndex;
  static int const u2SIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,u2Idx>::type::spaceIndex;
  static int const u3SIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,u3Idx>::type::spaceIndex;
/// \class DomainCache

public:
  class DomainCache : public CacheBase<SSTFunctional,DomainCache>
  {
  public:
    DomainCache(SSTFunctional<RType,VarSet> const& F_,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      F(F_),
      vars(vars_) 
    {}

    void moveTo(Cell const& entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& z, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>
                          xglob = e->geometry().global(z);

      u[0]  = at_c<u0Idx>(vars.data).value(at_c<u0SIdx>(evaluators));
      du[0] = at_c<u0Idx>(vars.data).gradient(at_c<u0SIdx>(evaluators))[0];
      u[1]  = at_c<u1Idx>(vars.data).value(at_c<u1SIdx>(evaluators));
      du[1] = at_c<u1Idx>(vars.data).gradient(at_c<u1SIdx>(evaluators))[0];
      u[2]  = at_c<u2Idx>(vars.data).value(at_c<u2SIdx>(evaluators));
      du[2] = at_c<u2Idx>(vars.data).gradient(at_c<u2SIdx>(evaluators))[0];
      u[3]  = at_c<u3Idx>(vars.data).value(at_c<u3SIdx>(evaluators));
      du[3] = at_c<u3Idx>(vars.data).gradient(at_c<u3SIdx>(evaluators))[0];
      x = xglob[0];  y=xglob[1];
      ka[0] = 0.5e-9; ka[1] = 0.5e-9; ka[2] = 0.5e-9; ka[3] = 0.5e-9; 
    sst1=360.0; small=1.0e-10;
    xleft =0.5-small; xright=0.6+small;
    if  ( (x >= xleft) && (x <= xright) &&  
      (y >= xleft) && (y <= xright) )
     {sst1=3250.0;};
    f[0] = 4.0e5 - 272.443800016*u[0] + 1.0e-4*u[1] + 0.007*u[3] - 3.67e-16*u[0]*u[1] - 4.13e-12*u[0]*u[3];
    f[1] = 272.4438*u[0] - 1.00016e-4*u[1] + 3.67e-16*u[0]*u[1] - 3.57e-15*u[1]*u[2];
    f[2] = -1.6e-8*u[2] + 0.007e0*u[3] + 4.1283e-12*u[0]*u[3] - 3.57e-15*u[1]*u[2]+800.0+sst1;
    f[3] = -7.000016e-3*u[3] + 3.57e-15*u[1]*u[2] - 4.1283e-12*u[0]*u[3] + 800.0;
    df[0][0] = -272.443800016-3.67e-16*u[1]-4.13e-12*u[3];
    df[1][0] = 272.4438  + 3.67e-16*u[1];
    df[2][0] = 4.1283e-12*u[3];
    df[3][0] = -4.1283e-12*u[3];
    df[0][1] = 1.0e-4 - 3.67e-16*u[0];
    df[1][1] = -1.00016e-4 + 3.67e-16*u[0] - 3.57e-15*u[2];
    df[2][1] = -3.57e-15*u[2];
    df[3][1] = 3.57e-15*u[2];
    df[0][2] = 0.0;
    df[1][2] = -3.57e-15*u[1];
    df[2][2] = -1.6e-8 - 3.57e-15*u[1];
    df[3][2] = 3.57e-15*u[1];
    df[0][3] = 0.007-4.13e-12*u[0];
    df[1][3] = 0.0;
    df[2][3] = 0.007+4.1283e-12*u[0];
    df[3][3] = -7.000016e-3-4.1283e-12*u[0];
      
    }
    
    Scalar
    d0() const
    {
      return (du[0]*du[0]+du[1]*du[1]+du[2]*du[2]+du[3]*du[3])/2 - f*(u[0]+u[1]+u[2]+u[3]);
    }

    template<int row> 
    Scalar d1_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const 
    {
      assert( row >=0 && row <=3 ); 
      return ka[row]*du[row]*arg.gradient[0] - f[row]*arg.value; 
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1,
                  VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg2) const   {
      assert( row >=0 && row <=3 && col >=0 && col <=3 ); 
      if (row==col)
         return ka[row]*arg1.gradient[0]*arg2.gradient[0]-df[row][col]*arg1.value*arg2.value;
      else
         return -df[row][col]*arg1.value*arg2.value;
    }

  private:
    SSTFunctional<RType,VarSet> const& F;
    typename AnsatzVars::VariableSet const& vars;
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Scalar u[4], f[4], x, y;
    Scalar sst1, small, xleft, xright;
    Scalar df[4][4], ka[4];
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> du[4];
  };

/// \class BoundaryCache
public:
  class BoundaryCache : public CacheBase<SSTFunctional,BoundaryCache>
  {
  public:
    static const bool hasInteriorFaces = false;
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(SSTFunctional<RType,VarSet> const&,
                  typename AnsatzVars::VariableSet const&,int flags=7)
    {}

    void moveTo(FaceIterator const& entity) {}

    template <class Position, class Evaluators>
    void evaluateAt(Position const&, Evaluators const& evaluators) { }

    Scalar
    d0() const 
    {
      return 0;
    }

    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim> const& arg) const
    {
      return 0; 
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim> const &arg1,
                    VariationalArg<Scalar,dim> const &arg2) const
    {
      return 0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    b2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return 0;
    }
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

  /**
   * Given an initial value, this is transfered to a properly scaled
   * finite element iterate.
   */
  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& us, typename AnsatzVars::VariableSet& u) const 
  {
    switch (row) { 
    case u0Idx: 
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(row==u0Idx,us,*this));
    break;
    case u1Idx: 
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(row==u1Idx,us,*this));
    break;
    case u2Idx: 
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(row==u2Idx,us,*this));
    break;
    case u3Idx: 
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(row==u3Idx,us,*this));
    break;
    default: assert(false);
    }
  }

  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static int const components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& us_,
                   SSTFunctional<RType,AnsatzVars> const& f_): doScaling(doScaling_), us(us_), f(f_) {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
    {
      ValueType u = us.value(cell,localCoordinate);
      return u;
    }

  private:
    bool doScaling;
    WeakFunctionView const& us;
    SSTFunctional<RType,AnsatzVars> const& f;
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
      return 2*shapeFunctionOrder;
  }

};

/// \example ht.cpp
/// show the usage of SSTFunctional describing a stationary heat transfer problem,
/// no adaptive grid refinement.
///
#endif
