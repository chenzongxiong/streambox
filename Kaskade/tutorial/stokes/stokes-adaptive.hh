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

#ifndef STOKESADAPTIVE_HH
#define STOKESADAPTIVE_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"

/**
 * Implements stokes' equation \f$ -\Delta u + \nabla p = 0,\ \mathrm{div}(u)=0.\f$
 */
template <class Scalar_, class VarSet>
class StokesFunctional : public FunctionalBase<VariationalFunctional>
{
public:
  using Scalar = Scalar_;
  using AnsatzVars = VarSet;
  using OriginVars = VarSet;
  using TestVars = VarSet;

  static constexpr int dim = AnsatzVars::Grid::dimension;

  static constexpr int uIdx = 0;
  static constexpr int pIdx = 1;
  static constexpr int uSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,uIdx>::type::spaceIndex;
  static constexpr int pSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,pIdx>::type::spaceIndex;

  // convenient template aliases
  template <int row> using TestComponents = std::integral_constant<size_t,TestVars::template Components<row>::m>;
  template <int row> using AnsatzComponents = std::integral_constant<size_t,AnsatzVars::template Components<row>::m>;

  class DomainCache : public CacheBase<StokesFunctional,DomainCache>
  {
  public:
    DomainCache(StokesFunctional const&, typename AnsatzVars::VariableSet const& vars_,int flags=7):
      vars(vars_) 
    {}

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      du = at_c<uIdx>(vars.data).gradient(at_c<uSpaceIdx>(evaluators));
      p = at_c<pIdx>(vars.data).value(at_c<pSpaceIdx>(evaluators));
    }

    Scalar d0() const 
    {
      return contraction(du,du)/2 - p*trace(du);
    }
    
    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if(row == uIdx)
      {
        Dune::FieldMatrix<Scalar,dim,dim> dv; dv << arg.gradient;

        return contraction(du,dv) - p*trace(dv);
      }
      if(row == pIdx)
      {
        Scalar dv = arg.value[0];

        return -dv*trace(du);
      }

      return 0.0;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      if(row == uIdx && col == uIdx)
      {
        Dune::FieldMatrix<Scalar,dim,dim> dv, dw; dv << arg1.gradient; dw << arg2.gradient;

        return contraction(dv,dw);
      }
      if(row == uIdx && col == pIdx)
      {
        Dune::FieldMatrix<Scalar,dim,dim> dv; dv << arg1.gradient;
        Scalar dw = arg2.value[0];

        return -dw*trace(dv);
      }

      if(row == pIdx && col == uIdx)
      {
        Scalar dv = arg1.value[0];
        Dune::FieldMatrix<Scalar,dim,dim> dw; dw << arg2.gradient;

        return -dv*trace(dw);
      }

      return 0.0;
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    Dune::FieldVector<Scalar,AnsatzComponents<pIdx>::value> p;
    Dune::FieldMatrix<Scalar,AnsatzComponents<uIdx>::value,dim> du;
  };

  class BoundaryCache : public CacheBase<StokesFunctional,BoundaryCache>
  {
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;
  public:
    BoundaryCache(StokesFunctional const&, typename AnsatzVars::VariableSet const& vars_,int flags=7)
     : vars(vars_), e(nullptr), gamma(1e9)
    {}

    void moveTo(FaceIterator const& entity)
    {
      e = &entity;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      Dune::FieldVector<Scalar,dim> xglob = (*e)->geometry().global(x);
      
      u = at_c<uIdx>(vars.data).value(at_c<uSpaceIdx>(evaluators));
      
      u0 = 0;
      if (xglob[1]>=.99999) u0[0] = 1;
    }

    Scalar d0() const 
    {
      return 0.5*gamma*(u-u0)*(u-u0);
    }
    
    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if(row == uIdx)
      {
        Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> dv; dv << arg.value;

        return gamma*(u-u0)*dv;
      }

      return 0.0;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      if(row == uIdx && col == uIdx)
      {
        Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> dv, dw; dv << arg1.value; dw << arg2.value;

        return gamma*dv*dw;
      }

      return 0.0;
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* e;
    Scalar gamma;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u, u0;
  };

public:
  template <int row, int col>
  struct D2
  {
    static constexpr bool present = !(row==pIdx && col==pIdx);
    static constexpr bool symmetric = col==row;
    static constexpr bool lumped = false;
    static constexpr int derivatives = 1;
  };

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder-2;
  }
};
#endif
