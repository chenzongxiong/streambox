/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ELASTOMECHANICS_HH
#define ELASTOMECHANICS_HH

#include <type_traits>
#include "fem/functional_aux.hh"
#include "fem/diffops/elastoVariationalFunctionals.hh"


// Deriving from FunctionalBase introduces default D1 and D2 structures.
template <class VarSet>
class ElasticityFunctional: public Kaskade::FunctionalBase<VariationalFunctional>
{
public:
  using Scalar = double;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using OriginVars = VarSet;
  static int const dim = AnsatzVars::Grid::dimension;
  using ElasticEnergy = Kaskade::Elastomechanics::LameNavier<dim,Scalar>;
  using Vector = Dune::FieldVector<Scalar, dim>;
  using Matrix = Dune::FieldMatrix<Scalar, dim, dim>;
  static int constexpr u_Idx = 0;
  static int constexpr u_Space_Idx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables, u_Idx>::type::spaceIndex;


  class DomainCache 
  {
  public:
    DomainCache(ElasticityFunctional const& functional, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_), energy(functional.moduli)
    {}

    template <class Entity>
    void moveTo(Entity const& entity) {}

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators)
    {
      energy.setLinearizationPoint( boost::fusion::at_c<u_Idx>(vars.data).derivative(boost::fusion::at_c<u_Space_Idx>(evaluators)) );
    }

    Scalar d0() const
    {
      return energy.d0();
    }

    template<int row>
    Vector d1 (VariationalArg<Scalar,dim> const& arg) const
    {
      return energy.d1(arg);
    }

    template<int row, int col>
    Matrix d2 (VariationalArg<Scalar,dim> const& argTest, VariationalArg<Scalar,dim> const& argAnsatz) const
    {
      return energy.d2(argTest,argAnsatz);
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    ElasticEnergy energy;
  };

  class BoundaryCache : public CacheBase<ElasticityFunctional,BoundaryCache>
  {
  public:
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(ElasticityFunctional const& f_, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_)
    {}


    void moveTo(FaceIterator const& face)
    {
      Vector up(0); up[1] = 1;                        // unit upwards pointing vector
      auto n = face->centerUnitOuterNormal();         // unit outer normal
      
      if (n*up > 0.5)          // top face: natural boundary conditions (homogeneous Neumann, zero normal stress)
      {
        alpha = 0;             // retardation force factor
        beta  = 0;             // absolute force vector
      }
      else if (n*up < -0.5)    // bottom face: force boundary condition (inhomogeneous Neumann, given normal stress)
      {
        alpha = 0;
        beta  = up;
      }
      else                     // side face: essential boundary conditions (homogeneous Dirichlet)
      {
        alpha = 1e14;          // requires large penalty for hard materials such as steel with a Young's modulus around 7e9
        beta  = 0;
      }
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,dim-1> const& x, Evaluators const& evaluators)
    {
      using namespace boost::fusion;

      u = at_c<u_Idx>(vars.data).value(at_c<u_Space_Idx>(evaluators));
    }

    Scalar
    d0() const
    {
      return alpha*(u*u) - beta*u;
    }

    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,dim> const& arg) const
    {
      return 2*alpha*(u*arg.value) - beta*arg.value;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,dim> const &arg1, VariationalArg<Scalar,dim,dim> const &arg2) const
    {
      return 2*alpha*(arg1.value*arg2.value);
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    Vector u, beta;
    Scalar alpha;
  };

  ElasticityFunctional(Kaskade::ElasticModulus const& moduli_): moduli(moduli_)
  {
  }


  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
  {
    if (boundary)
      return 2*shapeFunctionOrder;      // mass term u*u on boundary
    else
      return 2*(shapeFunctionOrder-1);  // energy term "u_x * u_x" in interior
  }

  Kaskade::ElasticModulus moduli;
};

#endif /* ELASTOMECHANICS_HH_ */
