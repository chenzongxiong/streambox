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

#ifndef NAVIERSTOKES_HH
#define NAVIERSTOKES_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"
#include "fem/variables.hh"

/// Example simple moving source equation
///

int printCount = 10; 

template <class RType, class VarSet>
class NavierStokesFunctional: public FunctionalBase<WeakFormulation>
{
  using Self = NavierStokesFunctional<RType,VarSet>;
public:
  using Scalar = RType;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using OriginVars = VarSet;

  static constexpr int dim = AnsatzVars::Grid::dimension;

  // phase 1
  static constexpr int uIdx = 0;
  static constexpr int uSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,uIdx>::type::spaceIndex;
  
  // pressure p
  static constexpr int pIdx = 1;
  static constexpr int pSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,pIdx>::type::spaceIndex;

  using Cell = typename AnsatzVars::Grid::template Codim<0>::Entity;
  using Position = Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension>;

  // convenient template aliases
  template <int row> using TestComponents = std::integral_constant<size_t,TestVars::template Components<row>::m>;
  template <int row> using AnsatzComponents = std::integral_constant<size_t,AnsatzVars::template Components<row>::m>;

/// \class DomainCache
///
  class DomainCache : public CacheBase<NavierStokesFunctional,DomainCache> 
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
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> xglob = e->geometry().global(x);

      u  = at_c<uIdx>(vars.data).value(at_c<uSpaceIdx>(evaluators));
      du = at_c<uIdx>(vars.data).derivative(at_c<uSpaceIdx>(evaluators));

      p  = at_c<pIdx>(vars.data).value(at_c<pSpaceIdx>(evaluators));
      dp = at_c<pIdx>(vars.data).derivative(at_c<pSpaceIdx>(evaluators));
      
      t = F.time();
      
      beta = 0.0;   // 
      betaC = 0.0;   // 
      h    = 1.0e-2;   // local meshwidth
      for (int i=0;i<dim;i++)
        for (int j=0;j<dim;j++)
        {
          // anisotrop streamline diffusion
//           double uCabs = 0.0, uAabs = 0.0;
//           for (int i=0;i<dim;i++) {uCabs += uC[i]*uC[i]; uAabs += uA[i]*uA[i];}
//           uCabs =sqrt(uCabs); uAabs =sqrt(uAabs);
//           if (uCabs < 1e-8) uCabs=1e-8; 
//           if (uAabs < 1e-8) uAabs=1e-8;
//           uCuC[i][j] = uC[i]*uC[j]/uCabs;
//           uAuA[i][j] = uA[i]*uA[j]/uAabs;
//           if(i==j) { uCuC[i][j] += 1.0; uAuA[i][j] += 1.0; }
          
          // isotrop
          uu[i][j] = 0.0;
          if(i==j) { uu[i][j] = 1.0;  }
        }
        
    }



    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if(row == 0)
      {
        Dune::FieldVector<Scalar,dim>     v;  v  << arg.value;
        Dune::FieldMatrix<Scalar,dim,dim> dv; dv << arg.derivative;

        return -F.mu*contraction(du,dv)  - p*trace(dv) - F.rho*v*(du*u);  
      }
      
      if(row == 1)
      {
        Scalar v = arg.value[0];

        return -v *trace(du);   
      }

      return 0.0;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      if(row == 0 && col == 0)
      {
        Dune::FieldVector<Scalar,dim>      v,  w;  v << arg1.value;       w << arg2.value;
        Dune::FieldMatrix<Scalar,dim,dim>  dv, dw; dv << arg1.derivative; dw << arg2.derivative;

        return -F.mu*contraction(dv,dw);		// uX*u explicit
      }


     
      if(row == 0 && col == 1)
      {
        Dune::FieldMatrix<Scalar,dim,dim> dv; dv << arg1.derivative;
        Scalar  w = arg2.value[0];
        return -w*trace(dv);
      }


      if(row == 1 && col == 0)
      {
        Scalar v = arg1.value[0]; 
        Dune::FieldVector<Scalar,dim> w; w << arg2.value;
        Dune::FieldMatrix<Scalar,dim,dim> dw; dw << arg2.derivative;

        return -v*trace(dw);
      }

      return 0.0;
    }

    template<int row, int col> 
    Scalar b2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const 
    {
      if (row==0 && col==0) 
      {
        Dune::FieldVector<Scalar,dim> v, w; v << arg1.value; w << arg2.value;
      	return F.rho*(v*w);
      }
            
      return 0.0; 

    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;

    Cell const* e;

    Scalar t, f;  
    Scalar beta, betaC, betaA, h;  
    Dune::FieldVector<Scalar,AnsatzComponents<pIdx>::value> p;
    Dune::FieldMatrix<Scalar,AnsatzComponents<pIdx>::value,dim> dp;
    
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u;
    Dune::FieldMatrix<Scalar,AnsatzComponents<uIdx>::value,dim> du;
    
    Dune::FieldMatrix<Scalar, dim,dim> uu;
  };

/// \class BoundaryCache
///
  class BoundaryCache : public CacheBase<NavierStokesFunctional,BoundaryCache> 
  {
  public:
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(Self const& F_,
                  typename AnsatzVars::VariableSet const& vars_,int flags=7):
      F(F_), vars(vars_), e(nullptr), gamma(1e12)
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
      
      t = F.time();      
      u = at_c<uIdx>(vars.data).value(at_c<uSpaceIdx>(evaluators));  
      u0 = 0; 
      if (xglob[1]>=.99999) 
      {
//         if (t < 1.0) u0[0] = -t;    // + 1e-6;  
//         else         u0[0] = -1.0;
        
        u0[0] = -1.0;
      }
          
    }
    
    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if(row == 0)
      {
        Dune::FieldVector<Scalar,dim> v; v << arg.value;

        return gamma*(u0-u)*v;
      }

      return 0.0;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      if(row == 0 && col == 0)
      {
        Dune::FieldVector<Scalar,dim> v, w; v << arg1.value; w << arg2.value;
        return -gamma*(v*w);
      }

      return 0.0;
    }

    template<int row, int col>
    Scalar b2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      return 0.0;
    }


  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* e;
    Scalar gamma, gamma1;
    
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u0;

    Dune::FieldVector<Scalar,AnsatzComponents<pIdx>::value> p;

    Scalar penalty, t;
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
  NavierStokesFunctional(Scalar tStart_,Scalar rho_,
                       Scalar lambda_,Scalar mu_):
                       tStart(tStart_), rho(rho_),
                       lambda(lambda_), mu(mu_) {}

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
    static constexpr bool present   = true;     //  !(row==pIdx && col==pIdx);
    static constexpr bool symmetric = false;    // col==row;
    static constexpr bool lumped    = false;
  };

  template <int row, int col>
  struct B2 
  {
    static bool const present   = ((row==0)&&(col==0));
    static bool const symmetric = true;
    static bool const constant  = false;
    static bool const lumped    = false;
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
  Scalar tStart, t, tau; 
  Scalar rho, lambda, mu;

  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static int const components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& u0_,
                   NavierStokesFunctional<RType,AnsatzVars> const& f_): doScaling(doScaling_), u0(u0_), f(f_) {}

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
    NavierStokesFunctional<RType,AnsatzVars> const& f;
  };

};
#endif
