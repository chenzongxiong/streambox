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

#ifndef SEMIEULER_HH
#define SEMIEULER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

namespace Kaskade
{
  /**
   * \ingroup timestepping
   * \brief Linearly implicit Euler method.
   *
   * A class that defines the weak formulation of the stationary
   * elliptic problem that results from the linearly implicit Euler
   * method for
   * \f[
   * B \dot u  = L(u) , \f]
   * which is, doing the linearization at \f$ u_0 \f$,  \f[
   * (B - \tau L'(u_0)) \delta u_k  = \tau  L(u_k), \quad u_{k+1} = u_k + \delta u_k.
   * \f]
   *
   * \f$ B \f$ may depend on \f$ u \f$, made explicit by setting \code
   * B2<row,col>::constant = false \endcode in the parabolic equation
   * class. In this case, \f$ B(u) \f$ has to be an invertable Nemyckii operator (a
   * diagonal operator, mapping \f$ u(x) \f$ to \f$ B(x,u(x)) u(x)
   * \f$. Then the linearly implicit Euler scheme is (see below for derivation)
   * \f[
   * (B(u_0) - \tau L'(u_0) + \tau_k B^{-1}(u_0)B'(u_0)\mathrm{diag}L(u_0)) \delta u_k  = \tau_k  L(u)+\frac{\tau_k}{\tau_{k-1}}(B(u_0)-B(u_k))\delta u_{k-1}.
   * \f]
   * The values \f$ u_0 \f$ (linearization point for the Jacobian), \f$ u_k \f$ (evaluation point for the right hand side), 
   * and \f$ \tau_k/\tau_{k-1}\, \delta u_{k-1}  \f$ (previous increment) have to be specified when constructing the DomainCache.
   * This is most conveniently done via \ref SemiLinearizationAt.
   * 
   * \tparam PE the weak formulation of the parabolic equation
   * 
   * \em Derivation: Lubich and Roche (Rosenbrock Methods for Differential-algebraic Systems with 
   * Solution-dependent Singular Matrix Multiplying the Derivative, Computing 43:325-342, 1990) suggest
   * to reformulate the parabolic system with non-constant \f$ B \f$ above as
   * \f[ \begin{aligned} \dot u &= z \\ 0 &= L(u)-B(u)z \end{aligned} \f]
   * and applying, e.g., Rosenbrock methods. Here, simply a linearly implicit Euler. Linearizing
   * the right hand side of the extended system at \f$ u_0, z_0 \f$ yields
   * \f[ \Bigg( \begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix}  
   *     - \tau_k \begin{bmatrix} 0 & I \\ L'(u_0)-B'(u_0)z_0 & -B(u_0) \end{bmatrix} \Bigg) \begin{bmatrix} \delta u_k \\ \delta z_k \end{bmatrix} 
   *     = \tau_k \begin{bmatrix} z_k \\ L(u_k)-B(u_k)z_k \end{bmatrix}. \f]
   * The first row immediately gives \f$ \delta u_k = \tau_k (z_k + \delta z_k) \f$, and the second 
   * one 
   * \f[ - \tau_k (L'(u_0)-B'(u_0)z_0) \delta u_k + \tau (B(u_0)\delta z_k + B(u_k)z_k) = \tau_k f(u_k). \f]
   * The second term on the left hand side is just 
   * \f[ \tau_k (B(u_0)(z_k+\delta z_k) + (B(u_k)-B(u_0))z_k) = B(u_0)\delta u_k +  (B(u_k)-B(u_0)) \tau_k z_k, \f]
   * hence we obtain
   * \f[ \left(B(u_0) - \tau_k(L'(u_0)-B'(u_0)z_0)\right)\delta u_k = \tau_k L(u_k)+(B(u_0)-B(u_k)) \frac{\tau_k}{\tau_{k-1}}\delta u_{k-1}. \f]
   * Choosing a consistent initial value \f$ z_0 = B(u_0)^{-1} L(u_0) \f$ finally results in the form given above.
   */
  template <class PE>
  class SemiImplicitEulerStep      
  {
  public:
    typedef PE ParabolicEquation;

    typedef typename ParabolicEquation::Scalar  Scalar;
    typedef typename ParabolicEquation::AnsatzVars AnsatzVars;
    typedef typename ParabolicEquation::TestVars TestVars;
    typedef typename ParabolicEquation::OriginVars OriginVars;
    typedef typename AnsatzVars::Grid Grid;

    static ProblemType const type = WeakFormulation;

    
     
   
    typedef PE  EvolutionEquation;

    typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,EvolutionEquation::TestVars::noOfVariables>::type CoefficientVectors;

    
    
    class DomainCache
    {
    public:
      /**
       * \brief Constructor.
       * \param f the weak formulation
       * \param vars evaluation point \f$ u_k \f$ for the right hand side
       * \param varsJ linearization point \f$ u_0 \f$ for the (extended) Jacobian
       * \param duVars previous increment \f$ \delta u_{k-1} \f$ (is irrelevant if \f$ B \f$ does not depend on \f$ u \f$)
       * \param tauRatio ratio of current and previous time step \f$ \tau_k / \tau_{k-1} \f$
       */
      DomainCache(SemiImplicitEulerStep<ParabolicEquation> const& f_, // TODO: implement tauRatio
                  typename AnsatzVars::VariableSet const& vars_,
                  typename AnsatzVars::VariableSet const& varsJ_,
                  typename AnsatzVars::VariableSet const& duVars_,
                  int flags):
        f(f_), pedc(*f_.eq,vars_,flags), pedcJ(*f_.eq,varsJ_,flags), duVars(duVars_) 
      {} 

      void moveTo(typename Grid::template Codim<0>::Entity const& entity) {
        pedc.moveTo(entity);
        pedcJ.moveTo(entity);
      }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators) 
      {
        pedc.evaluateAt(x, evaluators);
        pedcJ.evaluateAt(x, evaluators);
        du = evaluateVariables(duVars,evaluators,valueMethod);
      }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
      d1 (VariationalArg<Scalar,dim> const& arg) const 
      {
        Dune::FieldVector<Scalar, TestVars::template Components<row>::m> result( f.getTau() * pedc.template d1<row>(arg) );

        boost::fusion::for_each(typename AnsatzVars::Variables(),MultiplyDifference<row,dim>(pedc,pedcJ,du,arg,result));
        return result;
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const
      {
        Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);

        if (f.onlyJacobian || f.onlyStiffnessMatrix || f.onlyF_uMatrix)
        {
          if (ParabolicEquation::template D2<row,col>::present)
            return pedcJ.template d2<row,col>(arg1, arg2);
          else { result = 0; return result;}
        }

        if (f.onlyMassMatrix)
        {
          if (ParabolicEquation::template B2<row,col>::present)
            return pedc.template b2<row,col>(arg1, arg2);
          else { result = 0; return result;}
        }

        if (ParabolicEquation::template B2<row,col>::present)
          result += pedcJ.template b2<row,col>(arg1, arg2);

        if (ParabolicEquation::template D2<row,col>::present)
          result -= f.getTau()*pedcJ.template d2<row,col>(arg1, arg2);
        
        return result;
      }

    private:
      SemiImplicitEulerStep<ParabolicEquation> const& f;
      typename ParabolicEquation::DomainCache pedc;  // for evaluating rhs at u_k
      typename ParabolicEquation::DomainCache pedcJ; // for evaluating Jacobian at u_0
      typename AnsatzVars::VariableSet const& duVars;
      EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> du;


      template <int row, int dim>
      class MultiplyDifference
      {
      public:
        MultiplyDifference(typename ParabolicEquation::DomainCache pedc_,
                           typename ParabolicEquation::DomainCache pedcJ_,
                           EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> const& du_,
                           VariationalArg<Scalar,dim> const& arg_,
                           Dune::FieldVector<Scalar, TestVars::template Components<row>::m>& result_):
              pedc(pedc_), pedcJ(pedcJ_), du(du_), arg(arg_), result(result_)
        {}

        template <class VarDescription>
        void operator()(VarDescription const& vd) const {
          // compute <(B(u0)-B(uk))*du,arg>
          int const col = VarDescription::id;

          if (ParabolicEquation::template B2<row,col>::present && !ParabolicEquation::template B2<row,col>::constant) 
          {
            // @TODO: should work with vector-valued shape functions
            VariationalArg<Scalar,dim> duArg;
            duArg.value[0] = 1;
            result += (pedcJ.template b2<row,col>(arg,duArg)-pedc.template b2<row,col>(arg,duArg)) * boost::fusion::at_c<col>(du);
          }
        }

      private:
        typename ParabolicEquation::DomainCache const& pedc;
        typename ParabolicEquation::DomainCache const& pedcJ;
        EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> const& du;
        VariationalArg<Scalar,dim> const& arg;
        Dune::FieldVector<Scalar, TestVars::template Components<row>::m>& result;
      };

    };

    // TODO: currently  a dependence of the boundary part of B on u is not implemented
    class BoundaryCache
    {
    public:
//       typedef typename AnsatzVars::Grid::LeafIntersectionIterator FaceIterator;
      typedef typename AnsatzVars::GridView::IntersectionIterator FaceIterator;

      BoundaryCache(SemiImplicitEulerStep<ParabolicEquation> const& f_,
                    typename AnsatzVars::VariableSet const& vars_,
                    typename AnsatzVars::VariableSet const& varsJ_,
                    typename AnsatzVars::VariableSet const& du_,
                    int flags)
      : f(f_), pedc(*f_.eq,vars_,flags), pedcJ(*f_.eq,varsJ_,flags), du(du_)  {}

      void moveTo(FaceIterator const& entity) 
      {
        pedc.moveTo(entity);
        pedcJ.moveTo(entity);
      }

      template <class Evaluators>
      void evaluateAt(Dune::FieldVector<typename Grid::ctype,Grid::dimension-1> const& x,
                      Evaluators const& evaluators) 
      {
        pedc.evaluateAt(x, evaluators);
        pedcJ.evaluateAt(x, evaluators);
      }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m> d1 (VariationalArg<Scalar,dim> const& arg) const 
      {
        return f.getTau() * pedc.template d1<row>(arg);
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
      {
        Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m> result(0);

	if ( f.onlyStiffnessMatrix )
        {
          if (ParabolicEquation::template D2<row,col>::present)
	  {
            result = pedcJ.template d2<row,col>(arg1, arg2);
	    result *= -1;
	    return result;
	  }
          else { result = 0; return result;}
        }
	
        if (ParabolicEquation::template D2<row,col>::present)
          result -= f.getTau()*pedcJ.template d2<row,col>(arg1, arg2);

        return result;
      }

    private:
      SemiImplicitEulerStep<ParabolicEquation> const& f;
      typename ParabolicEquation::BoundaryCache pedc;
      typename ParabolicEquation::BoundaryCache pedcJ;
      typename AnsatzVars::VariableSet const& du;
    };



  public:
    template <int row, int col>
    struct D2: public ParabolicEquation::template D2<row,col>
    {
      static bool const present = (ParabolicEquation::template D2<row,col>::present ||
          ParabolicEquation::template B2<row,col>::present);
      static bool const symmetric = (!ParabolicEquation::template D2<row,col>::present ||
          ParabolicEquation::template D2<row,col>::symmetric) &&
          (!ParabolicEquation::template B2<row,col>::present ||
              ParabolicEquation::template B2<row,col>::symmetric);
      //     static bool const constant = (!ParabolicEquation::template D2<row,col>::present ||
      //                                   ParabolicEquation::template D2<row,col>::constant) &&
      //                                  (!ParabolicEquation::template B2<row,col>::present ||
      //                                   ParabolicEquation::template B2<row,col>::constant);
      static bool const lumped = ((ParabolicEquation::template D2<row,col>::lumped ||
          !ParabolicEquation::template D2<row,col>::present) &&
          (ParabolicEquation::template B2<row,col>::lumped||
              !ParabolicEquation::template B2<row,col>::present));

    };

    template <int row>
    struct D1: public ParabolicEquation::template D1<row>
    {
    };

    SemiImplicitEulerStep(ParabolicEquation *eq_, Scalar dt) : eq(eq_), tau(dt), onlyJacobian(false), 
                                                               onlyMassMatrix(false), 
                                                               onlyStiffnessMatrix(false), onlyF_uMatrix(false) {}

    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return eq->integrationOrder(cell,shapeFunctionOrder,boundary);
    }

    void setTau(Scalar dt) { tau = dt; }
    Scalar getTau() const { return tau; }


    void newTime(Scalar dt) { eq->newTime(dt); }

    double time() const { return eq->time(); }
    void time(double t) const { eq->time(t); }
    
    void setOnlyJacobian(bool p) { onlyJacobian = p; }
    void setOnlyMassMatrix(bool p) { onlyMassMatrix = p; }
    void setOnlyStiffnessMatrix(bool p) { onlyStiffnessMatrix = p; eq->onlyStiffnessMatrix=p;}
    void setOnlyF_uMatrix(bool p) { onlyF_uMatrix = p; eq->onlyF_uMatrix=p;}
    
    void setSecondRHS(bool p){ eq->secondRHS=p;}

    void setDirichletBd(bool p){ eq->dirichletBd=p;}
    
    void setR(CoefficientVectors v, CoefficientVectors w) 
    {
      CoefficientVectors u(v);
      u-=w;
      eq->r=u;
    }
    

    template <int row>
    void pointwiseEuler(typename AnsatzVars::VariableSet& u) const 
    { 
      eq->template pointwiseEuler<row>(u); 
    }

    Scalar f(Scalar u, Scalar v)  const {return eq->f(u,v);}
    Scalar fu(Scalar u, Scalar v) const {return eq->fu(u,v);}
    
    Scalar g(Scalar u, Scalar v)  const {return eq->g(u,v);}
    Scalar gv(Scalar u, Scalar v) const {return eq->gv(u,v);}



    ParabolicEquation& parabolicEquation() { return *eq; }

    /**
     * Returns the scaling for error estimation.
     */
    typename ParabolicEquation::Scaling scaling() const { return eq->scaling(); }

    void temporalEvaluationRange(double t0, double t1) { eq->temporalEvaluationRange(t0,t1); }


    private:
    ParabolicEquation *eq;
    Scalar tau;
    bool onlyJacobian, onlyMassMatrix, onlyStiffnessMatrix, onlyF_uMatrix;
  };
} // namespace Kaskade
#endif
