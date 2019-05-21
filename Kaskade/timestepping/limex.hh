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

#ifndef LIMEX_HH
#define LIMEX_HH

#include <boost/timer/timer.hpp>

#include "timestepping/extrapolation.hh"
#include "fem/iterate_grid.hh"
#include "timestepping/semieuler.hh"

namespace Kaskade
{
  /**
   * This class implements the extrapolated linearly implicit Euler
   * method for integrating time-dependent evolution problems. The
   * implementation follows Deuflhard/Bornemann Chapter 6.4.3.
   *
   * \todo Currently the matrix \f$ B \f$ has to be constant.
   */
  template <class Eq>
  class Limex
  {
  public:
    typedef Eq                                                     EvolutionEquation;
    typedef typename EvolutionEquation::AnsatzVars::VariableSet State;

  private:
    typedef SemiLinearizationAt<SemiImplicitEulerStep<EvolutionEquation> > Linearization;
    typedef VariationalFunctionalAssembler<Linearization> GOp;

  public:
    /**
     * Constructs an ODE integrator. The arguments eq and ansatzVars
     * have to exist during the lifetime of the integrator.
     */
    Limex(GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager,
        EvolutionEquation& eq_, typename EvolutionEquation::AnsatzVars const& ansatzVars_,
        std::vector<std::pair<double,double> > const& tolX):
          ansatzVars(ansatzVars_), eq(&eq_,0), gop(gridManager.signals,ansatzVars.spaces), extrap(0),
          rhsAssemblyTime(0.0), matrixAssemblyTime(0.0), factorizationTime(0.0), solutionTime(0.0) {}

    /**
     * Computes a state increment that advances the given state in
     * time. The time in the given evolution equation is increased by
     * dt.
     *
     * \param x the initial state to be evolved
     * \param dt the time step
     * \param order the extrapolation order
     * \return the state increment (references an internal variable that will be invalidated by a subsequent call of step)
     *
     * \todo (i) check for B constant, do not reassemble matrix in this case (ii) implement fixed point
     * iteration instead of new factorization in case B is not constant
     */
    State const& step(State const& x, double dt, int order)
    {
      boost::timer::cpu_timer timer;

      std::vector<double> stepFractions(order+1);
      for (int i=0; i<=order; ++i) stepFractions[i] = 1.0/(i+1);
      extrap.clear();

      int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
      int const neq = EvolutionEquation::TestVars::noOfVariables;
      size_t nnz = gop.nnz(0,neq,0,nvars,false);
      size_t size = ansatzVars.degreesOfFreedom(0,nvars);

      std::vector<int> ridx(nnz), cidx(nnz);
      std::vector<double> data(nnz), rhs(size), sol(size);

      State dx(x), dxsum(x), tmp(x);
      double const t = eq.time();

      eq.temporalEvaluationRange(t,t+dt);

      for (int i=0; i<=order; ++i) {
        double const tau = stepFractions[i]*dt;
        eq.setTau(tau);
        eq.time(t);

        // Evaluate and factorize matrix B(t)-tau*J
        dx *= 0;
        timer.start();
        gop.assemble(Linearization(eq,x,x,dx));
        matrixAssemblyTime += (double)(timer.elapsed().user)/1e9;

        timer.start();
        gop.toTriplet(0,neq,0,nvars,ridx.begin(),cidx.begin(),data.begin(),false);
        UMFFactorization<double> matrix(size,0,ridx,cidx,data);
        factorizationTime += timer.elapsed();

        // First right hand side (j=0) has been assembled together with matrix.
        timer.start();
        gop.toSequence(0,neq,rhs.begin());
        for (int k=0; k<rhs.size(); ++k) assert(finite(rhs[k]));
        matrix.solve(rhs,sol);
        for (int k=0; k<sol.size(); ++k) assert(finite(sol[k]));
        dx.read(sol.begin());
        dxsum = dx;
        solutionTime += (double)(timer.elapsed().user)/1e9;

        // propagate by linearly implicit Euler
        for (int j=1; j<=i; ++j) {
          // Assemble new right hand side tau*f(x_j)+(B(t)-B(t+j*tau))*dx_(j-1)
          eq.time(eq.time()+tau);
          tmp = x; tmp += dxsum;
          timer.start();
          gop.assemble(Linearization(eq,tmp,x,dx),GOp::RHS);
          rhsAssemblyTime += timer.elapsed();
          timer.start();
          gop.toSequence(0,neq,rhs.begin());
          for (int k=0; k<rhs.size(); ++k) assert(finite(rhs[k]));
          matrix.solve(rhs,sol);
          for (int k=0; k<sol.size(); ++k) assert(finite(sol[k]));
          dx.read(sol.begin());
          dxsum += dx;
          solutionTime += (double)(timer.elapsed().user)/1e9;
        }

        // insert into extrapolation tableau
        extrap.push_back(dxsum,stepFractions[i]);

        // restore initial time
        eq.time(t);
      }

      return extrap.back();
    }

    double estimateError(State const& x,int i, int j) const
    {
      assert(extrap.size()>1);

      std::vector<std::pair<double,double> > e(ansatzVars.noOfVariables);

      relativeError(typename EvolutionEquation::AnsatzVars::Variables(),extrap[i].data,
                    extrap[j].data,x.data,ansatzVars.spaces,eq.scaling(),
                    e.begin());

      return e[0].first/(0.1+e[0].second);
    }

    template <class OutStream>
    void reportTime(OutStream& out) const {
      out << "Limex time: " << matrixAssemblyTime << "s matrix assembly\n"
          << "            " << rhsAssemblyTime << "s rhs assembly\n"
          << "            " << factorizationTime << "s factorization\n"
          << "            " << solutionTime << "s solution\n";
    }

    void advanceTime(double dt) { eq.time(eq.time()+dt); }


  private:
    typename EvolutionEquation::AnsatzVars const& ansatzVars;
    SemiImplicitEulerStep<EvolutionEquation>  eq;
    GOp                                       gop;

  public:
    ExtrapolationTableau<State>  extrap;
    double rhsAssemblyTime, matrixAssemblyTime, factorizationTime, solutionTime;
  };
} // namespace Kaskade
#endif
