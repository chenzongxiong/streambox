#ifndef LIMEX_HH
#define LIMEX_HH

#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <boost/timer/timer.hpp>

#include "dune/istl/solvers.hh"
#include "dune/istl/preconditioners.hh"

#include "fem/embedded_errorest.hh"
#include "fem/iterate_grid.hh"
#include "timestepping/extrapolation.hh"
#include "timestepping/semieuler.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"

#include "linalg/factorization.hh"
#include "linalg/umfpack_solve.hh"
#include "linalg/mumps_solve.hh"
#include "linalg/superlu_solve.hh"

#include "linalg/trivialpreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/iluprecond.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/hyprecond.hh"
#include "linalg/jacobiPreconditioner.hh"

#include "utilities/enums.hh"

DirectType xyz = DirectType::UMFPACK;

namespace Kaskade
{
  /**
   * \ingroup timestepping
   * \brief Extrapolated linearly implicit Euler method.
   *
   * This class implements the extrapolated linearly implicit Euler
   * method for integrating time-dependent evolution problems. The
   * implementation follows Deuflhard/Bornemann Chapter 6.4.3.
   */
  template <class Eq>
  class Limex
  {
  public:
    typedef Eq                                                     EvolutionEquation;
    typedef typename EvolutionEquation::AnsatzVars::VariableSet State;

  private:
    typedef SemiLinearizationAt<SemiImplicitEulerStep<EvolutionEquation> > Linearization;
    typedef VariationalFunctionalAssembler<Linearization> Assembler;

  public:
    /**
     * Constructs an ODE integrator. The arguments eq and ansatzVars
     * have to exist during the lifetime of the integrator.
     */
    Limex(GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager_,
        EvolutionEquation& eq_, typename EvolutionEquation::AnsatzVars const& ansatzVars_, DirectType st=DirectType::MUMPS, PrecondType precondType_ = PrecondType::ILUK,
          int verbosity_ = 1):
          gridManager(gridManager_), ansatzVars(ansatzVars_), eq(&eq_,0),
          assembler(gridManager,ansatzVars.spaces), extrap(0),
          rhsAssemblyTime(0.0), matrixAssemblyTime(0.0), precAssemblyTime(0.0), factorizationTime(0.0), solutionTime(0.0),
          initSolutionTime(0.0), estimateTime(0), directType(st), precondType(precondType_), verbosity(verbosity_)
    {}


    /**
     * Computes a state increment that advances the given state in
     * time. The time in the given evolution equation is increased by
     * dt.
     *
     * \param x the initial state to be evolved
     * \param dt the time step
     * \param order the extrapolation order >= 0 (0 corresponds to the linearly implicit Euler)
     * \param tolX
     * \return the state increment (references an internal variable that will be
     *         invalidated by a subsequent call of step)
     *
     * \todo (i) check for B constant, do not reassemble matrix in this case
     *       (ii) implement fixed point iteration instead of new factorization
     *            in case B is not constant
     */
    State const& step(State const& x, double dt, int order,
        std::vector<std::pair<double,double> > const& tolX)
    {
      boost::timer::cpu_timer timer;

      std::vector<double> stepFractions(order+1);
      for (int i=0; i<=order; ++i) stepFractions[i] = 1.0/(i+1); // harmonic sequence
      extrap.clear();

      int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
      int const neq = EvolutionEquation::TestVars::noOfVariables;
      size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
      size_t  size = ansatzVars.degreesOfFreedom(0,nvars);

      MatrixAsTriplet<double> triplet(nnz);
      std::vector<double> rhs(size), sol(size);

      typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::type CoefficientVectors;

      typedef  AssembledGalerkinOperator<Assembler,0,neq,0,neq> Op;
      State dx(x), dxsum(x), tmp(x);
      double const t = eq.time();

      eq.temporalEvaluationRange(t,t+dt);

      bool const iterative = true;      // false: use a direct solver,      true: use an iterative solver

      for (int i=0; i<=order; ++i) {
        double const tau = stepFractions[i]*dt;
        eq.setTau(tau);
        eq.time(t);

        // mesh adaptation loop. Just once if i>0.
        gridManager.setVerbosity(verbosity);
        bool accurate = true;
        do {
          int dimension = gridManager.grid().dimension;

          // Evaluate and factorize matrix B(t)-tau*J
          dx *= 0;
          timer.start();
          assembler.assemble(Linearization(eq,x,x,dx));
          matrixAssemblyTime += (double)(timer.elapsed().user)/1e9;
          Op A(assembler);

          timer.start();
          if (iterative) {
            CoefficientVectors rhs(assembler.rhs());
            typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::type LinearSpaceX;
            //depricated: CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
            CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars.spaces));
            solution = 0;

            if ( verbosity>0 )
            {
              std::cout << " eq 0: dof = " << boost::fusion::at_c<0>(dx.data).space().degreesOfFreedom()
                              << " pts in mesh = " << gridManager.grid().size(dimension) << std::endl;
              if ( verbosity>1 ) std::cout << " dimension of grid = " << dimension << std::endl;
            };
            
            initSolutionTime += (double)(timer.elapsed().user)/1e9;
            elapsedTimeSinceReset = (double)(timer.elapsed().user)/1e9;

            // in der Folge finden sich vesrschiedene iterative L�ser
            // und Vorkonditionierer zur Auswahl per Ein-/Auskommentieren.
            // Der Mechanismus muss noch umgestaltet werden mittels
            // passender Parameter aus dem Main()-Programm via itegrate.hh !!!!

            int iteSteps = 5000;
            double iteEps = 1.0e-10;
            Dune::InverseOperatorResult res;

            switch (precondType)
            {
            case PrecondType::NONE:
            {
              TrivialPreconditioner<Op> trivial;
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,trivial,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::BOOMERAMG:
            {
              int steps = iteSteps;
              int coarsentype = 21;
              int interpoltype = 0;
              int cycleType = 1;
              int relaxType = 3;
              int variant = 0;
              int overlap = 1;
              int syseqn = 1;
              double tol = iteEps;
              double strongThreshold = (dimension==2)?0.25:0.6;
              BoomerAMG<Op> BoomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType, 
                            strongThreshold,variant,overlap,syseqn,verbosity);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,BoomerAMGPrecon,iteEps,iteSteps,verbosity-1);       // 1e-10
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::ILUK:
            {
              int fill_lev = 1;
              ILUKPreconditioner<Op> iluk(A,fill_lev,verbosity);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              //Dune::CGSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity-1);
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,iluk,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::JACOBI:
            default:
            {
              JacobiPreconditioner<Op> jacobi(A,1.0);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,jacobi,iteEps,iteSteps,verbosity-1);
              // alternative:  Dune::LoopSolver<LinearSpaceX> cg(A,jacobi,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            }

            if ( !(res.converged) || (res.iterations == 2001) ) {
              std::cout << "   no of iterations in cg = " << res.iterations << std::endl;
              std::cout << " convergence status of cg = " << res.converged  << std::endl;
              //std::cout << "A: "  << std::endl;
              //for (size_t ii=0; ii<nnz; ii++)
              //    std::cout   << A.getmat().ridx[ii] << ", "
              //                << A.getmat().cidx[ii] << ", "
              //                << A.getmat().data[ii] << std::endl;
              //std::cout << "A, end: "  << std::endl;
              assert(0);
            }

            dx.data = solution.data;
          }
          else {
            triplet = A.template get<MatrixAsTriplet<double> >();
            Factorization<double> *matrix = 0;
            switch (directType) {
            case DirectType::UMFPACK:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::UMFPACK3264:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::MUMPS:
              matrix = new MUMPSFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::SUPERLU:
              matrix = new SUPERLUFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            default:
              throw -321;
              break;
            }
            matrix->setVerbose(verbosity);
            factorizationTime += (double)(timer.elapsed().user)/1e9;

            // First right hand side (j=0) has been assembled together with matrix.
            timer.start();
            A.getAssembler().toSequence(0,neq,rhs.begin());
            for (int k=0; k<rhs.size(); ++k) assert(std::isfinite(rhs[k]));
            matrix->solve(rhs,sol);
            delete matrix;
            for (int k=0; k<sol.size(); ++k) assert(std::isfinite(sol[k]));
            dx.read(sol.begin());
          }
          solutionTime += (double)(timer.elapsed().user)/1e9;


          timer.start();
          // Mesh adaptation only if requested and only for the full
          // implicit Euler step (first stage).
          if (!tolX.empty() && i==0) {
            // Compute spatial error estimate
            tmp = dx;
            projectHierarchically(ansatzVars,tmp);
            tmp -= dx;

            // perform mesh adaptation
            accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,verbosity);
            if (!accurate) {
              nnz = assembler.nnz(0,neq,0,nvars,false);
              size = ansatzVars.degreesOfFreedom(0,nvars);
              rhs.resize(size);
              sol.resize(size);
            }
            estimateTime += (double)(timer.elapsed().user)/1e9;
          }
        } while (!accurate);

        dxsum = dx;

        // propagate by linearly implicit Euler
        if (order==0) 
        {
          // insert into extrapolation tableau
          extrap.push_back(dxsum,stepFractions[i]);
          // restore initial time
          eq.time(t);
          if (verbosity>1) std::cout << "   end of limex step" << std::endl;
          return extrap.back();
        }
        else
        for (int j=1; j<=i; ++j) {
          // Assemble new right hand side tau*f(x_j)
          eq.time(eq.time()+tau);
          tmp = x; tmp += dxsum;
          State zero(x); zero *= 0;
          timer.start();
          assembler.assemble(Linearization(eq,tmp,x,zero),Assembler::RHS|Assembler::MATRIX);
          rhsAssemblyTime += (double)(timer.elapsed().user)/1e9;

          //size_t dof = 2*boost::fusion::at_c<0>(dx.vars).space().degreesOfFreedom()
          //            - gridManager.grid().size(2);
          Op A(assembler);

          timer.start();
          if (iterative) {
            CoefficientVectors rhs(assembler.rhs());
            typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::type LinearSpaceX;
            //depricatedCoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
            CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars.spaces));
            solution = 0;

            initSolutionTime += (double)(timer.elapsed().user)/1e9;
            elapsedTimeSinceReset = (double)(timer.elapsed().user)/1e9;


            int iteSteps = 5000;
            double iteEps = 1.0e-10;
            Dune::InverseOperatorResult res;

            switch (precondType)
            {
            case PrecondType::NONE:
            {
              TrivialPreconditioner<Op> trivial;
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,trivial,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::BOOMERAMG:
            {
              int steps = iteSteps;
              int coarsentype = 21;
              int interpoltype = 0;
              int cycleType = 1;
              int relaxType = 3;
              int variant = 0;
              int overlap = 1;
              int syseqn = 1;
              double tol = iteEps;
              int dimension = gridManager.grid().dimension;
              double strongThreshold = (dimension==2)?0.25:0.6;
              BoomerAMG<Op> BoomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
                            strongThreshold,variant,overlap,syseqn,verbosity);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,BoomerAMGPrecon,iteEps,iteSteps,verbosity-1);       // 1e-10
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::ILUK:
            {
              int fill_lev = 1;
              ILUKPreconditioner<Op> iluk(A,fill_lev,verbosity);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              //Dune::CGSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity-1);
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,iluk,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            case PrecondType::JACOBI:
            default:
            {
              JacobiPreconditioner<Op> jacobi(A,1.0);
              precAssemblyTime += (double)(timer.elapsed().user)/1e9-elapsedTimeSinceReset;
              Dune::BiCGSTABSolver<LinearSpaceX> cg(A,jacobi,iteEps,iteSteps,verbosity-1);
              cg.apply(solution,rhs,res);
            }
            break;
            }

            if ( !(res.converged) || (res.iterations == 2001) ) {
              typedef typename Op::field_type field_type;
              MatrixAsTriplet<field_type> triplet = A.template get<MatrixAsTriplet<field_type> >();

              std::cout << "   no of iterations in cg = " << res.iterations << std::endl;
              std::cout << " convergence status of cg = " << res.converged  << std::endl;
              std::cout << "A: "  << std::endl;
              for (size_t ii=0; ii<nnz; ii++)
                std::cout   << triplet.ridx[ii] << ", " << triplet.cidx[ii] << ", " << triplet.data[ii] << std::endl;
              std::cout << "A, end: "  << std::endl;
              assert(0);
            }
            dx.data = solution.data;
          }
          else {
            A.getAssembler().toSequence(0,neq,rhs.begin());
            MatrixAsTriplet<double> triplet = A.template get<MatrixAsTriplet<double> >();
            Factorization<double> *matrix = 0;
            switch (directType) {
            case DirectType::UMFPACK:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::UMFPACK3264:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::MUMPS:
              matrix = new MUMPSFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::SUPERLU:
              matrix = new SUPERLUFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            default:
              throw -321;
              break;
            }
            matrix->setVerbose(verbosity);
            matrix->solve(rhs,sol);
            delete matrix;
            for (int k=0; k<rhs.size(); ++k) assert(std::isfinite(rhs[k]));
            for (int k=0; k<sol.size(); ++k) assert(std::isfinite(sol[k]));
            dx.read(sol.begin());
          }

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

    /**
     * Estimates the time discretization error of the previously
     * computed step by taking the difference between the diagonal and
     * subdiagonal extrapolation values of maximal order. This requires
     * that order>1 has been given for the last step.
     */
    std::vector<std::pair<double,double> > estimateError(State const& x,int i, int j) const
      {
      assert(extrap.size()>1);

      std::vector<std::pair<double,double> > e(ansatzVars.noOfVariables);

      relativeError(typename EvolutionEquation::AnsatzVars::Variables(),extrap[i].data,
          extrap[j].data,x.data,
          ansatzVars.spaces,eq.scaling(),e.begin());

      return e;
      }

    template <class OutStream>
    void reportTime(OutStream& out) const {
      out << "Limex time: " << matrixAssemblyTime << "s matrix assembly\n"
          << "            " << rhsAssemblyTime << "s rhs assembly\n"
          << "            " << precAssemblyTime << "s preconditioner assembly\n"
          << "            " << factorizationTime << "s factorization\n"
          << "            " << initSolutionTime << "s init solution\n"
          << "            " << solutionTime << "s solution\n"
          << "            " << estimateTime << "s estimate\n"
          << std::flush;
    }

    void advanceTime(double dt) { eq.time(eq.time()+dt); }


  private:
    GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager;
    typename EvolutionEquation::AnsatzVars const& ansatzVars;
    SemiImplicitEulerStep<EvolutionEquation>      eq;
    Assembler                                     assembler;

  public:
    ExtrapolationTableau<State>  extrap;
    double rhsAssemblyTime, matrixAssemblyTime, precAssemblyTime, factorizationTime, solutionTime;
    double elapsedTimeSinceReset, initSolutionTime, estimateTime;
    DirectType directType;
    PrecondType precondType;
    int verbosity;
  };
} // namespace Kaskade
#endif
