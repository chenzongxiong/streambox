#include <cmath>
#include <fstream>
#include <string>

#include <boost/timer/timer.hpp>

#include "optimization.hh"
#include "adaptive_refinement.hh"
#include "linalg/linearsystem.hh"
#include "linalg/triplet.hh"
#include "./opt_aux/src/include/Fmin.h"
#include "modelCreator.hh"
#include "modelFunctions.hh"

#ifndef Cygwin
using std::to_string;
#else
std::string to_string(size_t n)
{
  std::stringstream s; s << n; return s.str();
}
#endif

namespace Kaskade
{
  OptimizationParameters::OptimizationParameters(double desiredAccuracy_, int maxSteps_) : desiredAccuracy(desiredAccuracy_), maxSteps(maxSteps_) {}
  
  void OptimizationParameters::setThetaAim(double theta) { ThetaAim = theta; ensureAdmissibleThetas(); }
  void OptimizationParameters::setThetaNormal(double theta) { ThetaNormal = theta; ensureAdmissibleThetas(); }
  void OptimizationParameters::setThetaMax(double theta) { ThetaMax = theta; ensureAdmissibleThetas(); }
  void OptimizationParameters::setThetas(double thetaNormal, double thetaAim, double thetaMax)
  {
    ThetaNormal = thetaNormal;
    ThetaAim = thetaAim;
    ThetaMax = thetaMax;
    ensureAdmissibleThetas();
  }
  void OptimizationParameters::setEps(double eps_)
  {
    eps = eps_;
    sqrtEps = sqrt(eps);
    thirdSqrtEps = pow(eps,1./3.);
  }

  double OptimizationParameters::getThetaAim() const { return ThetaAim; }
  double OptimizationParameters::getThetaNormal() const { return ThetaNormal; }
  double OptimizationParameters::getThetaMax() const { return ThetaMax; }
  double OptimizationParameters::getEps() const { return eps; }
  double OptimizationParameters::getSqrtEps() const { return sqrtEps; }
  double OptimizationParameters::getThirdSqrtEps() const { return thirdSqrtEps; }

  void OptimizationParameters::ensureAdmissibleThetas()
  {
    if( !(ThetaAim < ThetaMax) ) 
    {
      std::cout << "OPTIMIZATION PARAMETERS: Warning: Inconsistent algorithmic parameters. (ThetaMax <= ThetaAIm)" << std::endl;
      std::cout << "OPTIMIZATION PARAMETERS: Adjusting ThetaAim from " << ThetaAim << " to ";
      ThetaAim = 0.9*ThetaMax;
      std::cout << ThetaAim << std::endl;
    } 
    if( !(ThetaNormal < ThetaAim) ) 
    {
      std::cout << "OPTIMIZATION PARAMETERS: Warning: Inconsistent algorithmic parameters. (ThetaAim <= ThetaNormal)" << std::endl;
      std::cout << "OPTIMIZATION PARAMETERS: Adjusting ThetaNormal from " << ThetaNormal << " to ";
      ThetaNormal = 0.9*ThetaAim;
      std::cout << ThetaNormal << std::endl;
    }     
  }
  
  Optimization::~Optimization(){}

  Optimization::Optimization(AbstractScalarProduct& nL, AbstractScalarProduct& nC, AbstractChart const& chart_, OptimizationParameters const& p_,
                             AbstractCompositeStepErrorEstimator* errorEstimator_, int verbose_)
    : normL(nL), normC(nC), chart(chart_.clone()), p(p_), errorEstimator(errorEstimator_),
      verbose(verbose_),L(verbose), C(verbose)
  {}

  Optimization::Optimization(AbstractNormalDirection& normalSolver,
                             AbstractTangentialSpace& tangentSpace_,
                             AbstractScalarProduct& nL,
                             OptimizationParameters const& p_, 
                             double omegaCinit, double omegaLinit, int verbose_):
    normL(nL),
    normC(nL),
    chart(new PrimalChart()),
    p(p_),
    normalDirection(&normalSolver),
    tangentSpace(&tangentSpace_),
    verbose(verbose_),
    L(verbose,omegaLinit), C(verbose,omegaCinit)
  {}


  Optimization::Optimization(AbstractNormalDirection& normalSolver,
                             AbstractScalarProduct& nL,
                             AbstractScalarProduct& nC,
                             AbstractChart const& chart_,
                             OptimizationParameters const& p_, int verbose_):
    normL(nL),
    normC(nC),
    chart(chart_.clone()),
    p(p_),
    normalDirection(&normalSolver),
    verbose(verbose_),
    L(verbose), C(verbose)
  {}

  /// Return true, if convergence is detected, false otherwise
  Convergence Optimization::convergenceTest(double nu, std::vector<double> const& tau, double normOfCorrection) const
  {
    if( verbose > 0 )
    {
      std::cout << csPre << "Desired accuracy: " << p.desiredAccuracy << ", ||x|| = " << normOfIterate << std::endl;
      std::cout << csPre << "||dx|| = " << normOfCorrection << std::endl;
    }

    //if( nu > 0.1 && normOfCorrection < p.getEps() ) return Convergence::Achieved;

    if( !tangentSpace->localConvergenceLikely() ) return Convergence::Missed;

    // no local convergence as long as at least one of the steps is damped
    if( !noDamping(nu) || !noDamping(tau) ) return Convergence::Missed;
    if( normOfCorrection < p.desiredAccuracy * normOfIterate ) return Convergence::Achieved;

    return Convergence::Missed;
  }


  bool Optimization::adaptiveMeshRefinement(LagrangeLinearization& lin, double nu, std::vector<double> tau, AbstractFunctionSpaceElement const& correction)
  {
    double safetyFactorForDesiredAccuracy = 0.9;

    if( !( noDamping(nu) && noDamping(tau) ) )
    {
      std::cout << csPre << "Damping factors != 1 -> Skipping error estimation." << std::endl;
      return false;
    }

    if(hbErrorEstimator->gridSize() > p.maxGridSize)
    {
      std::cout << csPre << "Grid size: " << hbErrorEstimator->gridSize() << " > " << p.maxGridSize << ". Skipping error estimation. " << std::endl;
      return false;
    }

      // estimate error
    (*hbErrorEstimator)(lin, *iterate, correction, step, *iterate);

    double absoluteError=hbErrorEstimator->estimatedAbsoluteError();
    double requiredAbsoluteError=std::max(p.requiredRelativeError*normOfLastCorrection,p.desiredAccuracy*safetyFactorForDesiredAccuracy*normL(*iterate));
    bool accuracyReached = absoluteError < requiredAbsoluteError;

    if(verbose > 0)
    {
      std::cout << csPre << "ERROR ESTIMATOR: err0=" << p.requiredRelativeError*normOfLastCorrection << ", err1=" << p.desiredAccuracy*safetyFactorForDesiredAccuracy << std::endl;
      std::cout << csPre << "ERROR ESTIMATOR: normOfAcceptedCorrection=" << normOfLastCorrection << ", p.desiredAccuracy=" << p.desiredAccuracy << std::endl;
      std::cout << csPre << "ERROR ESTIMATOR: absolute error: " << absoluteError << ", required error: " << requiredAbsoluteError << std::endl;
    }

    if(accuracyReached)
    {
      std::cout << "ERROR ESTIMATOR: Discretization accuracy reached!" << std::endl;
      return false;
    }

    hbErrorEstimator->refineGrid();
    return true;
  }



  int Optimization::runAlgorithm()
  {
    std::ofstream omegafile("omegas.log");
    std::ofstream coefffile("coeff.log");
    std::ofstream correctionfile("corrections.log");
    std::ofstream values("values.log");
    values.setf(std::ios::scientific,std::ios::floatfield);
    values.precision(16);
    std::ofstream steps("steps.log");

    double nu = 0;
    std::vector<double> tau(1,0);
    // Main loop
    for( step=1 ; step <= p.maxSteps; ++step )
    {
      if(verbose > 0) std::cout <<  "\n --------------- Optimization-Step: " << step << " ----------------------" << std::endl;

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * Normal Step * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      if( functionalN == nullptr ) std::cout << "no normal lineraization" << std::endl;
      normalLinearization = functionalN->getLinearization(*iterate);
      values << normalLinearization->eval() << std::endl;

      normC.setOrigin(*normalLinearization);
      normL.setOrigin(*normalLinearization);
      std::unique_ptr<AbstractFunctionSpaceElement> normalStepResidual(iterate->initZeroVector()), adjointResidual(iterate->initZeroVector()), undampedStep(iterate->initZeroVector());
      auto normalStepResult = computeNormalStep(normalStepResidual.get(),adjointResidual.get());

      std::unique_ptr<AbstractFunctionSpaceElement> normalStep(std::move(normalStepResult.first)), adjointCorrection(std::move(normalStepResult.second));
      undampedStep->axpy(1.0,*normalStep,"primal");
      double normOfUndampedCorrection = 0;

      double Lxdn_res = adjointResidual->applyAsDualTo(*normalStep,"primal")
                      - normalStepResidual->applyAsDualTo(*adjointCorrection,"primal")
                      + adjointResidual->applyAsDualTo(*normalStep,"dual");
      if( verbose > 1 ) std::cout << csPre << "Lxdn_res: " << Lxdn_res << std::endl;
      auto normNormal = normC( *normalStep );

      // Predict damping factor for normal step
      nu = (normNormal > 0) ? updateNormalStepDampingFactor(normNormal) : 1.;
      if( verbose > 0 ) printNormalStep(normNormal,nu);

      std::unique_ptr<AbstractFunctionSpaceElement> tmpIter(iterate->clone());
      tmpIter->axpy(nu,*normalStep,"primal");
      std::unique_ptr<AbstractLinearization> lin_x_nudn(functionalN->getLinearization(*tmpIter));

      std::unique_ptr<AbstractFunctionSpaceElement> dc_x0_dn ( iterate->initZeroVector() );
      normalLinearization->d2axpy(1.0,*dc_x0_dn, *normalStep, 2, 3, 0, 2);
      double p_dc_x0_dn = iterate->applyAsDualTo(*dc_x0_dn,"dual");
      double f_x1 = lin_x_nudn->eval();
      double f_x0 = normalLinearization->eval();

      double normalStepMonitor = ( ( std::fabs(f_x1 - f_x0) + std::fabs(p_dc_x0_dn) ) > 0 ) ? std::fabs( f_x1 - f_x0 + nu*p_dc_x0_dn - Lxdn_res ) / ( std::fabs(f_x1 - f_x0) + std::fabs(p_dc_x0_dn) ) : 0.;
      if( normNormal < p.getSqrtEps()*normOfIterate || step == 1 ) normalStepMonitor = 0;
      bool reliableQuadraticModel = true;//normalStepMonitor < 0.5;

      if( verbose > 1 )
      {
        std::cout << "f(x0) = " << f_x0 << ", f(x1) = " << f_x1 << ", nupc'(x0) = " << nu << "*" << p_dc_x0_dn << std::endl;
        std::cout << "f(x1) - f(x0) = " << f_x1 - f_x0 << std::endl;
      }
      if( verbose > 0 ) std::cout << csPre << "NORMAL STEP MONITOR: " << normalStepMonitor << ", is reliable: " << reliableQuadraticModel << std::endl;

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * Tangential Step * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      double normTangential = 0;
      normOfUndampedCorrection = normNormal;
      std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > tangentialBasis(1,iterate->clone());

      std::unique_ptr<LagrangeLinearization> lagrangeLinearization(nullptr);

      if(reliableQuadraticModel)
      {
        if(functionalT!=nullptr) tangentialLinearization = functionalT->getLinearization(*iterate);
        else tangentialLinearization = nullptr;

        lagrangeLinearization.reset(new LagrangeLinearization(normalLinearization.get(),tangentialLinearization.get(),*iterate));

        tangentialBasis = computeTangentialStep(*lagrangeLinearization,*normalStep, nu, tau);
        for(size_t i=0; i<tangentialBasis.size(); ++i) undampedStep->axpy(1.0,*tangentialBasis[i],"primal");
        normTangential = normL(*tangentialBasis[0]);
        normOfUndampedCorrection = normC(*undampedStep);
      }
      else
      {
        *tangentialBasis[0] *= 0;
        tau[0] = 0;
      }
      if( verbose > 0 ) std::cout << csPre << "Tangential step length: " << normTangential << std::endl;
      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * Lipschitz constants and damping * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      AcceptanceTest acceptanceTestResult = AcceptanceTest::Failed;
      std::unique_ptr<AbstractFunctionSpaceElement> trialIterate(iterate->initZeroVector()), secondOrderCorrected(iterate->clone()), correction(iterate->initZeroVector());
      double normCAtCorr = normNormal;
      while(acceptanceTestResult != AcceptanceTest::Passed)
      {
        double eta = 1;
        if( acceptanceTestResult == AcceptanceTest::LeftAdmissibleDomain ) nu *= 0.5;
        else nu = updateNormalStepDampingFactor(normNormal);
        if( verbose > 0 ) printNormalStep(normNormal,nu);
        if(!reliableQuadraticModel)
        {
          correction = createCorrection(nu,*normalStep, tau, tangentialBasis);
          normCAtCorr = normC(*correction);
          // Compute Trial Iterate: first candidate: trialIterate = chart_x(dx)
          *trialIterate *= 0;
          chart->addPerturbation(*trialIterate,*correction,*normalLinearization);

          std::unique_ptr<AbstractLinearization> lin_xplus(functionalN->getLinearization(*trialIterate));
          auto simplifiedNormalStep = computeSimplifiedNormalStep( *lin_xplus, normalStepResidual.get());
          simplifiedNormalStep->axpy( nu-1.0, *normalStep, "primal" );

          updateConstraintD1LipschitzConstant( normC( *simplifiedNormalStep ), normCAtCorr );
        }
        else
        {
          QuadraticModelCreator quadraticModelCreator(*normalStep,tangentialBasis,*lagrangeLinearization,Lxdn_res);
          NormModelCreator normModelCreatorL(*normalStep,tangentialBasis,normL);

          QuadraticFunction normLModel(normModelCreatorL.create(nu));
          QuadraticFunction quadraticModel(quadraticModelCreator.create(nu));
          L.setL_xx(quadraticModel.quadraticPart[0][0]);
          IsotropicCubicRegularization reg(normLModel,L.omega/6);
          RegularizedQuadraticFunction cubic(quadraticModel, reg);
          if( acceptanceTestResult == AcceptanceTest::LeftAdmissibleDomain ) tau[0] *= 0.5;
          else updateTangentialDampingFactor(nu, normNormal, normTangential, cubic, tau);

          // Compute Correction as a linear combination of normal step and basis vectors of tangent search space
          correction = createCorrection(nu,*normalStep,tau,tangentialBasis);
          auto normLAtCorr = sqrt(normLModel.d0(tau));
          normCAtCorr = normC(*correction);
          auto quadraticModelAtCorr = quadraticModel.d0(tau);
          // Compute Trial Iterate: first candidate: trialIterate = chart_x(dx)
          trialIterate = iterate->clone();
          chart->addPerturbation(*trialIterate,*correction,*lagrangeLinearization);
          if( !functionalN->inDomain(*trialIterate) || ( functionalT != nullptr && !functionalT->inDomain(*trialIterate) ) )
          {
            std::cout << csPre << "Iterate leaves admissible domain: rejecting step." << std::endl;
            acceptanceTestResult = AcceptanceTest::LeftAdmissibleDomain;
          }
          else
          {
            acceptanceTestResult = AcceptanceTest::Failed;
            std::unique_ptr<AbstractLinearization> lin_xplus(functionalN->getLinearization(*trialIterate));
            auto simplifiedNormalStep = computeSimplifiedNormalStep( *lin_xplus, normalStepResidual.get());
            simplifiedNormalStep->axpy( nu-1.0, *normalStep, "primal" );
            //for( size_t i = 0; i < tangentialBasis.size(); ++i ) simplifiedNormalStep->axpy( tau[i] - 1.0, *tangentialBasis[i], "primal" );

            updateConstraintD1LipschitzConstant( normC( *simplifiedNormalStep ), normCAtCorr );

            // Perform second order correction update secondOrderCorrected=chart_x(dx+ds)
            std::unique_ptr<AbstractFunctionSpaceElement> secondOrderCorrection(simplifiedNormalStep->clone());
            secondOrderCorrection->axpy(1.0,*correction,"primal");
            *secondOrderCorrected = *iterate->clone();
            chart->addPerturbation(*secondOrderCorrected,*secondOrderCorrection,*lagrangeLinearization);

            eta = updateLagrangianD2LipschitzConstant(*lagrangeLinearization, *secondOrderCorrected, normLAtCorr, quadraticModelAtCorr, cubic, nu, tau,*normalStep,*trialIterate,*simplifiedNormalStep,*correction,Lxdn_res);
          }
        }

        // Acceptance tests
        if( acceptanceTestResult != AcceptanceTest::LeftAdmissibleDomain ) acceptanceTestResult = acceptanceTest(eta,nu,tau,normCAtCorr);
        if(!regularityTest(nu,tau,reliableQuadraticModel)) return -1;

        if( reliableQuadraticModel && (acceptanceTestResult==AcceptanceTest::TangentialStepFailed) && L.omega < (1 + 0.25*(1 - p.etaMin))*L.oldOmega )
        {
          if( verbose > 0 ) std::cout << csPre << "Stagnating omegaL update. Accepting step." << std::endl;
          acceptanceTestResult = AcceptanceTest::Passed;
          if(eta < 0.1)
          {
            if( verbose > 0 ) std::cout << csPre << "Ignoring tangential step." << std::endl;
            for(size_t i=0; i<tangentialBasis.size(); ++i)
            {
              trialIterate->axpy(-tau[i],*tangentialBasis[i]);
              secondOrderCorrected->axpy(-tau[i],*tangentialBasis[i]);
            }
          }
        }

        if(acceptanceTestResult == AcceptanceTest::Passed)
        {
          normOfLastCorrection = normCAtCorr;
          normOfLastCorrection_Undamped = normOfUndampedCorrection;
          correctionfile << normOfLastCorrection << " " << normL(*normalStep) << " " << normL(*tangentialBasis[0]) << std::endl;
        }
      }

      // Error estimation
      bool refinedMesh = false;
      if( hbErrorEstimator != nullptr && tangentSpace->localConvergenceLikely()  && reliableQuadraticModel ) refinedMesh = adaptiveMeshRefinement(*lagrangeLinearization, nu,tau,*correction);

      // at this point, *trialIterate is the combined step
      // and *secondOrderCorrected is the second order corrected step
      // if Theta > 0.25, this might be worse w.r.t feasibility than the original step
      // use original step in this case
      if(normalDirection)
      {
        if(C.theta > 0.25 || !reliableQuadraticModel)
        {
          if(verbose > 1) std::cout << csPre << "Adding ordinary update" << std::endl;
        }
        else
        {
          *trialIterate = *secondOrderCorrected;
          if(verbose > 1) std::cout << csPre << "Adding second order correction" << std::endl;
        }
      }

      iterate->swap(*trialIterate);
      normOfIterate = normL(*iterate);
//      *trialIterate -= *iterate;

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * Write data  * * * * * * * * * * * * * * * * */
      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      if( verbose > 0)
      {
        std::cout << "\n" << csPre << "Accepted step" << std::endl;
        printNormalStep(normNormal, nu);
        printTangentialStep(normTangential, tau[0]);
      }

      if( verbose > 1 )
      {
        std::string name = std::string("iterate_") + to_string(step);
        iterate->writeToFile(name,false);
      }

      omegafile << " " << L.omega << " " << " " << C.omega << " " << std::endl;
      coefffile << nu << "  " << tau[0] << std::endl;

      if( !refinedMesh && convergenceTest(nu,tau,normL(*correction)) == Convergence::Achieved )
      {
        steps << step << std::endl;

        if( verbose > 1 )
        {
          std::cout << csPre << "overall normal step computation time: " << normalStepComputationTime << std::endl;
          std::cout << csPre << "overall tangential step computation time: " << tangentialStepComputationTime << std::endl;
        }
        return 1;
      }
    }
    return -2;
  }

  void Optimization::solve(AbstractFunctional &fN, AbstractFunctional& fT, AbstractFunctionSpaceElement &x) { solve(fN,&fT,x); }

  void Optimization::solve(AbstractFunctional& fN, AbstractFunctional* fT,AbstractFunctionSpaceElement& x)
  {
    functionalN = &fN;
    functionalT = fT;

    iterate = x.clone();
    normalStepComputationTime = tangentialStepComputationTime = 0;
    if( tangentSpace != nullptr ) tangentSpace->setEps(p.getEps());

    normOfLastCorrection = normOfLastCorrection_Undamped = -1;
    normOfIterate = 0;
    auto terminationFlag = algorithmWrapper();
    if(terminationFlag > 0 || terminationFlag == -2) x=*iterate;
  }

  void Optimization::solve(AbstractFunctional& fN, AbstractFunctional& fT, AbstractFunctionSpaceElement& x, AbstractHierarchicalErrorEstimator& hbErrorEstimator_)
  {
    hbErrorEstimator = &hbErrorEstimator_;
    solve(fN,fT,x);
  }


  void Optimization::solve(AbstractFunctional& f, AbstractFunctionSpaceElement& x, AbstractHierarchicalErrorEstimator& hbErrorEstimator_)
  {
    hbErrorEstimator = &hbErrorEstimator_;
    solve(f,nullptr,x);
  }

  bool Optimization::regularityTest(double nu, std::vector<double> const& tau, bool reliableQuadraticModel) const
  {
    bool passed = nu > p.minimalDampingFactor;
    if(reliableQuadraticModel) passed = passed && tau[0] > p.minimalDampingFactor;
    if(verbose > 0 && !passed) std::cout << csPre << "Regularity test failed!!!" << std::endl;
    return passed;
  }

  std::pair<std::unique_ptr<AbstractFunctionSpaceElement>,std::unique_ptr<AbstractFunctionSpaceElement> > Optimization::computeNormalStep(AbstractFunctionSpaceElement* normalStepResidual, AbstractFunctionSpaceElement* adjointResidual)
  {
    std::pair<std::unique_ptr<AbstractFunctionSpaceElement>,std::unique_ptr<AbstractFunctionSpaceElement> > result(iterate->clone(),iterate->clone());
    *(result.first) *= 0;
    *(result.second) *= 0;

    if(normalDirection) // constrained problem
    {
      // Compute normal step, compute update for Lagrange multiplier
      if( step > 1 && std::fabs(normOfLastCorrection_Undamped-normOfLastCorrection) < p.getEps() * normOfLastCorrection_Undamped) normalDirection->setRelativeAccuracy(std::max(std::min(p.minimalAccuracy, C.omega*normOfLastCorrection),p.desiredAccuracy));
      else normalDirection->setRelativeAccuracy(p.minimalAccuracy);
      boost::timer::cpu_timer timer;
      normalDirection->ordinaryAndAdjoint(*(result.first), *(result.second), *normalLinearization, normalStepResidual, adjointResidual);
      normalStepComputationTime += (double)timer.elapsed().wall*1e-9;

      if( verbose > 1 ) std::cout << "normal step computation time: " << boost::timer::format(timer.elapsed()) << std::endl;

      iterate->axpy(1.0,*(result.second),"dual");
    }
    else // unconstrained problem
      if(verbose > 0) std::cout << csPre << "Unconstrained Problem" << std::endl;
    
    return result;
  }

  std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > Optimization::computeTangentialStep(LagrangeLinearization& lagrangeLinearization, AbstractFunctionSpaceElement const& normalStep, double nu, std::vector<double> const& tau)
  {
    std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > tangentialBasis;
    tangentialBasis.resize(1,std::shared_ptr<AbstractFunctionSpaceElement>(iterate->clone())); *(tangentialBasis[0]) *= 0;

    if(!tangentSpace) return tangentialBasis;

    double relativeAccuracy = std::max(p.desiredAccuracy,std::min(p.minimalAccuracy, L.omega * std::fabs(normOfLastCorrection)));
    if(step > 1 && std::fabs(normOfLastCorrection_Undamped-normOfLastCorrection) < p.getEps() * normOfLastCorrection_Undamped && noDamping(nu))
    {
      if(tangentSpace) tangentSpace->setRelativeAccuracy(relativeAccuracy);
      if(verbose > 1)
      {
        std::cout << csPre << "desired contraction: " << p.getThetaAim() << std::endl;
      }
    }
    else tangentSpace->setRelativeAccuracy(p.minimalAccuracy);// * std::min(1.,normOfLastCorrection_Undamped));

    boost::timer::cpu_timer timer;
    tangentSpace->setLipschitzConstant(L.omega);
    tangentSpace->basis(tangentialBasis, lagrangeLinearization, normalStep, nu, nullptr);
    tangentialStepComputationTime += (double) timer.elapsed().wall * 1e-9;

    if( verbose > 1 ) std::cout << "tangential step computation time: " << boost::timer::format(timer.elapsed()) << std::endl;

    return tangentialBasis;
  }

  double Optimization::updateNormalStepDampingFactor(double normNormal) const
  {
    double nu = 1;
    if(normNormal > p.getEps() && std::fabs(normNormal*C.omega) > p.getEps()) nu = std::min(1.0,/*2.0*/p.getThetaNormal()/(C.omega*normNormal));
    if(noDamping(nu)) nu = 1;
    return nu;
  }

  void Optimization::updateTangentialDampingFactor(double nu, double normNormal, double normTangential, RegularizedQuadraticFunction const& cubic, std::vector<double> &tau) const
  {
    if( tau.size()==1 )
    {
      tau[0] = 1.0;

      if( normTangential < p.getSqrtEps() ) return;

      double maxTauMax = 1;
      double taumax = maxTauMax;
      if( C.omega != 0.0  )
      {
        if( ( pow(p.getThetaAim()/C.omega,2)-pow(nu*normNormal,2) ) > 0 && normTangential > 0) taumax = sqrt(pow(p.getThetaAim()/C.omega,2)/*-pow(nu*normNormal,2)*/)/normTangential;//std::min(1.,normTangential);

        if( taumax > maxTauMax && tangentSpace->localConvergenceLikely() )
        {
          if( verbose > 1 ) std::cout << csPre << "reducing taumax: " << taumax << " to " << maxTauMax << std::endl;
          taumax = maxTauMax;
        }
      }
      if( verbose > 1 ) std::cout << csPre << "taumax: " << taumax << std::endl;

      if( taumax > 0.0 )
      {
        CubicModel1dForFmin cubicmodel(cubic);
        double toleranceFmin = p.getThirdSqrtEps()*std::min(taumax,1.);
        double taumin = 0;
        tau[0] = Fmin(taumin,taumax,cubicmodel,toleranceFmin);
      }

      if( noDamping(tau) ) tau[0]=1.0;
    }
    else
    {
      if( tau.size() > 1 )
      {
        std::cout << csPre << "Multidimensional search space not implemented" << std::endl;
        exit(-1);
      }
    }

  }

  std::unique_ptr<AbstractFunctionSpaceElement> Optimization::computeSimplifiedNormalStep(AbstractLinearization const& lin_xplus, AbstractFunctionSpaceElement* normalStepResidual)
  {
    std::unique_ptr<AbstractFunctionSpaceElement> simplifiedNormalStep(iterate->initZeroVector());
    if(!normalDirection) return simplifiedNormalStep;

    normalDirection->simplified(*simplifiedNormalStep,lin_xplus,nullptr);
    return simplifiedNormalStep;
  }

  void Optimization::updateConstraintD1LipschitzConstant(double normSimplifiedNormal, double normCAtCorr)
  {
    if( verbose > 1 ) std::cout << csPre << "normCAtCorr = " << normCAtCorr << ", normSimplifiedNormal = " << normSimplifiedNormal << std::endl;
//    if( normCAtCorr > p.getEps() )
    {
      C.theta = normSimplifiedNormal/normCAtCorr; // contraction rate 
      bool lock = (C.theta < 0.25 && (normCAtCorr < p.getSqrtEps()*normOfIterate || normSimplifiedNormal < p.getEps()*normOfIterate) );
      C.update(2*C.theta/normCAtCorr, lock);
    }
  }

  double Optimization::updateLagrangianD2LipschitzConstant(AbstractLinearization const& lin_x0, AbstractFunctionSpaceElement const& secondOrderCorrected, double normLAtCorr, double quadraticModelAtCorr, RegularizedQuadraticFunction const& cubic, double nu, std::vector<double> const& tau,
                                                           AbstractFunctionSpaceElement const& normalStep, AbstractFunctionSpaceElement const& trialIterate, AbstractFunctionSpaceElement const& sNormalStep, AbstractFunctionSpaceElement const& correction, double Lxdn_res)
  {
    double eta = 1;
    if( tangentSpace == nullptr ) return eta;

    std::unique_ptr<AbstractLinearization> lin_xplus(functionalN->getLinearization(trialIterate));
    auto f_x0 = lin_x0.eval();
    // compute: f(x_)-q(x)dx or alternative if high round-off is expected
    double deltaf(0.0);
    auto zero = tau;
    for(double& d : zero) d = 0;

//    if(normalDirection)
//    {
//      std::unique_ptr<AbstractFunctionSpaceElement> dc_x0_dn ( iterate->initZeroVector() ), lxxdn( iterate->initZeroVector() );
//      lin_x0.d2axpy(1.0, *lxxdn, normalStep, 0, 2, 0, 2);
//      lin_x0.d2axpy(1.0, *dc_x0_dn, normalStep, 2, 3, 0, 2);
//      double p_dc_x0_dn = iterate->applyAsDualTo(*dc_x0_dn,"dual");
//      double lxxdndn = normalStep.applyAsDualTo(*lxxdn);
      std::unique_ptr<AbstractLinearization> lin_xbar=functionalN->getLinearization(secondOrderCorrected);
      deltaf=lin_xbar->eval()-f_x0;
      double deltaf_cor = deltaf - cubic.d0(zero);
      if(std::fabs(cubic.d0(tau) - cubic.d0(zero)) > p.getSqrtEps()*normOfIterate/* && std::fabs(deltaf_cor) > p.getSqrtEps()*/) eta = deltaf_cor/( cubic.d0(tau) - cubic.d0(zero) ) ;
      else eta = 1;
      if( verbose > 1 )
      {
        std::cout << csPre << "predicted decrease: " << (cubic.d0(tau) - cubic.d0(zero)) << std::endl;
        std::cout << csPre << "actual decrease: " << deltaf_cor << std::endl;
        std::cout << csPre << " eta = " << eta << std::endl;
      }
      L.setFirstOrder( normLAtCorr, C.theta, deltaf - quadraticModelAtCorr );


//    }
//    else { assert("not implemented"); }

    if(L.highRoundOffError(true,p.getEps()*f_x0))
    {
      std::unique_ptr<AbstractFunctionSpaceElement> tmp1=iterate->clone();
      // tmp1 = chart_x(0)-chart_x(dx+ds)
      *tmp1 -= secondOrderCorrected;

      lin_xbar->evald(*tmp1);
      double p0_c_xbar=tmp1->applyAsDualTo(*iterate,"dual");  // p0 is dual component of *iterate

      double p0_c_x0(0.0);
      if(!noDamping(nu))
      {
        lin_x0.evald(*tmp1);
        p0_c_x0=tmp1->applyAsDualTo(*iterate,"dual");     // p0 is dual component of *iterate
      }

      // Lx_xplus_ds = L_x(x_+,p_0) ds
      double Lx_xplus_ds(0.0);
      lin_xplus->evald(*tmp1);
      Lx_xplus_ds=tmp1->applyAsDualTo(sNormalStep,"primal"); // delta s is primal component of *snormalStep

      // secondorderestimate [ (L_x(x_+,p0)-L_x(x0,p0)-L_xx(x0,p0)dx)dx ] + [ L_x(x_+,p0)ds ]-[ p0 c(xbar)-(1-nu)p0 c(x0) ]
      // each term in [...] is third order.
      // the summands in the first (...) are first order (close to solution), the difference is second order -> first order round-off error effects
      // close to a solution where nu = 1, the very last term vanishes -> then no other round-off error effects
      // if nu is very small, and c(x_0) approx c(xbar), then round-off error is expected (try alternative?)

      // tmp1 = L_x(x0,p0)+L_xx(x0,p0)dx
      lin_x0.evald(*tmp1);
      lin_x0.ddxpy(*tmp1,correction);

      // tmp2=L_x(x_+,p0)-L_x(x0,p0)-L_xx(x0,p0)dx
      std::unique_ptr<AbstractFunctionSpaceElement> tmp2=iterate->clone();
      lin_xplus->evald(*tmp2);
      (*tmp2)-=(*tmp1);

      L.setSecondOrder(0.5*tmp2->applyAsDualTo(correction,"primal")+0.5*Lx_xplus_ds-(p0_c_xbar-(1-nu)*(p0_c_x0-Lxdn_res)));
    }

    L.update( std::fabs( eta - 1. ) < (1.-p.etaLock) );

    if(verbose > 1)
    {
      std::cout << csPre << "f(x_0): " << f_x0 << " f(x): "  << f_x0+deltaf << "  |  Decrease: " << -deltaf << std::endl;
      std::cout << csPre << "Quadratic model q(x):" << f_x0+quadraticModelAtCorr << std::endl;
      std::cout << csPre << "q(x)-f(x_0): " << quadraticModelAtCorr << std::endl;
      std::cout << csPre << "q(x)-f(x): " << quadraticModelAtCorr-deltaf << std::endl;
    }

    return eta;
  }


  AcceptanceTest Optimization::acceptanceTest(double eta, double nu, std::vector<double> const& tau, double normOfCorrection) const
  {
    if( nu > 0.1 && normOfCorrection < p.getEps() ) return AcceptanceTest::Passed;
    if( noDamping(nu) && noDamping(tau) && (normOfCorrection) < p.getSqrtEps() ) return AcceptanceTest::Passed;
    
    if( eta < p.etaMin )
    {
      if( verbose > 0 ) std::cout << csPre << "Rejecting due to small eta: " << eta << " < " << p.etaMin << std::endl;
      return AcceptanceTest::TangentialStepFailed;
    }
    
    if(C.theta > p.getThetaMax())
    {
      if(verbose > 0) std::cout << csPre << "Insufficient contraction:" << C.theta << " > " << p.getThetaMax() << std::endl;
      return AcceptanceTest::NormalStepFailed;
    }

    return AcceptanceTest::Passed;
  }

  bool Optimization::noDamping(double d) const
  {
    return std::fabs(d-1.) < dampingFactorTolerance;
  }

  bool Optimization::noDamping(std::vector<double> const& tau) const
  {
    for(double d : tau) if(!noDamping(d)) return false;
    return true;
  }


  void Optimization::terminationMessage(int flag)
  {
    switch(flag)
    {
      case 1  : std::cout << csPre << "Desired accuracy reached!" << std::endl; break;
      case -1 : std::cout << csPre << "Regularity test failed!" << std::endl;; break;
      case -2 : std::cout << csPre << "Maximum number of iterations reached!" << std::endl; break;
      default : Algorithm::terminationMessage(flag);
    }
  }

  void Optimization::printNormalStep(double normNormal, double nu) const
  {
    std::cout << csPre << "Normal step: " << normNormal << std::endl;
    std::cout << csPre << "Damping factor: " << nu << std::endl << std::endl;
  }

  void Optimization::printTangentialStep(double normTangential, double tau) const
  {
    std::cout << csPre << "Tangential step: " << normTangential << std::endl;
    std::cout << csPre << "Damping factor: " << tau << std::endl << std::endl;

  }

  std::unique_ptr<AbstractFunctionSpaceElement> Optimization::createCorrection(double nu, const AbstractFunctionSpaceElement &normalStep, const std::vector<double> &tau, const std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > &tangentialBasis) const
  {
    std::unique_ptr<AbstractFunctionSpaceElement> correction(normalStep.clone());
    *correction *= 0.0;
    correction->axpy(nu,normalStep,"primal");
    for(size_t i=0; i<tangentialBasis.size();++i) correction->axpy(tau[i],*tangentialBasis[i],"primal");
    return correction;
  }

  void Optimization::addCorrection(AbstractFunctionSpaceElement &v, double nu, const AbstractFunctionSpaceElement &normalStep, const std::vector<double> &tau, const std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > &tangentialBasis) const
  {
    v.axpy(nu,normalStep,"primal");
    for(size_t i=0; i<tangentialBasis.size(); ++i) v.axpy(tau[i],*tangentialBasis[i],"primal");
  }
}  // namespace Kaskade
