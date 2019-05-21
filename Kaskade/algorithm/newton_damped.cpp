#include "newton_damped.hh"
#include <iostream>
#include <cmath>


namespace Kaskade
{

//------------------------------------------------------------------------------------------------------

void DampedCovariantNewtonMethod::initialize()
{
  p.reset();
  NewtonsMethod::initialize();
  startIterate= iterate->clone();  
  auxvector = iterate->clone();  
  scorrection = iterate->clone();  
  p.Theta = 0.0;
  p.dampingFactor = p.initDampingFactor;
//  linearSolver.setAbsoluteAccuracy(p.desiredAccuracy/2);
}

void DampedCovariantNewtonMethod::initNewtonStep()
{
  double accLinSlv=0.0;
  if(p.totalCorrection.isValid())
    accLinSlv = p.relativeAccuracy*p.totalCorrection;
//  accLinSlv = std::max(accLinSlv,p.desiredAccuracy);
  linearSolver.setAbsoluteAccuracy(accLinSlv);
}


void DampedCovariantNewtonMethod::predictNextDampingFactor(AbstractFunctionSpaceElement& correction)
{
  p.lastrejected.rejected=false;
  if(p.normCorr.isValid()) 
  {
    p.absoluteAccuracyLast=p.normCorr*linearSolver.getRelativeAccuracy();
    p.normCorrLast=p.normCorr;
    p.lengthOfStep=norm(*auxvector);
  }
  p.normCorr=norm(correction);
  if(p.normSCorr.isValid() && p.normCorrLast.isValid() && p.absoluteAccuracyLast.isValid())
  {
    *scorrection -= correction;
    p.omega=norm(*scorrection)/(p.lengthOfStep*(p.normSCorr+p.absoluteAccuracyLast));
    double lambdas = 1.0/(p.omega*p.normCorr);
    p.dampingFactor = std::min(1.0,lambdas);
    if(report) std::cout << "Damping factor predicted:" << p.dampingFactor << std::endl;
  } 
  else 
  {
    p.dampingFactor = p.initDampingFactor;
    if(report) std::cout << "Initial Damping factor:" << p.dampingFactor << std::endl;
  }
}

RegularityTest DampedCovariantNewtonMethod::regularityTest(double scalingFactor)
{
  if(scalingFactor < p.minDampingFactor) return RegularityTest::Failed;
  return NewtonsMethod::regularityTest(scalingFactor);
}

void DampedCovariantNewtonMethod::computeTrialIterate
(AbstractFunctionSpaceElement& trialIterate, AbstractFunctionSpaceElement const& direction, AbstractLinearization const& lin)
  {
    *auxvector = direction;
    *auxvector *= p.dampingFactor;
    chart.addPerturbation(trialIterate, *auxvector, lin);
  }

AcceptanceTest DampedCovariantNewtonMethod:: evaluateTrialIterate(AbstractFunctionSpaceElement const& trialIterate,
                                                                                 AbstractFunctionSpaceElement const& correction, 
                                                                                 AbstractLinearization const& lin)
{ 
  if(!functional->inDomain(trialIterate)) 
  {
    p.dampingFactor *= p.reduceOnOutsideDomain;
    if(report) std:: cout << "Iterate outside domain of definition. New Damping factor:" 
                          << p.dampingFactor << std::endl;
    p.dampingForDomain=p.dampingFactor;
    return AcceptanceTest::Failed;
  }
  simplifiedLinearization=functional->getLinearization(trialIterate);
  linearSolver.simplified(*scorrection, *simplifiedLinearization);

  p.normSCorr=norm(*scorrection);

  *auxvector=*scorrection;

  auxvector->axpy(p.dampingFactor-1.0,correction);

  double normDSCorr = norm(*auxvector);

  double normDampedCorr=p.normCorr*p.dampingFactor;

  p.Theta = normDSCorr/ std::max(normDampedCorr,1e-16);

  p.SCorrByCorr = p.normSCorr/p.normCorr;

  p.omega=2.0*p.Theta/normDampedCorr;

  if(p.dampingFactor==1.0 && !p.omega0.isValid()) p.omega0=p.omega;

  if(!p.omega0.isValid()) 
    p.omega0=p.omega;

  
  if(report) std::cout << "|Dx|:" << p.normCorr << " Theta:" << p.Theta << " |Dx_|/|Dx|:" << p.SCorrByCorr << " Omega:" << p.omega << std::endl;

  double lambdas=1.0/(p.normCorr*p.omega);
  
  if(p.SCorrByCorr >= 1.0)  
  { // New damping factor had to be reduced, test if this was sufficient

    p.lastrejected.rejected=true;
    p.lastrejected.dampingFactor=p.dampingFactor;
    p.lastrejected.omega=p.omega;
    p.lastrejected.Theta=p.Theta;

    p.dampingFactor=std::min(lambdas,0.5*p.dampingFactor);  // A least halve damping factor
    if(report) std::cout << "Damping factor reduced:" <<  p.dampingFactor << std::endl;
    return AcceptanceTest::Failed;
  }

  double dampingFactorNew = std::min(1.0,lambdas);
  if(dampingFactorNew >= p.dampingFactor && p.lastrejected.rejected==true)
  {
// If new damping factor is increased, do not increas it beyond the last rejected daming factor
// rather compute the optimal damping factor under the assumption that omega interpolates linearly
// the omegas measured at the current damping factor and at the last rejected damping factor

//  omega(lambda) = a+lambda b so that omega(lambda_k)=p.omega and omega(lambda_rej)= p.lastrejected.omega
    double dlambda=p.lastrejected.dampingFactor-p.dampingFactor;
    assert(dlambda > 0);
    std::cout << "Omegas:" << p.omega << " " << p.lastrejected.omega << std::endl;
    double b=std::max((p.lastrejected.omega-p.omega)/dlambda,0.0);
    double a=p.omega-b*p.dampingFactor;
//    assert(b > 0);
// 0 = -2/|Dx|+2a lambda+3b lambda^2
    dampingFactorNew=std::min(1.0,2.0/(p.normCorr*(a+std::sqrt(a*a+6*b/p.normCorr))));
    
  }


  if(dampingFactorNew >= 4.0*p.dampingFactor || (dampingFactorNew==1.0 && p.dampingFactor < 1.0))  
  { // New damping factor has increased dramatically, better test again!
    p.dampingFactor = dampingFactorNew;
    if(report) std::cout << "Damping factor increased:" <<  p.dampingFactor << std::endl;
    return AcceptanceTest::Failed;
  }  

  // damping factor has successfully passed the test
  if(report) std::cout << "Damping Factor accepted:" <<  p.dampingFactor << std::endl;
  return AcceptanceTest::Passed;
}

void DampedCovariantNewtonMethod::updateIterate(AbstractFunctionSpaceElement& iterate, 
                                                AbstractFunctionSpaceElement& trialIterate,
                                                AbstractLinearization const& lin)
{
  if(p.Theta <= 0.25 && dampingFactor() == 1.0)
  {
    *auxvector=*scorrection;
    auxvector->axpy(dampingFactor(),*correction);

    chart.addPerturbation(trialIterate,*auxvector,lin);

    std::cout << "Simplified correction added!" << std::endl;
  }
  *auxvector = iterate;
  *auxvector -= trialIterate;
  iterate.swap(trialIterate);
  p.dampingForDomain=1.0;
}

Convergence DampedCovariantNewtonMethod::convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate)
{
  *trialIterate = iterate;
  *trialIterate -= *startIterate;
  p.totalCorrection = norm(*trialIterate);
  if(report) std::cout << "|dx|: " << p.normCorr << " Required:" << p.desiredAccuracy << std::endl; 
  if(p.Theta < 1e-3 && dampingFactor() == 1.0 && p.normSCorr+p.normCorr*linearSolver.getRelativeAccuracy() <= p.desiredAccuracy) return Convergence::Achieved;
  if(p.normCorr <= p.desiredAccuracy && dampingFactor() == 1.0 && p.Theta <= 0.5) return Convergence::Achieved;
  return Convergence::Missed;
}


////-----------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------


void StateConstraintsNewtonMethod::initialize()
{
  NewtonsMethod::initialize();
  startIterate= iterate->clone();  
  scorrection = iterate->clone();  
  p.Theta = 0.0;
  p.dampingFactor = p.initDampingFactor;
  linearSolver.setAbsoluteAccuracy(0.0);
}

void StateConstraintsNewtonMethod::initNewtonStep()
{
  failedbydamping = false;
  double accLinSlv=0.0;
  if(p.totalCorrection.isValid())
    accLinSlv = p.relativeAccuracy*p.totalCorrection;
  accLinSlv = std::max(accLinSlv,p.desiredAccuracy);
  linearSolver.setAbsoluteAccuracy(accLinSlv);
}


RegularityTest StateConstraintsNewtonMethod::regularityTest(double scalingFactor)
{
  if(p.Theta > p.ThetaMaxAllowed || failedbydamping) return  RegularityTest::Failed;
  return NewtonsMethod::regularityTest(scalingFactor);
}

AcceptanceTest StateConstraintsNewtonMethod::evaluateTrialIterate(
  AbstractFunctionSpaceElement const& trialIterate, 
  AbstractFunctionSpaceElement const& correction, 
  AbstractLinearization const& lin) 
{
// Norm of pointwise damped Correction
  *scorrection=lin.getOrigin();
  *scorrection-=trialIterate;
  p.normModCorr=norm(*scorrection);
  if(report) std::cout << "|Dx|:" << p.normCorr<< " |Dx_C|:" << p.normModCorr << std::endl;

  std::cout << "Damping: " << (p.normModCorr/p.normCorr) << std::endl;

// Norm of simplified Correction
  simplifiedLinearization=functional->getLinearization(trialIterate);
  linearSolver.simplified(*scorrection, *simplifiedLinearization);

  p.normSCorr=norm(*scorrection);
  
  failedbydamping=(p.normCorr/p.normModCorr < 0.05);

  p.Theta = p.normSCorr/ p.normCorr;
  if(!p.Theta0.isValid())   p.Theta0 = p.normSCorr/ p.normCorr;
  if(!p.omega0.isValid())  p.omega0= p.Theta0/p.normModCorr;
  
  if(report) std::cout << "Contraction: " << p.Theta << std::endl;
  return AcceptanceTest::Passed;
}

void StateConstraintsNewtonMethod::updateIterate(AbstractFunctionSpaceElement& iterate, 
                                                   AbstractFunctionSpaceElement& trialIterate,
                                                   AbstractLinearization const& lin)
{
  nsc=p.normModCorr*p.Theta/(1-p.Theta);
  if(p.Theta > 1) nsc = 1e300;
  std::cout << "|x_h-x_*h|: " << nsc << " " << std::flush;
  if(p.Theta <= 0.25)
  {
    if(report) std::cout << " + dx_: " << p.normModSCorr << " " << std::flush;
    scorrection->axpy(dampingFactor(),*correction);

    chart.addPerturbation(iterate,*scorrection,lin);
  }
  else
    iterate.swap(trialIterate);
  if(linearSolver.improvementPossible()) nsc += linearSolver.getAbsoluteAccuracy();
  std::cout << std::endl;
}


Convergence StateConstraintsNewtonMethod::convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate)
{
  p.accuracyReached = nsc;
  *scorrection = iterate;
  *scorrection -= *startIterate;
  p.totalCorrection = norm(*scorrection);
  std::cout << "|x_h-x_*|:" << nsc << " Required:" << std::max(p.relativeAccuracy*p.totalCorrection,p.desiredAccuracy) <<  std::endl;
  if(p.normModCorr/p.normCorr  < 0.5 || p.Theta > 0.5) return Convergence::Missed;
  if(nsc < p.relativeAccuracy*p.totalCorrection || nsc < p.desiredAccuracy) return Convergence::Achieved;
  return Convergence::Missed;
}

}  // namespace Kaskade
