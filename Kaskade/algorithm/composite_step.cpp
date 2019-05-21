#include "composite_step.hh"
#include <cmath>
#include <fstream>

namespace Kaskade
{

void TangentialStepArmijo::getSearchDirection(AbstractFunctionSpaceElement& scorrection, AbstractFunctionSpaceElement const& correction)
  {  
    assert(linPtr);

//    *iterate += correction;
    performArmijoLoop();
    scorrection=*trialIterate;
    scorrection -= linPtr->getOrigin();
  }


int TangentialStepArmijo::performArmijoLoop()
{
  dampingFactor = 1.0;
  armijopar = 0.001;
  sp.setLinearization(*linPtr,linSolver);
  sp.getSearchDirection(*correction);
  functionalAtIterate=linPtr->eval();
  do
  {
    sp.getTrialIterate(*trialIterate, *correction , *iterate, dampingFactor);
    if(std::fabs(dampingFactor)<1e-9) return -1;
  }
  while(evaluateTrialIterate(*trialIterate, *correction, *linPtr)==false);
  return 0;      
}

bool TangentialStepArmijo::evaluateTrialIterate(
  AbstractFunctionSpaceElement const& trialIterate, 
  AbstractFunctionSpaceElement const& correction, 
  AbstractLinearization const& lin)
{ 
  simplifiedLinearization=fPtr->getLinearization(trialIterate);
  double fx = functionalAtIterate;
  double fxt = simplifiedLinearization->eval();
  scorrection = iterate->clone();
  *scorrection -= trialIterate;
  std::cout << "|x-tx|=" << norm(*scorrection) << std::endl;
  std::vector<double> scaling(3,1.0), scorig;
  lin.getScalePars(scorig);
  lin.scaleRHS(scaling);
  lin.evald(*gradient);
  lin.scaleRHS(scorig);
  double dfxdx = -gradient->applyAsDualTo(correction);
  std::cout << "f(x): " << fx << " f(x)-f(x+tx): " << fx-fxt << " <df,dx> " << dfxdx << " " << (fx-fxt)/dampingFactor/dfxdx << std::endl;
  if(fxt <= std::min(0.0,fx-dampingFactor*armijopar*dfxdx)) return true;
  else dampingFactor *= 0.5;
  if(dfxdx <= 0) 
  {
    std::cout << "Warning: maybe ascent direction!" << std::endl;
    if(dampingFactor > 0) 
    {
      dampingFactor *= -1.0;
      return false;
    }
  }
  std::cout << "DampingFactor: " << dampingFactor << std::endl;
  return false;
}


void CompositeStep::terminationMessage(int flag)
{
  std::cout << "Composite Step: "; 
  switch((int)flag)
  {
  case 1  : std::cout << "Desired accuracy reached!" << std::endl; break;
  case -1 : std::cout << "Regularity test failed!" << std::endl;; break;
  case -2 : std::cout << "Maximum number of iterations reached!" << std::endl; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

  /// Return true, if convergence is detected, false otherwise
  NewtonsMethod::Convergence CompositeStep::convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate)
  {
    double nrm=norm(correction);
    std::cout << "|dx|" << nrm << std::endl;
    if(nrm < p.desiredAccuracy) return NewtonsMethod::ConvergenceAchieved;
    
    return NewtonsMethod::ConvergenceMissed;
  }


int CompositeStep::runAlgorithm() 
{
  // Main loop
  for(step=1; step <= maxSteps(); step++)
  { 
    newtonLinearization = functional->getLinearization(*iterate);
    if(report) std::cout <<  "\n Optimization-Step: " << step << std::endl;
    normalStep.setLinearization(functional, newtonLinearization.get());
    tangentialStep.setLinearization(functional, newtonLinearization.get());
    norm.setOrigin(*newtonLinearization);
    do
    {
    if(report >= 2) std::cout <<  "\n Normal-Step: " << std::endl;
      normalStep.getSearchDirection(*correction);
      normalStep.getTrialIterate(*trialIterate,*correction,*iterate,1.0);
    if(report >= 2) std::cout <<  "\n Tangential-Step: " << std::endl;
      tangentialStep.getSearchDirection(*scorrection,*correction);
      tangentialStep.getTrialIterate(*trialIterate2,*scorrection,*trialIterate,1.0);
      *correction=*iterate;
      *correction -= *trialIterate2;
    }
    while(evaluateTrialIterate(*trialIterate2, *correction ,  *newtonLinearization)==NewtonsMethod::AcceptanceTestFailed);

    updateIterate(*iterate,*trialIterate2, *newtonLinearization);
    newtonLinearization.reset();
//    logQuantities();
    if(convergenceTest(*correction, *iterate)==NewtonsMethod::ConvergenceAchieved) return 1;
  }
  return -2;
}

void CompositeStep::solve(AbstractFunctional* f,AbstractFunctionSpaceElement& x)
{
  functional=f;
  iterate = x.clone();
  trialIterate=x.clone(); 
  trialIterate2=x.clone(); 
  correction=x.clone();
  scorrection=x.clone(); 
  p.termination=algorithmWrapper();
  if(p.termination > 0 || p.termination == -2) x=*iterate;
}

} // namespace Kaskade