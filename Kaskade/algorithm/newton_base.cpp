#include "newton_base.hh"
#include <iostream>
#include <limits>
#include <cmath>

namespace Kaskade
{

void NewtonsMethod::solve(AbstractFunctional* f,AbstractFunctionSpaceElement& x)
{
  doOneStep=false;
  functional=f;
  iterate = x.clone();
  p.termination=algorithmWrapper();
  if(p.termination > 0 || p.termination == -2) x=*iterate;
}

void NewtonsMethod::oneStep(AbstractFunctional* f, AbstractLinearization* l,AbstractFunctionSpaceElement& x)
{
  doOneStep=true;
  functional=f;
  newtonPtr=l;
  iterate = x.clone();
  p.termination=oneStepWrapper();
  x=*iterate;
}

void NewtonsMethod::updateIterate(AbstractFunctionSpaceElement& iterate, 
                                  AbstractFunctionSpaceElement& trialIterate,
                                  AbstractLinearization const& lin)
{
  iterate=trialIterate;
}

void NewtonsMethod::terminationMessage(int flag)
{
  std::cout << "Newtons Method: "; 
  switch((int)flag)
  {
  case 1  : std::cout << "Desired accuracy reached!" << std::endl; break;
  case -1 : std::cout << "Regularity test failed!" << std::endl; break;
  case -2 : std::cout << "Maximum number of iterations reached!" << std::endl; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

void NewtonsMethod::initialize()
{
  trialIterate=iterate->clone();
  correction=iterate->clone();
}

void NewtonsMethod::finalize(int flag)
{
}


int NewtonsMethod::runAlgorithm() 
{
  // Main loop
  gridHasChanged = false;

  for(step=1; step <= maxSteps(); step++)
  { 
    if(!doOneStep)
    {
      newtonLinearization=functional->getLinearization(*iterate);
      if(report) std::cout <<  "\n Newton-Step: " << step << std::endl;
      newtonPtr=newtonLinearization.get();
    }
//    sp.setLinearization(*newtonPtr,linearSolver);
    initNewtonStep();
    getSearchDirection(*correction);

    gridHasChanged |= linearSolver.changedGrid();

    norm.setOrigin(*newtonPtr);
    predictNextDampingFactor(*correction);
 
    // Damping Loop
    CycleDetection<double> cd;
    do
    {
      if(cd.detectCycling(dampingFactor())) break;
      computeTrialIterate(*trialIterate,*correction,*newtonLinearization);
      //sp.getTrialIterate(*trialIterate, *correction , *iterate, dampingFactor());
      if(regularityTest(dampingFactor())==RegularityTest::Failed) return -1;
    }
    while(evaluateTrialIterate(*trialIterate, *correction ,  *newtonPtr)==AcceptanceTest::Failed);
    if(regularityTest(dampingFactor())==RegularityTest::Failed) return -1;

    updateIterate(*iterate,*trialIterate, *newtonPtr);
    if(doOneStep) return 1;
    if(!p.reuseFactorization) newtonLinearization.reset();
    logQuantities();
    if(convergenceTest(*correction, *iterate)==Convergence::Achieved) return 1;
  }
  return -2;
}

}  // namespace Kaskade
