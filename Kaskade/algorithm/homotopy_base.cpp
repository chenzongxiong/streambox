#include "homotopy_base.hh"
#include "newton_base.hh"
//#include "newton_bridge.hh"


namespace Kaskade
{

void HomotopyBase::solve(AbstractParameterFunctional* f, AbstractFunctionSpaceElement& x)
{
  functional=f;
  trialIterate = x.clone();
  p.termination=algorithmWrapper();
  finalize();
  x= *trialIterate;

}

void HomotopyBase::logQuantities()
{
  p.logStep();
}

void HomotopyBase::terminationMessage(int flag)
{
  std::cout << "Interior Point Method: "; 
  switch((int)flag)
  {
  case  1: std::cout << " Desired accuracy reached!" << std::endl; break;
  case -1: std::cout << " Did not converge in \"maxSteps\" steps!"; break;
  case -2: std::cout << " Stepsize too small!"; break;
  case -3: std::cout << " Too many trial steps!"; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

int HomotopyBase::convergenceTest()  
{
  if(p.stepSuccessful)
  {
    p.mu=p.muTrial;
  }
  std::cout << "||x_mu-x_opt||=" <<lengthOfPath() << std::endl;
  if(p.muTrial <= muFinal() || (lengthOfPath() <= p.accfinal && p.stepSuccessful))
    return 1;       //Method converged
 
   if(step >=2 && 1.0-p.sigma < (1.0-p.reductionFactorWorst)/2.0)
     return -2;       // Stepsize too small
    
  return 0;    // No special event
}

void HomotopyBase::computeGapParameter()
{
  // Step size strategy in case of failure of corrector
  if(!p.stepSuccessful)
  {
    if(step==1) // Bad case: initial mu was too small
    {
      p.mu *= p.firstStepFailureFactor;
      p.muTrial = p.mu;
    }
    else
      p.muTrial=muOnFailure();
  }
  else
  {
    switch(step) {
    case 1  : p.muTrial=p.mu; break;
    default : p.muTrial=muOnSuccess(step);
    }
  }
  p.muTrial= std::max(p.muTrial.value(), 0.99*muFinal());
  p.muTrial= std::max(p.muTrial.value(), p.reductionFactorBest*p.mu);
  p.sigma=p.muTrial/p.mu;
  p.deltaMu = p.mu-p.muTrial;
  std::cout << "sigma:" << p.sigma << " ";
}


  
bool HomotopyBase::stepSuccessful()
{
  // Corrector was successful:
  if(p.correctorTermination==1 || p.correctorTermination==-2)
  {
    p.stepSuccessful=true;
    updateModelOfHomotopy(step);
    updateIterate();
    return true;
  } 
  // Corrector failed:
  else
  {
    p.stepSuccessful=false;
    if(step > 1) recoverIterate();
    return false;
  }
}

void HomotopyBase::computePredictor(int step)
{
}

  
void HomotopyBase::computeCorrector()
{
  corrector.resetParameters();
  initializeCorrector();
  std::unique_ptr<AbstractFunctional> fuPtr(functional->getFunctional(makePars(p.muTrial.value())));
  std::cout << "mu: " << p.muTrial;
  corrector.reportOnIteration(true);
  corrector.solve(fuPtr.get(), *trialIterate);
  p.correctorTermination=corrector.getParameters().termination;
  finalizeCorrector();
}


int HomotopyBase::runAlgorithm() 
{
  // Main loop
  int convFlag=0;
  for(step=1;step<=p.maxSteps; step++)
  {
    if(report) std::cout << "------- Interior Point Step: " << step << " -------: ";
    int trialstep=0;
    do{
      trialstep++;
      if(trialstep > this->p.maxTrialSteps) 
        return -3;         // Too many trial steps

      computeGapParameter();
      computePredictor(step);
      computeCorrector();


    } while(!stepSuccessful());
    logQuantities();

    if((convFlag=convergenceTest()))
    {
      finalizeHomotopy();
      return convFlag;
    }
  }
  return -1;  // Not converged after MaxNoOfSteps
}

double InteriorPointSqrtModel::muOnFailure()
{ p.sigma=0.5+p.sigma/2.0; 
  if(p.sigma >= 0.99) p.sigma = 1.0;
  return p.mu*p.sigma; 
};


double InteriorPointSqrtModel::muOnSuccess(int step) {
  double sigma;
  if(step==2) 
    sigma=p.reductionStart;
  else {
    double b = pp.omega*pp.accuracyCorrector;
    double c = pp.eta*pp.omega*p.mu;
    double a = c+p.desiredContraction;
    sigma = (b+sqrt(b*b+4*a*c))/(2*a);
    sigma = sigma*sigma;
    sigma = std::min(sigma,std::max(2*p.sigma,std::sqrt(p.sigma)));
    sigma = std::min(sigma,0.995);
    sigma = std::max(sigma,p.reductionFactorBest);
  }
  return sigma*p.mu;
};

void InteriorPointSqrtModel::updateModelOfHomotopy(int step) 
{
  if(step > 1)
  {
    *iterate -= *trialIterate;
    double n(norm(*iterate));
    double nP;
    if(normPlain)
      nP=(*normPlain)(*iterate);
    else
      nP=n;
    *iterate *= 0.0;
    std::cout << p.mu.value() << " " << p.muTrial.value() << std::endl;
    pp.eta=n/p.deltaMu;
    pp.slope=nP/p.deltaMu;
    GuardedCovariantNewtonParameters const& cnp(dynamic_cast<GuardedCovariantNewtonParameters const&>(corrector.getParameters()));
    pp.omega=cnp.omega0.value();
    pp.accuracyCorrector=cnp.accuracyReached;
    std::cout << "Norm:" << n << " ||x'||_2:" << nP << " Dmu:" << p.deltaMu << " ||x'||_x:"<< pp.eta << " omega:" << pp.omega << std::endl;
  }
};

void InteriorPointSqrtModel::initializeCorrector() 
{ 
  if(p.muTrial < p.muFinal) corrector.setDesiredAccuracy(p.accfinal);
  else
  {
    if(pp.eta.isValid()) corrector.setDesiredAccuracy(std::max(p.accfinal*0.25,lengthOfPath()*p.relDistanceToPath));
    else corrector.setDesiredAccuracy(sqrt(p.muTrial)*p.relDistanceToPath);
    corrector.setDesiredRelativeAccuracy(p.relDistanceToPath);//std::sqrt(p.muTrial/p.mu));
  }
}

void InteriorPointSqrtModel::updateIterate() { 
  *iterate=*trialIterate;
  p.j=functional->getFunctional(makePars(p.muTrial.value()))->getLinearization(*iterate)->eval();
  jModel.update(p.muTrial,p.j);
  jModelL.update(p.muTrial,p.j);
  jModel.fixUpdate();
  jModelL.fixUpdate();
  p.jp=jModel.getValue(p.muTrial);
  pp.jpl=jModelL.getValue(0.0);
  std::cout << "Est. Error in Functional:" << p.jp << std::endl;
  pp.timemeasured=(double)(overalltime.elapsed().user)/1e9;
}


void InteriorPointSqrtModel::printDiagnosis()
{
  std::ofstream eout("eta.log"); pp.eta.print(eout); 
  std::ofstream epout("slope.log"); pp.slope.print(epout); 
  std::ofstream nout("newtonSum.log"); pp.newtonSum.print(nout); 
  std::ofstream oout("omega.log"); pp.omega.print(oout); 
  std::ofstream sout("sigma.log"); p.sigma.print(sout); 
  std::ofstream outmu("mu.log"); p.muTrial.print(outmu); 
  std::ofstream outj("j.log"); p.j.print(outj); 
  std::ofstream outjp("jp.log"); p.jp.print(outjp); 
  std::ofstream outjpl("jplin.log"); pp.jpl.print(outjpl); 
  std::ofstream outac("ac.log"); pp.accuracyCorrector.print(outac); 
  std::ofstream outtm("time.log"); pp.timemeasured.print(outtm); 
  
  std::string filename;
  std::stringstream ss;
  ss << "iter_";
  ss << step << "_";
  ss << p.muTrial;
  ss >> filename;
  iterate->writeToFile(filename,false);
}

struct StepEquationSlope
{
  double err; // estimated iteration error of last step
  double eta; // slope; 
  double mu;  // homotopy parameter
  double Theta; // desired contraction
  double omega; // aff. cov. Lipschitz constant for Newton's method

  double operator()(double sigma) const 
  //  { return -(err+eta*mu*(1-sigma)/std::sqrt(sigma))+Theta/omega*std::sqrt(sigma); }
  { return -(err+eta*2*mu*(1-std::sqrt(sigma)))+Theta/omega*std::sqrt(sigma); }
};


double InteriorPointSlopeEstimate::muOnSuccess(int step)
{
  if(pp.eta.isValid() && pp.omega.isValid())
  {
    double sigma;
    
    StepEquationSlope f;
    f.err = pp.accuracyCorrector;
    f.eta = pp.eta;
    f.mu = p.mu;
    f.Theta = p.desiredContraction;
    f.omega = pp.omega;
    int iter;
    double sig = bisection(1e-20,1.0-1e-15,f,1e-6,iter);
    sigma = sig; 
    sigma = std::min(sigma,p.reductionFactorWorst);
    sigma = std::max(sigma,p.reductionFactorBest);
    return sigma*p.mu;
  }
  else
  {
    return p.reductionStart*p.mu;
  }

}

void InteriorPointSlopeEstimate::updateModelOfHomotopy(int step) 
{
  double h;
  if(corrector.changedGrid() || step == 1 || p.deltaMu == 0.0 || (!p.deltaMu.isValid()))
  {

    std::unique_ptr<AbstractFunctional> fuPtr(functional->getParameterLinFunctional(makePars(p.muTrial.value())));
    std::unique_ptr<AbstractLinearization> paraLinearization(fuPtr->getLinearization(*trialIterate));
    solver.doSolve(*tangent,*paraLinearization);
    h = 1.0;
  }
  else
  {
    *tangent = *trialIterate;
    *tangent -= *iterate;
    h= p.deltaMu;
  }
  double nD=norm(*tangent)/h;
  std::unique_ptr<AbstractFunctional> fuVPtr(functional->getLinFunctionValue(makePars(p.muTrial.value())));
  std::unique_ptr<AbstractLinearization> paraLinFu(fuVPtr->getLinearization(*trialIterate));
  paraLinFu->evald(*iterate);
  double slopeJ = iterate->applyAsDualTo(*tangent)/h;
  double nDP=(*normPlain)(*tangent)/h;
  pp.jpl=p.mu.value()*slopeJ;

  pp.slope=nDP;
  pp.eta=nD;

  GuardedCovariantNewtonParameters const& cnp(dynamic_cast<GuardedCovariantNewtonParameters const&>(corrector.getParameters()));
  pp.omega=cnp.omega0.value();
  pp.accuracyCorrector=cnp.accuracyReached;
  std::cout  << " ||x'||_2:" << nDP << " Dmu:" << p.deltaMu << " ||x'||_x:"<< pp.eta << " omega:" << pp.omega << std::endl;
}

void  InteriorPointSlopeEstimate::initializeCorrector() 
{ 
  if(p.muTrial < p.muFinal) corrector.setDesiredAccuracy(p.accfinal);
  else
  {
    if(pp.eta.isValid()) corrector.setDesiredAccuracy(std::max(p.accfinal*0.5,lengthOfPath()*p.relDistanceToPath));
    else corrector.setDesiredAccuracy(sqrt(p.muTrial)*p.relDistanceToPath);
    corrector.setDesiredRelativeAccuracy(p.relDistanceToPath);//std::sqrt(p.muTrial/p.mu));
  }
}

void InteriorPointSlopeEstimate::finalizeHomotopy()
{
  if(finalsolver)
  {
    std::cout << "----------------------- Final Correction Step ---------------------" << std::endl;
    finalsolver->resetParameters();
    finalsolver->setDesiredAccuracy(p.accfinal); 
    finalsolver->setDesiredRelativeAccuracy(0.0); 
    std::unique_ptr<AbstractFunctional> fuPtr(functional->getFunctional(makePars(p.muTrial.value())));
    std::cout << "muTrial: " << p.muTrial << std::endl;
    finalsolver->reportOnIteration(true);
    finalsolver->solve(fuPtr.get(), *trialIterate);
    p.correctorTermination=corrector.getParameters().termination;
    updateIterate();
    pp.timemeasured=(double)(overalltime.elapsed().user)/1e9;
    logQuantities();
  }
}
  
void InteriorPointSlopeEstimate::updateIterate()
{ 
  *iterate=*trialIterate;
  p.j=functional->getFunctional(makePars(p.muTrial.value()))->getLinearization(*iterate)->eval();
  jModelL.update(p.muTrial,p.j);
  jModelL.fixUpdate();
  p.jp=jModelL.getValue(p.muTrial);
  std::cout << "Est. Error in Functional:" << p.jp << std::endl;
  pp.timemeasured=(double)(overalltime.elapsed().user)/1e9;
}

void  InteriorPointSlopeEstimate::printDiagnosis()
{
  std::ofstream eout("eta.log"); pp.eta.print(eout); 
  std::ofstream epout("slope.log"); pp.slope.print(epout); 
  std::ofstream cpout("curvature.log"); pp.curvature.print(cpout); 
  std::ofstream nout("newtonSum.log"); pp.newtonSum.print(nout); 
  std::ofstream oout("omega.log"); pp.omega.print(oout); 
  std::ofstream sout("sigma.log"); p.sigma.print(sout); 
  std::ofstream outmu("mu.log"); p.muTrial.print(outmu); 
  std::ofstream outj("j.log"); p.j.print(outj); 
  std::ofstream outjp("jp.log"); p.jp.print(outjp); 
  std::ofstream outjpl("jplin.log"); pp.jpl.print(outjpl); 
  std::ofstream outac("ac.log"); pp.accuracyCorrector.print(outac); 
  std::ofstream outtm("time.log"); pp.timemeasured.print(outtm); 

  std::string filename;
  std::stringstream ss;
  ss << "iter_";
  ss << step << "_";
  ss << p.muTrial;
  ss >> filename;
  iterate->writeToFile(filename,false);
}

}  // namespace Kaskade
