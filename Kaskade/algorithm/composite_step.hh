#ifndef COMPOSITE_STEP_HH
#define COMPOSITE_STEP_HH

#include <memory>
#include "newton_base.hh"
#include "newton_damped.hh"

namespace Kaskade
{

class NormalStep
{
public:
  virtual ~NormalStep() {};
  virtual void setLinearization(AbstractFunctional* f, AbstractLinearization* lin) = 0;
  virtual void getTrialIterate(
    AbstractFunctionSpaceElement& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractFunctionSpaceElement const& iterate,
    double damping) = 0;

  virtual void getSearchDirection(AbstractFunctionSpaceElement& correction) = 0;
};

class TangentialStep
{
public:
  virtual ~TangentialStep() {};
  virtual void setLinearization(AbstractFunctional* f, AbstractLinearization* lin)=0;
  virtual void getTrialIterate(
    AbstractFunctionSpaceElement& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractFunctionSpaceElement const& iterate,
    double damping) =0;

  virtual void getSearchDirection(AbstractFunctionSpaceElement& scorrection, AbstractFunctionSpaceElement const& correction) =0;
};


template<class Newton>
class TangentialStepNewton : public TangentialStep
{
public:
  TangentialStepNewton(StepPolicy& s, AbstractNorm& n, typename Newton::Parameters& p,  AbstractLinearSolver& linSolver_) 
    : sp(s), norm(n), pars(p), linSolver(linSolver_)
  {}

  void setLinearization(AbstractFunctional* f, AbstractLinearization* lin)
  {
    pars.reset();
    newton.reset(new Newton(linSolver,sp,norm,pars));
    linPtr = lin;  
    fPtr = f;
  }
  virtual void getTrialIterate(
    AbstractFunctionSpaceElement& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractFunctionSpaceElement const& iterate,
    double damping)
  {
    trialIterate = iterate;
    trialIterate.axpy(damping,correction);
  }

void getSearchDirection(AbstractFunctionSpaceElement& scorrection, AbstractFunctionSpaceElement const& correction)
  {  
    assert(linPtr);
    assert(newton.get());
    scorrection=linPtr->getOrigin();
    scorrection += correction;
    newton->reportOnIteration(report);
    newton->oneStep(fPtr, linPtr, scorrection);
    scorrection -= linPtr->getOrigin();
    scorrection -= correction;
  }

private:
  std::unique_ptr<Newton> newton;
  StepPolicy& sp;
  AbstractNorm& norm;
  AbstractLinearization* linPtr;
  AbstractFunctional* fPtr;
  typename Newton::Parameters& pars;
  AbstractLinearSolver& linSolver;
public:
  int report;
};

class TangentialStepArmijo : public TangentialStep
{
public:
  TangentialStepArmijo(StepPolicy& s, AbstractNorm& n,  AbstractLinearSolver& linSolver_) 
    : sp(s), norm(n), linSolver(linSolver_)
  {}

  void setLinearization(AbstractFunctional* f, AbstractLinearization* lin)
  {
    linPtr = lin;  
    fPtr = f;
    iterate=linPtr->getOrigin().clone();
    correction=linPtr->getOrigin().clone();
    trialIterate=linPtr->getOrigin().clone();
    gradient=linPtr->getOrigin().clone();
  }
  virtual void getTrialIterate(
    AbstractFunctionSpaceElement& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractFunctionSpaceElement const& iterate,
    double damping)
  {
    trialIterate = iterate;
//    trialIterate.axpy(damping,correction);
    trialIterate+=correction;
    sp.getTrialIterate(trialIterate,correction, iterate,damping);
  }

  void getSearchDirection(AbstractFunctionSpaceElement& scorrection, AbstractFunctionSpaceElement const& correction);

private:
  int performArmijoLoop();

bool evaluateTrialIterate(
  AbstractFunctionSpaceElement const& trialIterate, 
  AbstractFunctionSpaceElement const& correction, 
  AbstractLinearization const& lin);


  std::unique_ptr<AbstractFunctionSpaceElement> iterate,correction, trialIterate, scorrection, gradient;
  StepPolicy& sp;
  AbstractNorm& norm;
  AbstractLinearization* linPtr; 
  std::unique_ptr<AbstractLinearization> simplifiedLinearization;
  AbstractFunctional* fPtr;
  AbstractLinearSolver& linSolver;
  
  double dampingFactor, functionalAtIterate, armijopar;
public:
  int report;
};



template<class UpdatePolicy>
class StepPolicyProjectedRHS : public UpdatePolicy
{
public:
  StepPolicyProjectedRHS(int rhsstart_, int rhsend_) : rhsstart(rhsstart_), rhsend(rhsend_) {}

  virtual void getSearchDirection(AbstractFunctionSpaceElement& correction)
  {  
    assert(UpdatePolicy::linPtr);
    assert(UpdatePolicy::linearSolver);
    std::vector<double> sorig, sp;
    for(int i=0; i<rhsstart; ++i) sp.push_back(0.0);
    for(int i=rhsstart; i<rhsend; ++i) sp.push_back(1.0);
    for(int i=rhsend; i <correction.nComponents(); ++i) sp.push_back(0.0);
    UpdatePolicy::linPtr->getScalePars(sorig);
    UpdatePolicy::linPtr->scaleRHS(sp);
    UpdatePolicy::linearSolver->solve(correction, *UpdatePolicy::linPtr);
    UpdatePolicy::linPtr->scaleRHS(sorig);
    correction *= -1.0;
  }

  virtual void getSimplifiedSearchDirection(AbstractFunctionSpaceElement& correction, 
                                    AbstractLinearization& linRHS)
  {
    assert(UpdatePolicy::linPtr);
    assert(UpdatePolicy::linearSolver);

    std::vector<double> sorig, sp;
    for(int i=0; i<rhsstart; ++i) sp.push_back(0.0);
    for(int i=rhsstart; i<rhsend; ++i) sp.push_back(1.0);
    for(int i=rhsend ; i <correction.nComponents(); ++i) sp.push_back(0.0);
    linRHS.getScalePars(sorig);
    linRHS.scaleRHS(sp);
    UpdatePolicy::linearSolver->resolve(correction, linRHS, *UpdatePolicy::linPtr);
    linRHS.scaleRHS(sorig);
    correction *= -1.0;
  }
  int rhsstart, rhsend;
};

template <class Newton>
class NormalStepNewton : public NormalStep
{
public:
  NormalStepNewton(StepPolicy& s, AbstractNorm& n, typename Newton::Parameters& p, AbstractLinearSolver& linSolver_) 
    : sp(s), norm(n), pars(p), linSolver(linSolver_)
  {
    
  }

  void setLinearization(AbstractFunctional* f, AbstractLinearization* lin)
  {
    pars.reset();
    newton.reset(new Newton(linSolver,sp,norm,pars));
    linPtr = lin;  
    fPtr = f;
  }
  virtual void getTrialIterate(
    AbstractFunctionSpaceElement& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractFunctionSpaceElement const& iterate,
    double damping)
  {
    trialIterate = iterate;
    trialIterate.axpy(damping,correction);
//    sp.getTrialIterate(trialIterate,correction, iterate,damping);
  }

  virtual void getSearchDirection(AbstractFunctionSpaceElement& correction)
  {  
    assert(linPtr);
    assert(newton.get());
    correction=linPtr->getOrigin();
    newton->reportOnIteration(report);
    newton->oneStep(fPtr, linPtr, correction);
    correction -= linPtr->getOrigin();
  }
private:
  std::unique_ptr<Newton> newton;
  StepPolicy& sp;
  AbstractNorm& norm;
  AbstractLinearization* linPtr;
  AbstractFunctional* fPtr;
  typename Newton::Parameters& pars;
  AbstractLinearSolver& linSolver;
public:
  int report;
};


class CompositeStep : public Algorithm
{
public:

  virtual ~CompositeStep() {}

/// Create Newton's Method, providing a solver, a norm and algorithmic parameters
  CompositeStep(NormalStep& normalStep_, TangentialStep& tangentialStep_, AbstractNorm& n, IterationParameters& p_): 
    normalStep(normalStep_), tangentialStep(tangentialStep_), norm(n), p(p_) {}

/// Solve the system f=0 with starting value x. On (successful) exit, the solution is x, otherwise it is left unmodified.
  void solve(AbstractFunctional* f, AbstractFunctionSpaceElement& x);

  virtual IterationParameters const& getParameters() {return p;}
/// set the desired accuracy
  void setDesiredAccuracy(double da) { assert(da>0); p.desiredAccuracy=da;}
/// set the desired accuracy
  virtual void setDesiredRelativeAccuracy(double ra) { assert(false); }
/// Reset all algorithmic parameters to their default values
  void resetParameters() {p.reset();}

  int stepsPerformed() {return step;}
  int maxSteps() {return p.maxSteps;}
protected:
  virtual NewtonsMethod::AcceptanceTest evaluateTrialIterate(
    AbstractFunctionSpaceElement const& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractLinearization const& lin) { return NewtonsMethod::AcceptanceTestPassed;}


  /// Return true, if convergence is detected, false otherwise
  virtual NewtonsMethod::Convergence convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate);

  virtual void updateIterate(AbstractFunctionSpaceElement& iterate, 
                             AbstractFunctionSpaceElement& trialIterate,
                             AbstractLinearization const& lin)
  {
    iterate.swap(trialIterate);
  }


  virtual NewtonsMethod::RegularityTest regularityTest(double scalingFactor)
  {
    return NewtonsMethod::RegularityTestPassed;
  }

  void terminationMessage(int flag);

private:
  int runAlgorithm();
  AbstractFunctional* functional;
  NormalStep& normalStep;
  TangentialStep& tangentialStep;
  AbstractNorm& norm;
  IterationParameters& p;
  int step;

  std::unique_ptr<AbstractFunctionSpaceElement> iterate, trialIterate, correction, scorrection, trialIterate2;

  std::unique_ptr<AbstractLinearization> newtonLinearization;
};

}  // namespace Kaskade
#endif
