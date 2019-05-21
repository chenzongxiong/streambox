#ifndef NEWTON_DAMPED_HH
#define NEWTON_DAMPED_HH
  
#include <memory>
#include "newton_base.hh"

namespace Kaskade
{

/// Paramters that are used and Logged by GuardedCovariantNewtonMethod
class GuardedCovariantNewtonParameters : public NewtonParameters
{
public:
  GuardedCovariantNewtonParameters(double desiredAccuracy_, int maxSteps_) :
    NewtonParameters(desiredAccuracy_, maxSteps_),
    maxContraction(1.0)
  {
    reset();
  }


// Parameters with predefined values that can be changed by client
  double maxContraction;          // default = 1.0
  LoggedQuantity<double> omega0;
  LoggedQuantity<double> omega;
  double accuracyReached;

protected:
  friend class GuardedCovariantNewtonMethod;

  virtual void doForAll(LQAction::ToDo td)
  {
    NewtonParameters::doForAll(td);    
    Theta.doAction(td,"Theta");
    omega.doAction(td,"omega");
    normSCorr.doAction(td,"|SCorr|");
    normCorrLast.doAction(td,"|corrLast|");
    omega0.doAction(td,"omega0");
  }

  LoggedQuantity<double> Theta;
  LoggedQuantity<double> normSCorr;
  LoggedQuantity<double> normCorrLast;

};



//------------------------------------------------------------------------------

/// Paramters that are used and Logged by GuardedCovariantNewtonMethod
class DampedCovariantNewtonParameters : public GuardedCovariantNewtonParameters
{
public:
  DampedCovariantNewtonParameters(double desiredAccuracy_, int maxSteps_, double minDampingFactor_=1e-12) :
    GuardedCovariantNewtonParameters(desiredAccuracy_, maxSteps_),
    initDampingFactor(1.0),
    minDampingFactor(minDampingFactor_),
    reduceOnOutsideDomain(0.5),
    relativeAccuracy(0.0)
  {
    lastrejected.rejected=false;
    reset();
  }
public:
  friend class DampedCovariantNewtonMethod;

  void doForAll(LQAction::ToDo td) { 
    GuardedCovariantNewtonParameters::doForAll(td);
    dampingFactor.doAction(td,"dampingFactor");
    SCorrByCorr.doAction(td,"SCorrByCorr");
    totalCorrection.doAction(td,"totalCorrection");
    lengthOfStep.doAction(td,"lengthOfStep");
    absoluteAccuracyLast.doAction(td,"absoluteAccuracyLast");
    dampingForDomain=1.0;
    lastrejected.rejected=false;
  }
// Initial damping factor
  double initDampingFactor;
// Minimal damping factor, below which the Newton iteration terminates unsuccessfully
  double minDampingFactor;
// How to reduce the damping factor, if iterate is outside the domain of definition of f
  double reduceOnOutsideDomain;

  double dampingForDomain;

  double relativeAccuracy;

  LoggedQuantity<double> dampingFactor;
  LoggedQuantity<double> SCorrByCorr;
  LoggedQuantity<double> absoluteAccuracyLast;
  LoggedQuantity<double> totalCorrection;
  LoggedQuantity<double> lengthOfStep;

  struct LastRejected
  {
    double dampingFactor, omega, Theta;
    bool rejected;
  };

  LastRejected lastrejected;  
};


/// Damped Newton's method that measures the Newton contraction via a simplified Newton step
/** Affine covariant damping strategy after Deuflhard '04   
 */
class DampedCovariantNewtonMethod : public NewtonsMethod
{
public:
  typedef DampedCovariantNewtonParameters Parameters; 

  DampedCovariantNewtonMethod(AbstractNewtonDirection& l, 
                              AbstractChart& chart_,
                              AbstractNorm& n, 
                              DampedCovariantNewtonParameters& p_): 
    NewtonsMethod(l, chart_,n,p_), p(p_) {}

  virtual DampedCovariantNewtonParameters& getParameters() {return p;}

  void setDampingFactorStart(double df) { p.initDampingFactor=df; }

protected:
  virtual RegularityTest regularityTest(double scalingFactor);
  virtual void computeTrialIterate(AbstractFunctionSpaceElement& trialIterate, AbstractFunctionSpaceElement const& direction, AbstractLinearization const& lin);

  virtual void initialize();
  virtual void initNewtonStep();
  virtual double& dampingFactor() {return p.dampingFactor.value_nonconst(); }


  virtual void predictNextDampingFactor(AbstractFunctionSpaceElement& correction);

  void setDesiredRelativeAccuracy(double ra) { assert(ra>0); p.relativeAccuracy=ra; }

  virtual AcceptanceTest evaluateTrialIterate(AbstractFunctionSpaceElement const& trialIterate, 
                                              AbstractFunctionSpaceElement const& correction, 
                                              AbstractLinearization const& lin);

  virtual void updateIterate(AbstractFunctionSpaceElement& iterate, 
                             AbstractFunctionSpaceElement& trialIterate,
                             AbstractLinearization const& lin);

  virtual Convergence convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate);

  DampedCovariantNewtonParameters &p;
  std::unique_ptr<AbstractFunctionSpaceElement> scorrection, startIterate;
  std::unique_ptr<AbstractFunctionSpaceElement> auxvector;

};

//-----------------------------------------------------------------------------------------

class ModifiedCovariantNewtonParameters : public DampedCovariantNewtonParameters
{
public:
  ModifiedCovariantNewtonParameters(double desiredAccuracy_, int maxSteps_, double ThetaMaxAllowed_, double mdf_=1e-12) :
    DampedCovariantNewtonParameters(desiredAccuracy_, maxSteps_, mdf_) , ThetaMaxAllowed(ThetaMaxAllowed_) {}

 
protected:
  friend class ModifiedCovariantNewtonMethod;
  friend class StateConstraintsNewtonMethod;

  virtual void doForAll(LQAction::ToDo td)
  {
    DampedCovariantNewtonParameters::doForAll(td);    
    Theta0.doAction(td,"|ModCorr|");
    normModCorr.doAction(td,"|ModCorr|");
    normModSCorr.doAction(td,"|ModSCorr|");
    totalCorrection.doAction(td,"|totCorr|");
  }
  LoggedQuantity<double> Theta0;
  LoggedQuantity<double> normModCorr;
  LoggedQuantity<double> normModSCorr;
  LoggedQuantity<double> totalCorrection;
  double relativeAccuracy, ThetaMaxAllowed;
};

///// A Newton method that uses a pointwise damping strategy to cope with bounds
//class ModifiedCovariantNewtonMethod : public NewtonsMethod
//{
//
//public:
//  typedef ModifiedCovariantNewtonParameters Parameters;
//
//
//  ModifiedCovariantNewtonMethod(AbstractLinearSolver& l,
//                                AbstractChart& chart_,
//                                AbstractNorm& n,
//                                ModifiedCovariantNewtonParameters& p_):
//    NewtonsMethod(l, chart_,n,p_), p(p_) { failedbydamping=false; }
//
//protected:
//
//  void setDesiredRelativeAccuracy(double ra) { assert(ra>0); p.relativeAccuracy=ra; }
//
//  void initialize();
//
//  void predictNextDampingFactor(AbstractFunctionSpaceElement& correction) { p.normCorr=norm(correction); }
//
//  void initNewtonStep();
//
//  virtual AcceptanceTest evaluateTrialIterate(
//    AbstractFunctionSpaceElement const& trialIterate,
//    AbstractFunctionSpaceElement const& correction,
//    AbstractLinearization const& lin);
//
//  void updateIterate(AbstractFunctionSpaceElement& iterate,
//                     AbstractFunctionSpaceElement& trialIterate,
//                     AbstractLinearization const& lin);
//
//  NewtonsMethod::Convergence convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate);
//
//  NewtonsMethod::RegularityTest regularityTest(double scalingFactor);
//
//private:
//  double mu;
//  double nsc;
//  bool failedbydamping;
//  ModifiedCovariantNewtonParameters &p;
//  std::unique_ptr<AbstractFunctionSpaceElement> scorrection;
//  std::unique_ptr<AbstractFunctionSpaceElement> startIterate;
//};



/// A Newton method that uses a pointwise damping strategy to cope with bounds
class StateConstraintsNewtonMethod : public NewtonsMethod
{

public:
  typedef ModifiedCovariantNewtonParameters Parameters;


  StateConstraintsNewtonMethod(AbstractNewtonDirection& l,
                                AbstractChart& chart_,
                                AbstractNorm& n,
                                Parameters& p_):
    NewtonsMethod(l, chart_,n,p_), p(p_) {}

protected:

  void setDesiredRelativeAccuracy(double ra) { assert(ra>0); p.relativeAccuracy=ra; }

  void initialize();

  void predictNextDampingFactor(AbstractFunctionSpaceElement& correction) { p.normCorr=norm(correction); }

  void initNewtonStep();

  virtual AcceptanceTest evaluateTrialIterate(
    AbstractFunctionSpaceElement const& trialIterate,
    AbstractFunctionSpaceElement const& correction,
    AbstractLinearization const& lin);

  void updateIterate(AbstractFunctionSpaceElement& iterate,
                     AbstractFunctionSpaceElement& trialIterate,
                     AbstractLinearization const& lin);

  Convergence convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate);

  RegularityTest regularityTest(double scalingFactor);

private:
  double mu;
  double nsc;
  bool failedbydamping;
  Parameters &p;
  std::unique_ptr<AbstractFunctionSpaceElement> scorrection;
  std::unique_ptr<AbstractFunctionSpaceElement> startIterate;
};

}  // namespace Kaskade
#endif
