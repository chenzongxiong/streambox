#ifndef HOMOTOPY_BASE_HH
#define HOMOTOPY_BASE_HH

#include <vector>
#include <memory> // unique_ptr
#include <cmath>
#include <fstream>

#include "dune/grid/config.h"

#include "abstract_interface.hh"
#include "algorithm_base.hh"
#include "newton_damped.hh"
#include "modelfunction.hh"
#include <boost/timer/timer.hpp>

namespace Kaskade
{

class NewtonsMethod;

  /** Abstract algorithm of interior point method.
  */
class InteriorPointSqrtModel;
template<class IPF, class DomainVector> class InteriorPointTangentialPredictor;
class InteriorPointSlopeEstimate;


class InteriorPointParameters : IterationParameters
{
public:
  InteriorPointParameters(double mu_, double desiredAccuracy_, double desiredContraction_, int maxSteps_, double muFinal_, double relDist, double accfinal_,
                          double reductionStart_,double reductionFactorBest_, double redWorst_=1-1e-4)
    : IterationParameters(desiredAccuracy_, maxSteps_),
      desiredContraction(desiredContraction_), 
      muFinal(muFinal_),
      reductionFactorWorst(redWorst_),
      firstStepFailureFactor(2.0),
      reductionStart(reductionStart_),
      relDistanceToPath(relDist),
      accfinal(accfinal_),
      reductionFactorBest(reductionFactorBest_),
      maxTrialSteps(10)
  {
    reset();
    mu=mu_;
  }

//Parameters with values that must be supplied by client
  double desiredContraction, muFinal, reductionFactorWorst, firstStepFailureFactor, reductionStart, relDistanceToPath, accfinal, reductionFactorBest;
  int maxTrialSteps;

  virtual ~InteriorPointParameters() {}

  int termination;

  void reset()
  {
    IterationParameters::reset();
    stepSuccessful=true;
    correctorTermination=-1;
  }

  LoggedQuantity<double> mu, sigma;
  
protected:

  bool stepSuccessful; 

  int correctorTermination;


  LoggedQuantity<double> deltaMu, j, jp, muTrial;

  friend class HomotopyBase;
  friend class InteriorPointSimple;
  friend class InteriorPointSqrtModel;
  template<class A, class B> friend class InteriorPointTangentialPredictor;
  friend class InteriorPointSlopeEstimate;

  virtual void doForAll(LQAction::ToDo td)
  {
    IterationParameters::doForAll(td);    
    mu.doAction(td,"mu");
    muTrial.doAction(td,"muTrial");
    j.doAction(td,"jp");
    jp.doAction(td,"j");
    sigma.doAction(td,"sigma");
    deltaMu.doAction(td,"deltaMu");
  }


};

/// Base class for homotopy methods. Here, the main algorithm is programmed
class HomotopyBase : public Algorithm
  {
  public:
    HomotopyBase(NewtonsMethod& n_, InteriorPointParameters& p_) :
      corrector(n_), p(p_)
    {}
    virtual ~HomotopyBase()
    { }

    void solve(AbstractParameterFunctional* f,AbstractFunctionSpaceElement& x);
  protected:

/// Optional overload by derived class
    virtual void initialize() {};

/// Optional overload by derived class
    virtual void finalize() {};

    virtual void logQuantities();

/// Default: classical predictor
    virtual void computePredictor(int step);

/// New mu if corrector was successful
    virtual double muOnSuccess(int step) = 0;

/// New mu if corrector failed
    virtual double muOnFailure() = 0;

/// mu for successful termination
    virtual double muFinal() = 0;

/// update of path parameters, e.g. eta and omega
    virtual void updateModelOfHomotopy(int step) = 0;

/// what should be done before corrector
    virtual void initializeCorrector() = 0;

/// what should be done after corrector
    virtual void finalizeCorrector() = 0;

/// make trial iterate to current iterate (on successful corrector)
    virtual void updateIterate() = 0;

/// make old iterate to current iterate (on failed corrector)
    virtual void recoverIterate() = 0;

/// estimate of length of homopoty path
    virtual double lengthOfPath() = 0;

/// convergence test. returns true if method converged and solution is found.
    virtual int convergenceTest();
    
/// to be performed, if convergence of path has occured
    virtual void finalizeHomotopy(){};

/// output of a termination message
    virtual void terminationMessage(int errorFlag);

    NewtonsMethod& corrector;
    InteriorPointParameters& p;
    AbstractParameterFunctional* functional;

    std::unique_ptr<AbstractFunctionSpaceElement> trialIterate;
    int step;

  private:

    bool stepSuccessful();

    void computeCorrector();
    void computeGapParameter();
    int runAlgorithm();
  };

/// Very simple implementation of homotopy method: fixed stepsize
class InteriorPointSimple : public HomotopyBase
{
public:
  InteriorPointSimple(NewtonsMethod& n_, InteriorPointParameters& p_) : 
    HomotopyBase(n_, p_) {}

  virtual double muOnSuccess(int step) { return p.mu*p.reductionStart;};
  virtual double muOnFailure() { return p.mu*p.reductionStart;};
  virtual double muFinal() { return 1e-7;};

  virtual void updateModelOfHomotopy(int step) {};

  virtual void initializeCorrector() {};
  virtual void finalizeCorrector() {};

  virtual void updateIterate() {};
  virtual void recoverIterate() {};

  virtual double lengthOfPath() { return std::sqrt(p.mu); };
};


class InteriorPointParametersSqrt : IterationParameters
{
public:
  InteriorPointParametersSqrt()
    : IterationParameters(0.0,0)
  {
    reset();
  }
  friend class InteriorPointSqrtModel;
  template<class A, class B> friend class InteriorPointTangentialPredictor;
  friend class InteriorPointSlopeEstimate;

  virtual ~InteriorPointParametersSqrt() {}
  
protected:

  LoggedQuantity<double> omega;
  LoggedQuantity<double> eta;
  LoggedQuantity<double> slope;
  LoggedQuantity<double> curvature;  
  LoggedQuantity<double> jpl;
  LoggedQuantity<double> accuracyCorrector;
  LoggedQuantity<int> newtonSum;
  LoggedQuantity<double> timemeasured;

  virtual void doForAll(LQAction::ToDo td)
  {
    omega.doAction(td,"omega");
    eta.doAction(td,"eta");
    slope.doAction(td,"slope");
    curvature.doAction(td,"curvature");
    jpl.doAction(td,"jpl");
    accuracyCorrector.doAction(td,"accuracyCorrector");
    newtonSum.doAction(td,"newtonSum");
    timemeasured.doAction(td,"timemeasured");
  }
};


/// InteriorPointMethod with a sqrt-model of the path: eta ~ mu^{-1/2}, omega ~mu^{-1/2}
class InteriorPointSqrtModel : public HomotopyBase
{
public:
  InteriorPointSqrtModel(NewtonsMethod& n_, AbstractNorm const& norm_, 
                         InteriorPointParameters& p_, AbstractNewtonDirection& solver_, 
                         AbstractNorm const* normPlain_=0) : 
    HomotopyBase(n_, p_),
    norm(norm_),
    solver(solver_),
    normPlain(normPlain_)
  {
    if(!normPlain) normPlain = &norm;
  }

  virtual void initialize() { iterate=trialIterate->clone(); nNewton=0;};
  virtual double muOnSuccess(int step);
  virtual double muOnFailure();
  virtual double muFinal() { return p.muFinal;};
  virtual void updateModelOfHomotopy(int step); 
  virtual void initializeCorrector(); 
  virtual void finalizeCorrector() { nNewton += corrector.stepsPerformed(); pp.newtonSum=nNewton;};
  virtual void updateIterate();
  virtual void recoverIterate() { *trialIterate=*iterate; }
  virtual double lengthOfPath() {if(pp.eta.isValid()) return 2*p.mu*pp.slope; else return 1e300;};
  void printDiagnosis();

  virtual void logQuantities()
  {
    HomotopyBase::logQuantities();
    pp.logStep();
    printDiagnosis();
  }

private:
  std::unique_ptr<AbstractFunctionSpaceElement> iterate;
  AbstractNorm const& norm;
  AbstractNewtonDirection& solver;
  InteriorPointParametersSqrt pp;
  JModel jModel;
  JModelLin jModelL;
  AbstractNorm const* normPlain;
  int nNewton;
  boost::timer::cpu_timer overalltime;
};

/// Performs bisection algorithm to find a zero of a 
/**************************************************
 * Assumptions
 *
 *  a < b
 *
 *  Equation f with has double operator(double ) monotonically INCREASING
 *
 ************************************************************/

template<class Equation>
double bisection(double a, double b, Equation const& f, double accuracy, int& iterations)
{
  iterations = 0;
  double u;
  do {
    ++iterations;
    u=(a+b)/2;
    if(u == a) return b;
    if(f(u) > 0) b=u; else a=u;
  } 
  while(b-a >= accuracy && (a+b)/2 != a && (a+b)/2 != b && iterations<50);
  if(iterations > 48) std::cout << "Warning: bisection algorithm performed > 48 iterations" << std::endl;
  return u;  
}

class InteriorPointSlopeEstimate : public HomotopyBase
{

public:
  InteriorPointSlopeEstimate(NewtonsMethod& n_, AbstractNorm const& norm_, InteriorPointParameters& p_, 
                          AbstractNewtonDirection& solver_, AbstractNorm const* normPlain_=0,
                          NewtonsMethod* finalsolver_=0) : 
    HomotopyBase(n_, p_),
    norm(norm_),
    solver(solver_),
    stp(0),
    normPlain(normPlain_),
    finalsolver(finalsolver_)
  {
    if(!normPlain_) normPlain=&norm;
  }

  virtual void initialize() { 
    iterate=trialIterate->clone(); 
    tangent=trialIterate->clone(); 
    nNewton=0;
  };
  virtual void finalize() {};

  virtual double muOnSuccess(int step);

  virtual double muOnFailure() { 
    p.sigma=0.5+p.sigma/2.0;    
    return p.mu*p.sigma;
  };

  virtual double muFinal() { return p.muFinal;};

  virtual void updateModelOfHomotopy(int step);

  virtual void computePredictor(int step){}


  virtual void initializeCorrector(); 
  virtual void finalizeCorrector() {    nNewton += corrector.stepsPerformed(); pp.newtonSum=nNewton;};

  virtual void finalizeHomotopy();

  virtual void updateIterate(); 

  virtual void recoverIterate() { *trialIterate=*iterate; }


  virtual double lengthOfPath() 
  { 
    if(pp.eta.isValid())
      return 2*p.mu*pp.slope; 
    else
      return 1e300;
  };

  virtual void logQuantities()
  {
    HomotopyBase::logQuantities();
    pp.logStep();
    printDiagnosis();
  }

  void printDiagnosis();

private:
  std::unique_ptr<AbstractFunctionSpaceElement> iterate, tangent;

  AbstractNorm const& norm;
  AbstractNewtonDirection& solver;
  InteriorPointParametersSqrt pp;
  JModelLin jModelL;
  int stp;
  AbstractNorm const* normPlain;
  NewtonsMethod* finalsolver;
  int nNewton;
  boost::timer::cpu_timer overalltime;
};

}  // namespace Kaskade

#endif
