#ifndef NEWTON_BASE_HH
#define NEWTON_BASE_HH
  
#include <vector>
#include <memory>

#include "abstract_interface.hh"
#include "algorithm_base.hh"

namespace Kaskade
{

/** \addtogroup alg */
/** @{ */


//------------------------------------------------------------------------------------------

class NewtonParameters : public IterationParameters
{
public:
  NewtonParameters(double desiredAccuracy_, int maxSteps_)
    : IterationParameters(desiredAccuracy_, maxSteps_), reuseFactorization(false)
  {
    reset();
  }


//Parameters with values that must be supplied by client
//  double desiredAccuracy;
//  int maxSteps;

  virtual ~NewtonParameters() {}
  
  virtual void reset() {
    doForAll(LQAction::Reset);
  }

  bool reuseFactorization;

protected:
  friend class NewtonsMethod;
  friend class SimpleNewtonMethod;

  void logStep() { 
    doForAll(LQAction::LogValue);
  }
  
  virtual void doForAll(LQAction::ToDo td)
  {
    normCorr.doAction(td,"|Corr|");
  }

  LoggedQuantity<double> normCorr;
};

class IdentityChart : public AbstractChart
{
public:
  void addPerturbation(AbstractFunctionSpaceElement& newIterate, 
                       AbstractFunctionSpaceElement const& perturbation, 
                       AbstractLinearization const& lin,
                       std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > basis = std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >()) const
  {
    newIterate = lin.getOrigin();
    newIterate += perturbation;
  }
};

//class TrivialAbstractChart : public AbstractChart
//{
//public:
//  void addPerturbation(AbstractFunctionSpaceElement& newIterate,
//                       AbstractFunctionSpaceElement const& perturbation,
//                       AbstractLinearization const& lin,
//                       std::vector<AbstractFunctionSpaceElement* > basis = std::vector<AbstractFunctionSpaceElement* >()) const
//  {
//    newIterate = lin.getOrigin();
//    newIterate += perturbation;
//  }
//};

/// Base class for Newton's method. Defines the main algorithm, potentially using damping. 
/** This algorithm is defined in terms of abstract quantities. This is for the sake of reusability.
It is design via virtual methods should provide flexibility to implement variants of this algorithm.
As often with virtual functions: the base class specifies, when they are called, the derived class specifies what they do
 */
class NewtonsMethod : public Algorithm
{
public:
  typedef NewtonParameters Parameters; 

  virtual ~NewtonsMethod() {}

/// Create Newton's Method, providing a solver, a norm and algorithmic parameters
  NewtonsMethod(AbstractNewtonDirection& l, AbstractChart& chart_, AbstractNorm& n, Parameters& p_): 
    linearSolver(l), chart(chart_), norm(n), p(p_) {}

/// Solve the system f=0 with starting value x. On (successful) exit, the solution is x, otherwise it is left unmodified.
  void solve(AbstractFunctional* f, AbstractFunctionSpaceElement& x);

/// Perform one Newton step (which may be several trial steps, if damping is applied)
  void oneStep(AbstractFunctional* f, AbstractLinearization* l, AbstractFunctionSpaceElement& x);

  void resolve(AbstractFunctionSpaceElement& x, AbstractLinearization const&  l) 
  { 
    if(p.reuseFactorization)
      linearSolver.simplified(x,l); 
    else
    {
      std::cout << "No Factorization available: set NewtonParameters::reuseFactorization=true" << std::endl;
    }
}

  virtual Parameters& getParameters() {return p;}
/// set the desired accuracy
  void setDesiredAccuracy(double da) { assert(da>0); p.desiredAccuracy=da;}
/// set the desired accuracy
  virtual void setDesiredRelativeAccuracy(double ra) { assert(false); }
/// Reset all algorithmic parameters to their default values
  void resetParameters() {p.reset();}

  int stepsPerformed() {return step;}

  bool changedGrid() { return gridHasChanged; }

  virtual void getSearchDirection(AbstractFunctionSpaceElement& direction) 
  {
    linearSolver.ordinary(direction, *newtonLinearization);
    direction *= -1.0;    
  }


  virtual void computeTrialIterate(AbstractFunctionSpaceElement& trialIterate, AbstractFunctionSpaceElement const& direction, AbstractLinearization const& lin)
  {
          chart.addPerturbation(trialIterate, direction, lin);
  }

  AbstractLinearization const& getLastLinearization() { assert(newtonLinearization.get()); return *newtonLinearization; }

  AbstractLinearization const& getLastSimplifiedLinearization() { assert(simplifiedLinearization.get()); return *simplifiedLinearization; }

protected:

/// Called before Iteration
  virtual void initialize();
/// Called after Iteration
  virtual void finalize(int);
/// Performs logging of quantities
  virtual void logQuantities() { p.logStep();}

  virtual void initNewtonStep() {};

  /// Return true, if convergence is detected, false otherwise
  virtual Convergence convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate) = 0;

/// Should compute a damping factor
  virtual void predictNextDampingFactor(AbstractFunctionSpaceElement& correction) {}

/// Return a damping factor
  virtual double& dampingFactor() { dummy=1.0; return dummy;}
/// return the maximal number of steps
  virtual double maxSteps() {return p.maxSteps;}

/// update an accepted iterate, default: iterate=trialIterate
  virtual void updateIterate(AbstractFunctionSpaceElement& iterate, 
                             AbstractFunctionSpaceElement& trialIterate,
                             AbstractLinearization const& lin);

  /** trial iterate = iterate + scalingFactor * correction
   * if this was successful, return RegularityTest::Passed, otherwise
   * (e.g. if scalingFactor is too small) RegularityTest::Failed
   * false will lead to an unsuccessful exit of algorithm
   */
  virtual RegularityTest regularityTest(double scalingFactor)
  {
    return RegularityTest::Passed;
  }

  /// function to evaluate the quality of the trial iterate.
  /// returns true if trial iterate is acceptable as new iterate
  /// i.e. if some damping factor is of appropriate size
  virtual AcceptanceTest evaluateTrialIterate(
    AbstractFunctionSpaceElement const& trialIterate, 
    AbstractFunctionSpaceElement const& correction, 
    AbstractLinearization const& lin) = 0;


  template<class T>
  class CycleDetection
  {
  public:
    CycleDetection() : hasIncreased(false) {};

    bool detectCycling(T& quantity)
    {
      if(quantityCycle.size()!=0)
        if(*std::min_element(quantityCycle.begin(),quantityCycle.end())<quantity) hasIncreased=true;

      if(hasIncreased && quantity < *std::min_element(quantityCycle.begin(),quantityCycle.end()))
      {
        quantity = *std::min_element(quantityCycle.begin(),quantityCycle.end());
        std::cout << "Detected Cycling: decreased below already accepted step size:" << quantity << std::endl;
        return true;
      }

      for(int i = 0; i< quantityCycle.size(); i++)
        if(std::abs(quantityCycle[i]-quantity)<std::abs(quantity)*std::numeric_limits<double>::epsilon()){
          quantity = quantityCycle[i+1];
          std::cout << "Detected Cycling: two identical step sizes:" << quantity << std::endl;
          return true;
        }

      if(quantityCycle.size() > 20){
        quantity=quantity/4.0;
        std::cout << "Detected Cycling: too many trials:" << quantity << std::endl;
        return true;
      }
      quantityCycle.push_back(quantity);
      return false;
    }

  private:
    std::vector<T> quantityCycle;
    bool hasIncreased;
  };

  virtual void terminationMessage(int flag);

  AbstractFunctional* functional;
  AbstractNewtonDirection& linearSolver;
  AbstractChart& chart;
  AbstractNorm& norm;
  Parameters& p;

  std::unique_ptr<AbstractFunctionSpaceElement> iterate, trialIterate, correction;

  std::unique_ptr<AbstractLinearization> newtonLinearization;
  std::unique_ptr<AbstractLinearization> simplifiedLinearization;
  AbstractLinearization* newtonPtr;

  int step;

private:
  int runAlgorithm();

  NewtonsMethod(NewtonsMethod&);

  double dummy;


  bool doOneStep;
  bool gridHasChanged;
};



/** @}*/

}  // namespace Kaskade
#endif
