#ifndef HYP_IP_HH
#define HYP_IP_HH

#include <memory>
#include "newton_base.hh"
#include "newton_damped.hh"
#include "searchspace.hh"

namespace Kaskade
{

class HypIP;

class HypIPParameters : public IterationParameters
{
public:
  HypIPParameters(double desiredAccuracy_, int maxSteps_) :
    IterationParameters(desiredAccuracy_, maxSteps_),
    gamma(1.0),
    alpha(1.0),
    sigma(1.0/3.0),
    Thetad(0.25),
    ThetaAim(0.9),
    ThetaNormal(0.5),
    ThetaMax(1.0),
    fracBdr(0.25),
    safetyFactor(1.0),
    nu0(1.0),
    tau0(0.0),
    firstmu(true)
  {
    reset();
  }


// Parameters with predefined values that can be changed by client
  double gamma;      // weighing between optimality and feasibility
  double alpha;      // weighing between optimality and feasibility
  double sigma;      // Bit-counting Lemma 
  double Thetad;     // desired contraction close to constraints
  double ThetaAim;     // required contraction 
  double ThetaNormal;     // required contraction 
  double ThetaMax;     // required contraction 
  double fracBdr;     // fraction, how much slacks may approach the boundary
  double safetyFactor;
  double nu0,tau0;   // starting values for damping     
  bool firstmu;
  LoggedQuantity<double> omegaC, omegaL;

protected:
  friend class HypIP;

  virtual void doForAll(LQAction::ToDo td)
  {
    omegaC.doAction(td,"omegaC");
    omegaL.doAction(td,"omegaL");
  }

};


class HypIP : public Algorithm
{
public:

  virtual ~HypIP() {}

/// Create Newton's Method, providing a solver, a norm and algorithmic parameters
  HypIP(SearchSpaceCreator& searchSpace_, AbstractScalarProduct& nL, AbstractScalarProduct& nC, AbstractChart& chart_,  HypIPParameters& p_,int dimx_,
        AbstractAdaptiveGrid* grid_=0,  AbstractCompositeStepErrorEstimator* errorEstimator_=0): 
    RemainderTerm(0.0), nu0(0.0), searchSpace(searchSpace_), normL(nL), normC(nC), chart(chart_), grid(grid_), p(p_), dimx(dimx_), errorEstimator(errorEstimator_)
  {
  }


/// Solve the system f=0 with starting value x. On (successful) exit, the solution is x, otherwise it is left unmodified.
void solve(AbstractFunctional* fN, AbstractFunctional* fT, AbstractFunctionSpaceElement& x);

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

  void dampingLoop();
  int refinementLoop();



  virtual AcceptanceTest evaluateCorrection(AbstractFunctionSpaceElement& correction, AbstractLinearization& lin, CUBThetaModelFunction& mF, double f0, std::vector<double>& coeff);

//   virtual NewtonsMethod::AcceptanceTest evaluateCorrectionOld(
//     AbstractFunctionSpaceElement& correction, 
//     AbstractLinearization& lin,
//     CUBModelFunction& mF,
//     double f0,
//     std::vector<double>& coeff);



  /// Return true, if convergence is detected, false otherwise
  virtual Convergence
  convergenceTest(AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate, std::vector<double>& coeff);

  virtual void updateIterate(AbstractFunctionSpaceElement& iterate, 
                             AbstractFunctionSpaceElement& trialIterate,
                             AbstractLinearization const& lin)
  {
//    normC.add(trialIterate,iterate,*scorrection,1.0,p.fracBdr); 
//    iterate.swap(trialIterate);
  }


  virtual RegularityTest regularityTest(double scalingFactor)
  {
    return RegularityTest::Passed;
  }

  void terminationMessage(int flag);

private:
  int runAlgorithm();
  double omegaL,omegaC, Theta, RelNormalSteps, RemainderTerm, normalstepnorm, snormalstepnorm, sbarstepnorm, theta2, nu0;
  SearchSpaceCreator& searchSpace;
  AbstractFunctional* functionalN;
  AbstractFunctional* functionalT;
  AbstractNormalSolver* normalSolver;
  AbstractScalarProduct& normL;
  AbstractScalarProduct& normC;
  AbstractChart& chart;
  AbstractAdaptiveGrid* grid;
  HypIPParameters& p;
  int step;
  int dimx;
  std::vector<double> coeff;

  std::unique_ptr<AbstractFunctionSpaceElement> iterate, trialIterate, correction, scorrection;
  std::vector<AbstractFunctionSpaceElement* > basisVectors;
  std::vector<std::tr1::shared_ptr<AbstractFunctionSpaceElement> > bV;
  std::vector<double> coeffs;

  std::unique_ptr<AbstractLinearization> normalLinearization;
  std::unique_ptr<AbstractLinearization> tangentialLinearization;
  std::unique_ptr<AbstractErrorEstimate> estimate;
  AbstractCompositeStepErrorEstimator* errorEstimator;
};

}  // namespace Kaskade
#endif
