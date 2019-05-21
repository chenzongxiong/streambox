#ifndef ADAPTIVE_REFINEMENT_HH
#define ADAPTIVE_REFINEMENT_HH

#include "abstract_interface.hh"
#include "algorithm_base.hh"
#include <limits>
#include <memory>
#include <utility>

namespace Kaskade
{

/** \addtogroup alg */
/** @{ */

/// Parameters that are used by SimpleAdaptiveRefinement
class AdaptiveParameters
{
public:
  AdaptiveParameters(double desiredAccuracy_, int maxSteps_, double bulkmin_, double bulkmax_,
                     double desAbsAcc_=0.0, int maxpatches_=std::numeric_limits<int>::max(), double minRelAcc=10.0)
    : desiredAccuracy(desiredAccuracy_), 
      maxSteps(maxSteps_),
      bulkmin(bulkmin_),
      bulkmax(bulkmax_),
      desiredAbsAccuracy(desAbsAcc_),
      maxpatches(maxpatches_),
      maxSizeReached(false),
      minRelAccuracy(minRelAcc),
      gridHasChanged(false),
      relativeAccuracy("relativeAccuracy"),
      absoluteAccuracy("absoluteAccuracy")
  {}

//Parameters with values that must be supplied by client
  double desiredAccuracy;
  int maxSteps;
  double bulkmin, bulkmax;
  double desiredAbsAccuracy;
  int maxpatches;
  bool maxSizeReached;
  double minRelAccuracy;
  bool gridHasChanged;

  virtual ~AdaptiveParameters() {}
  
protected:
  friend class SimpleAdaptiveRefinement;
  friend class TestErrorEstimator;
  friend class FinalAdaptiveRefinement;
  friend class FixedSolverWithErrorEstimate;

  virtual void logStep() { relativeAccuracy.logValue(); absoluteAccuracy.logValue(); }
  LoggedQuantity<double> relativeAccuracy;
  LoggedQuantity<double> absoluteAccuracy;
};


/// Adaptive refinement algorithm
/** Implementation uses the strategy pattern, and should be flexible to use. In order
 * to make it work one has to provide implementations of the following interfaces:
 *
 * - AbstractLinearSolver : responsible for the solution of the linear discrete problems
 * - AbstractGrid : coarse grained interface to a grid, mostly mark/refine
 */

class SimpleAdaptiveRefinement : public Algorithm, public AbstractNewtonDirection
{
// Public interface as Algorithm
public:
/// Construction
  SimpleAdaptiveRefinement(AbstractNewtonDirection& als,
                           AbstractErrorEstimator& aee,
                           AbstractAdaptiveGrid& ag,
                           AbstractNorm& an,
                           AdaptiveParameters& p_)
    : fixedSolver(als), errorEstimator(aee), norm(an), grid(ag), p(p_), alwaysestimate(true)
  {}

/// To be used as an algorithm 
  virtual int runAlgorithm();
  virtual ~SimpleAdaptiveRefinement() {};

//Protected interface as Algorithm for inheritance  
protected:
  virtual bool convergenceTest(AbstractErrorEstimate const& estimate, AbstractAdaptiveGrid const&);

  virtual double getBulk(int step) const;

public:
/// set relative accuracy
  virtual void setRelativeAccuracy(double accuracy)  {p.desiredAccuracy=accuracy;}
/// set absolute accuracy
  virtual void setAbsoluteAccuracy(double accuracy) { p.desiredAbsAccuracy=accuracy;}
/// get achieved relative accuracy
  double getRelativeAccuracy() { return p.relativeAccuracy; }
/// get achieved absolute accuracy
  double getAbsoluteAccuracy() { return p.absoluteAccuracy; }

  bool improvementPossible() { return !p.maxSizeReached; }

/// solve problem
  virtual void doSolve(AbstractFunctionSpaceElement& corr, AbstractLinearization& lin);

  virtual void doResolve(AbstractFunctionSpaceElement& correction, 
                         AbstractLinearization const& lin,
                         AbstractLinearization const& olin) const
  { fixedSolver.doResolve(correction, lin, olin); }

  virtual void doResolve(AbstractFunctionSpaceElement& correction, 
                         AbstractLinearization const& lin) const
  { fixedSolver.doResolve(correction, lin); }


  virtual bool changedGrid() { return p.gridHasChanged; }


protected:
  std::shared_ptr<AbstractErrorEstimate> estimate;
  AbstractNewtonDirection& fixedSolver;
  AbstractErrorEstimator& errorEstimator;
  AbstractNorm& norm;
  AbstractAdaptiveGrid& grid;
  AbstractLinearization* linearization;
  AbstractFunctionSpaceElement* correction;
  std::unique_ptr<AbstractFunctionSpaceElement> oldcorrection;
  double noldold;
  virtual void initialize();
  virtual void finalize(int flag);
  virtual void terminationMessage(int flag);
public:
  AdaptiveParameters& p;
  bool alwaysestimate;
};


/// Performs error estimation and creates error indicators, but does not refine. Good for testing error estimators
class FinalAdaptiveRefinement : public SimpleAdaptiveRefinement
{
public:
  FinalAdaptiveRefinement(AbstractNewtonDirection& als,
                          AbstractErrorEstimator& aee,
                          AbstractAdaptiveGrid& ag,
                          AbstractNorm& an,
                          AdaptiveParameters& p_) : SimpleAdaptiveRefinement(als,aee,ag,an,p_)
  {}
  virtual double getBulk(int step) const;
/// set relative accuracy
  virtual void setRelativeAccuracy(double accuracy)  {p.desiredAccuracy=accuracy*0.9;}
/// set absolute accuracy
  virtual void setAbsoluteAccuracy(double accuracy) { p.desiredAbsAccuracy=accuracy*0.9;}
};

/// Performs error estimation and creates error indicators, but does not refine. Good for testing error estimators
class FixedSolverWithErrorEstimate : public SimpleAdaptiveRefinement
{
public:
  FixedSolverWithErrorEstimate(AbstractNewtonDirection& als,
                               AbstractErrorEstimator& aee,
                               AbstractAdaptiveGrid& ag,
                               AbstractNorm& an,
                               AdaptiveParameters& p_) : SimpleAdaptiveRefinement(als,aee,ag,an,p_)
  {}
  virtual int runAlgorithm();
};

class TestErrorEstimator : public Algorithm
{
// Public interface as Algorithm
public:
/// Construction
  TestErrorEstimator(AbstractNewtonDirection& als,
                     AbstractErrorEstimator& aee,
                     AbstractAdaptiveGrid& ag,
                     AbstractNorm& an,
                     AdaptiveParameters& p_)
    : fixedSolver(als), errorEstimator(aee), norm(an), grid(ag), p(p_), alwaysestimate(true)
  {}

/// To be used as an algorithm 
  virtual int runAlgorithm();
//  virtual ~SimpleAdaptiveRefinement() {};

//Protected interface as Algorithm for inheritance  
protected:
  virtual bool convergenceTest(AbstractErrorEstimate const& estimate, AbstractAdaptiveGrid const&);

  virtual double getBulk(int step) const;

public:
/// set relative accuracy
  void setRelativeAccuracy(double accuracy)  {p.desiredAccuracy=accuracy;}
/// set absolute accuracy
  void setAbsoluteAccuracy(double accuracy) { p.desiredAbsAccuracy=accuracy;}
/// get achieved relative accuracy
  double getRelativeAccuracy() { return p.relativeAccuracy; }
/// get achieved absolute accuracy
  double getAbsoluteAccuracy() { return p.absoluteAccuracy; }

  bool improvementPossible() { return !p.maxSizeReached; }

/// solve problem
  void solve(AbstractFunctionSpaceElement& corr, AbstractLinearization& lin);

protected:
  AbstractNewtonDirection& fixedSolver;
  AbstractErrorEstimator& errorEstimator;
  AbstractNorm& norm;
  AbstractAdaptiveGrid& grid;
  AbstractLinearization* linearization;
  AbstractFunctionSpaceElement* correction;
  std::unique_ptr<AbstractFunctionSpaceElement> oldcorrection;
  double noldold;
  virtual void initialize();
  virtual void finalize(int flag);
  virtual void terminationMessage(int flag);
public:
  AdaptiveParameters& p;
  bool alwaysestimate;
};


/** @}*/
} // namespace Kaskade
#endif
