#include "adaptive_refinement.hh"
#include <cmath>

namespace Kaskade
{

void SimpleAdaptiveRefinement::doSolve(AbstractFunctionSpaceElement& corr, AbstractLinearization& lin)
{
  p.gridHasChanged = false;
  correction = &corr;  
  *correction *= 0.0;
  if(report >= 2)
  {
    oldcorrection = corr.clone();  
    *oldcorrection *= 0.0;
  }
  linearization=&lin;
  algorithmWrapper();
}

void SimpleAdaptiveRefinement::terminationMessage(int flag)
{
  std::cout << "Adaptive Refinement: "; 
  switch((int)flag)
  {
  case 1  : std::cout << "Desired accuracy reached!" << std::endl; break;
  case 2 : std::cout << "Maximum number of steps reached!" << std::endl; break;
  case 3  : std::cout << "Maximum number of cells reached!" << std::endl; p.relativeAccuracy=0.0; break;
  case 4  : std::cout << "Number of refinement steps: 0" << std::endl; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

void SimpleAdaptiveRefinement::initialize()
{
  // WARNING: No connectToSignalForFlush in AbstractLinearization. Ws 2014-12-11
  std::cerr << "couldn't connect in file " << __FILE__ << ":" << __LINE__ << "\n";
  abort();
//   linearization->connectToSignalForFlush(grid.gridWillChange);
  noldold=0.0;
}

void SimpleAdaptiveRefinement::finalize(int flag)
{
  if(flag < 0) grid.flushMarks();
}


bool SimpleAdaptiveRefinement::convergenceTest(AbstractErrorEstimate const& estimate, AbstractAdaptiveGrid const&)
{
  double ncorr = norm(*correction);
  if(report >= 2) {
    std::cout << "|Corr|: "<< ncorr << " ";
     *oldcorrection -= *correction;     
     double noldnew=norm(*oldcorrection);
     std::cout << "|DCorr|: "<< noldnew << " |DC|/|DC_| " << noldnew/noldold;
     if(p.absoluteAccuracy.isValid()) std::cout << " Eff:" << p.absoluteAccuracy/noldnew; 
     std::cout << std::endl;
     noldold=noldnew;
     *oldcorrection = *correction;
  }
  double abslast=1.0;
  if(p.absoluteAccuracy.isValid()) abslast=p.absoluteAccuracy;
  p.absoluteAccuracy = estimate.absoluteError();
  p.relativeAccuracy.value_nonconst() = p.absoluteAccuracy/ncorr;
  if(report) 
  {
    if(p.desiredAbsAccuracy==0)
      std::cout << "   |Err|/|Corr|:" << p.relativeAccuracy.value();
    else
      std::cout << "   |Err|/|Corr|:" << p.relativeAccuracy.value() << " |Err|/Tol:" << p.absoluteAccuracy/p.desiredAbsAccuracy << " Tol:" << p.desiredAbsAccuracy;
    if(p.absoluteAccuracy.isValid() && report >=2)
      std::cout << " |Err|:" << p.absoluteAccuracy << " |Err|/|Err_-|: " << p.absoluteAccuracy/abslast << std::endl;
    else
      std::cout << " |Err|:" << p.absoluteAccuracy << std::endl;
  }
  if(p.relativeAccuracy.value() < p.desiredAccuracy || p.absoluteAccuracy < p.desiredAbsAccuracy) return true;
  if(grid.size() >= p.maxpatches)
  {
    p.maxSizeReached = true;
    return true;
  }
  return false;
}

double SimpleAdaptiveRefinement::getBulk(int step) const
{
  double cbulk = 1.0;
//  if(p.relativeAccuracy.isValid() && step!=1) 
  cbulk = std::sqrt(1-std::pow(p.desiredAccuracy/p.relativeAccuracy*0.9,2));
//  cbulk = std::min(cbulk,std::sqrt(1-std::pow(p.desiredAbsAccuracy/p.absoluteAccuracy*0.9,2)));
  double bulk = std::max(p.bulkmin, std::min(cbulk,p.bulkmax));
  return bulk;
}

int SimpleAdaptiveRefinement::runAlgorithm()
{
//  std::cout << "Function value: " << linearization->eval() << std::endl;
  p.gridHasChanged = false;

  for(int i=1; i<= std::max(p.maxSteps,1); i++)
  {
    if(grid.getNMarked() != 0) p.gridHasChanged = true;
    grid.adapt();
    if(report) 
    {
      std::cout << "#Patches: " << grid.size() << " ";
    }
    fixedSolver.doSolve(*correction,*linearization);
    norm.setOrigin(*linearization);
    if(grid.size() >= p.maxpatches) p.maxSizeReached=true; else p.maxSizeReached=false;
    if(p.maxSteps==0) { p.relativeAccuracy=0.0; return 4; }
    if(p.maxSizeReached && !alwaysestimate) return 3;
    estimate.reset(errorEstimator.createEstimate(*correction,*linearization).release());
    bool accuracyReached = convergenceTest(*estimate,grid);;
    if(estimate.get() && grid.size() < p.maxpatches && p.maxSteps!=0 && getBulk(i) >= 1e-5)
    {
      std::cout << "Marking..." << getBulk(i) << std::flush;
      grid.mark(*estimate,getBulk(i));
      std::cout << " #marked:" << grid.getNMarked() << std::endl;
    }

    if(p.relativeAccuracy.value() > p.minRelAccuracy &&!p.maxSizeReached && (p.maxSteps == i))
    {
      if(grid.getNMarked() != 0) p.gridHasChanged = true;
      grid.adapt();
      fixedSolver.doSolve(*correction,*linearization);
      grid.flushMarks();
    }
    if(accuracyReached) return 1;
  }
  return 2;
}


int FixedSolverWithErrorEstimate::runAlgorithm()
{
  p.relativeAccuracy = 0;
  fixedSolver.doSolve(*correction,*linearization);
  norm.setOrigin(*linearization);
  // ?? estimate is not used
  std::unique_ptr<AbstractErrorEstimate> estimate(errorEstimator.createEstimate(*correction,*linearization));
  return 1;
}

double FinalAdaptiveRefinement::getBulk(int step) const
{
    double cbulk = std::sqrt(1-std::pow(p.desiredAbsAccuracy/p.absoluteAccuracy*0.9,2));
    double bulk = std::max(p.bulkmin, std::min(cbulk,p.bulkmax));
    return bulk;
}

void TestErrorEstimator::solve(AbstractFunctionSpaceElement& corr, AbstractLinearization& lin)
{
  p.gridHasChanged = false;
  correction = &corr;  
  *correction *= 0.0;
  if(report >= 2)
  {
    oldcorrection = corr.clone();  
    *oldcorrection *= 0.0;
  }
  linearization=&lin;
  algorithmWrapper();
}

void TestErrorEstimator::terminationMessage(int flag)
{
  std::cout << "Adaptive Refinement: "; 
  switch((int)flag)
  {
  case 1  : std::cout << "Desired accuracy reached!" << std::endl; break;
  case 2 : std::cout << "Maximum number of steps reached!" << std::endl;; break;
  case 3  : std::cout << "Maximum number of cells reached!" << std::endl; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

void TestErrorEstimator::initialize()
{
  // WARNING: No connectToSignalForFlush in AbstractLinearization. Ws 2014-12-11
  std::cerr << "couldn't connect in file " << __FILE__ << ":" << __LINE__ << "\n";
  abort();
//  linearization->connectToSignalForFlush(grid.gridWillChange);
  noldold=0.0;
}

void TestErrorEstimator::finalize(int flag)
{
  if(flag < 0) grid.flushMarks();
}


bool TestErrorEstimator::convergenceTest(AbstractErrorEstimate const& estimate, AbstractAdaptiveGrid const&)
{
  double ncorr = norm(*correction);
  if(report >= 2) {
    std::cout << "|Corr|: "<< ncorr << " ";
     *oldcorrection -= *correction;     
     double noldnew=norm(*oldcorrection);
     std::cout << "|DCorr|: "<< noldnew << " |DC|/|DC_| " << noldnew/noldold;
     if(p.absoluteAccuracy.isValid()) std::cout << " Eff:" << p.absoluteAccuracy/noldnew; 
     std::cout << std::endl;
     noldold=noldnew;
     *oldcorrection = *correction;
  }
  double abslast=1.0;
  if(p.absoluteAccuracy.isValid()) abslast=p.absoluteAccuracy;
  p.absoluteAccuracy = estimate.absoluteError();
  p.relativeAccuracy = p.absoluteAccuracy/ncorr;
  if(report) 
  {
    if(p.desiredAbsAccuracy==0)
      std::cout << "   |Err|/|Corr|: " << p.relativeAccuracy;
    else
      std::cout << "   |Err|/|Corr|: " << p.relativeAccuracy << " |Err|/Tol: " << p.absoluteAccuracy/p.desiredAbsAccuracy;
    if(p.absoluteAccuracy.isValid() && report >=2)
      std::cout << " |Err|: " << p.absoluteAccuracy << " |Err|/|Err_-|: " << p.absoluteAccuracy/abslast << std::endl;
    else
      std::cout << " |Err|: " << p.absoluteAccuracy << std::endl;
  }
  if(p.relativeAccuracy < p.desiredAccuracy || p.absoluteAccuracy < p.desiredAbsAccuracy) return true;
  if(grid.size() >= p.maxpatches)
  {
    p.maxSizeReached = true;
    return true;
  }
  return false;
}

double TestErrorEstimator::getBulk(int step) const
{
    return p.bulkmin;
}

int TestErrorEstimator::runAlgorithm()
{
  typedef std::shared_ptr<AbstractFunctionSpaceElement> CorrPtr;
  std::vector<CorrPtr> corrections; 
  std::vector<double> estimates;

  p.gridHasChanged = false;

  std::cout << "Function value: " << linearization->eval() << std::endl;
  for(int i=1; i<= std::max(p.maxSteps,1); i++)
  {
    if(grid.getNMarked() != 0) p.gridHasChanged = true;
    grid.adapt();
    if(report) 
    {
      std::cout << "#Patches: " << grid.size() << " ";
    }
    fixedSolver.doSolve(*correction,*linearization);
    CorrPtr corrptr((correction->clone()).release());
    corrections.push_back(corrptr);
    norm.setOrigin(*linearization);
    if(grid.size() >= p.maxpatches) p.maxSizeReached=true; else p.maxSizeReached=false;
    if((p.maxSteps==0 || p.maxSizeReached) && !alwaysestimate) return 1;
    std::unique_ptr<AbstractErrorEstimate> estimate(errorEstimator.createEstimate(*correction,*linearization));
    bool accuracyReached = convergenceTest(*estimate,grid);;

    estimates.push_back(estimate->absoluteError());

    if(grid.size() < p.maxpatches && p.maxSteps!=0) grid.mark(*estimate,getBulk(i));
    if(accuracyReached)
    {
      for(int k=0; k<corrections.size(); ++k)
      {
        *corrections[k] -= *corrections[corrections.size()-1];
        double normc=norm(*corrections[k]);
        std::cout << "|c_1-c_end|:" << normc << "  est:" << estimates[k] << "  ets/dc:" <<  estimates[k]/normc << std::endl;
      }
      return 1;
    }
  }
  for(int k=0; k<corrections.size(); ++k)
  {
    *corrections[k] -= *corrections[corrections.size()-1];
    double normc=norm(*corrections[k]);
    std::cout << "|c_1-c_end|:" << normc << "  est:" << estimates[k] << "  est/dc:" <<  estimates[k]/normc << std::endl;
  }
  if(p.maxSizeReached) return 3;
  return 2;
}

}  // namespace Kaskade
