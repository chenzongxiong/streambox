#include "abstract_interface.hh"
#include "algorithm_base.hh"
#include <boost/timer.hpp>
#include <time.h>


namespace Kaskade
{

int Algorithm::algorithmWrapper()
{
  this->initialize();
  int terminationFlag;

  if(measureTime)
  {
    clock_t t0=time(0);
    boost::timer refTimer;
    terminationFlag=runAlgorithm();
    if(refTimer.elapsed() > 0.2)
      std::cout << "CPU Time: " << refTimer.elapsed() << " sec. "<< "Wallclock-Time: " << time(0)-t0 << " sec."<< std::endl;
  } else
  {
    terminationFlag=runAlgorithm();
  }
  if(report) terminationMessage(terminationFlag);
  this->finalize(terminationFlag);
  return terminationFlag;
}

int Algorithm::oneStepWrapper()
{
  this->initialize();
  int terminationFlag;

  if(measureTime)
  {
    boost::timer refTimer;
    terminationFlag=runAlgorithm();
    if(refTimer.elapsed() > 0.2)
      std::cout << "Time Measured: " << refTimer.elapsed() << " sec."<< std::endl;
  } else
  {
    terminationFlag=runAlgorithm();
  }
  return terminationFlag;
}


void Algorithm::terminationMessage(int flag)
{
  std::cout << "Termination with (unknown) exit code: " << flag << "!" << std::endl;
}

}  // namespace Kaskade