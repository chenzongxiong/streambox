#ifndef ALGORITHM_BASE_HH
#define ALGORITHM_BASE_HH

#include <vector>
#include <string>
#include <iostream>

namespace LQAction 
{

enum toDo {LogValue, Reset, Devalidate};

typedef enum toDo ToDo;
}

namespace Kaskade
{
  enum class RegularityTest {Failed, Passed};
  enum class AcceptanceTest {InitialGuess, FailedRecompute, Failed, NormalStepFailed, TangentialStepFailed, Passed, LeftAdmissibleDomain};
  enum class Convergence {Missed, Achieved};

/// Class that represents a quantity that can be logged during the course of an algortihm
/**
Useful as a member of IterationParameters or a derived class
Supports a vaildity check and an assertion which returns the name of the quantity
Examples of use: e.g. GuardedCovariantNewtonParameters
 */
template<typename T> 
class LoggedQuantity
{
public:
  void doAction(LQAction::ToDo td, std::string const& name_="noName")
  {
    switch(td)
    {
    case LQAction::LogValue : logValue(); break;
    case LQAction::Reset : reset(name_); break;
    case LQAction::Devalidate : devalidate(); break;
    default : assert(!"strange error in doAction");
    }
  }

/// Delete any information contained in this class, on exit: state like freshly constructed
  void reset(std::string const& name_) 
  { 
    logged=false; 
    valid=false; 
    name=name_; 
    logger.resize(0); 
    zero=0.0; 
  }


  LoggedQuantity(std::string const& name_) { reset(name_); }
  LoggedQuantity() { reset("noName");}

///Insert current value into log buffer
  void logValue() { 
    if(valid && !logged) logger.push_back(currentValue); 
    if(!valid) logger.resize(logger.size()+1);
    logged=true; 
  }

///Insert current value into log buffer at a certain number i (used seldomly)
  void logValue(int i) { if(valid) logger[i]=currentValue; logged=true; }

/// Set current value
  LoggedQuantity<T>& operator=(T const& rvalue) { currentValue=rvalue; logged=false; valid=true; return *this; }
 
/// Get size of log buffer
  int size() {int sz= logger.size(); if(valid && !logged) sz++; return sz;}

/// Returns true, iff there is a valid current value
  bool isValid() {return valid;}

/// devalidate current value
  void devalidate() {valid=false;}

/// returns true, iff current value is already logged
  bool isLogged() {return logged;}

/// Return value from log buffer (used seldomly)
  T const& operator[](int i) 
  { 
    if(i < logger.size() && i >= 0) return logger[i]; 
    if(i==logger.size() && valid && !logged) return currentValue; 
    namedAssertion(false,"Index failure!");
    return currentValue;
  }

/// Directly access current value (may fail in special cases due to type non-uniqueness issues)
  operator T&() { namedAssertion(valid,"TypeCast: operator T&()"); return currentValue;}

/// access current value via a function member: Use this if operator
/// T&() fails

  T& value_nonconst() { valid = true;  namedAssertion(valid,"TypeCast: value_nonconst()"); return currentValue;}

  T const& value() const { namedAssertion(valid,"TypeCast: value()"); return currentValue;}

/// print log-buffer into a stream, to be used for analysis of an algorithm
  void print(std::ostream& s)
  {
    s.setf(std::ios::scientific,std::ios::floatfield);
    s.precision(16);
    for(int i=0; i!=size(); i++)
      s << (*this)[i] << std::endl;
  }

private:
  std::vector<T> logger;
  T currentValue;
  bool logged;
  bool valid;
  std::string name;
  void namedAssertion(bool val,std::string caller) const
  { 
    if(!val) 
    {
      std::cout << name << " from " << caller << ":" << std::endl; 
      assert(0); 
    }
  }
  double zero;
};

//--------------------------------------------------
/// Base class for algorithmic parameters. 
class IterationParameters
{
public:
  IterationParameters(double desiredAccuracy_, int maxSteps_)
    : desiredAccuracy(desiredAccuracy_), 
      maxSteps(maxSteps_)
  {
  }

//Parameters with values that must be supplied by client
  double desiredAccuracy;
  int maxSteps;
  int termination;

  virtual ~IterationParameters() {}

/// Log all quantities in this class
  void logStep() { 
    doForAll(LQAction::LogValue);
  }
  
/// Reset all quantities in this class
  virtual void reset() {
    doForAll(LQAction::Reset);
  }
  
protected:

/// To be overloaded by derived class
/** For any LoggedQuantity quant, declared in the derived class, insert quant.doAction(td) 
 */
  virtual void doForAll(LQAction::ToDo td)
  {
  }


  friend class NewtonsMethod;
};


//-------------------------------------------------------------------------------------------

class AbstractFunctionSpaceElement;

/// Base class for algorithms. Provides a unified interface and some simple wrapper routines, which perform optional timing, etc.
class Algorithm
{
public:
    
  Algorithm() : report(0), measureTime(false) {}
  virtual ~Algorithm(){}

  void performTiming(bool doit) {measureTime=doit; }
  void reportOnIteration(int level) {report=level; }

protected:
  virtual void initialize(){};
  virtual void finalize(int){};
  virtual void terminationMessage(int);
/// Run algorithm, completely with initialization and finalization
  int algorithmWrapper();
/// Run one step of algorithm
  int oneStepWrapper();

  int report;

private:
  virtual int runAlgorithm()=0;

  bool measureTime;
};

}  // namespace Kaskade

#endif
