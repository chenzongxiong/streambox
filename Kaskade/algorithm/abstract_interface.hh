#ifndef ABSTRACT_INTERFACE_HH
#define ABSTRACT_INTERFACE_HH

#include <memory>
#include <iostream>
#include <boost/signals2.hpp>
#include <string>
#include <cmath>
#include <vector>

#include <dune/grid/config.h>
#include "linalg/triplet.hh"

namespace Kaskade
{
  /**
   * @file
   * @brief  Interfaces for function space oriented algorithms
   * @author Anton Schiela
   */
  /** \addtogroup abstract*/
  /**@{*/

  /// Abstract Vector for function space algorithms
  class AbstractFunctionSpaceElement
  {
  public:
    // No copy constructor.
    AbstractFunctionSpaceElement(AbstractFunctionSpaceElement const &) = delete;

    virtual ~AbstractFunctionSpaceElement() {}

    /// *this += alpha*l
    AbstractFunctionSpaceElement& axpy(double alpha, AbstractFunctionSpaceElement const& l, int component)
    {
      if(this==&l)
      {
        std::vector<double> alphaV(nComponents(),alpha+1.0);
        return doscale(alphaV);
      }
      else
        return doaxpy(alpha, l, component);
    }

    AbstractFunctionSpaceElement& axpy(double alpha, AbstractFunctionSpaceElement const&l, std::string const role)
    {
      return axpy_role(alpha,l,role);
    }

    AbstractFunctionSpaceElement& axpy_role(double alpha, AbstractFunctionSpaceElement const&l, std::string const role)
    {
      for(int i=0; i<nComponents(); ++i)
      {
        if(getRole(i)==role) axpy(alpha,l,i);
      }
      return *this;
    }

    /// *this += alpha*l
    AbstractFunctionSpaceElement& axpy(double alpha, AbstractFunctionSpaceElement const& l)
    {
      if(this==&l)
      {
        std::vector<double> alphaV(nComponents(),alpha+1.0);
        return doscale(alphaV);
      }
      else
        return doaxpy_(alpha, l);
    }


    /// Assignment
    AbstractFunctionSpaceElement& operator=(AbstractFunctionSpaceElement const& v)
    {
      if(this!=&v) return doassign(v);
      return *this;
    }

    /// Basic vector arithmetic
    AbstractFunctionSpaceElement& operator+=(AbstractFunctionSpaceElement const& v) { return doaxpy_(1.0,v); }
    /// Basic vector arithmetic
    AbstractFunctionSpaceElement& operator-=(AbstractFunctionSpaceElement const& v) { return doaxpy_(-1.0,v); }
    /// Basic vector arithmetic
    AbstractFunctionSpaceElement& operator*=(double lambda)
    {
      std::vector<double> lambdaV(nComponents(),lambda);
      return doscale(lambdaV);
    }

    /// Scaling each component of the vector separately
    AbstractFunctionSpaceElement& operator*=(std::vector<double>const& lambda)
    {
      return doscale(lambda);
    }

    /// Interpret *this as a dual vector, and apply it to v
    /** Duality is currently not represented by types, hence users have to make
  sure that the dual pairing <*this,v> is mathematically meaningful. The standard 
 implementation of the dual pairing in Bridge::Vector is
 the dot-product between the coordinate representations of *this and v. 
 For that *this should be the result of AbstractLinearization::evald
     */
    double applyAsDualTo(AbstractFunctionSpaceElement const& v, int component) const
    {
      return doapplyAsDualTo(v,component,component+1);
    }

    double applyAsDualTo(AbstractFunctionSpaceElement const& v) const
    {
      return doapplyAsDualTo(v,0,nComponents());
    }

    double applyAsDualTo(AbstractFunctionSpaceElement const& v, std::string const& role) const
    {
      return applyAsDualTo_role(v,role);
    }

    double applyAsDualTo_role(AbstractFunctionSpaceElement const&v, std::string const role) const
    {
      double res(0.0);
      for(int i=0; i<nComponents(); ++i)
      {
        if(getRole(i)==role) res+=applyAsDualTo(v,i);
      }
      return res;
    }

    virtual std::string getRole(int component) const = 0;

    /// Shallow swap
    void swap(AbstractFunctionSpaceElement& v) { doswap(v); }

    virtual int nComponents() const = 0;

    /// Optional output
    virtual void writeToFile(std::string const& file, bool append, int order=1) const {}

    /// Optional output
    virtual void print(std::string const& message="") const { std::cout << "AbstractFunctionSpaceElement: No printing available!" << std::endl;}

    /// Construction of a vector of the same type
    virtual std::unique_ptr<AbstractFunctionSpaceElement> clone() const = 0;

    /// Construction of a vector of the same type
    virtual std::unique_ptr<AbstractFunctionSpaceElement> initZeroVector() const = 0;

  private:
    /// Interface for implementation
    virtual AbstractFunctionSpaceElement& doaxpy(double alpha, AbstractFunctionSpaceElement const& l, int component) = 0;
    virtual AbstractFunctionSpaceElement& doaxpy_(double alpha, AbstractFunctionSpaceElement const& l) = 0;
    /// Interface for implementation
    virtual AbstractFunctionSpaceElement& doscale(std::vector<double> const& lambda) = 0;
    /// Interface for implementation
    virtual AbstractFunctionSpaceElement& doassign(AbstractFunctionSpaceElement const& l) = 0;
    /// Interface for implementation
    virtual void doswap(AbstractFunctionSpaceElement& l) = 0;
    virtual double doapplyAsDualTo(AbstractFunctionSpaceElement const& v,int vbegin, int vend) const = 0;

  protected:
    AbstractFunctionSpaceElement() {}
  };
  
  /// Abstract linearization.
  /** Given a nonlinear functional, this class represents second order information at the point of linearization
   */
  class AbstractLinearization
  {
  public:
    /// Evaluate functional
    virtual double eval() const = 0;
    /// Evaluate L1 norm of integrand
    virtual double evalL1norm() const = 0;
    /// Evaluate derivative
    virtual void evald(AbstractFunctionSpaceElement &g, int rbegin=0, int rend=-1) const = 0;
    /// Evaluate hessian times second argument: y = y+ddf*x
    void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
    {
      d2axpy(1.0,y,x,rbegin,rend,cbegin,cend);
    }
    /// Evaluate hessian times second argument: y = y+ddf*x
    void ddtxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
    {
      d2taxpy(1.0,y,x,rbegin,rend,cbegin,cend);
    }
    /// Evaluate hessian times second argument: y = y+a*ddf*x
    virtual void d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const = 0;
    /// Evaluate hessian times second argument: y = y+a*ddf*x
    virtual void d2taxpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const = 0;
    /// Access matrix representation if available.
    virtual void getMatrixBlocks(MatrixAsTriplet<double>& mat, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const = 0;//{ assert(!"not implemented"); }
    /// Get point of linearization
    virtual AbstractFunctionSpaceElement const& getOrigin() const = 0;

    virtual void precompute() = 0;

    virtual void flush() = 0;

    virtual ~AbstractLinearization() {}
  };

  class AbstractPreconditioner
  {
  public:
    virtual void setLinearization(AbstractLinearization&) = 0;
    virtual void apply(AbstractFunctionSpaceElement const& x, AbstractFunctionSpaceElement& Px) = 0;
  };

  /// Abstract connection.
  /**
   * Connect signal to flush() (in general for invalidating assembled data on grid refinement)
   */
  class AbstractFlushConnection
  {
  public:
    virtual ~AbstractFlushConnection(){}
    virtual void connectToSignalForFlush(boost::signals2::signal<void()>&) = 0;
  };

  class AbstractConnectedLinearization : public AbstractLinearization, public AbstractFlushConnection
  {};

  /// Representation of a nonlinear functional
  class AbstractFunctional
  {
  public:
    virtual double evaluate(AbstractFunctionSpaceElement const& x) const { assert(!"not implemented"); return 0; }
    virtual std::unique_ptr<AbstractLinearization> getLinearization(AbstractFunctionSpaceElement const& x) const = 0;
    virtual std::unique_ptr<AbstractFunctionSpaceElement> getImageVector(AbstractFunctionSpaceElement const& x) const = 0;
    virtual bool inDomain(AbstractFunctionSpaceElement const& x) const { return true; }
    virtual ~AbstractFunctional() {}
  };

  /// Representation of parameters
  class AbstractParameters
  {
  public:
    virtual ~AbstractParameters() {}
  };


  /// ...for parameter dependent functionals, implements AbstractParameters
  template<typename ParameterType>
  class Parameters : public AbstractParameters
  {
  public:
    explicit Parameters(ParameterType const& p) : pref(p) {}
    ParameterType const& getPars() const {return pref;}
  private:
    ParameterType const& pref;
  };

  ///
  template<class ParameterType>
  Parameters<ParameterType> makePars(ParameterType const& p) { return Parameters<ParameterType>(p); }


  /// Creates a functional from a homotopy of functionals by inserting a parameter
  class AbstractParameterFunctional
  {
  public:
    virtual std::unique_ptr<AbstractFunctional> getFunctional(AbstractParameters const&) const = 0;
    virtual std::unique_ptr<AbstractFunctional> getParameterLinFunctional(AbstractParameters const& p) const
          { assert(!"ParameterFunctional: Linearization not Implemented"); return getFunctional(p);}
    virtual std::unique_ptr<AbstractFunctional> getLinFunctionValue(AbstractParameters const& p) const
          { assert(!"ParameterFunctional: Linearization not Implemented"); return getFunctional(p);}
    virtual ~AbstractParameterFunctional() {}
  };

  /// Class that models the functionality of a (possibly inexact) linear solver
  class AbstractNewtonDirection
  {
  public:
    /** This method is not const, since the internal state of the solver may be modified */
    /** Computes an undamped ordinary Newton step, i.e., \f$correction=-F'(x)^{-1}F(x)\f$*/
    void ordinary(AbstractFunctionSpaceElement& correction, AbstractLinearization& linearization)
    {  doSolve(correction, linearization);
    correction *= -1.0;
    }
    /** Computes an undamped simplified Newton step, i.e., correction=-F'(xold)^{-1}F(x)*/
    void simplified(AbstractFunctionSpaceElement& correction, AbstractLinearization const& linearization, AbstractLinearization const& oldlinearization) const
    {  doResolve(correction, linearization, oldlinearization);
    correction *= -1.0;
    }

    /** Computes an undamped simplified Newton step, i.e., correction=-F'(xold)^{-1}F(x)*/
    /** where F'(xold)^{-1} has been kept from an earlier ordinary() call */
    void simplified(AbstractFunctionSpaceElement& correction, AbstractLinearization const& linearization) const
    {  doResolve(correction, linearization);
    correction *= -1.0;
    }

    /// Specify accuracy that should be achieved
    virtual void setRelativeAccuracy(double accuracy) = 0;
    /// Specify accuracy that should be achieved
    virtual void setAbsoluteAccuracy(double) {}
    /// Get accuracy that was actually achieved
    virtual double getRelativeAccuracy() = 0;
    /// Get accuracy that was actually achieved
    virtual double getAbsoluteAccuracy() = 0;

    virtual bool improvementPossible() = 0;

    virtual bool changedGrid() { return false; }

    virtual ~AbstractNewtonDirection() {}

    /// Solve Newton System: correction = +F'(iterate)^{-1}F(iterate)
    virtual void doSolve(AbstractFunctionSpaceElement& correction,
        AbstractLinearization& lin) = 0;
    /// Solve simplified Newton system: scorrection = +F'(iterate)^{-1}F(trialIterate)
    virtual void doResolve(AbstractFunctionSpaceElement& correction,
        AbstractLinearization const& lin,
        AbstractLinearization const& olin) const = 0;

    /// Solve simplified Newton system: scorrection = +F'(iterate)^{-1}F(trialIterate)
    virtual void doResolve(AbstractFunctionSpaceElement& correction,
        AbstractLinearization const& lin) const = 0;

    mutable boost::signals2::signal<void()> changed;
  };


  class AbstractChart
  {
  public:
    virtual void addPerturbation
    (AbstractFunctionSpaceElement& newIterate,
        AbstractFunctionSpaceElement const& perturbation,
        AbstractLinearization const& lin,
        std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > basis = std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >()) const = 0;

    virtual std::unique_ptr<AbstractChart> clone() const = 0;
  };


  class AbstractNorm
  {
  public:
    virtual ~AbstractNorm(){}
    virtual double operator()(AbstractFunctionSpaceElement const&) const = 0;
    virtual void setOrigin(AbstractLinearization const&) {}
  };

  class AbstractScalarProduct : public AbstractNorm
  {
  public:
    virtual double operator()(AbstractFunctionSpaceElement const&, AbstractFunctionSpaceElement const&) const = 0;
    double operator()(AbstractFunctionSpaceElement const& v) const
    {
      double sc = this->operator()(v,v);
      if(sc < 0) std::cout << "Warning: scalar product not positive definite! Taking absolute value" << std::endl;
      return sqrt(std::fabs(sc));
    }
  };

  ///Representation of an error estimate, i.e. the output of an error estimator
  class AbstractErrorEstimate
  {
  public:
    virtual double absoluteError() const = 0;
    virtual ~AbstractErrorEstimate() {}
  };

  /// Representation of an error estimator
  class AbstractErrorEstimator
  {
  public:
    virtual std::unique_ptr<AbstractErrorEstimate> createEstimate(AbstractFunctionSpaceElement const& correction,
        AbstractLinearization const& lin) const = 0;
    virtual ~AbstractErrorEstimator() {}
  };

  ///Representation of an adaptive grid and a simple set of operations thereon
  class AbstractAdaptiveGrid
  {
  public:
    /// Change grid according to marked elements
    virtual void adapt() = 0;
    /// Mark elements, using an error estimate that hold error indicators
    virtual void mark(AbstractErrorEstimate&, double portion) = 0;
    /// Remove all marks in the grid
    virtual void flushMarks() = 0;
    /// Number of patches
    virtual int size() = 0;
    virtual ~AbstractAdaptiveGrid() {}

    virtual int getNMarked() = 0;

    /// Inform others that the grid will change
    mutable boost::signals2::signal<void()> gridWillChange;
  };

  class AbstractHierarchicalErrorEstimator
  {
  public:
    virtual ~AbstractHierarchicalErrorEstimator(){}

    virtual void operator()(AbstractLinearization const& lin, AbstractFunctionSpaceElement const& x, AbstractFunctionSpaceElement const& dx, int, AbstractFunctionSpaceElement const& rhs) = 0;

    virtual void refineGrid() = 0;

    virtual double estimatedAbsoluteError() const = 0;

    virtual size_t gridSize() const = 0;

//    virtual void setNorm(AbstractNorm const&) = 0;
  };

  class ContinuousScalarFunction
  {
  public:
    virtual ~ContinuousScalarFunction() {}
    virtual double d0(std::vector<double> const&) const = 0;
  };

  class DifferentiableScalarFunction : public ContinuousScalarFunction
  {
  public:
    virtual ~DifferentiableScalarFunction() {}

    virtual std::vector<double> d1(std::vector<double> const&) const = 0;
  };

/** @}*/

}  // namespace Kaskade
#endif
