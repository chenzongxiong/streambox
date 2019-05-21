#ifndef OPT_INTERFACE_HH
#define OPT_INTERFACE_HH

#include <memory>
#include <vector>
#include "dune/common/fvector.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/matrix.hh"

namespace Kaskade
{

  class LagrangeLinearization;
  class AbstractFunctionSpaceElement;
  class AbstractLinearization;
  class AbstractErrorEstimate;

  /// Class that models the functionality of a (possibly inexact) linear solver
  class AbstractTangentialSpace
  {
  public:
    typedef Dune::FieldVector<double,1> field_type;

    virtual ~AbstractTangentialSpace() {}

    /** Solve with possibly multiple solutions, return value: number of computed solutions*/
    int basis(std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >& corrections,
              LagrangeLinearization& linearization,
	      AbstractFunctionSpaceElement const& normalStep,
	      double nu0, AbstractFunctionSpaceElement* residual=nullptr)
    {   return computeBasis(corrections, linearization, normalStep, nu0, residual); }

    /// Specify accuracy that should be achieved
    virtual void setRelativeAccuracy(double accuracy) = 0;

    /// The maximal number of solution vectors, returned by basis
    virtual int nSolutionVectors() const = 0;
    
    virtual AbstractFunctionSpaceElement& getCorrectRhs() = 0;

    /// Returns true, if some information on the norm is available
    virtual bool getNorms(Dune::Matrix<field_type>& M) const { return false; }

    virtual bool localConvergenceLikely() { return true; }

    virtual void regularize(bool) {}

    virtual bool regularizationEnabled() const { return false; }

    virtual void setEps(double) {}

    virtual void setLipschitzConstant(double) {}

  private:
    virtual int computeBasis(std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >& correction,
                        LagrangeLinearization& lin,
			AbstractFunctionSpaceElement const& normalStep, double nu0, AbstractFunctionSpaceElement* residual) = 0;
  };

  class AbstractNormalDirection
  {
  public:
    /// compute min 1/2 <dn,dn> s.t. c'(x_0)dn+c(x_0)=0
    /// compute Lagrangemultiplier for: min 1/2 <w,w>+f'(x_0)  s.t. c'(x_0)w=0
    /// performs factorization
    /// uses normal linearization at x_0
    void ordinaryAndAdjoint(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& linearization, AbstractFunctionSpaceElement* correctionResidual=nullptr, AbstractFunctionSpaceElement* adjointResidual=nullptr)
    {  return computeCorrectionAndAdjointCorrection(correction, adjointCorrection, linearization, correctionResidual, adjointResidual); }

    /// compute min 1/2 <dn,dn> s.t. c'(x_0)dn+c(x)=0
    /// reuses factorization from ordinary(...) or ordinaryAndAdjoint(...)
    /// linearization at x
    void simplified(AbstractFunctionSpaceElement& correction, AbstractLinearization const& linearization, AbstractFunctionSpaceElement *residual=nullptr) const
    {  computeSimplifiedCorrection(correction, linearization, residual);  }

    virtual void setRelativeAccuracy(double accuracy) {}

    virtual void setEps(double) {}

  private:
    virtual void computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& linearization, AbstractFunctionSpaceElement* correctionResidual, AbstractFunctionSpaceElement* adjointResidual) = 0;

    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction, 
					     AbstractLinearization const& lin, AbstractFunctionSpaceElement* residual) const = 0;

  };

  class SearchSpaceCreator;

  /// Representation of an error estimator
  class AbstractCompositeStepErrorEstimator
  {
  public:
    virtual std::unique_ptr<AbstractErrorEstimate> createEstimate(
								  AbstractFunctionSpaceElement& iterate,
								  SearchSpaceCreator const& searchSpace, 
								  std::vector<AbstractFunctionSpaceElement*> const& basisVectors, 
								  std::vector<double>const & coeff,
								  AbstractLinearization const& linN,
								  AbstractLinearization const& linT) const = 0;
    virtual ~AbstractCompositeStepErrorEstimator() {};
  };

} // namespace Kaskade
#endif
