#ifndef OPTIMIZATION_HH
#define OPTIMIZATION_HH

#include <memory>
#include <vector>
#include "dune/grid/config.h"
#include "opt_interface.hh"
#include "algorithm_base.hh"
#include "abstract_interface.hh"
#include "newton_bridge.hh"
#include "lipschitzConstants.hh"
#include "linalg/triplet.hh"
#include "lagrangeLinearization.hh"

namespace Kaskade
{
  class QuadraticFunction;
  class RegularizedQuadraticFunction;

  template <class X>
  struct EnergyScalarProduct : public AbstractScalarProduct
  {
    virtual ~EnergyScalarProduct(){}

    virtual void setOrigin(AbstractLinearization const& linearization)
    {
      SparseLinearSystem const& m = dynamic_cast<SparseLinearSystem const&>(linearization);
      m.getMatrixBlocks(Hu,0,1,0,1);
      m.getMatrixBlocks(Hy,1,2,1,2);
    }

    virtual double operator()(AbstractFunctionSpaceElement const& v, AbstractFunctionSpaceElement const& w) const
    {
      X const& x = Bridge::getImpl<X>(v);
      X const& y = Bridge::getImpl<X>(w);

      using namespace boost::fusion;

      int sizeu=x.descriptions.degreesOfFreedom(0,1);
      int sizey=x.descriptions.degreesOfFreedom(1,2);

      double sumu(0.0);
      {
        std::vector<double> Hx(sizeu,0.0), argument1(sizeu,0.0), argument2(sizeu,0.0);
        vectorToSequence(at_c<0>(x.data).coefficients(),argument1.begin());
        vectorToSequence(at_c<0>(y.data).coefficients(),argument2.begin());
        Hu.axpy(Hx, argument1);
        for(int i=0; i<sizeu; ++i) sumu += Hx[i]*argument2[i];
      }


      double sumy(0.0);
      {
        std::vector<double> Hx(sizey,0.0), argument1(sizey,0.0),argument2(sizey,0.0);
        vectorToSequence(at_c<1>(x.data).coefficients(),argument1.begin());
        vectorToSequence(at_c<1>(y.data).coefficients(),argument2.begin());
        Hy.axpy(Hx, argument1);
        for(int i=0; i<sizey; ++i) sumy += Hx[i]*argument2[i];
      }

      return sumy + sumu;
    }

  private:
    MatrixAsTriplet<double> Hy, Hu;
  };


  class OptimizationParameters
  {
  public:
    OptimizationParameters(double desiredAccuracy_, int maxSteps_);
    
    OptimizationParameters(OptimizationParameters const&) = default;
    OptimizationParameters& operator=(OptimizationParameters const&) = default;

    // Parameters with predefined values that can be changed by client
    double desiredAccuracy = 1e-12;
    double desiredEstimatorAccuracy = 1e-12;
    double minimalDampingFactor = 1e-12;
    double etaMin = 0.5;
    double etaLock = 0.99;
    double minimalAccuracy = 1e-3;
    size_t maxSteps = 1000;
    size_t maxGridSize = 10000;
    double requiredRelativeError = 0.1;

    void setThetaAim(double theta);
    void setThetaNormal(double theta);
    void setThetaMax(double theta);
    void setThetas(double thetaNormal, double thetaAim, double thetaMax);
    void setEps(double eps_);

    double getThetaAim() const;
    double getThetaNormal() const;
    double getThetaMax() const;
    double getEps() const;
    double getSqrtEps() const;
    double getThirdSqrtEps() const;

  private:
    void ensureAdmissibleThetas();
  
    double ThetaAim = 0.25;     // required contraction 
    double ThetaNormal = 0.1;   // required contraction
    double ThetaMax = 0.5;     // 
    double eps = 1e-15;
    double sqrtEps = sqrt(1e-15);
    double thirdSqrtEps = 1e-5;
  };


  class PrimalChart : public AbstractChart
  {
  public:
    void addPerturbation(AbstractFunctionSpaceElement& newIterate,
        AbstractFunctionSpaceElement const& perturbation,
        AbstractLinearization const& lin,
        std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > basis = std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >()) const
    {
      newIterate = lin.getOrigin();
      newIterate.axpy_role(1.0,perturbation,"primal");
    }

    std::unique_ptr<AbstractChart> clone() const { return std::unique_ptr<PrimalChart>(new PrimalChart(*this)); }
  };

  class Optimization : public Algorithm
  {
  public:
    /// Create Newton's Method, providing a solver, a norm and algorithmic parameters
    Optimization(AbstractNormalDirection& normalSolver,
        AbstractTangentialSpace& tangentSpace_,
        AbstractScalarProduct& nL,
        OptimizationParameters const& p_,
        double omegaCinit=1,
        double omegaLinit=1,
        int verbose=1);

    Optimization(AbstractNormalDirection& normalSolver,
        AbstractScalarProduct& nL,
        AbstractScalarProduct& nC,
        AbstractChart const& chart_,
        OptimizationParameters const& p_,
        int verbose=1);

    /// Create Newton's Method, providing a solver, a norm and algorithmic parameters
    Optimization(AbstractScalarProduct& nL,
        AbstractScalarProduct& nC,
        AbstractChart const& chart_,
        OptimizationParameters const& p_,
        AbstractCompositeStepErrorEstimator* errorEstimator_,
        int verbose=1);

    virtual ~Optimization();

    /// Solve the system f=0 with starting value x. On (successful) exit, the solution is x, otherwise it is left unmodified.
    void solve(AbstractFunctional& fN, AbstractFunctional& fT, AbstractFunctionSpaceElement& x);

    /// Solve the system f=0 with starting value x. On (successful) exit, the solution is x, otherwise it is left unmodified.
    /**
     * Do adaptive mesh refinement as soon as local convergence sets in
     */
    void solve(AbstractFunctional& fN, AbstractFunctional& fT, AbstractFunctionSpaceElement& x, AbstractHierarchicalErrorEstimator& hbErrorEstimator_);

    void solve(AbstractFunctional& f, AbstractFunctionSpaceElement& x, AbstractHierarchicalErrorEstimator& errorEstimator);

private:
    void solve(AbstractFunctional& fN, AbstractFunctional* fT, AbstractFunctionSpaceElement& x);
    int runAlgorithm();

    std::pair<std::unique_ptr<AbstractFunctionSpaceElement>,std::unique_ptr<AbstractFunctionSpaceElement> > computeNormalStep(AbstractFunctionSpaceElement* normalStepResidual, AbstractFunctionSpaceElement* adjointResidual);
//    std::unique_ptr<AbstractFunctionSpaceElement> computeAdjoint();
    std::unique_ptr<AbstractFunctionSpaceElement> computeSimplifiedNormalStep(AbstractLinearization const& lin_xplus, AbstractFunctionSpaceElement* normalStepResidual = nullptr);
    std::unique_ptr<AbstractFunctionSpaceElement> computeSecondOrderCorrection(AbstractLinearization const& lin_xplus, AbstractFunctionSpaceElement const& normalStep, double nu);
    std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > computeTangentialStep(LagrangeLinearization&, AbstractFunctionSpaceElement const& normalStep, double nu, std::vector<double> const& tau);

    void updateConstraintD1LipschitzConstant(double normSimplifiedNormal, double normCAtCorr);
    double updateLagrangianD2LipschitzConstant(AbstractLinearization const& lin_x0, AbstractFunctionSpaceElement const& secondOrderCorrected, double normCAtCorr, double quadraticModelAtCorr, RegularizedQuadraticFunction const& cubic, double nu, std::vector<double> const& tau, AbstractFunctionSpaceElement const& normalStep, AbstractFunctionSpaceElement const& trialIterate, AbstractFunctionSpaceElement const& sNormalStep, AbstractFunctionSpaceElement const& correction, double Lxdn_res);

    double updateNormalStepDampingFactor(double normNormal) const;
    void updateTangentialDampingFactor(double nu, double normNormal, double normTangential, RegularizedQuadraticFunction const& cubic, std::vector<double>& tau) const;
    bool adaptiveMeshRefinement(LagrangeLinearization& lin, double nu, std::vector<double> tau, AbstractFunctionSpaceElement const& correction);

    AcceptanceTest acceptanceTest(double eta, double nu, std::vector<double> const& tau, double normOfCorrection) const;
    bool regularityTest(double nu, std::vector<double> const& tau, bool reliableQuadraticModel) const;
    Convergence convergenceTest(double nu, std::vector<double> const& tau, double normOfCorrection) const;

    void terminationMessage(int flag);
    void printNormalStep(double normNormal, double nu) const;
    void printTangentialStep(double normTangential, double tau) const;
    
    bool noDamping(double d) const;
    bool noDamping(std::vector<double> const& tau) const;

    std::unique_ptr<AbstractFunctionSpaceElement> createCorrection(double nu, AbstractFunctionSpaceElement const& normalStep, std::vector<double> const& tau, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > const& tangentialBasis) const;
    void addCorrection(AbstractFunctionSpaceElement& v, double nu, AbstractFunctionSpaceElement const& normalStep, std::vector<double> const& tau, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > const& tangentialBasis) const;

    AbstractFunctional* functionalN = nullptr, * functionalT = nullptr;
    AbstractScalarProduct &normL, &normC;
    std::unique_ptr<AbstractChart const> chart = std::unique_ptr<AbstractChart const>(new PrimalChart());
    OptimizationParameters p;

    std::shared_ptr<AbstractLinearization> normalLinearization = nullptr, tangentialLinearization = nullptr;
    AbstractCompositeStepErrorEstimator* errorEstimator = nullptr;
    AbstractHierarchicalErrorEstimator* hbErrorEstimator = nullptr;
    AbstractNormalDirection* normalDirection = nullptr;
    AbstractTangentialSpace* tangentSpace = nullptr;
    int verbose = 0;
    std::string csPre = std::string("COMPOSITE STEP: ");
    std::unique_ptr<AbstractFunctionSpaceElement> iterate;

    // algorithmic parameters
    double dampingFactorTolerance = 1e-2;
    double normalStepComputationTime = 0, tangentialStepComputationTime = 0;
    
    // internal variables
    double normOfLastCorrection = -1., normOfLastCorrection_Undamped = -1., normOfIterate = 0;

    LagrangianD2LipschitzConstant L;
    ConstraintD1LipschitzConstant C;
  public:
    size_t step = 0;
  };

}  // namespace Kaskade
#endif
