#ifndef LIPSCHITZCONSTANTS_HH
#define LIPSCHITZCONSTANTS_HH

namespace Dune
{
  template <class,int> class FieldVector;
  template <class,class> class Matrix;
}

namespace Kaskade
{
  struct LagrangianD2LipschitzConstant
  {
    explicit LagrangianD2LipschitzConstant(int verbose_, double initialOmega = 1e-3);

    LagrangianD2LipschitzConstant(LagrangianD2LipschitzConstant const&) = default;
    LagrangianD2LipschitzConstant& operator=(LagrangianD2LipschitzConstant const&) = default;

    void setFirstOrder(double norm_dx_, double thetaC_, double modelError_);

    void setSecondOrder(double secondOrderEstimate_=0.0);

    /// compute robust value for omegaL
    void update(bool doLock);

    /// decide whether their might be high round off error
    bool highRoundOffError(bool hasNormalDirection, double epsilon);

    void setL_xx(double m);

    void print() const;

    bool isPositiveDefinite() const;

    double omega = 1e-3; // new robust value for omega
    double oldOmega = 1e-3; // last value for omega
    
  private:
    double lowFactor = 1e-4;  // factor of maximal decrease of omega between two steps
    double highFactor = 1e3; // factor of maximal increase of omega between two steps
    double modelError = 0;    // Model Error: difference between second order model and function
    double norm_dx = 1;    // norm of correction
    double thetaC = 0;     // contraction for normal step
    bool lock = false;
    int verbose = 0;
    double secondOrderEstimate = 0;
    bool secondOrderUsed = false; // false: third order estimate (round-off errors may be large) true: second order estimate
    double L_xx=1;
    bool firstEstimate = true;
  };

  struct ConstraintD1LipschitzConstant
  {
    explicit ConstraintD1LipschitzConstant(int verbose_, double initialOmega = 1e-3);

    void update(double newOmega, bool lock = false);

    double omega = 1e-3, oldOmega = 1e-3, theta = 0;

  private:
    double lowFactor = 1e-3, highFactor = 1e3;
    int verbose = 0;
    bool firstEstimate = true;
  };
}

#endif // LIPSCHITZCONSTANTS_HH
