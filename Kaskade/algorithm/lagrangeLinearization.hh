#ifndef LAGRANGELINEARIZATION_HH
#define LAGRANGELINEARIZATION_HH

#include "linalg/triplet.hh"

namespace Kaskade
{
  class AbstractLinearization;
  class AbstractFunctionSpaceElement;

  class LagrangeLinearization : public AbstractLinearization
  {
  public:
    LagrangeLinearization(AbstractLinearization* Nlin_, AbstractLinearization* Tlin_, AbstractFunctionSpaceElement const& adjointCorrection, int stateId_=1, int controlId_=0, bool hasTlin_uu_=false);

    virtual ~LagrangeLinearization();

    /// Evaluate f(origin)
    virtual double eval() const;

    virtual double evalL1norm() const;

    /// Evaluate f'(origin)(.) + <Lag. Multiplier, c'(origin)(.)  > this is dual element
    virtual void evald(AbstractFunctionSpaceElement &g, int rbegin=0, int rend=-1) const;

    /// Evaluate scaled hessian (of Lagrangian) times second argument
    virtual void d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const;

    virtual void d2taxpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const;

    /// Get point of linearization
    virtual AbstractFunctionSpaceElement const& getOrigin() const;

    /// Precompute data
    // be careful, are Lagrange-Multipliers (needed for Tlin) up to date?
    virtual void precompute();

    virtual void flush();

  //    virtual void touch();
  //
  //    virtual void connectToSignalForFlush(boost::signals2::signal0<void>& sig);

    AbstractLinearization& getTangentialLinearization();
    AbstractLinearization const& getTangentialLinearization() const;

    AbstractLinearization& getNormalLinearization();
    AbstractLinearization const& getNormalLinearization() const;

    virtual void getMatrixBlocks(MatrixAsTriplet<double>& mat, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const;

  private:
    AbstractLinearization* Nlin;
    AbstractLinearization* Tlin;
    std::unique_ptr<AbstractFunctionSpaceElement> origin;
    std::unique_ptr<AbstractFunctionSpaceElement> lagrangeCorrection;
    int stateId, controlId;
    bool hasTlin_uu = false;
  //    bool existsLagrangian;
  };
}

#endif // LAGRANGELINEARIZATION_HH
