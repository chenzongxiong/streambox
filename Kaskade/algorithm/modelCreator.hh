#ifndef MODELCREATOR_HH
#define MODELCREATOR_HH

#include <memory>
#include <vector>

#include <dune/grid/config.h>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>

namespace Kaskade
{
  class AbstractScalarProduct;
  class AbstractFunctionSpaceElement;
  class LagrangeLinearization;
  class QuadraticFunction;

  /// Creates a quadratic model, parametrized by nu with argument(vector) tau, such that
  /// \f$ q_{\nu}(\tau) = f'(x_0)\delta x+\frac{1}{2}L_{xx}(x_0,p_0)\delta x^2 \f$
  /// where \f$ \delta x = \nu \delta n + \tau \delta t \f$
  /// Uses the identity \f$ f'(x_0)\delta x = L_x(x_0,p_0)\delta t+\nu p_0 c(x_0) \f$
  struct QuadraticModelCreator
  {
    QuadraticModelCreator(AbstractFunctionSpaceElement const& normalStep, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >const& tangentialBasis, LagrangeLinearization const& lagrangeLinearization, double residualCorrection = 0);

    QuadraticFunction create(double nu) const;

  private:
    double p0c0, Lxxdndn;
    Dune::BlockVector<Dune::FieldVector<double,1>> Lxdt, Lxxdndt;
    Dune::Matrix<Dune::FieldVector<double,1>> Lxxdtdt;
  };

  struct NormModelCreator
  {
    /**
     * \param tangentialBasis basis of subspace of tangential space
     * \param sp ?
     */
    NormModelCreator(AbstractFunctionSpaceElement const& normalStep, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > const& tangentialBasis, AbstractScalarProduct const& sp);

    // Create scalar product in tangential subspace.
    QuadraticFunction create(double nu) const;

  private:
    double dndn;
    Dune::BlockVector<Dune::FieldVector<double,1>> dndt;
    Dune::Matrix<Dune::FieldVector<double,1>> dtdt;
  };
}

#endif // MODELCREATOR_HH
