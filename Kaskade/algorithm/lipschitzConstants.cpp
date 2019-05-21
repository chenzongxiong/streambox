#include <dune/grid/config.h>
#include <dune/common/fvector.hh>
#include <dune/istl/matrix.hh>

#include "lipschitzConstants.hh"

namespace Kaskade
{
  LagrangianD2LipschitzConstant::LagrangianD2LipschitzConstant(int verbose_, double initialOmega)
    : omega(initialOmega), oldOmega(initialOmega), verbose(verbose_), L_xx(0.)
  {}

  void LagrangianD2LipschitzConstant::setFirstOrder(double norm_dx_, double thetaC_, double modelError_)
  {
    norm_dx = norm_dx_;
    thetaC = thetaC_;
    modelError = modelError_;
    secondOrderUsed = false;
  }

  void LagrangianD2LipschitzConstant::setSecondOrder(double secondOrderEstimate_)
  {
    secondOrderEstimate = secondOrderEstimate_;
    secondOrderUsed = true;
  }

  void LagrangianD2LipschitzConstant::update(bool doLock)
  {
//    lock = doLock;
    double thirdOrderOmega = modelError/(norm_dx*norm_dx*norm_dx);
    double secondOrderOmega(thirdOrderOmega);
    if(secondOrderUsed) secondOrderOmega = secondOrderEstimate/(norm_dx*norm_dx);

    oldOmega = omega;
    omega = std::min(secondOrderOmega,thirdOrderOmega)*6.0;
  /*  if( omega < 0 ) omega = (lock || doLock) ? oldOmega : 0.1*oldOmega;

    if( (lock) && omega > oldOmega ) omega = oldOmega;
*/
    if( !doLock )
    {
      omega = std::max(oldOmega*lowFactor,omega);
      omega = std::min(oldOmega*highFactor,omega);
    }
    else
    {
      if(omega > oldOmega)
      {
        std::cout << "LIPSCHITZ CONSTANT Lagrangian d2 LOCKED." << std::endl;
        omega = oldOmega;
      }
      else omega = std::max(omega,oldOmega*lowFactor);
    }

    if( verbose > 0 ) std::cout << "LIPSCHITZ CONSTANT Lagrangian d2: " << oldOmega << " -> " << omega << std::endl;
  }


  bool LagrangianD2LipschitzConstant::highRoundOffError(bool hasNormalDirection, double epsilon)
  {
    bool highRoundOffError(false);

    if(secondOrderUsed && norm_dx * norm_dx < epsilon) highRoundOffError = true;
    if(!secondOrderUsed && norm_dx * norm_dx * norm_dx < epsilon) highRoundOffError = true;
    if( highRoundOffError && verbose > 0) std::cout << "LIPSCHITZ CONSTANT: Possibly high round-off errors in computation of Lipschitz constant for L_xx!" << std::endl;

    if( verbose > 1 )
    {
      std::cout << "LIPSCHITZ CONSTANT: constraint contraction: " << thetaC << std::endl;
      std::cout << "LIPSCHITZ CONSTANT: norm_dx: " << norm_dx << ", cubed: " << norm_dx*norm_dx*norm_dx<< std::endl;
      std::cout << "LIPSCHITZ CONSTANT: Model error (me): " << modelError << ", Second order estimate (soe): " << secondOrderEstimate << ", me/soe: " << modelError/secondOrderEstimate << std::endl;
      std::cout << "LIPSCHITZ CONSTANT: Lower bound (eps): " << epsilon << ", eps/me: " << epsilon/modelError << std::endl;
    }

//    lock = highRoundOffError;
    return highRoundOffError;
  }

  void LagrangianD2LipschitzConstant::setL_xx(double m) { L_xx = m; }

  void LagrangianD2LipschitzConstant::print() const
  {
    if(secondOrderUsed) std::cout << "LIPSCHITZ CONSTANT: omegaL second order:" << secondOrderEstimate/(norm_dx*norm_dx*norm_dx)*6.0 << " first order:" << modelError/(norm_dx*norm_dx*norm_dx)*6.0 << std::endl;
    else std::cout << "LIPSCHITZ CONSTANT: omegaL:" << modelError/(norm_dx*norm_dx*norm_dx)*6.0 << std::endl;
  }

  bool LagrangianD2LipschitzConstant::isPositiveDefinite() const { return L_xx > 0; }

  /**********************************************************************************************************/
  // Lipschitz constant for constraint derivative.
  /**********************************************************************************************************/
  ConstraintD1LipschitzConstant::ConstraintD1LipschitzConstant(int verbose_, double initialOmega)
    : omega(initialOmega), oldOmega(initialOmega), verbose(verbose_)
  {}

  void ConstraintD1LipschitzConstant::update(double newOmega, bool lock)
  {
    oldOmega = omega;
    if( !lock || newOmega <= omega ) omega = newOmega;

    omega = std::max(oldOmega*lowFactor,omega);
    omega = std::min(oldOmega*highFactor,omega);

    if( verbose > 0 ) std::cout << "LIPSCHITZ CONSTANT Constrained d1: " << oldOmega << " -> " << omega << std::endl;
  }
}
