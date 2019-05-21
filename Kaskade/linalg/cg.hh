#ifndef CG_HH
#define CG_HH

#include "cgImplementation.hh"

namespace Kaskade
{
  /**
   * Standard preconditioned conjugate gradient method.
   *
   * @tparam X type of domain space element
   * @tparam Xstar type of image space element
   * @tparam TimerPolicy enable timers for dual pairings, axpy-operations, matrix-vector products and preconditioners with struct MeasureTime
   */
  template <class X, class Xstar, class TimerPolicy = DoNotMeasureTime>
  using CG = CGBase<X,Xstar,CGImplementationType::STANDARD,TimerPolicy>;

  /**
   * Truncated preconditioned conjugate gradient method for nonconvex problems. Stops iteration if a direction of negative curvature is encountered.
   *
   * @tparam X type of domain space element
   * @tparam Xstar type of image space element
   * @tparam TimerPolicy enable timers for dual pairings, axpy-operations, matrix-vector products and preconditioners with struct MeasureTime
   */
  template <class X, class Xstar, class TimerPolicy = DoNotMeasureTime>
  using TCG = CGBase<X,Xstar,CGImplementationType::TRUNCATED,TimerPolicy>;

  /**
   * Regularized preconditioned conjugate gradient method for nonconvex problems. Denote the used operator by \f$A\f$ and the preconditioner by \f$P\f$.
   * Then if a direction of negative curvature is encountered \f$A\f$ is implicitly replaced by the regularization $\f$A+\thetaP\f$. Then the IterateType::CG method is
   * restarted for the regularized problem. The necessary quantities are available during the standard cg implementation, thus the costs for computing the
   * regularization are neglishible.
   *
   * @tparam X type of domain space element
   * @tparam Xstar type of image space element
   * @tparam TimerPolicy enable timers for dual pairings, axpy-operations, matrix-vector products and preconditioners with struct MeasureTime
   */
  template <class X, class Xstar, class TimerPolicy = DoNotMeasureTime>
  using RCG = CGBase<X,Xstar,CGImplementationType::REGULARIZED,TimerPolicy>;

  /**
   * Hybrid preconditioned conjugate gradient method for nonconvex problems, mixing the truncated with the regularized conjugate gradient method. If a direction
   * of negative curvature is encounted and the termination criterion indicates sufficient decrease in the used norm the iteration is stopped. Else, denoting
   * the used operator by \f$A\f$ and the preconditioner by \f$P\f, \f$A\f$ is implicitly replaced by the regularization $\f$A+\thetaP\f$. Then the IterateType::CG method is
   * restarted for the regularized problem. The necessary quantities are available during the standard cg implementation, thus the costs for computing the
   * regularization are neglishible.
   *
   * @tparam X type of domain space element
   * @tparam Xstar type of image space element
   * @tparam TimerPolicy enable timers for dual pairings, axpy-operations, matrix-vector products and preconditioners with struct MeasureTime
   */
  template <class X, class Xstar, class TimerPolicy = DoNotMeasureTime>
  using HCG = CGBase<X,Xstar,CGImplementationType::HYBRID,TimerPolicy>;
}

#endif // CG_HH
