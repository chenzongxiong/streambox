/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef Cygwin
// availability of std::to_string
#define _GLIBCXX_USE_C99 1
#endif

#include <cassert>
#include <cmath>
#include <numeric>
#include <dune/grid/config.h>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "linalg/apcg.hh"

namespace Kaskade
{

  template <class R>
  PCGEnergyErrorTerminationCriterion<R>::PCGEnergyErrorTerminationCriterion(Real tol_, int maxit_)
  : tol(tol_), maxit(maxit_), lookah(6)
  {}

  template <class R>
  void PCGEnergyErrorTerminationCriterion<R>::clear() { gammas2.clear(); }

  template <class R>
  PCGEnergyErrorTerminationCriterion<R>& PCGEnergyErrorTerminationCriterion<R>::tolerance(Real tol_) 
  { 
    assert(tol>=0); 
    tol = tol_; 
    return *this; 
  }

  template <class R>
  PCGEnergyErrorTerminationCriterion<R>& PCGEnergyErrorTerminationCriterion<R>::lookahead(int lah) 
  { 
    assert(lah>=0); 
    lookah = lah; 
    return *this;
  }


  template <class R>
  void PCGEnergyErrorTerminationCriterion<R>::step(Real gamma2) { gammas2.push_back(gamma2); }

  template <class R>
  PCGEnergyErrorTerminationCriterion<R>::operator bool() const {
    return (gammas2.size()>=5 && error() <= tol) || gammas2.size()>=maxit;
  }

  template <class R>
  R PCGEnergyErrorTerminationCriterion<R>::error() const {


    //
    // estimate parameters a2 and b2 of smoothing convergence model.
    // First build up the normal equations
    //
    Dune::FieldMatrix<R,2,2> AtA(0);
    Dune::FieldVector<R,2> f(0);

    int const n = gammas2.size();

    // Check whether there is enough data for robust lookahead error estimation
    if (n<lookah)
      return std::numeric_limits<R>::max();

    // energy error estimator [e_k]_n^2
    R eps2 = 0;

    // step backwards to do partial summation on the fly
    for (int k=n-1; k>=0; --k) {
      eps2 += gammas2[k];

      R wk = 1.0/(k+1); // discount weight
      R pk = 1.0/(std::sqrt(static_cast<R>(k))+1);
      pk *= pk;

      AtA[0][0] += wk*wk*pk*pk;
      AtA[0][1] -= wk*wk*pk;
      AtA[1][1] += wk*wk;

      f[0] += eps2*wk*wk*pk;
      f[1] -= eps2*wk*wk;
    }
    AtA[1][0] = AtA[0][1];

    // solve normal equations
    Dune::FieldVector<R,2> ab2(0);
    AtA.solve(ab2,f);

    R const a2 = ab2[0];

    //
    // estimate linear contraction rate
    //
    R const eps0n2 = eps2; // [eps_0]_n^2


    eps2 = 0;
    R theta = 0;
    int kThetaMax = -1;
    // step backwards to do partial summation on the fly
    for (int k=n-1; k>=1; --k) {
      eps2 += gammas2[k];

      R theta_k = 0;
      for (int i=0; i<9; ++i) // fixed point iteration
        theta_k = std::pow( (eps2+pow(theta_k,2*n)*(eps0n2-eps2))/eps0n2, 0.5/k );
      if (theta_k>=theta) {
        theta = theta_k;
        kThetaMax = k;
      }
    }

    //
    // compute crossover of both models
    //
    int k = 1;
    for ( ; k<=n; ++k)
      if (pow(theta,2*k) <= 1 + a2 * (1./pow(std::sqrt(k)+1,2)-1)*(1-pow(theta,2*n))/eps0n2 )
        break;

    // Check which convergence phase is active
    R epsn2;
    if (kThetaMax>0.8*n || kThetaMax>=n-6) {
      // first phase: smoothing
      epsn2 = a2/(std::sqrt(n)+1);
      epsn2 *= epsn2;
    } else {
      // second phase: linear convergence
      epsn2 = pow(theta,2*n)*eps0n2/(1-pow(theta,2*n));
    }

    return std::sqrt(std::accumulate(gammas2.begin()+(n-lookah),gammas2.end(),epsn2));
  }


  // Explicit instantiation for float and double

  /**
   * The first explicit instantation (for float) does not work with Dune2.1 as there is a bug in Dune::DenseMatrix<float>
   * yielding calls to std::max(double,float) and std::max(float,double), both not existent.
   */
  template class PCGEnergyErrorTerminationCriterion<float>;
  template class PCGEnergyErrorTerminationCriterion<double>;

} // namespace Kaskade7
