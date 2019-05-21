/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BSPLINES_HH
#define BSPLINES_HH

#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

namespace Kaskade {

  /**
   * \brief computes the differentiation matrix for bspline interpolation of order p on t 
   * 
   * B-Splines of order p define an interpolation \f$ u(t) \f$ for all \f$ t \f$ between 
   * the \f$ n+1 \f$ timepoints \f$ t_i \f$ with 
   * given values \f$ u(t_i) \f$. The derivative's values \f$ u'(t_i) \f$ at the \f$ m+1 \f$ 
   * timepoints \f$ s_i \f$ can then be computed as a linear combination
   * of the values \f$ u(t_i) \f$: \f[ u'(s_i) = \sum_{k=0}^n d_{ik} u(t_k) \f]
   * 
   * For \f$ p=n \f$, the interpolation is classical polynomial interpolation.
   * 
   * \param[in] t the time points from which to interpolate. Entries have to be strictly increasing.
   * \param[in] s the time points from which to interpolate. Entries have to be strictly increasing.
   * \param[in] p the order of B-spline interpolation (nonnegative)
   * \param[out] d the \f$ m+1 \times n+1 \f$ differentiation matrix
   */
  template <class Scalar>
  void bsplineDifferentiationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, int p, 
				    Dune::DynamicMatrix<Scalar>& d);
  
  /**
   * \brief computes the evaluation on s matrix for bspline interpolation of order p on t 
   * 
   * B-Splines of order p define an interpolation \f$ u(t) \f$ for all \f$ t \f$ between 
   * the \f$ n+1 \f$ timepoints \f$ t_i \f$ with 
   * given values \f$ u(t_i) \f$. The values \f$ u(s_i) \f$ can then be computed as a linear combination
   * of the values \f$ u(t_i) \f$: \f[ u(s_i) = \sum_{k=0}^n e_{ik} u(t_k) \f]
   * 
   * For \f$ p=n \f$, the interpolation is classical polynomial interpolation.
   * 
   * \param[in] t the time points from which to interpolate. Entries have to be strictly increasing.
   * \param[in] s the time points from which to interpolate. Entries have to be strictly increasing.
   * \param[in] p the order of B-spline interpolation (nonnegative)
   * \param[out] e the \f$ m+1 \times n+1 \f$ evaluation matrix
   */
  template <class Scalar>
  void bsplineEvaluationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, int p,
			       Dune::DynamicMatrix<Scalar>& e);
}

#endif 

