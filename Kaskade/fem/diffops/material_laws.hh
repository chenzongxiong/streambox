/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * material_laws.hh
 * Collection of hyperelastic material laws in weak formulation.
 * Currently implemented: Hooke, St.Venant-Kirchhoff, compressible Mooney-Rivlin
 *  Created on: Oct 10, 2012
 *      Author: Lars Lubkoll
 */

#ifndef MATERIAL_LAWS_HH_
#define MATERIAL_LAWS_HH_

#include <string>
#include <type_traits>

#include <dune/common/fmatrix.hh>

#include "fem/fixdune.hh"
#include "linalg/determinant.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/elementary_functions.hh"
#include "utilities/linalg/scalarproducts.hh"
#include "utilities/straintensors.hh"

namespace Kaskade{
  namespace HyperelasticMaterialLaw
  {
    /// Implementation details are hidden here
    namespace ImplementationDetail{



      template <typename CompressibleFunction>
      struct GammaSignWrapper{ static constexpr int value = -1; };

      template <>
      struct GammaSignWrapper<LN>{ static constexpr int value = -1; };

      template <>
      struct GammaSignWrapper<PenaltyLN> { static constexpr int value = -1; };

      /// Gamma function for compressible material laws
      /**
       * \param Scalar scalar type
       * \param dim dimension
       * \param monomial_exponent exponent of the first part (>0)
       * \param CompressibleFunction type of the compressible part (NaturalLogarithm, Monomial<k>, with k<0)
       */
      template <typename Scalar, int dim, typename CompressibleFunction>
      class Gamma{
      public:
        typedef Gamma<Scalar,dim,CompressibleFunction> This;
        typedef Monomial<2> Monom;
        typedef Dune::FieldVector<Scalar,dim> Vector;
        typedef Dune::FieldMatrix<Scalar,dim,dim> Matrix;
        enum{ SignOfCompressibleFunction = ImplementationDetail::GammaSignWrapper<CompressibleFunction>::value };

        /// Constructor
        /**
         * \param strain_ strain tensor
         * \param lambda first Lam\'e constant
         * \param mu second Lam\'e constant
         */
        Gamma(Matrix &deformationGradient, Scalar const weight1, Scalar const weight2, Scalar const det_min_val_=1.0e-9) :
          det(deformationGradient), det_value(det.d0()), det_min_val(det_min_val_), monom(weight1), cFunction( SignOfCompressibleFunction * weight2)
        {}

        /// Function value
        /**
         * \return function value
         */
        Scalar d0() const
        {
          det_value = det.d0();
          if(det_value < det_min_val) det_value = fabs(det_min_val);
          return monom.d0(det.d0())+cFunction.d0(det_value);
        }

        /// First derivative
        /**
         * \param dA increment
         * \return first derivative
         */
        Scalar d1(Matrix const& dA) const
        {
          det_value = det.d0();
          if(det_value < det_min_val) det_value = fabs(det_min_val);
          return monom.d1(det.d0())*det.d1(dA)+cFunction.d1(det_value)*det.d1(dA);
        }

        /// Second derivative
        /**
         * \param dA increment
         * \param dB increment
         * \return second derivative
         */
        Scalar d2(Matrix const& dA, Matrix const& dB) const
        {
          det_value = det.d0();
          if(det_value < det_min_val) det_value = fabs(det_min_val);
          return monom.d2(det.d0())*det.d1(dA)*det.d1(dB) + monom.d1(det.d0())*det.d2(dA,dB) +
              cFunction.d2(det_value)*det.d1(dA)*det.d1(dB) + cFunction.d1(det_value)*det.d2(dA,dB);
        }

        /// Third derivative
        /**
         * \param du1 increment (first derivative)
         * \param du2 increment (second derivative)
         * \param du3 increment (third derivative)
         * \return third derivative
         */
        Scalar d3(Matrix const& du1, Matrix const& du2, Matrix const& du3) const
        {
          det_value = det.d0();
          if(det_value < det_min_val) det_value = fabs(det_min_val);
          return /*monom.d3(det.d0())*det.d1(du1)*det.d1(du2)*det.d1(du3) +*/
          monom.d2(det.d0())*det.d2(du1,du3)*det.d1(du2) +
          monom.d2(det.d0())*det.d1(du1)*det.d2(du2,du3) +
          monom.d2(det.d0())*det.d2(du1,du2)*det.d1(du3) +
          monom.d1(det.d0())*det.d3(du1,du2,du3) +
          cFunction.d3(det_value)*det.d1(du1)*det.d1(du2)*det.d1(du3) +
          cFunction.d2(det_value)*det.d2(du1,du3)*det.d1(du2) +
          cFunction.d2(det_value)*det.d1(du1)*det.d2(du2,du3) +
          cFunction.d2(det_value)*det.d2(du1,du2)*det.d1(du3) +
          cFunction.d1(det_value)*det.d3(du1,du2,du3);
        }

        /// Set weights
        /**
         * \param weight1 weight for the monomial part
         * \param weight2 weight for the "compressible" penalty part
         */
        void setWeights(Scalar const weight1, Scalar const weight2)
        {
          monom.setA(weight1);    cFunction.setA(SignOfCompressibleFunction * weight2);
        }

      private:
        Determinant<dim> det;
        mutable Scalar det_value;
        Scalar const det_min_val;
        Monom monom;
        CompressibleFunction cFunction;
      };
    }

  }


}
#endif /* MATERIAL_LAWS_HH_ */
