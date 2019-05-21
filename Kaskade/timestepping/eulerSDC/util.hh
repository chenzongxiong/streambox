/*
 * util.hh
 *
 *  Created on: May 3, 2015
 *      Author: sunayana_ghosh
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//In this file we implement some of the utility functions necessary for
//the classes WorkModel, ErrorEstimate and the main integrate function of the class EulerSDC.

#ifndef KASKADE_TIMESTEPPING_EULERSDC_UTIL_HH_
#define KASKADE_TIMESTEPPING_EULERSDC_UTIL_HH_

#include "norm.hh"

#include <cmath>          //std::abs
#include <functional>     //std::function
#include <vector>

#include <boost/timer/timer.hpp>

#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

#include "fem/fixdune.hh"
#include "timestepping/sdc.hh"

namespace Kaskade {
  /**
   * \ingroup eulerSDC
   * \brief Some utility functions.
   */
//=================================================================================================
//  Utility function to find norm of a difference of vector of two Vectors
//=================================================================================================

/**
 * \ingroup eulerSDC
 * \brief Function to compute the norm of the difference of two vector of Vectors.
 *
 * @param y0  : A vector of Vectors, usually of type std::vector<Vector>
 * @param y1  : A vector of Vectors, usually of type std::vector<Vector>
 * @param norm  : A type of norm, in this case an object of type Kaskade::Norm
 * @return An object of type Vector::value_type is returned, it is usually a double.
 */
template<class Vector, class Norm>
typename Vector::value_type normVecDiff(std::vector<Vector> const& y0, std::vector<Vector> const& y1, Norm norm);


//=================================================================================================
//  Implementation for single sdc iteration step based on explicit Euler method.
//=================================================================================================
 /**
  * \ingroup eulerSDC
  * This function performs one spectral defect correction (SDC) iteration for a system of ODE's of the form:
  * \f$ y'(t) = f(y(t)) \f$, \f$ y(0) = y_0\f$ and \f$ t \in [0, T] \f$,
  * using the explicit Euler method on a given time grid. Here SDC is interpreted as fixed point iteration.
  * Given an approximate solution \f$ y^{[j]} \in \mathbb{P}_N \f$, the error function
  * \f[ \delta^{[j]} = y - y^{[j]} \f]
  * satisfies the defect equation
  * \f[ \delta^{[j]'}(t) = y'(t) - y^{[j]'}(t) = f(y(t)) - y^{[j]'}(t). \f]
  * The equivalent Picard equation is
  * \f[ \delta^{[j]}(t) = \int_{\tau=0}^t \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - y^{[j]'}(\tau) \,d\tau\right)\f]
  * Evaluated at the grid nodes \f$ t_i\f$ we obtain
  * \f{eqnarray*}{
  * \delta^{[j]}_i & = & \delta^{[j]}_{i-1} + \int_{\tau=t_{i-1}}^{t_i} \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - {y^{[j]}}'(\tau) \right) \,d\tau \\
  *                & = & \delta^{[j]}_{i-1} + \int_{\tau=t_{i-1}}^{t_i} \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - f(y^{[j]}(\tau))\right) \,d\tau
  *               + \int_{\tau=t_{i-1}}^{t_i} \hspace{-0.5em}f(y^{[j]}(\tau)) \,d\tau - ( y^{[j]}_i-y^{[j]}_{i-1})
  * \f}
  * starting at \f$ \delta_0^{[j]} = 0. \f$ Using left looking rectangular rule for approximating the first integral, and the canonical quadrature by
  * polynomial interpolation on the nodes \f$ t_1,\ldots t_N \f$ for second integral, approximate values \f$ \hat{\delta}_i^{[j]} \f$ for
  * \f$ \delta_i^{[j]} \f$ can be evaluated. With the left looking rectangular rule, we obtain the explicit scheme
  * \f{eqnarray*}{
  *  \hat\delta^{[j]}_i & = &  \hat\delta^{[j]}_{i-1}
  * + (t_i-t_{i-1}) \left( f(y^{[j]}_{i-1}+\hat{\delta}^{[j]}_{i-1}) - f(y^{[j]}_{i-1})\right)
  * +  \sum_{k=1}^N S_{ik} f(y^{[j]}_k) - y^{[j]}_i+y^{[j]}_{i-1},
  * \f}
  * where the entries of the spectral quadrature matrix \f$ S \in \mathbb{R}^{N \times N} \f$ are defined in terms of the Lagrange polynomials
  * \f$ L_k \in \mathbb{P}_N \f$ satisfying \f$ L_k(t_i) = \delta_{ik} \f$ as
  * \f{eqnarray*}{
  *   S_{ik} & = & \int_{\tau=t_{i-1}}^{t_i} L_k(\tau) \,d\tau, \quad i,k =1,\dots, N.
  * \f}
  * An improved approximation \f$ y^{[j+1]} \f$ is then obtained by polynomial interpolation of \f$ \hat \delta^{[j]}_i \f$,
  * \f{eqnarray*}{
  * y^{[j+1]} = y^{[j]} + \hat\delta^{[j]}.
  * \f}
  *
  *\tparam Vector  a vector type, usually Dune::DenseVector.
  *
  * \param[in] grid the collocation time grid
  * \param[in] rhsFunc is the function which represents the right hand side \f$ f \f$ and is used in the computation of \f$ \delta y \f$.
  * \param[out] yi stores the current iterate of the solution \f$ y \f$.
  * \param[out] dyi is the approximate correction \f$ \delta y \f$
  * \param[in] verbose a boolean used to print out the iterates of \f$ y \f$.
  *
  */

  template<class Vector>
  void sdcExplicitEulerIterationStep(Kaskade::SDCTimeGrid const& grid,
                                     std::function<Vector(typename Vector::value_type, Vector const&)> rhsFunc,
                                     std::vector<Vector> & yi, std::vector<Vector> & dyi,
                                     bool verbose = false)
  {
    //extract time points from the grid
    auto const& pts = grid.points();
    //extracting values of the integration matrix corresponding to the time grid.
    auto const& integ = grid.integrationMatrix();
    //number of intervals
    int const n = pts.size() - 1;

    //including starting point for yi and dyi
    assert(yi.size() == n+1);
    assert(dyi.size() == n+1);

    //initialize all elements of dyi to zero.
    for (int i = 0; i <= n; ++i)
    {
      dyi[i] = Vector(0.0);
    }

    //declare Vector total to be used inside the for loop
    Vector total(0.0);

    //the SDC step for explicit Euler
    //perform n Euler steps
    boost::timer::cpu_timer timer;
    for (auto i = 1; i <= n; ++i)
    {
      //Initialize all the entries of total to zero
      total = Vector(0.0);
      //computation of the Lagrange interpolation
      for (int k = 0; k <= n; ++k)
      {
        total += integ[i-1][k] * rhsFunc(pts[k], yi[k]);
      }
      //approximate correction computation
      dyi[i] = dyi[i-1] + (pts[i] - pts[i-1]) *(rhsFunc(pts[i-1],yi[i-1]+dyi[i-1]) - rhsFunc(pts[i-1],yi[i-1])) + total - yi[i] + yi[i-1];
    }

    //Compute the next iteration of yi.
    for (auto i = 0; i <= n; ++i)
    {
      yi[i] += dyi[i];
      if (verbose)
        std::cout << "yi[" << i << "] = " << yi[i] << "\n";
    }

  }

  //=================================================================================================
   //  Implementation for single inexact sdc iteration step based on explicit Euler method.
   //=================================================================================================


   template<class Vector, class RealVector>
   void inexactSDCExplicitEulerIterationStep(Kaskade::SDCTimeGrid const& grid,
                                             std::function<Vector(typename Vector::value_type, Vector const&)> rhsFunc,
                                             RealVector const& toleranceVector, typename Vector::value_type rho,
                                             Kaskade::NormType norm_t,
                                             std::vector<Vector> & yi, std::vector<Vector> & dyi,
                                             bool verbose = false)
     {
       //extract time points from the grid
       auto const& tpts = grid.points();
       //extracting values of the integration matrix corresponding to the time grid.
       auto const& integ = grid.integrationMatrix();
       //number of intervals
       unsigned int const n = tpts.size() - 1;

       //including starting point for yi and dyi
       assert(yi.size() == n+1);
       assert(dyi.size() == n+1);

       //initialize all elements of dyi to zero.
       for (auto i = 0u; i <= n; ++i)
       {
         dyi[i] = Vector(0.0);
       }

       //declare Vector total to be used inside the for loop
       Vector total(0.0);

       //convert the toleranceVector which is a vector of doubles to a vector of Vectors.
       std::vector<Vector> errVector(n+1, Vector(0.0));
       //depending on the normType fill in err vector. This step is taken since the tolerance vector is a vector of scalars
       switch(norm_t)
       {
       case Kaskade::NormType::ONE_NORM:
       {
         auto sz = dyi[0].size();
         for (auto i = 1u; i <= n; ++i)
           errVector[i] = Vector(toleranceVector[i]/sz);
         break;
       }
       case Kaskade::NormType::MAX_NORM:
       {
         for (auto i = 1u; i <= n; ++i)
           errVector[i] = Vector(toleranceVector[i]);
         break;
        }
       }

       //the perturbed SDC step for explicit Euler
       //perform n Euler steps
       boost::timer::cpu_timer timer;
       for (auto i = 1u; i <= n; ++i)
       {
         //Initialize all the entries of total to zero
         total = Vector(0.0);

         //computation of the Lagrange interpolation
         for (auto k = 0u; k <= n; ++k)
         {
           total += integ[i-1][k] * (rhsFunc(tpts[k], yi[k]) + errVector[k]);
         }
         //approximate correction computation
         dyi[i] = dyi[i-1] + (tpts[i] - tpts[i-1]) *((rhsFunc(tpts[i-1],yi[i-1]+dyi[i-1]) - rhsFunc(tpts[i-1],yi[i-1]))- (1 - rho)*errVector[i-1]) + total - yi[i] + yi[i-1];
       }

       //Compute the next iteration of yi.
       for (auto i = 0u; i <= n; ++i)
       {
         yi[i] += dyi[i];
         if (verbose)
           std::cout << "yi[" << i << "] = " << yi[i] << "\n";
       }
     }


  //===================================================================================================
  //   Implementation of the abstract base class SDCUtil for all utility methods depending on the norms.
  //===================================================================================================

  template <class Vector, class Norm>
  class SDCUtil
  {
  public:
    using field_type = typename Vector::value_type;
    using RealVector = Dune::DynamicVector<double>;
    using RealMatrix = Dune::DynamicMatrix<double>;

    //virtual function
    /**
     * \ingroup eulerSDC
     * This function computes an estimate of the SDC contraction factor (\f$ \rho \f$) given three consecutive iterations of \f$ y \f$ and
     * depends on the given norm.
     * The estimate of \f$ \rho \f$ is given by:
     * \f{eqnarray*}{
     *  \rho & = & \frac{\|y^{[j+1]} - y^{[j]}\|}{\|y^{[j]} - y^{[j-1]}\| + \|y^{[j+1]} - y^{[j]}\|}
     * \f}
     *
     * \param[in] yPrev a Vector denoting previous iterate of \f$ y \f$.
     * \param[in] yCurrent a Vector denoting current iterate of \f$ y \f$.
     * \param[in] yNext a Vector denoting next iterate of \f$ y \f$   as obtained from SDC iteration step.
     * \param[in] norm an object of abstract base class type Norm, where we can provide the specific norm we consider for a problem.
     *
     * \return \f$ \rho \f$ an estimate of the SDC contraction factor. Return type is Vector::field_type
     */
    field_type sdcContractionFactor(std::vector<Vector> const& yPrev,
                                    std::vector<Vector> const& yCurrent,
                                    std::vector<Vector> const& yNext,
                                    Norm& norm);

    //pure virtual functions
    /**
     * \ingroup eulerSDC
     * The function definition changes with the associated norm. For derivation details for the vectors \f$ \alpha \f$
     * and \f$ \Gamma \f$ see M. Weiser, S. Ghosh: Adaptive inexact SDC Methods.
     * The vector \f$ \alpha \f$ is given by \f$  \alpha = \{\alpha_1, \ldots , \alpha_N \}\f$, where
     * \f{eqnarray*}{
     * \alpha_i & = & \sum_{k=1}^N |S_{ki}| + (1 - \delta_{i,N})\,(t_{i+1} - t_i)\,(1 + \rho),
     * \f}
     * this is computed in case of \f$ 1\f$-norm and in the case of \f$ \max\f$-norm we have
     * \f$ \Gamma = \{\Gamma_1, \ldots, \Gamma_N \} \f$, where
     *  \f{eqnarray*}{
     *    \Gamma_i & = & \max_{1\leq n \leq N} \alpha_{n,i}, \textrm{ such that } \\
     *    \alpha_{n,i} & = & |S_{ni}| + \delta_{i,n-1}\,(t_i - t_{i-1})\,(1+\rho)
     *  \f}
     * @param[in] timePoints : Time points generated using the points() method of SDCTimeGrid.
     * @param[in] integrationMatrix : Integration matrix generated using the integrationMatrix() method of SDCTimeGrid.
     *                            Computed once for a given time grid.
     * @param[in] yPrev : Estimate of exact solution \f$ y \f$ using SDC iteration step
     * @param[in] yCurrent : The following estimate of \f$ y \f$ using SDC iteration step
     * @param[in] yNext : The next estimate of \f$ y \f$ after yCurrent using SDC iteration step
     *
     * \return Reference to a Vector, where the vector is \f$ \alpha \f$ or \f$ \Gamma \f$ depending on the norm.
     */

    virtual RealVector const& computeAlphaVec(RealVector const& timePoints,
                                              RealMatrix const& integrationMatrix,
                                              std::vector<Vector> const& yPrev,
                                              std::vector<Vector> const& yCurrent,
                                              std::vector<Vector> const& yNext) = 0;


    //virtual destructor
    virtual ~SDCUtil(){}
  };

  //===================================================================================================
  //   Implementation of the derived class SDCUtilOneNorm w.r.t. 1-norm.
  //===================================================================================================
  template<class Vector, class Norm>
  class SDCUtilOneNorm : public SDCUtil<Vector, Norm>
  {
  public:
    using field_type = typename Vector::value_type;
    using RealVector = Dune::DynamicVector<double>;
    using RealMatrix = Dune::DynamicMatrix<double>;

    /**
     *
     * @param timePts : Represents a vector of time points
     * @param integrationMatrix
     * @param yPrev
     * @param yCurrent
     * @param yNext
     * @return
     */

    virtual RealVector const& computeAlphaVec(RealVector const& timePts,
                                              RealMatrix const& integrationMatrix,
                                              std::vector<Vector> const& yPrev,
                                              std::vector<Vector> const& yCurrent,
                                              std::vector<Vector> const& yNext)
    {
      //number of collocation points
      auto nCollocationPts = timePts.size();
      alphaVec = RealVector(nCollocationPts-1, 0.0);

      //create one norm object
      Kaskade::OneNorm<Vector> on;
      //compute sdc contraction rate
      auto rho = Kaskade::SDCUtil<Vector, Norm>::sdcContractionFactor(yPrev, yCurrent, yNext, on);
      for (auto i = 0u; i < nCollocationPts-1; ++i)
      {
        field_type total = 0.0;
        for (auto k = 0u; k < nCollocationPts-1; ++k)
          total += std::abs(integrationMatrix[k][i]);
        if (i != nCollocationPts-2)
          alphaVec[i] = (total + (timePts[i+1] - timePts[i]) * (1 + rho));
        else
          alphaVec[i] = total;
      }
      return alphaVec;
    }


    virtual ~SDCUtilOneNorm(){}

  private:
    RealVector alphaVec;
    RealMatrix alphaMat;
  };

  //===================================================================================================
  //   Implementation of the derived class SDCUtilMaxNorm w.r.t. max-norm.
  //===================================================================================================
  template<class Vector, class Norm>
    class SDCUtilMaxNorm : public SDCUtil<Vector, Norm>
    {
    public:
      using field_type = typename Vector::value_type;
      using RealVector = Dune::DynamicVector<double>;
      using RealMatrix = Dune::DynamicMatrix<double>;


      //TODO: Doc me
      /**
       *
       * @param timePts
       * @param integrationMatrix
       * @param yPrev
       * @param yCurrent
       * @param yNext
       * @return
       */
      virtual RealVector const& computeAlphaVec(RealVector const& timePts,
                                                RealMatrix const& integrationMatrix,
                                                std::vector<Vector> const& yPrev,
                                                std::vector<Vector> const& yCurrent,
                                                std::vector<Vector> const& yNext)
      {
        //number of collocation points.
        auto nCollocationPts = timePts.size();
        gamma = RealVector(nCollocationPts, 0.0);

        //initialize the matrix alpha with size (nCollocationPts - 1) x nCollocationPts
        auto alphaMat = computeAlphaMat(timePts, integrationMatrix, yPrev, yCurrent, yNext, nCollocationPts);
        //Compute the vector Gamma.
        for (int i = 0; i < nCollocationPts -1; ++i)
        {
          for (int j = 0; j < nCollocationPts; ++j)
          {
            if (alphaMat[i][j] > gamma[j])
              gamma[j] = alphaMat[i][j];
          }
        }
        return gamma;
       }



      virtual ~SDCUtilMaxNorm() {};
    private:
      RealVector gamma;
      RealMatrix alphaMat;

      /**
       *
       * @param timePts
       * @param integrationMatrix
       * @param yPrev
       * @param yCurrent
       * @param yNext
       * @param nCollocationPts
       * @return
       */
      RealMatrix const& computeAlphaMat(RealVector const& timePts,
                                        RealMatrix const& integrationMatrix,
                                        std::vector<Vector> const& yPrev,
                                        std::vector<Vector> const& yCurrent,
                                        std::vector<Vector> const& yNext,
                                        int nCollocationPts);


    };


  //=================================================================================================
  //   FUNCTION DEFINITION BEGINS FROM HERE
  //=================================================================================================

  //=================================================================================================
  //  Utility function to find norm of a difference of vector of two Vectors
  //=================================================================================================

  template<class Vector, class Norm>
  typename Vector::value_type normVecDiff(std::vector<Vector> const& y0, std::vector<Vector> const& y1, Norm norm)
   {
     //Initialize the all the entries of dy to 0.0
     std::vector<Vector> dy(y0.size(), Vector(0.0));
     for (auto i = 0u; i < y0.size(); ++i)
     {
       dy[i] = y1[i] - y0[i];
     }
     return norm.value(dy);
   }

  //===================================================================================================
  //   Implementation of the abstract base class SDCUtil for all utility methods depending on the norms.
  //===================================================================================================

  template <class Vector, class Norm>
  typename Vector::value_type SDCUtil<Vector, Norm>::sdcContractionFactor(std::vector<Vector> const& yPrev,
                                                                          std::vector<Vector> const& yCurrent,
                                                                          std::vector<Vector> const& yNext,
                                                                          Norm& norm)
   {
    size_t n = yPrev.size();
    //Compute the differences: yPrev - yCurrent and yCurrent - yNext
    std::vector<Vector> y01(n, Vector(0.0));
    std::vector<Vector> y12(n, Vector(0.0));
    for (auto i = 0u; i < n; ++i)
    {
      y01[i] = yCurrent[i] - yPrev[i] ;
      y12[i] = yNext[i] - yCurrent[i];
    }
    auto val01 = norm.value(y01);
    auto val12 = norm.value(y12);
    auto beta = val12/(val01 + val12);
    auto val = val12/val01;
    if (val > 1 || val != val) //check if val > 1 or val is NaN.
    {
      if (beta != beta) //check if beta = NaN
        val = 0.70;
      else
        val = beta;
    }
    return val;
   }


  //===================================================================================================
  //   Implementation of the private method computeAlphaMat for the derived class SDCUtilMaxNorm
  //===================================================================================================
  using RealVector = Dune::DynamicVector<double>;
  using RealMatrix = Dune::DynamicMatrix<double>;

  template<class Vector, class Norm>
  RealMatrix const& SDCUtilMaxNorm<Vector, Norm>::computeAlphaMat(RealVector const& timePts,
                                                                  RealMatrix const& integrationMatrix,
                                                                  std::vector<Vector> const& yPrev,
                                                                  std::vector<Vector> const& yCurrent,
                                                                  std::vector<Vector> const& yNext,
                                                                  int nCollocationPts)
    {
    //create max norm object
    Kaskade::MaxNorm<Vector> mn;
    //compute sdc contraction rate
    auto rho = Kaskade::SDCUtil<Vector, Norm>::sdcContractionFactor(yPrev, yCurrent, yNext, mn);
    alphaMat = RealMatrix(nCollocationPts-1, nCollocationPts, 0.0);
    for (auto i = 0; i < nCollocationPts-1; ++i)
    {
      for (auto k = 0; k < nCollocationPts; ++k)
      {
        if (i == k)
          alphaMat[i][k] = std::abs(integrationMatrix[i][k]) + (timePts[i+1] - timePts[i]) * (1 + rho);
        else
          alphaMat[i][k] = std::abs(integrationMatrix[i][k]);
      }
    }
    return alphaMat;
    }


}//end namespace Kaskade




#endif /* KASKADE_TIMESTEPPING_EULERSDC_UTIL_HH_ */
