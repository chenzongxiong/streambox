/*
 * workModel.hh
 *
 *  Created on: May 4, 2015
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
#ifndef KASKADE_TIMESTEPPING_EULERSDC_WORKMODEL_HH_
#define KASKADE_TIMESTEPPING_EULERSDC_WORKMODEL_HH_

//includes from current project
#include "norm.hh"
#include "util.hh"

//includes from c, c++
#include <cmath>
#include <iostream>
#include <vector>

//includes from boost
#include <boost/math/special_functions/round.hpp>
#include <boost/math/tools/minima.hpp>   //used for Brent's method to compute minima of a function

//includes from DUNE
#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

namespace Kaskade {
  //========================================================================================
  //     Implementation for enum class WorkModelType
  //========================================================================================

  /**
   * \ingroup eulerSDC
   * \brief Enum class to define the norm type is class Norm.
   */

  //--C++11 strongly typed enums.
  enum class WorkModelType {
    ITERATION,      //!< ITERATION
    FED             //!< FED
  };

  /**
   * \ingroup eulerSDC
   * \brief Abstract base class for different work models
   *
   * This class represents a generic base class for an arbitrary work model
   */

  //================================================================================================
  //      Implementation of the abstract base class WorkModel to compute the local tolerances.
  //================================================================================================

  /**
   * \ingroup eulerSDC
   * \brief Abstract base class WorkModel which provides an interface for different kinds of work models for
   * computation of the cost function, maximum number of iterations and local tolerances.
   */

   template<class Vector, class Norm, class Utils>
   class WorkModel
   {
   public:
     using field_type = typename Vector::value_type;
     using RealVector = Dune::DynamicVector<double>;
     using RealMatrix = Dune::DynamicMatrix<double>;



     //pure virtual functions

       /**
        * \ingroup eulerSDC
        * \brief Pure virtual function implemented in derived classes for different work models. Computes the
        *        local tolerances for every time point and every iteration step given global tolerance.
        *
        * @param yPrev      An (i-1)-th iterative approximation to the solution using sdc Iteration step.
        * @param yCurrent   An i-th iterative approximation to the solution using sdc Iteration step.
        * @param yNext      An (i+1)-th iterative approximation to the solution using sdc Iteration step.
        * @return Returns the matrix of local tolerances.
        */
       virtual RealMatrix const& computeLocalTolerances(std::vector<Vector> const& yPrev,
                                                        std::vector<Vector> const& yCurrent,
                                                        std::vector<Vector> const& yNext) = 0;

       /**
        * \ingroup eulerSDC
        *
        * @param yPrev
        * @param yCurrent
        * @param rho
        *
        * @return
        */
       virtual field_type lowerBoundIterJ(std::vector<Vector> const& yPrev,
                                          std::vector<Vector> const& yCurrent,
                                          field_type rho) = 0;


       //virtual destructor
       virtual ~WorkModel() {}
   };

   //========================================================================================================
   //             Implementation of the derived class Iteration representing the iteration work model.
   //========================================================================================================

   template<class Vector, class Norm, class Utils>
   class Iteration : public WorkModel<Vector, Norm, Utils>
   {
   public:
     using field_type = typename Vector::value_type;
     using RealVector = Dune::DynamicVector<double>;
     using RealMatrix = Dune::DynamicMatrix<double>;

     //constructor for the derived class Iteration
     Iteration(field_type tol_,
               std::vector<Vector> const& y0_, std::vector<Vector> const& y1_,
               Norm& norm_, Utils& util_,
               RealVector const& timePts_, RealMatrix const& integrationMatrix_,
               field_type tolNewton_, int maxIterNewton_,
               field_type rhoit_, int mmin_);

     // function to compute cost
     /**
      * \ingroup eulerSDC
      * @param yPrev
      * @param yCurrent
      * @param yNext
      * @param tolMat
      *
      * @return
      */

     field_type computeCost(std::vector<Vector> const& yPrev,
                            std::vector<Vector> const& yCurrent,
                            std::vector<Vector> const& yNext,
                            RealMatrix const& tolMat);


     /**
      * \ingroup eulerSDC
      * @param yPrev
      * @param yCurrent
      * @param yNext
      * @param rho
      *
      * @return
      */
     field_type getIterJ(std::vector<Vector> const& yPrev,
                         std::vector<Vector> const& yCurrent,
                         std::vector<Vector> const& yNext,
                         field_type rho);


     virtual RealMatrix const& computeLocalTolerances(std::vector<Vector> const& yPrev,
                                                      std::vector<Vector> const& yCurrent,
                                                      std::vector<Vector> const& yNext)
     {
       //To compute optimal local tolerances we first need the value of maxSDCIterations
       auto rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
       //compute maxSDCIterations
       auto maxSDCIterations = getIterJ(yPrev, yCurrent, yNext, rho);
       //then we need to compute the value of mu which depends on maxSDCIterations
       auto mu = computeMu(yPrev, yCurrent, rho, maxSDCIterations, nIntervals, tol, norm);
       //compute the value of alphaVec
       auto alphaVec = util.computeAlphaVec(timePts, integrationMatrix, yPrev, yCurrent, yNext);
       //initialize the matrix
       int itJ = (int) maxSDCIterations;
       xijMat = RealMatrix(nIntervals, itJ, 0.0);
       for (auto i = 0; i < nIntervals; ++i)
       {
         for (auto j = 0; j < itJ; ++j)
         {
           xijMat[i][j] = -1/(mu * std::pow(rho, itJ-1-j) * alphaVec[i]);
         }
       }
       return xijMat;
     }


     //function computes the lower bound of maxSDCIterations
     virtual field_type lowerBoundIterJ(std::vector<Vector> const& yPrev,
                                        std::vector<Vector> const& yCurrent,
                                        field_type rho)
     {
       auto val = Kaskade::normVecDiff(yPrev,yCurrent,norm);
       auto lboundJ = std::log((1-rho) * tol / val) / std::log(rho);
       return boost::math::round(lboundJ);
     }


     //delete later
     field_type static f (std::vector<Vector> const& yPrev,
         std::vector<Vector> const& yCurrent,
         std::vector<Vector> const& yNext,
         double rho)
       {
         auto yTotal = yPrev[0] + yCurrent[0] + yNext[0];
         return (yTotal[0] + yTotal[1]) * rho;
       }

     //destructor
     virtual ~Iteration(){}

   private:
     field_type tol;
     std::vector<Vector> y0;
     std::vector<Vector> y1;
     Norm norm;
     Utils util;
     field_type normy01;
     RealVector timePts;
     RealMatrix integrationMatrix;
     int nIntervals;
     field_type maxPrecisionBrent = 20;  //computation with 20 bit precision
     field_type rhoit;
     int mmin;
     RealMatrix xijMat;

     //private cost function

     field_type static costFunction(std::vector<Vector> const& yPrev,
                                    std::vector<Vector> const& yCurrent,
                                    std::vector<Vector> const& yNext,
                                    field_type maxSDCIterations,
                                    Norm& norm,
                                    Utils& util,
                                    field_type normy01,
                                    RealVector const& timePts,
                                    RealMatrix const& integrationMatrix,
                                    int nIntervals, field_type tol);


     field_type static computeMu(std::vector<Vector> const& yPrev,
                          std::vector<Vector> const& yCurrent,
                          field_type rho,
                          field_type maxSDCIterations,
                          int nIntervals,
                          field_type tol,
                          Norm& norm);

     field_type computeDerivativeMu(std::vector<Vector> const& yPrev,
                                    std::vector<Vector> const& yCurrent,
                                    field_type rho,
                                    field_type maxSDCIterations);

     field_type derivativeCostFunction(std::vector<Vector> const& yPrev,
                                       std::vector<Vector> const& yCurrent,
                                       field_type rho,
                                       field_type alphaProd,
                                       field_type c,
                                       field_type maxSDCIterations);

     //compute product of entries of alphaVec
     field_type static computeProduct(RealVector const& alpha);


     //This method is based on the ***assumption*** that cost function is monotonically decreasing
     field_type computeIterJ(std::vector<Vector> const& yPrev,
                       std::vector<Vector> const& yCurrent,
                       std::vector<Vector> const& yNext,
                       field_type rho);


   };


   //=======================================================================================================
   //   Implementation of the derived class FiniteElementDiscretization representing the FED work model.
   //======================================================================================================

   template<class Vector, class Norm, class Utils>
   class FiniteElementDiscretization : public WorkModel<Vector, Norm, Utils>
   {
   public:
     using field_type = typename Vector::value_type;
     using RealVector = Dune::DynamicVector<double>;
     using RealMatrix = Dune::DynamicMatrix<double>;


     //constructor for the derived class FiniteElementDiscretization
     FiniteElementDiscretization(field_type tol_,
                                 std::vector<Vector> const& y0_, std::vector<Vector> const& y1_,
                                 Norm& norm_, Utils& util_,
                                 RealVector const& timePts_, RealMatrix const& integrationMatrix_,
                                 int dim_);


     //function to compute cost
     //(cannot be a pure virtual function since the parameters are different for different work models.)
     field_type computeCost(RealMatrix const& tolMat);

     field_type getIterJ(std::vector<Vector> const& yPrev,
                         std::vector<Vector> const& yCurrent,
                         field_type rho);

     virtual RealMatrix const& computeLocalTolerances(std::vector<Vector> const& yPrev,
                                                      std::vector<Vector> const& yCurrent,
                                                      std::vector<Vector> const& yNext)
     {
       field_type d1 = 0.0;
       field_type d2 = -1.0 / ((double) dim + 1.0);
       auto rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
       //std::cout << "rho = " <<  rho << std::endl;
       //compute max sdc iterations
       maxSDCIterations = getIterJ(yPrev, yCurrent, rho);
       auto alphaVec = util.computeAlphaVec(timePts, integrationMatrix, yPrev, yCurrent, yNext);
       auto lambda = computeLambda(yPrev, yCurrent, yNext, maxSDCIterations);
       //std::cout << "lambda = " << lambda << std::endl;
       //initialize the matrix
       int itJ = (int) maxSDCIterations;
       xijMat = RealMatrix(nIntervals, itJ, 0.0);
       for (auto i = 0; i < nIntervals; ++i)
       {
         auto val = std::pow(alphaVec[i], d2);
         for (auto j = 0; j < itJ; ++j)
         {
           d1 = (maxSDCIterations-1-j)/(d2);
           xijMat[i][j] = lambda * std::pow(rho, d1) * val;
         }
       }
       return xijMat;
     }


     //function computes the lower bound of maxSDCIterations
     virtual field_type lowerBoundIterJ(std::vector<Vector> const& yPrev,
                                        std::vector<Vector> const& yCurrent,
                                        field_type rho)
     {
       auto val = Kaskade::normVecDiff(yPrev, yCurrent, norm);
       auto lboundJ = std::log((1-rho)*tol/val)/std::log(rho);
       return boost::math::round(lboundJ);
     }


     //destructor
     virtual ~FiniteElementDiscretization(){}

   private:
     field_type tol;
     std::vector<Vector> y0;
     std::vector<Vector> y1;
     Norm norm;
     field_type normy01;
     Utils util;
     RealVector timePts;
     RealMatrix integrationMatrix;
     int nIntervals;
     int dim;
     field_type maxSDCIterations;
     RealMatrix xijMat;

     //private methods
     //computes A or B depending on the norm
     field_type computeA(std::vector<Vector> const& yPrev,
                         std::vector<Vector> const& yCurrent,
                         std::vector<Vector> const& yNext);

     field_type computeLambda(std::vector<Vector> const& yPrev,
                              std::vector<Vector> const& yCurrent,
                              std::vector<Vector> const& yNext,
                              field_type maxSDCIterations);
   };

   //=================================================================================================
   //   FUNCTION DEFINITION BEGINS FROM HERE
   //=================================================================================================

   //=================================================================================================
   //   ITERATION  WORKMODEL
   //=================================================================================================

   //=================================================================================================
   //   Constructor for the class Iteration
   //=================================================================================================

   template<class Vector, class Norm, class Utils>
   Iteration<Vector, Norm, Utils>::Iteration(field_type tol_,
               std::vector<Vector> const& y0_, std::vector<Vector> const& y1_,
               Norm& norm_, Utils& util_,
               RealVector const& timePts_, RealMatrix const& integrationMatrix_,
               field_type tolNewton_, int maxPrecisionBrent_,
               field_type rhoit_, int mmin_)
               : tol(tol_), y0(y0_), y1(y1_),
                 norm(norm_), util(util_), normy01(Kaskade::normVecDiff(y0,y1,norm)),
                 timePts(timePts_), integrationMatrix(integrationMatrix_), nIntervals(timePts.size()-1),
                 maxPrecisionBrent(maxPrecisionBrent_),
                 rhoit(rhoit_), mmin(mmin_){}

   //=================================================================================================
   //   Public member function for Iteration WorkModel to compute total cost.
   //=================================================================================================

   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeCost(std::vector<Vector> const& yPrev,
                                                                           std::vector<Vector> const& yCurrent,
                                                                           std::vector<Vector> const& yNext,
                                                                           RealMatrix const& tolMat)
   {
     //initialize the total cost
     field_type totalCost = 0;
     //compute SDC contraction factor
     auto rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
     //compute c = log(|y1-y0|)
     auto c = std::log(normy01);
     //compute the constant k = -mmin * log(rhoit)
     auto k = -mmin * std::log(rhoit);
     //dimensions of the tolerance matrix
     auto rows = tolMat.rows();
     auto cols = tolMat.cols();
     for (auto i = 0u; i < rows; ++i)
     {
       for (auto j = 0u; j < cols; ++j)
       {
         totalCost += std::max(k, c - std::log(tolMat[i][j]/std::pow(rho,j)));
       }
     }
     return totalCost;
   }

   //=================================================================================================
   //   Public member function for Iteration WorkModel to get the maxSDCIterations.
   //=================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::getIterJ(std::vector<Vector> const& yPrev,
                                                                        std::vector<Vector> const& yCurrent,
                                                                        std::vector<Vector> const& yNext,
                                                                        field_type rho)
   {
     auto maxSDCIterations = computeIterJ(yPrev, yCurrent, yNext, rho);
     return boost::math::round(maxSDCIterations);
   }


   //=================================================================================================
   //   Private member function for Iteration WorkModel for the costFunction.
   //=================================================================================================


   //make it a static function, and do not use any private member variables.
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::costFunction(std::vector<Vector> const& yPrev,
                                                                            std::vector<Vector> const& yCurrent,
                                                                            std::vector<Vector> const& yNext,
                                                                            field_type maxSDCIterations,
                                                                            Norm& norm,
                                                                            Utils& util,
                                                                            field_type normy01,
                                                                            RealVector const& timePts,
                                                                            RealMatrix const& integrationMatrix,
                                                                            int nIntervals,
                                                                            field_type tol)
   {
     //compute SDC contraction factor
     auto rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
     //compute c = log(|y1-y0|)
     auto c = std::log(normy01);
     //compute mu
     auto mu = computeMu(yPrev, yCurrent, rho, maxSDCIterations, nIntervals, tol, norm);
     auto alphaVec = util.computeAlphaVec(timePts, integrationMatrix, yPrev, yCurrent, yNext);
     auto alphaProd = computeProduct(alphaVec);
     auto cost = nIntervals * c + nIntervals * maxSDCIterations * std::log(rho) + nIntervals * std::log(-mu) + std::log(alphaProd);
     return cost;
   }


   //==========================================================================================================
   //   Private member function for Iteration WorkModel for computing the mu function as described in the paper
   //==========================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeMu(std::vector<Vector> const& yPrev,
                                                                         std::vector<Vector> const& yCurrent,
                                                                         field_type rho,
                                                                         field_type maxSDCIterations,
                                                                         int nIntervals,
                                                                         field_type tol,
                                                                         Norm& norm)
   {
     auto mu = -((1-rho) * (maxSDCIterations - 1) * nIntervals)/((1-rho) * tol - std::pow(rho, maxSDCIterations) * Kaskade::normVecDiff(yPrev,yCurrent,norm));
     return mu;
   }

   //========================================================================================================================
   //   Private member function for Iteration WorkModel for computing the derivative of mu function as described in the paper
   //========================================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeDerivativeMu(std::vector<Vector> const& yPrev,
                                                                                   std::vector<Vector> const& yCurrent,
                                                                                   field_type rho,
                                                                                   field_type maxSDCIterations)
   {
     auto val = Kaskade::normVecDiff(yPrev,yCurrent,norm);
     auto rhoJ = std::pow(rho, maxSDCIterations);
     auto dervMu = nIntervals * (1-rho) *
         (-(1-rho) * tol + val * rhoJ * (1 - (maxSDCIterations - 1) * std::log(rho)))/
         (std::pow((rhoJ * val - (1-rho) * tol), 2));
     return dervMu;
   }

   //========================================================================================================================
   //   Private member function for Iteration WorkModel for computing the derivative of cost Function
   //========================================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::derivativeCostFunction(std::vector<Vector> const& yPrev,
                                                                                      std::vector<Vector> const& yCurrent,
                                                                                      field_type rho,
                                                                                      field_type alphaProd,
                                                                                      field_type c,
                                                                                      field_type maxSDCIterations)
   {

     auto mu = computeMu(yPrev, yCurrent, rho, maxSDCIterations, nIntervals, tol, norm);
     auto dervMu = computeDerivativeMu(yPrev, yCurrent, rho, maxSDCIterations);
     auto dervCost = nIntervals * c + nIntervals * (2*maxSDCIterations - 1) * std::log(rho) + std::log(alphaProd)
                     + nIntervals * std::log(-mu) + nIntervals * (maxSDCIterations - 1) * dervMu/mu;
     return dervCost;
   }

   //==========================================================================================================================
   //   Private member function for Iteration WorkModel for computing the product of alphaVector where all entries are positive
   //==========================================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeProduct(RealVector const& alpha)
   {
     field_type alphaProd = 1;
     for(auto i = 0u; i < alpha.size(); ++i)
       alphaProd *= alpha[i];

     return alphaProd;
   }

   //============================================================================================
   //   Private member function for Iteration WorkModel for computing maxSDCIterations
   //============================================================================================

//   //Use Brent's method to compute the minimum of the cost function
//   template<class Vector, class Norm, class Utils>
//   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeIterJ(std::vector<Vector> const& yPrev,
//                                                                            std::vector<Vector> const& yCurrent,
//                                                                            std::vector<Vector> const& yNext,
//                                                                            field_type rho)
//   {
//     using Result = std::pair<double, double>;
//     //compute lower bound of maxSDCIterations to give as starting value for to Brent's method
//     field_type lbmaxSDCIterations = lowerBoundIterJ(yPrev, yCurrent, rho);
//     //initialize interval where the minima is searched
//     field_type left = lbmaxSDCIterations;
//     //initialize the lower bound on maxSDCIterations as input guess.
//     //auto func = std::bind(f, yPrev, yCurrent, yNext, std::placeholders::_1);
//     auto func = std::bind(costFunction, yPrev, yCurrent, yNext, std::placeholders::_1, norm, util, normy01,
//                           timePts, integrationMatrix, nIntervals, tol);
//     Result soln = boost::math::tools::brent_find_minima(func, left, left+6, maxPrecisionBrent);
//
//     return soln.first;
//   }

   //Numeric method to find minimia of the function
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type Iteration<Vector, Norm, Utils>::computeIterJ(std::vector<Vector> const& yPrev,
                                                                            std::vector<Vector> const& yCurrent,
                                                                            std::vector<Vector> const& yNext,
                                                                            field_type rho)
     {
     //compute lower bound of maxSDCIterations
     field_type lbmaxSDCIterations = lowerBoundIterJ(yPrev, yCurrent, rho);
     field_type maxSDCIterations = lbmaxSDCIterations + 1;
     //compute the cost at lbmaxSDCIterations
     field_type minCost = costFunction(yPrev, yCurrent, yNext, lbmaxSDCIterations, norm, util, normy01,
                                       timePts, integrationMatrix, nIntervals, tol);
     field_type val = costFunction(yPrev, yCurrent, yNext, maxSDCIterations, norm, util, normy01,
                                   timePts, integrationMatrix, nIntervals, tol);
     while(val > minCost)
     {
       minCost = val;
       maxSDCIterations++;
       val = costFunction(yPrev, yCurrent, yNext, maxSDCIterations, norm, util, normy01,
                          timePts, integrationMatrix, nIntervals, tol);
     }
     return maxSDCIterations;
     }



   //=================================================================================================
   //   FINITE ELEMENT DISCRETIZATION  WORKMODEL
   //=================================================================================================

   //=================================================================================================
   //   Constructor for the class FiniteElementDiscretization
   //=================================================================================================
   template<class Vector, class Norm, class Utils>
   FiniteElementDiscretization<Vector, Norm, Utils>::FiniteElementDiscretization(field_type tol_,
                                    std::vector<Vector> const& y0_, std::vector<Vector> const& y1_,
                                    Norm& norm_, Utils& util_,
                                    RealVector const& timePts_, RealMatrix const& integrationMatrix_,
                                    int dim_)
                                    : tol(tol_),
                                      y0(y0_), y1(y1_), norm(norm_), normy01(Kaskade::normVecDiff(y0, y1, norm)),
                                      util(util_),
                                      timePts(timePts_), integrationMatrix(integrationMatrix_), nIntervals(timePts.size()-1),
                                      dim(dim_){}

   //========================================================================================================
   //   Public member function computeCost for the class FiniteElementDiscretization, computes the total cost
   //========================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type FiniteElementDiscretization<Vector, Norm, Utils>::computeCost(RealMatrix const& tolMat)
   {
     //initialize the total cost
     field_type totalCost = 0.0;
     //dimensions of the tolerance matrix
     auto rows = tolMat.rows();
     auto cols = tolMat.cols();
     for (auto i = 0u; i < rows; ++i)
     {
       for (auto j = 0u; j < cols; ++j)
         totalCost += 1/std::pow(tolMat[i][j], dim);
     }
     return totalCost;
   }


   //===========================================================================================================
   //   Public member function getIterJ for the class FiniteElementDiscretization, returns the maxSDCIterations
   //===========================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type FiniteElementDiscretization<Vector, Norm, Utils>::getIterJ(std::vector<Vector> const& yPrev,
                                                                                          std::vector<Vector> const& yCurrent,
                                                                                          field_type rho)
   {
     auto val = Kaskade::normVecDiff(yPrev, yCurrent, norm);
     maxSDCIterations = (dim + 1)/(std::log(rho)) * std::log((1-rho)*tol / val);

     return boost::math::round(maxSDCIterations);
   }

   //===============================================================================================================
   //   Private member function for FiniteElementDiscretization WorkModel for computing A or B depending on the norm
   //===============================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type FiniteElementDiscretization<Vector, Norm, Utils>::computeA(std::vector<Vector> const& yPrev,
                                                                                          std::vector<Vector> const& yCurrent,
                                                                                          std::vector<Vector> const& yNext)
   {
     //compute alphaVec
     auto alphaVec = util.computeAlphaVec(timePts, integrationMatrix, yPrev, yCurrent, yNext);
     //initialize alphaTildeVec
     auto alphaTildeVec = RealVector(alphaVec.size(), 0.0);
     //casting dim/(dim + 1) to double
     double d = ((double) dim)/((double) dim + 1);
     field_type sumAlphaTilde = 0.0;
     //compute alphaTildeVec and A the sum of alphaTilde
     for (auto i = 0u; i < alphaVec.size(); ++i)
     {
       alphaTildeVec[i] = std::pow(alphaVec[i], d);
       sumAlphaTilde += alphaTildeVec[i];
     }
     return sumAlphaTilde;
   }

   //====================================================================================================================
   //   Private member function for FiniteElementDiscretization WorkModel for computing lambda as described in the paper.
   //====================================================================================================================
   template<class Vector, class Norm, class Utils>
   typename Vector::value_type FiniteElementDiscretization<Vector, Norm, Utils>::computeLambda
                                                                           (std::vector<Vector> const& yPrev,
                                                                            std::vector<Vector> const& yCurrent,
                                                                            std::vector<Vector> const& yNext,
                                                                            field_type maxSDCIterations)
   {
     //compute rho, sdc-contraction factor
     auto rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
     //casting
     double d = ((double) dim)/((double) dim + 1);
     field_type tildeRho = std::pow(rho, d);
     auto sumAlphaTilde = computeA(yPrev, yCurrent, yNext);
     auto val = Kaskade::normVecDiff(yPrev, yCurrent, norm);
     auto lambda = ((1-rho) * tol - val * std::pow(rho, maxSDCIterations)) * (1 - tildeRho) / ((1 - rho) * (1 - std::pow(tildeRho, maxSDCIterations)) * sumAlphaTilde);

     return lambda;
   }

} //end namespace Kaskade

#endif /* KASKADE_TIMESTEPPING_EULERSDC_WORKMODEL_HH_ */
