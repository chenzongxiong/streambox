/*
 * sdcEuler.hh
 *
 *  Created on: May 17, 2015
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
#ifndef KASKADE_TIMESTEPPING_EULERSDC_SDCEULER_HH_
#define KASKADE_TIMESTEPPING_EULERSDC_SDCEULER_HH_

//includes from current project
#include "norm.hh"
#include "util.hh"
#include "workModel.hh"


//includes from c, c++
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

//includes from DUNE
#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

//includes from Kaskade
#include "timestepping/sdc.hh"

namespace Kaskade {
//=======================================================================================================
//  Implementation for enum class AlgorithmType
//=======================================================================================================
/**
 * \ingroup eulerSDC
 *
 * \brief enum class to steer the main integrate algorithm in EulerSDC class in a certain way.
 */
enum class AlgorithmType {
  ADAPTIVE,       //!< ADAPTIVE
  NAIVE           //!< NAIVE
};

//=======================================================================================================
//  Implementation for template class EulerSDC
//=======================================================================================================

/**
 * \ingroup eulerSDC
 * \brief Base class to perform SDC iterations based on the forward Euler Method.
 * Total iterations performed depends on the algorithm in the paper Adaptive Inexact SDC Methods by Ghosh and Weiser.
 */

template<class Vector, class Norm, class Utils, class TimeGrid=LobattoTimeGrid>
class EulerSDC
{
public:
  using field_type = typename Vector::value_type;

  using RealVector = Dune::DynamicVector<double>;
  using RealMatrix = Dune::DynamicMatrix<double>;

  using Iter = Kaskade::Iteration<Vector, Norm, Utils>;
  using Fed  = Kaskade::FiniteElementDiscretization<Vector, Norm, Utils>;

  //constructor for class EulerSDC
  //The constructor consists of member initializer lists which are relevant for all work models
  EulerSDC(field_type t0, field_type t1, size_t nIntervals_,
      field_type tol_, field_type maxSDCIterations_,
      Norm& norm_, Utils& util_,
      Kaskade::NormType norm_type, Kaskade::WorkModelType work_type,
      bool verbose_ = false) :
        tbegin(t0), tend(t1), nIntervals(nIntervals_), grid(nIntervals, tbegin, tend),
        tol(tol_), maxSDCIterations(maxSDCIterations_),
        norm(norm_), util(util_), norm_t(norm_type), work_t(work_type), verbose(verbose_) {}


  //private members necessary for different work models are set using the set methods.
  //set tolerance for Newton method used for the Iteration work model
  void setNewtonTolerance(field_type newtTol)
  {
    tolNewton = newtTol;
  }

  //set max iterations for Newton method, used for the Iteration work model
  void setNewtonIterations(int newtIter)
  {
    maxIterNewton = newtIter;
  }

  //set rhoit: contraction factor for the iteration work model
  void setRhoIT(field_type rit)
  {
    rhoit = rit;
  }

  //set mmin: minimum number of iterations necessary for the iteration work model.
  void setMmin(int minnum)
  {
    mmin = minnum;
  }

  //set dimensions, used in case of the FiniteElementDiscretization work model.
  void setDimFED(int d)
  {
    dim = d;
    //std::cout << "dim = " << dim << std::endl;
  }

  //set localTolerance for NAIVE ALGORITHM
  void setNaiveLocalTol(field_type locTol)
  {
    localTol = locTol;
  }



  /**
   * \ingroup eulerSDC
   * @param initialValue
   * @param rightHandSide
   *
   * @return
   */

  Vector integrate(Vector const& initialValue,
                   std::function<Vector (typename Vector::value_type, Vector const&)> rightHandSide,
                   Kaskade::AlgorithmType integrate_t)
  {
    //initialize the solution vector yi and the error vector dyi
    std::vector<Vector> yi(grid.points().size(), initialValue), dyi(grid.points().size(), Vector(0.0));
    std::vector<Vector> y0 = yi;
    std::vector<Vector> yPrev = yi;

    //TODO: Change to dump output to a log file
    //print the entries of the vector yi
    if (verbose)
    {
      std::cout << "print entries to the vector yi" << std::endl;
      for (auto i = 0u; i < yi.size(); ++i)
        std::cout << "[ " << yi[i] << " ]" << std::endl;
    }

    //extract time points and integration matrix from timegrid
    auto const& tPts = grid.points();
    auto const& integMat = grid.integrationMatrix();

    //initialize the tolerance vector to zero entries
    //usually the tolerance vector is the last column of the tolerance matrix.
    //We can compute an actual tolerance vector only after the first two computes of the inexact SDC iteration step.
    RealVector tolVector(nIntervals+1, 0.0);
    field_type rho = 0.5; //is there a need to set the value of rho?

    //first inexactSDC step, where the tolVector is set to a zero vector and rho = 0.5.
    Kaskade::inexactSDCExplicitEulerIterationStep(grid, rightHandSide, tolVector, rho, norm_t, yi, dyi, verbose);
    //create vectors y1, yCurrent
    std::vector<Vector> y1 = yi;
    std::vector<Vector> yCurrent = yi;
    std::vector<Vector> yNext = yi;
    int totalJ = 2;


    while (maxSDCIterations > 0)
    {
      //print dyi
      if (verbose)
      {
        std::cout << "dyi = ";
        for (int i = 0; i < dyi.size(); ++i)
        {
          std::cout << "dyi[" << i << "] = " << dyi[i] << std::endl;
        }
      }

      //carry out an inexact SDC step
      Kaskade::inexactSDCExplicitEulerIterationStep(grid, rightHandSide, tolVector, rho, norm_t, yi, dyi, verbose);
      yNext = yi;
      //compute the sdc contraction factor
      rho = util.sdcContractionFactor(yPrev, yCurrent, yNext, norm);
      //TODO: change to setw() instead of tab
      if(totalJ == 2)
        std::cout << "\t \t \t \t \t \t \t \t \t \t \t " << rho << std::endl;
      else
        std::cout << "\t \t \t " << rho << std::endl;

      //compute alphaVec or gammaVec depending on the norm
      auto alphaVec = util.computeAlphaVec(tPts, integMat, yPrev, yCurrent, yNext);

      //Create the work model object depending on the work model type
      switch(work_t)
      {
      //=====================================================================================
      //      Iteration Work Model
      //=====================================================================================
      case Kaskade::WorkModelType::ITERATION:
      {
        //initialize the iteration work model object
        Iter workModel(tol, y0, y1, norm, util, tPts, integMat, tolNewton, maxIterNewton, rhoit, mmin);
        //compute maxSDCIterations
        maxSDCIterations = workModel.getIterJ(yPrev, yCurrent, yNext, rho);
        std::cout << "\t \t" << maxSDCIterations;
        if(maxSDCIterations > 0)
        {
          switch(integrate_t)
          {
          case Kaskade::AlgorithmType::ADAPTIVE:
          {
            tolMat = workModel.computeLocalTolerances(yPrev, yCurrent, yNext);
            break;
          }
          case Kaskade::AlgorithmType::NAIVE:
          {
            //create the matrix of local tolerances
            tolMat = RealMatrix(nIntervals, maxSDCIterations, localTol);
            break;
          }
         }
         //extract the tolVector
         tolVector = extractLastColumn(tolMat);
         totalCost = workModel.computeCost(yPrev, yCurrent, yNext, tolMat);
         val = std::log(totalCost);
         std::cout << "\t\t" << val;
        }
        break;
      }
      //=====================================================================================
      //      FiniteElementDiscretization Work Model
      //=====================================================================================
      case Kaskade::WorkModelType::FED:
      {
        //initialize the finite element discretization work model object
        Fed workModel(tol, y0, y1, norm, util, tPts, integMat, dim);
        //compute maxSDCIterations
        maxSDCIterations = workModel.getIterJ(yPrev, yCurrent, rho);
        std::cout << "\t \t" << maxSDCIterations;
        if (maxSDCIterations > 0)
        {
          switch(integrate_t)
          {
          case Kaskade::AlgorithmType::ADAPTIVE:
          {
            tolMat = workModel.computeLocalTolerances(yPrev, yCurrent, yNext);
            break;
          }
          case Kaskade::AlgorithmType::NAIVE:
          {
            //create the matrix of local tolerances
            tolMat = RealMatrix(nIntervals, maxSDCIterations, localTol);
            break;
          }
         } //end switch for type of integration algorithm
         //extract the tolVector
         tolVector = extractLastColumn(tolMat);
         totalCost = workModel.computeCost(tolMat);
         val = std::log(totalCost);
         std::cout << "\t \t" << val;
        }
        break;
      }
     }

      ++totalJ;
      if (maxSDCIterations > 0)
        std::cout << "\t \t \t \t " << totalJ;
      else
        std::cout << "\t \t \t \t \t \t " << totalJ;
      yPrev = yCurrent;
      yCurrent = yNext;

    }
    //reset maxSDCIterations so as to be able to use another algorithm (NAIVE or ADAPTIVE)
    maxSDCIterations = 1;
    return yi.back();
  } //end integrate function

private:
  field_type tbegin;
  field_type tend;
  size_t nIntervals;
  TimeGrid grid;
  field_type tol;
  field_type maxSDCIterations;
  Norm norm;
  Utils util;
  Kaskade::NormType norm_t;
  Kaskade::WorkModelType work_t;
  bool verbose;
  field_type tolNewton = 1e-2;
  int maxIterNewton = 1000;
  int dim = 2;
  field_type rhoit = 0.5;
  int mmin = 3;
  RealMatrix tolMat;
  RealVector tolVec;
  field_type totalCost;
  field_type localTol = 1e-3;
  field_type val;



  //private function to extract the last column of a given matrix
  RealVector const& extractLastColumn(RealMatrix const& toleranceMatrix)
  {
    auto rows = toleranceMatrix.rows();
    auto cols = toleranceMatrix.cols();
    tolVec = RealVector(rows, 0.0);
    for (auto i = 0; i < rows; ++i)
    {
      tolVec[i] = toleranceMatrix[i][cols-1];
    }
    return tolVec;
  }
};

} //end namespace Kaskade


#endif /* KASKADE_TIMESTEPPING_EULERSDC_SDCEULER_HH_ */
