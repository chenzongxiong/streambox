/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>


#include "fem/errorest.hh"
#include "utilities/power.hh"

namespace Kaskade
{
  
  std::vector<double> RefinementCriterion::threshold(boost::multi_array<double,2> const& normalizedErrors) const
  {
    std::vector<double> thr(normalizedErrors.shape()[1]);
    size_t nCells = normalizedErrors.shape()[0];
    
    std::vector<double> tmp(nCells);
    
    bool refine = false;

    // step through all variables
    for (int j=0; j<thr.size(); ++j)
    {
      // extract the normalized errors of variable j
      double totalError = 0;
      for (size_t i=0; i<nCells; ++i) 
      {
        tmp[i] = normalizedErrors[i][j];
        totalError += tmp[i];
      }
      
      if (totalError <= 1)
        thr[j] = 2; // definitely larger than any error contribution
      else
      {
        refine = true; // remember to refine
        thr[j] = computeThreshold(tmp,totalError,j);
      }
    }

    // If all variables are accurate enough, signal no refinement necessary
    if (!refine)
      thr.clear();
    
    return thr;
  }
  

  FixedFractionCriterion::FixedFractionCriterion(double fraction_): fraction(fraction_) 
  {
    assert(0<fraction && fraction<1);
  }
    
  double FixedFractionCriterion::computeThreshold(std::vector<double>& normalizedErrors, double, int) const
  {
    size_t offset = std::min(normalizedErrors.size()-1,static_cast<size_t>((1-fraction)*normalizedErrors.size())); // obtain position of fraction median
    
    // obtain the fraction median value
    std::nth_element(normalizedErrors.begin(),normalizedErrors.begin()+offset,normalizedErrors.end());
    return normalizedErrors[offset];
  }
  
  BulkCriterion::BulkCriterion(double fraction_): fraction(fraction_)
  {
    assert(0<fraction && fraction<1);
  }
  
  double BulkCriterion::computeThreshold(std::vector<double>& normalizedErrors, double totalError, int) const
  {
    // Naive implementation: sort and sum up. TODO: look for threshold by bisection, this leads to O(N) 
    // instead of O(N log(N)) runtime. Profile runtime on significantly large grids.
    std::sort(normalizedErrors.begin(),normalizedErrors.end());
    double error = 0;
    size_t i;
    for (i=0; i<normalizedErrors.size()-1 && error<(1-fraction)*totalError; ++i)
      error += normalizedErrors[i];
    return normalizedErrors[i];
  }
  
  MaxValueCriterion::MaxValueCriterion(double fraction_): fraction(fraction_)
  {
    assert(0<fraction && fraction<1);
  }
  
  double MaxValueCriterion::computeThreshold(std::vector<double>& normalizedErrors, double , int) const
  {
    return fraction * *std::max_element(normalizedErrors.begin(),normalizedErrors.end());
  }
  
  BabuskaRheinboldtCriterion::BabuskaRheinboldtCriterion(std::vector<int> const& order_): order(order_) { }
  
  double BabuskaRheinboldtCriterion::computeThreshold(std::vector<double>& normalizedErrors, double , int j) const
  {
    assert(j<order.size());
    return power(2.0,-order[j]) * *std::max_element(normalizedErrors.begin(),normalizedErrors.end());
  }
  
}
