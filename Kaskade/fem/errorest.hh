/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ERROREST_HH
#define ERROREST_HH

#include <vector>

#include <boost/fusion/include/for_each.hpp>
#include <boost/utility.hpp>
#include <boost/multi_array.hpp>

namespace Kaskade 
{
  /**
   * \ingroup refcrit
   * \brief Base class for refinement criteria.
   */
  class RefinementCriterion
  {
  public:
    /**
     * \brief Computes the thresholds for refinement.
     * \param normalizedErrors a twodimensional array containing the normalized error contribution for each cell and
     *                         each variable, i.e. normalizedErrors[i][j] contains the error contribution of cell i
     *                         to the error in variable j, divided by the tolerance, such that the sum over i
     *                         should not exceed one (for acceptance)
     * \returns an array containing a threshold value such that any cell i for which normalizedErrors[i][j] > return[j]
     *          for any j is to be refined, or an empty array in case all variables are accurate enough
     */
    std::vector<double> threshold(boost::multi_array<double,2> const& normalizedErrors) const;
    
  private:
    virtual double computeThreshold(std::vector<double>& normalizedErrors, double totalError, int j) const = 0;
  };
  
  /**
   * \ingroup refcrit
   * \brief Refines the grid according to the normalized errors and the given thresholds.
   */
  template <class GridManager, class GridView>
  void markAndRefine(GridManager& gridManager, GridView const& gridView, boost::multi_array<double,2> const& normalizedErrors, 
                     std::vector<double> const& threshold)
  {
    auto end = gridView.template end<0>();
    for (auto ci=gridView.template begin<0>(); ci!=end; ++ci)
    {
      auto idx = gridView.indexSet().index(*ci);
      for (int j=0; j<threshold.size(); ++j)
        if (normalizedErrors[idx][j] >= threshold[j])
        {
          gridManager.mark(1,*ci);
          break;
        }
    }
    gridManager.adaptAtOnce();
  }
  
  /**
   * \ingroup refcrit
   * \brief Fixed fraction refinement criterion.
   * Determines the refinement thresholds such that at least the specified fraction of the cells are
   * marked for refinement. This ensures that not too few cells are marked, which could lead to very many 
   * refinement iterations, but may lead to erratic refinement if there are indeed only very local refinements
   * necessary.
   * 
   * In very rare circumstances (namely when there are many cells with exactly the same local error estimate),
   * much more cells than the required fraction can be marked for refinement.
   */
  class FixedFractionCriterion: public RefinementCriterion
  {
  public:
    FixedFractionCriterion(double fraction = 0.2);
    
  private:
    double fraction;
    virtual double computeThreshold(std::vector<double>& normalizedErrors, double totalError, int j) const;
  };
  
  /**
   * \ingroup refcrit
   * \brief Bulk refinement criterion.
   * Determines the refinement thresholds such that approximately the specified fraction of the total error 
   * is removed by the refinement (under the unrealistically optimistic assumption that refinement eliminates
   * the error completely in that cell - in reality, it's only reduced by a certain factor, so the total 
   * error is reduced somewhat less).
   */
  class BulkCriterion: public RefinementCriterion
  {
  public:
    BulkCriterion(double fraction = 0.2);
    
  private:
    double fraction;
    virtual double computeThreshold(std::vector<double>& normalizedErrors, double totalError, int j) const;
  };
  
  /**
   * \ingroup refcrit
   * \brief Max value refinement criterion.
   * Determines the refinement thresholds such that all cells with an error contribution exceeding a certain 
   * fraction of the maximum error contribution are refined.
   */
  class MaxValueCriterion: public RefinementCriterion
  {
  public:
    MaxValueCriterion(double fraction = 0.2);
    
  private:
    double fraction;
    virtual double computeThreshold(std::vector<double>& normalizedErrors, double totalError, int j) const;
  };
  
  /**
   * \ingroup refcrit
   * \brief Babuska-Rheinboldt refinement criterion.
   * Determines the refinement thresholds such that all cells with an error contribution exceeding the 
   * expected error contribution of the worst cell \em after refinement will be marked. This tends to 
   * equilibrate the error contributions quite fast.
   * 
   * Note that the local convergence order should be estimated (important in the vicinity of singularities)
   * but is here assumed to be the fixed specified order.
   */
  class BabuskaRheinboldtCriterion: public MaxValueCriterion
  {
  public:
    BabuskaRheinboldtCriterion(std::vector<int> const& order);
    
  private:
    std::vector<int> order;
    virtual double computeThreshold(std::vector<double>& normalizedErrors, double totalError, int j) const;
  };
  
  //---------------------------------------------------------------------------------------------
  
  namespace ErrorestDetail
  {
    template <class GroupByCell>
    struct GroupedSummationCollector
    {
      GroupedSummationCollector(GroupByCell const& group_):
        group(group_)
      {}
      
      template <class Cell>
      int integrationOrder(Cell const& , int shapeFunctionOrder) const
      {
        return shapeFunctionOrder;
      }
      
      // need to define this because boost::multi_array does not perform assignment if the arrays have
      // different shape.
      GroupedSummationCollector<GroupByCell>& operator=(GroupedSummationCollector<GroupByCell> const& c)
      {
        auto shape = c.sums.shape();

        sums.resize(boost::extents[shape[0]][shape[1]]);
        sums = c.sums;
        group = c.group;

        return *this;
      }
      
      template <class CellPointer, class Index, class Sequence>
      void operator()(CellPointer const& ci, Index idx, double weight, Sequence const& x)
      {
        if (sums.num_elements()==0)
          sums.resize(boost::extents[group.nGroups][boost::fusion::size(x)]);
        int i = 0; // for_each needs constant functor...
        boost::fusion::for_each(x,Add(sums,group[idx],i,weight));
      }
      
      void join(GroupedSummationCollector<GroupByCell> const& c)
      {
        auto shape = c.sums.shape();
        auto myshape = sums.shape();
        
        if(c.sums.num_elements()==0) return;
        if (sums.num_elements()==0 || myshape[0]!=shape[0] || myshape[1]!=shape[1]) // not yet initialized
        {
          sums.resize(boost::extents[shape[0]][shape[1]]);
          sums = c.sums;                                  // -> do a simple copy
        }
        else
          for (size_t i=0; i<shape[0]; ++i)               // otherwise add up
            for (size_t j=0; j<shape[1]; ++j)
              sums[i][j] += c.sums[i][j];
      }

      boost::multi_array<double,2> sums;

    private:
      GroupByCell group;

      struct Add
      {
        Add(boost::multi_array<double,2>& sums_, int idx_, int& i_, double w_):
          sums(sums_), idx(idx_), i(i_) , w(w_)
        {}

        template <class T>
        void operator()(T const& t) const { sums[idx][i++] += w*t; }

      private:
        boost::multi_array<double,2>& sums;
        int idx;
        int& i;
        double w;
      };
    };
    
  }
}

#endif
