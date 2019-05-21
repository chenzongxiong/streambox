/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ADAPTATION_STRATEGY_HH
#define ADAPTATION_STRATEGY_HH

#include <algorithm>
#include <utility>
#include <vector>

namespace Kaskade
{
  namespace AdaptationStrategy_Detail
  {
    struct BiggerThanAbs
    {
      template <class Scalar, class IndexInt>
      bool operator()(const std::pair<Scalar, IndexInt>& x1, const std::pair<Scalar,IndexInt>& x2)
      {
        return fabs(x1.first) > fabs(x2.first);
      }
    };

    struct Add
    {
      template <class Scalar, class IndexInt>
      Scalar operator()(Scalar x1, const std::pair<Scalar,IndexInt>& x2)
      {
        return x1 + x2.first;
      }
    };
  }

  namespace Adaptivity
  {
    template <class Grid>
    class FixedCellFraction
    {
      static constexpr int dim = Grid::dimension;
      typedef size_t IndexInt;
    public:
      template <typename... Args>
      explicit FixedCellFraction(GridManagerBase<Grid>& gridManager_, double fraction_, Args...)
       : gridManager(gridManager_), fraction(fraction_)
      {}

      template <class Err, class ErrorRepresentation, class Scalar>
      void refineGrid_impl(Err const& err, ErrorRepresentation& errorDistribution, Scalar tol)
      {
        std::sort(errorDistribution.begin(), errorDistribution.end(), AdaptationStrategy_Detail::BiggerThanAbs());
        std::cout << "max error2: " << std::max_element(errorDistribution.begin(), errorDistribution.end())->first << std::endl;
        std::cout << "min error2: " << std::min_element(errorDistribution.begin(), errorDistribution.end())->first << std::endl;

        std::vector<std::pair<Scalar,IndexInt> > bulkErrorDistribution(errorDistribution.begin(), errorDistribution.begin() + static_cast<size_t>(errorDistribution.size() * fraction));
	Scalar bulkErrorSquared = std::accumulate(bulkErrorDistribution.begin(),bulkErrorDistribution.end(),0.0,AdaptationStrategy_Detail::Add());

        std::cout << "bulkErrorDistribution.size()=" << bulkErrorDistribution.size() << std::endl;
        std::cout << "bulkErrorSquared=" << bulkErrorSquared << std::endl;
        std::cout << "remaining error contribution (squared): " << std::accumulate(errorDistribution.begin()+bulkErrorDistribution.size(),errorDistribution.end(), 0.0, AdaptationStrategy_Detail::Add()) << std::endl;

        // do refinement
        size_t counter = 0;
        std::vector<int> visitedCells;
        auto cend = err.descriptions.gridView.template end<0>();
        for (auto ci=gridManager.grid().leafView().template begin<0>(); ci!=cend; ++ci) // iterate over cells
        {
          for(int i=0; i<bulkErrorDistribution.size(); ++i) // iterate over chosen part of the error distribution
            if(gridManager.grid().leafIndexSet().index(*ci) == bulkErrorDistribution[i].second)
            {

              visitedCells.push_back(gridManager.grid().leafIndexSet().index(*ci));

              ++counter;
              gridManager.mark(1,*ci);
            }
        }

        std::cout << "FIXED CELL FRACTION: marked " << counter << " of " << gridManager.grid().size(0) << " cells for refinement." << std::endl;
        gridManager.adaptAtOnce();
      }

    private:
      GridManagerBase<Grid>& gridManager;
      double fraction;
    };

    template <class Grid>
    class FixedErrorFraction
    {
      static constexpr int dim = Grid::dimension;
      typedef size_t IndexInt;
    public:
      template <typename... Args>
      explicit FixedErrorFraction(GridManagerBase<Grid>& gridManager_, double fraction_, Args...)
       : gridManager(gridManager_), fraction(fraction_)
      {}

      template <class Err, class ErrorRepresentation, class Scalar>
      void refineGrid_impl(Err const& err, ErrorRepresentation& errorDistribution, Scalar tol)
      {
        tol *= fraction*fraction;
        std::cout << "used tolerance: " << tol << std::endl;
        std::sort(errorDistribution.begin(), errorDistribution.end(), AdaptationStrategy_Detail::BiggerThanAbs());
        std::cout << "max error2: " << std::max_element(errorDistribution.begin(), errorDistribution.end())->first << std::endl;
        std::cout << "min error2: " << std::min_element(errorDistribution.begin(), errorDistribution.end())->first << std::endl;

        Scalar bulkErrorSquared = 0;
        std::vector<std::pair<Scalar,IndexInt> > bulkErrorDistribution;

        while(bulkErrorSquared < tol)
        {
          bulkErrorDistribution.push_back(errorDistribution[bulkErrorDistribution.size()]);
          bulkErrorSquared += bulkErrorDistribution.back().first;
        }

        std::cout << "bulkErrorDistribution.size()=" << bulkErrorDistribution.size() << std::endl;
        std::cout << "bulkErrorSquared=" << bulkErrorSquared << std::endl;
        std::cout << "remaining error contribution (squared): " << std::accumulate(errorDistribution.begin()+bulkErrorDistribution.size(),errorDistribution.end(), 0.0, AdaptationStrategy_Detail::Add()) << std::endl;

        // do refinement
        size_t counter = 0;
        std::vector<int> visitedCells;
        auto cend = err.descriptions.gridView.template end<0>();
        for (auto ci=gridManager.grid().leafView().template begin<0>(); ci!=cend; ++ci) // iterate over cells
        {
          for(int i=0; i<bulkErrorDistribution.size(); ++i) // iterate over chosen part of the error distribution
            if(gridManager.grid().leafIndexSet().index(*ci) == bulkErrorDistribution[i].second)
            {

              visitedCells.push_back(gridManager.grid().leafIndexSet().index(*ci));

              ++counter;
              gridManager.mark(1,*ci);
            }
        }

        std::cout << "FIXED ERROR FRACTION: marked " << counter << " of " << gridManager.grid().size(0) << " cells for refinement." << std::endl;
        gridManager.adaptAtOnce();
      }

    private:
      GridManagerBase<Grid>& gridManager;
      double fraction;
    };

    template <class Grid>
    class ErrorEquilibration2
    {
    public:
      template <typename... Args>
      explicit ErrorEquilibration2(GridManagerBase<Grid>& gridManager_, Args...) : gridManager(gridManager_)
      {}

      template <class Err, class ErrorRepresentation, class Scalar>
      void refineGrid_impl(Err const& err, ErrorRepresentation& errorDistribution, Scalar tol)
      {
        std::cout << "ERROR EQUILIBRATION" << std::endl;
        std::sort(errorDistribution.begin(), errorDistribution.end(), AdaptationStrategy_Detail::BiggerThanAbs());

        std::cout << "err distr" << std::endl;
//        for(auto a : errorDistribution) std::cout << errorDistribution.first << ", " <<
//
//        size_t lastIndex = errorDistribution.size()-1;
//        size_t firstIndex = 0;
//
//        while(firstIndex+3 < lastIndex)
//        {
//          if(bulkErrorDistribution.size() > 4)
//          {
//            Scalar firstContribution = 15.0/256.0 * errorDistribution[firstIndex].first;
//
//            Scalar lastContributions = bulkErrorDistribution[lastIndex--].first;
//            lastContributions += bulkErrorDistribution[lastIndex--].first;
//            lastContributions += bulkErrorDistribution[lastIndex].first;
//            lastContributions *= 15.0/16.0;
//
//            if(lastContributions < firstContribution)
//            {
//              ++lastIndexForDoubleRefinement;
//              ++firstIndex;
//              --lastIndex;
//              bulkErrorDistribution.erase(bulkErrorDistribution.end()-3, bulkErrorDistribution.end());
//            }
//            else break;
//          }
//        }
        Scalar N2 = gridManager.grid().size(0);
//        N2 *= N2;
        Scalar relTol = tol/N2;
        // do refinement
        auto cend = gridManager.grid().leafView().template end<0>();
        size_t counter = 0;
        Scalar minErr = errorDistribution.back().first, maxErr = errorDistribution[0].first;
        //constexpr int dim = ErrorRepresentation::Descriptions::Grid::dimension;
        //Dune::FieldVector<Scalar,dim> x0(0.2);

        std::cout << "grr" << std::endl;
        std::cout << "marking cells..." << std::flush;
        for (auto ci=gridManager.grid().leafView().template begin<0>(); ci!=cend; ++ci) // iterate over cells
        {
//          if(counter > N3) break;
//          std::cout << "counter: " << counter << std::endl;
//          std::cout << "id: " << is.index(*ci) << std::endl;
//          std::cout << "corners: " << std::endl;
//          for(size_t i=0; i<ci->geometry().corners(); ++i) std::cout << "i: " << ci->geometry().corner(i) << std::endl;
          for(size_t i=0; i<errorDistribution.size(); ++i) // iterate over chosen part of the error distribution
            if(gridManager.grid().leafIndexSet().index(*ci) == errorDistribution[i].second)
            {
//          auto const& tmp = boost::fusion::at_c<0>(err.data);
//          std::cout << "gsg" << std::endl;
//          Scalar val = tmp.value(ci->geometry().global(x0));
//          val = std::fabs(val);
//          //Scalar val = boost::fusion::at_c<0>(err.data).value(*ci,x0);
//          std::cout << "val=" << val << std::endl;
//          if(minErr > val) minErr = val;
//          if(maxErr < val) maxErr = val;
//          if(val > relTol){
//            ++counter;
//            gridManager.mark(1,*ci);
//          }
//              if(minErr > errorDistribution[i].first) minErr = errorDistribution[i].first;
//                  if(maxErr < errorDistribution[i].first) maxErr = errorDistribution[i].first;
              if(errorDistribution[i].first > relTol){
                std::cout << "marking cell at position " << i << std::endl;
                ++counter;
                gridManager.mark(1,*ci);
              }
            }
        }
        std::cout << "done." << std::endl;
        std::cout << "totalErrorSquared: " << tol << " -> tolerance: " << relTol << std::endl;
        std::cout << "maxErr: " << maxErr << ", minErr: " << minErr << std::endl;
        std::cout << "refining " << counter << " of " << gridManager.grid().size(0) << " cells" << std::endl;

        std::cout << "adapting..." << std::flush;
        gridManager.adaptAtOnce();
        std::cout << "done." << std::endl;
        // adjust error function
      }

    private:
      GridManagerBase<Grid>& gridManager;
    };

    template <class Grid>
    class ErrorEquilibration
    {
    public:
      template <typename... Args>
      explicit ErrorEquilibration(GridManagerBase<Grid>& gridManager_, Args...) : gridManager(gridManager_)
      {}

      template <class Err, class ErrorRepresentation, class Scalar>
      void refineGrid_impl(Err const& err, ErrorRepresentation& errorDistribution, Scalar tol)
      {
        std::cout << "ERROR EQUILIBRATION" << std::endl;
        std::sort(errorDistribution.begin(), errorDistribution.end(), AdaptationStrategy_Detail::BiggerThanAbs());

        Scalar n2 = gridManager.grid().size(0);
        Scalar relTol = tol/n2;

        // do refinement
        auto cend = gridManager.grid().leafView().template end<0>();
        size_t counter = 0;
        Scalar minErr = errorDistribution.back().first, maxErr = errorDistribution[0].first;
        std::vector<size_t> mainContributors;
        //constexpr int dim = ErrorRepresentation::Descriptions::Grid::dimension;
        //Dune::FieldVector<Scalar,dim> x0(0.2);

        std::cout << "marking cells..." << std::flush;
        for (auto ci=gridManager.grid().leafView().template begin<0>(); ci!=cend; ++ci) // iterate over cells
        {
          for(size_t i=0; i<errorDistribution.size(); ++i) // iterate over chosen part of the error distribution
            if(gridManager.grid().leafIndexSet().index(*ci) == errorDistribution[i].second)
            {
              if(errorDistribution[i].first > relTol){
                if(errorDistribution[i].first > 4 * relTol) mainContributors.push_back(errorDistribution[i].second);
                ++counter;
                gridManager.mark(1,*ci);
              }
            }
        }
        std::cout << "done." << std::endl;
        std::cout << "totalErrorSquared: " << tol << " -> tolerance: " << relTol << std::endl;
        std::cout << "maxErr: " << maxErr << ", minErr: " << minErr << std::endl;
        std::cout << "refining " << counter << " of " << gridManager.grid().size(0) << " cells" << std::endl;

        std::cout << "adapting..." << std::flush;
        gridManager.adaptAtOnce();
        std::cout << "done." << std::endl;

	/*        std::cout << "second refinement" << std::endl;
        auto cend2 = gridManager.grid().leafView().template end<0>();
        for(auto ci=gridManager.grid().leafView().template begin<0>(); ci!=cend2; ++ci)
        {
          if(ci->level() > 0)
          {
            auto father = ci->father();
            if(ci->level() > father->level())
              for(size_t i : mainContributors) if(i == gridManager.grid().levelIndexSet(gridManager.grid().maxLevel()-1).index(*father)) { gridManager.mark(1,*ci); return; }
          }
        }

        std::cout << "adapting..." << std::flush;
        gridManager.adaptAtOnce();
        std::cout << "done." << std::endl;*/
      }

    private:
      GridManagerBase<Grid>& gridManager;
    };
  } // end namespace Adaptivity
} // end namespace Kaskade

#endif
