/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef EMBEDDED_ERROREST_HH
#define EMBEDDED_ERROREST_HH

#include <type_traits>
#include "boost/multi_array.hpp"

#include "fem/errorest.hh"
#include "fem/fixfusion.hh"
#include "fem/iterate_grid.hh"
#include "utilities/power.hh"

namespace Kaskade
{
  /// \cond internals
  namespace EmbeddedErrorestimatorDetail {

    using namespace boost::fusion;

    struct GetProjections
    {
      template <class Evaluator>
      typename Evaluator::Space::Mapper::ShapeFunctionSet::Matrix operator()(Evaluator const& evaluator) const
      {
        typename Evaluator::Space::Mapper::ShapeFunctionSet::Matrix p = evaluator.shapeFunctions().hierarchicProjection();
        evaluator.combiner().rightTransform(p);
        evaluator.combiner().leftPseudoInverse(p);
        return p;
      }
    };


    template <class Pairs>
    struct ApplyProjections
    {
      ApplyProjections(Pairs const& data_): data(data_) {}

      template <class Pair>
      void operator()(Pair const& var) const
      {
        int const sid = std::remove_reference<typename boost::fusion::result_of::value_at_c<Pair,0>::type>::type::spaceIndex;

        typedef typename boost::fusion::result_of::value_at_c<Pairs,sid>::type DataPair;
        typename boost::fusion::result_of::value_at_c<DataPair,0>::type evaluator = at_c<0>(at_c<sid>(data));
        typename boost::fusion::result_of::value_at_c<DataPair,1>::type p = at_c<1>(at_c<sid>(data));

        typedef typename std::remove_reference<typename boost::fusion::result_of::value_at_c<Pair,1>::type>::type Function;
        Function& f = at_c<1>(var);
        typename Function::StorageType x(p.N());

        for (int i=0; i<p.N(); ++i) {
          x[i] = 0;
          for (int j=0; j<p.M(); ++j)
            x[i].axpy(p[i][j][0][0],f.coefficients()[evaluator.globalIndices()[j]]);
        }
        for (int i=0; i<p.N(); ++i) {
          f.coefficients()[evaluator.globalIndices()[i]] = x[i];
        }
      }

    private:
      Pairs const& data;
    };

    template <class Sequence>
    ApplyProjections<Sequence> applyProjections(Sequence const& data)
    {
      return ApplyProjections<Sequence>(data);
    }



  }; // end of namespace EmbeddedErrorestimatorDetail
  /// \endcond

  /** \addtogroup adapt */
  /** @{ */
  /**
   * Projects the given FE function onto the the polynomial ansatz
   * subspace of one order lower. This relies on the hierarchic
   * projection defined by the shape function sets on each cell.
   *
   * For continuous spaces and order=1 the result is probably undefined.
   */
  template <class VariableSetDescription>
  void projectHierarchically(VariableSetDescription const& varDesc,
                             typename VariableSetDescription::VariableSet& f)
  {
    using namespace boost::fusion;
    using namespace EmbeddedErrorestimatorDetail;

    // Step through all cells.
    typedef typename VariableSetDescription::Spaces    Spaces;
    typedef typename SpaceType<Spaces,0>::type::Scalar Scalar;
    typedef typename VariableSetDescription::Grid      Grid;
    typedef typename VariableSetDescription::GridView  GridView;

    GridView const& gridView = varDesc.gridView ;

    // Evaluators for shape functions of all FE spaces. Remember that
    // every thread has to use its own collection of evaluators!
    typedef ShapeFunctionCache<Grid,Scalar> SfCache;
    auto evaluators = getEvaluators(varDesc.spaces,static_cast<SfCache*>(nullptr));
    using Evaluators = decltype(evaluators);
    

    // A sequence of projection matrices.
    typedef typename boost::fusion::result_of::as_vector<
                  typename boost::fusion::result_of::transform<Evaluators,GetProjections>::type
             >::type Projections;
    Projections projections;

    // Step through all cells. On each cell, compute K^+ P K a, where a
    // is the coefficient vector, K the combiner (often the identity),
    // and P the projection matrix that projects to a subspace of lower
    // order shape functions.
    for (auto ci=gridView.template begin<0,Dune::All_Partition>();
         ci!=gridView.template end<0,Dune::All_Partition>(); ++ci) {
      // for all spaces involved, compute the shape functions and
      // their global indices, which are needed for evaluating the
      // functional's derivative.
      moveEvaluatorsToCell(evaluators,*ci);

      // Compute projections.
      projections = transform(evaluators,GetProjections());

      // Apply projections.
      typename VariableSetDescription::Variables vars;
      for_each2(vars,f.data,
                applyProjections(zip(evaluators,projections)));
    }
  }



  //---------------------------------------------------------------------
  // Embedded error estimation.
  /// \cond internals
  namespace EmbeddedErrorestDetail {

    class CellByCell
    {
    public:
      template <class IndexSet>
      CellByCell(IndexSet const& is): nGroups(is.size(0)) {}
      
      size_t operator[](size_t idx) const { assert(idx<nGroups); return idx; }
      
      size_t nGroups;
    };

  } // End of namespace EmbeddedErrorestDetail
  /// \endcond
  
  /**
   * \brief Embedded error estimation and mesh refinement. 
   * 
   * Given an approximate solution \f$ u \f$ and an explicit error approximation representation \f$ \delta u \f$
   * from the same finite element spaces, this error estimator computes
   * \f[ \epsilon_i = \int_\Omega (s_u(u)\delta u)^2 \, dx, \quad \eta_i = \int_\Omega s(u)^2 \, dx. \f]
   * If \f$ \epsilon_i > \mathrm{aTOL}_i^2 + \mathrm{rTOL}_i^2 \eta_i \f$ for any \f$ i \f$, the mesh is refined.
   * 
   * The cells with largest contribution to \f$ \epsilon_i \f$ are marked for refinement.
   */
  template <class VariableSetDescription, class Scaling=IdentityScaling>
  class EmbeddedErrorEstimator
  {
  public:
    typedef GridManager<typename VariableSetDescription::Grid> GManager;
    typedef EmbeddedErrorEstimator<VariableSetDescription,Scaling> Self;
    
    /**
     * \brief Constructor
     * 
     * Both the grid manager and the variable set description objects need to exist until the 
     * error estimator is destructed.
     */
    EmbeddedErrorEstimator(GManager& gridManager_, VariableSetDescription const& varDesc_, Scaling const& scaling_=Scaling())
    : gridManager(gridManager_), varDesc(varDesc_), tol(varDesc.noOfVariables,std::make_pair(0.0,0.01)), scaling(scaling_) {}
    
    /**
     * \brief perform error estimation
     */
    bool estimate(typename VariableSetDescription::VariableSet const& err,
                  typename VariableSetDescription::VariableSet const& sol) const
    {
      using namespace boost::fusion;
      using namespace EmbeddedErrorestDetail;
      
      assert(size(err.data)==size(sol.data));

// std::cerr << "emb. error. est\n";      
      ErrorestDetail::GroupedSummationCollector<CellByCell> sum(CellByCell(varDesc.indexSet));
      scaledTwoNormSquared(join(typename VariableSetDescription::Variables(),typename VariableSetDescription::Variables()),
                           join(err.data,sol.data),varDesc.spaces,scaling,sum); 

// std::cerr << "Exit in " << __FILE__ << ":" << __LINE__ << "\n"; abort();
      
      
      int const s = varDesc.noOfVariables;
      int const n = sum.sums.shape()[0];
      
      assert(sum.sums.shape()[1]==2*s);
      
      // Compute scaled L2 norms of solution and errors.
      std::vector<double> norm2(s), error2(s);
      for (int idx=0; idx<n; ++idx)
        for (int j=0; j<s; ++j) {
          error2[j] += sum.sums[idx][j];
          norm2[j]  += sum.sums[idx][j+s];
        }
        
        
      if (verbosity>0)
      {
        std::cout << "   ||solution||^2 = " << std::scientific << std::setprecision(3); std::copy(norm2.begin(),norm2.end(),std::ostream_iterator<double>(std::cout," "));
        std::cout << std::endl;
        std::cout << "     ||errors||^2 = "; std::copy(error2.begin(),error2.end(),std::ostream_iterator<double>(std::cout," "));
        std::cout << std::resetiosflags( ::std::ios::scientific ) << std::setprecision(6) << std::endl;
      };
        
      // Compute error limit and normalized errors
      std::vector<double> errorLimit(s);
      boost::multi_array<double,2> normalizedErrors(boost::extents[n][s]);
        
      for (int j=0; j<s; ++j) 
      {
        errorLimit[j] = power<2>(tol[j].first) + power<2>(tol[j].second)*norm2[j];
        for (int i=0; i<n; ++i) normalizedErrors[i][j] = sum.sums[i][j] / errorLimit[j];
      }
        
      // Get refinement thresholds
      FixedFractionCriterion marker(0.1);
      auto threshold = marker.threshold(normalizedErrors);
        
      // Check for refinement necessity
      if (threshold.empty()) // all variables accurate enough
        return true;
      
      // Refine the mesh
      markAndRefine(gridManager,varDesc.gridView,normalizedErrors,threshold);
      return false;
    }
    
    /**
     * \brief set the tolerances
     * 
     * The tolerance for each variable consists of an absolute part (aTOL) and a relative part (rTOL).
     * Sufficient accuracy is assumed, if the error \f$ \epsilon \f$ of the approximate solution \f$ u \f$ 
     * satisfies \f[ \|\epsilon\|^2 \le \mathrm{aTOL}^2 + \mathrm{rTOL}^2 \|u\|^2. \f]
     * 
     * Defaults to aTOL=0 and rTOL=0.01 for all variables.
     */
    Self& setTolerances(std::vector<std::pair<double,double> > const& tol_)
    {
      tol = tol_;
      return *this;
    }
    
    /**
     * \brief set the scaling
     */
    Self& setScaling(Scaling const& scaling_)
    {
      scaling = scaling_;
      return *this;
    }
    
    /**
     * \brief set the refinement criterion
     */
    Self& setCriterion(RefinementCriterion const& criterion_)
    {
      criterion = &criterion_;
      return *this;
    }
    
    /**
     * \brief set the verbosity level 
     * 
     * Verbosity levels:
     * - 0: no output
     * - 1: squared norms of error and solution are written to standard output
     */
    Self& setVerbosity(int verbosity_)
    {
      verbosity = verbosity_;
      return *this;
    }
    
    
  private:
    GManager&                              gridManager;
    VariableSetDescription const&          varDesc;
    std::vector<std::pair<double,double> > tol;
    Scaling                                scaling;
    RefinementCriterion const*             criterion;
    int                                    verbosity;
  };

  /**
   * \brief Embedded error estimation and mesh refinement. 
   * 
   * Returns true if the error is below the tolerance, in which case no mesh refinement is
   * performed.
   * 
   * \param err an explicit representation of the estimated error
   * \param sol the (approximate) solution
   * \param scaling A pointwise nonlinear scaling to be applied before error estimation.
   */
  template <class VariableSetDescription, class Scaling>
  bool embeddedErrorEstimator(VariableSetDescription const& varDesc,
                              typename VariableSetDescription::VariableSet const& err,
                              typename VariableSetDescription::VariableSet const& sol,
                              Scaling const& scaling,
                              std::vector<std::pair<double,double> > const& tol,
                              GridManager<typename VariableSetDescription::Grid>& gridManager,
                              int verbosity=1)
  {
    EmbeddedErrorEstimator<VariableSetDescription,Scaling> estimator(gridManager,varDesc,scaling);
    FixedFractionCriterion marker(0.1);
    estimator.setTolerances(tol).setCriterion(marker).setVerbosity(verbosity);
    return estimator.estimate(err,sol);
  }
} /* end of namespace Kaskade */

/** @} */



#endif
