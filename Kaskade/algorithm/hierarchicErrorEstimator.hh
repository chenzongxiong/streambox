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

#ifndef ALGORITHM_HIERARCHIC_ERROR_ESTIMATOR_HH
#define ALGORITHM_HIERARCHIC_ERROR_ESTIMATOR_HH

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/timer/timer.hpp>
#include <boost/fusion/include/at_c.hpp>

#include "algorithm/newton_bridge.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/variables.hh"
#include "linalg/direct.hh"
#include "linalg/tcg.hh"
#include "utilities/enums.hh"

// forward declarations
namespace Dune
{
  struct InverseOperatorResult;
  template <class,int> class FieldVector;
}

namespace Kaskade
{
  namespace HierarchicErrorEstimator_Detail{
    bool biggerThanAbs(const std::pair<double,int>& x1, const std::pair<double,int>&  x2)
    {
      return fabs(x1.first)>fabs(x2.first);
    }

    double add(const double& x1, const std::pair<double,int>&  x2)
    {
      return x1+x2.first;
    }
  }

  template <class CorrectionVector>
  struct LeaveRHS
  {
    explicit LeaveRHS(CorrectionVector const&) {}

    template <class RhsVector>
    void operator()(RhsVector& rhs)
    {}    
  };

  template <class CorrectionVector>
  struct SubstractCorrection
  {
    explicit SubstractCorrection(CorrectionVector const& cor_) 
    : cor(cor_)
    {}

    template <class RhsVector>
    void operator()(RhsVector& rhs)
    {
      rhs -= cor;
    }

  private:
    CorrectionVector const& cor;
  };


  template <class Functional, /*class LF,*/ class ExtensionVariableSetDescription, class ExtensionSpace, class NormFunctional=Functional, template <class> class AdjustRHS=LeaveRHS>
  class HierarchicalBasisErrorEstimator : public AbstractHierarchicalErrorEstimator
  {
    static constexpr int yId = 1, uId = 0, pId = 2;
    static constexpr int noOfVariables = ExtensionVariableSetDescription::noOfVariables;
    //    typedef typename LF::AnsatzVars LVars;
  public:
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::AnsatzVars::VariableSet VariableSet;
    static constexpr int dim = VariableSet::Descriptions::Grid::dimension;
    typedef HierarchicErrorEstimator<LinearizationAt<Functional>,ExtensionVariableSetDescription,ExtensionVariableSetDescription,HierarchicErrorEstimatorDetail::TakeAllD2<LinearizationAt<Functional> > > ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> Assembler;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<0,2>::type CoefficientVector02;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<0,noOfVariables>::type CoefficientVector;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;

    //    HierarchicalBasisErrorEstimator(Functional& f_, LF& lf_, LVars& lvars, ExtensionVariableSetDescription& extensionVariableSetDescription_, ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    //    : f(f_), lf(lf_), normFunctional(f_), extensionVariableSetDescription(extensionVariableSetDescription_), extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
    //    squaredFraction(fraction*fraction), verbose(verbose_)
    //    {}

    HierarchicalBasisErrorEstimator(Functional& f_, /*LF& lf_,*/ /*LVars& lvars_,*/ NormFunctional& normFunctional_, ExtensionVariableSetDescription& extensionVariableSetDescription_, 
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    : f(f_), /*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), extensionVariableSetDescription(extensionVariableSetDescription_), extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
      squaredFraction(fraction*fraction), verbose(verbose_)
    {}

    virtual ~HierarchicalBasisErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const& lowerOrderRhs)
    {
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<VariableSet> const& x = dynamic_cast<const Bridge::Vector<VariableSet>&>(x_);
      Bridge::Vector<VariableSet> const& dx = dynamic_cast<const Bridge::Vector<VariableSet>&>(dx_);
      //    Bridge::Vector<VariableSet> const& loRhs = dynamic_cast<const Bridge::Vector<VariableSet>&>(lowerOrderRhs);


      boost::timer::cpu_timer atimer;
      assembler.assemble(ErrorEstimator(LinearizationAt<Functional>(f,x.get()),dx.get()));
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << "s" << std::endl;

      typedef AssembledGalerkinOperator<Assembler> Operator;
      typedef AssembledGalerkinOperator<Assembler,0,2,0,2> CorrectionOperator;
      typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<0,2>::type CorrectionVector;

      Operator op(assembler);
      //      CorrectionOperator(assembler);

      CoefficientVector estRhs(assembler.rhs());
      std::vector<Scalar> tmpRepVec;
      IstlInterfaceDetail::toVector(estRhs,tmpRepVec);
      //estRhs = Bridge::getImpl<VariableSet>(lowerOrderRhs);
      std::cout << "rhs entries:" << std::endl;
      Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
      for(size_t i=0; i<tmpRepVec.size(); ++i){
        if(std::isnan(tmpRepVec[i])) std::cout << "nan at " << i << std::endl;

        Scalar val = tmpRepVec[i];
        if(val < minRhs) minRhs = val;
        if(val > maxRhs) maxRhs = val;
        rhsSum += val*val;
      }
      std::cout << "max rhs: " << maxRhs << std::endl;
      std::cout << "min rhs: " << minRhs << std::endl;
      std::cout << "sum: " << rhsSum << std::endl;
      // possible adjustment of rhs
      //adjustRHS(estRhs);

      boost::timer::cpu_timer timer;
      CoefficientVector estSol(ExtensionVariableSetDescription::template CoefficientVectorRepresentation<0,noOfVariables>::init(extensionVariableSetDescription));

      directInverseOperator(op,DirectType::UMFPACK,MatrixProperties::GENERAL).apply(estRhs,estSol);


      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(timer.elapsed(), 2, "%t") << "s" << std::endl;
      typename ExtensionVariableSetDescription::VariableSet mySol(extensionVariableSetDescription);
      mySol = estSol;
      using namespace boost::fusion;

      Scalar maxVal = std::numeric_limits<Scalar>::min(),
          minVal = std::numeric_limits<Scalar>::max(),
          sum = 0;

      if(verbose)
      {
        for(int i=0; i<at_c<yId>(mySol.data).size();  ++i)
        {
          for(int j=0; j<at_c<yId>(mySol.data).coefficients()[i].size; ++j)
          {
            Scalar val = at_c<yId>(mySol.data).coefficients()[i][j];
            sum += val*val;
            if(val > maxVal) maxVal = val;
            if(val < minVal) minVal = val;
          }
        }
        std::cout << "y: minimal value: " << minVal << std::endl;
        std::cout << "y: maximal value: " << maxVal << std::endl;
        std::cout << "y: l2-sum: " << sum << std::endl;
        Scalar tmpMax = maxVal, tmpMin = minVal, tmpSum = sum;
        maxVal = std::numeric_limits<Scalar>::min();
        minVal = std::numeric_limits<Scalar>::max();
        sum = 0;

        for(int i=0; i<at_c<uId>(mySol.data).size();  ++i)
        {
          for(int j=0; j<at_c<uId>(mySol.data).coefficients()[i].size; ++j)
          {
            Scalar val = at_c<uId>(mySol.data).coefficients()[i][0];
            sum += val*val;
            if(val > maxVal) maxVal = val;
            if(val < minVal) minVal = val;
          }
        }
        std::cout << "u: minimal value: " << minVal << std::endl;
        std::cout << "u: maximal value: " << maxVal << std::endl;
        std::cout << "u: l2-sum: " << sum << std::endl;
        tmpMax = std::max(tmpMax,maxVal), tmpMin = std::min(tmpMin,minVal), tmpSum += sum;
        maxVal = std::numeric_limits<Scalar>::min();
        minVal = std::numeric_limits<Scalar>::max();
        sum = 0;

        for(int i=0; i<at_c<pId>(mySol.data).size();  ++i)
        {
          for(int j=0; j<at_c<pId>(mySol.data).coefficients()[i].size; ++j)
          {
            Scalar val = at_c<pId>(mySol.data).coefficients()[i][j];
            sum += val*val;
            if(val > maxVal) maxVal = val;
            if(val < minVal) minVal = val;
          }
        }
        std::cout << "p: minimal value: " << minVal << std::endl;
        std::cout << "p: maximal value: " << maxVal << std::endl;
        std::cout << "p: l2-sum: " << sum << std::endl;

        tmpSum += sum;
        std::cout << "minimal value: " << std::min(tmpMin,minVal) << std::endl;
        std::cout << "maximal value: " << std::max(tmpMax,maxVal) << std::endl;
        std::cout << "l2-sum: " << tmpSum << std::endl;
      }
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      EnergyError energyError(normFunctional,x.get(),mySol);
      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());

      eeAssembler.assemble(linearization(energyError,x.get()), EnergyErrorAssembler::RHS);
      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
      typename EnergyError::AnsatzVars::VariableSet mde( energyError.getVariableSetDescription() );
      mde = distError;

      if(verbose)
      {
        std::string name = "errorDistribution_";
        name += boost::lexical_cast<std::string>(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);
    }

    void refineGrid()
    {
      //      int lastIndexForDoubleRefinement = -1;
      std::sort(errorDistribution.begin(), errorDistribution.end(), HierarchicErrorEstimator_Detail::biggerThanAbs);

      Scalar bulkCriterionTolerance = squaredFraction*totalErrorSquared;
      if(verbose)
      {
        std::cout << "ERROR ESTIMATOR: totalErrorSquared: " << totalErrorSquared << std::endl;
        std::cout << "ERROR ESTIMATOR: bulkCriterionTolerance: " << bulkCriterionTolerance << std::endl;
      }
      Scalar bulkErrorSquared = 0;

      std::vector<std::pair<Scalar,size_t> > bulkErrorDistribution;
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();

      while(bulkErrorSquared < bulkCriterionTolerance)
      {
        bulkErrorDistribution.push_back(errorDistribution[bulkErrorDistribution.size()]);
        bulkErrorSquared += bulkErrorDistribution.back().first;
      }

      if(verbose) std::cout << "ERROR ESTIMATOR: number of candidates for refinement: " << bulkErrorDistribution.size() << std::endl;

      //      size_t lastIndex = bulkErrorDistribution.size()-1;
      //      size_t firstIndex = 0;

      /*while(firstIndex+3 < lastIndex)
      {
        if(bulkErrorDistribution.size() > 4)
        {
          Scalar firstContribution = bulkErrorDistribution[firstIndex].first;
          firstContribution *= 15.0/256.0;

          Scalar lastContributions = bulkErrorDistribution[lastIndex--].first;
          lastContributions += bulkErrorDistribution[lastIndex--].first;
          lastContributions += bulkErrorDistribution[lastIndex].first;
          lastContributions *= 15.0/16.0;

          if(lastContributions < firstContribution)
          {
            ++lastIndexForDoubleRefinement;
            ++firstIndex;
            --lastIndex;
            bulkErrorDistribution.erase(bulkErrorDistribution.end()-3, bulkErrorDistribution.end());
          }
          else break;
        }
      }*/


      // Refine mesh.
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci) // iterate over cells
      {
        for(int i=0; i<bulkErrorDistribution.size(); ++i) // iterate over chosen part of the error distribution
          if(is.index(*ci) == bulkErrorDistribution[i].second)
          {
            //            if(i <= lastIndexForDoubleRefinement) {
            //              extensionSpace.gridManager().mark(2,*ci);
            //            }
            //            else{
            extensionSpace.gridManager().mark(1,*ci);
            //            }
          }
      }

      extensionSpace.gridManager().adaptAtOnce();
    }

    double estimatedAbsoluteError() const final
        {
      return sqrt(fabs(totalErrorSquared));
        }

    size_t gridSize() const final
        {
      return extensionSpace.gridManager().grid().size(0);
        }

  private:
    Functional& f;
    //    LF& lf;
    //    LVars& lvars;
    NormFunctional& normFunctional;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Assembler assembler;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose;
    // AdjustRHS<CoefficientVector> adjustRHS;
  };
}

#endif
