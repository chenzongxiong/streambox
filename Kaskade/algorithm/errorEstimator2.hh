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

#ifndef ALGORITHM_3_ERROR_ESTIMATOR_HH
#define ALGORITHM_3_ERROR_ESTIMATOR_HH

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
#include "linalg/blockDiagonalSchurPreconditioner.hh"
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
  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
            class NormFunctional, class LinearSolverLA, class LinearSolverHA, class LinearSolverLM=LinearSolverLA, class LinearSolverHM=LinearSolverHA,
            bool lumpM=false>
  class HierarchicalBasisErrorEstimator2 : public AbstractHierarchicalErrorEstimator
  {
    static constexpr int noOfVariables = ExtensionVariableSetDescription::noOfVariables;
//    typedef typename LF::AnsatzVars LVars;
    static constexpr int yId = 1, uId = 0, pId = 2;
  public:
    typedef Functional<VariableSetDescription,VariableSetDescription,VariableSetDescription,lumpM>                             FLL;
    typedef Functional<ExtensionVariableSetDescription,VariableSetDescription,ExtensionVariableSetDescription,lumpM>         FHL;
    typedef Functional<VariableSetDescription,ExtensionVariableSetDescription,VariableSetDescription,lumpM>                    FLH;
    typedef Functional<ExtensionVariableSetDescription,ExtensionVariableSetDescription,ExtensionVariableSetDescription,lumpM>  FHH;

    typedef HierarchicErrorEstimator<LinearizationAt<FLL>,ExtensionVariableSetDescription,ExtensionVariableSetDescription,HierarchicErrorEstimatorDetail::TakeAllD2<LinearizationAt<FLL> > > ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> Assembler;

    typedef VariationalFunctionalAssembler<LinearizationAt<FLL> >                                 Assll;
    typedef VariationalFunctionalAssembler<LinearizationAt<FLH> >                                 Asslh;
    typedef VariationalFunctionalAssembler<LinearizationAt<FHL> >                                 Asshl;
    typedef VariationalFunctionalAssembler<LinearizationAt<FHH> >                                 Asshh;

    typedef AssembledGalerkinOperator<Asshh,pId,pId+1,yId,yId+1>                                              Ahh;
    typedef AssembledGalerkinOperator<Assll,pId,pId+1,yId,yId+1>                                              All;
    typedef AssembledGalerkinOperator<Asslh,pId,pId+1,yId,yId+1>                                              Alh;
    typedef AssembledGalerkinOperator<Asshl,pId,pId+1,yId,yId+1>                                              Ahl;
    typedef AssembledGalerkinOperator<Asshh,yId,yId+1,yId,yId+1>                                              Myyhh;

    typedef AssembledGalerkinOperator<Assll,pId,pId+1,uId,uId+1>                                              Bll;
    typedef AssembledGalerkinOperator<Asslh,pId,pId+1,uId,uId+1>                                              Blh;
    typedef AssembledGalerkinOperator<Asshl,pId,pId+1,uId,uId+1>                                              Bhl;
    typedef AssembledGalerkinOperator<Asshh,pId,pId+1,uId,uId+1>                                              Bhh;

    typedef AssembledGalerkinOperator<Assll,uId,uId+1,uId,uId+1>                                             Muull;
    typedef AssembledGalerkinOperator<Asslh,uId,uId+1,uId,uId+1>                                             Muulh;
    typedef AssembledGalerkinOperator<Asshl,uId,uId+1,uId,uId+1>                                             Muuhl;
    typedef AssembledGalerkinOperator<Asshh,uId,uId+1,uId,uId+1>                                             Muuhh;

    typedef typename FLL::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;
    typedef VariableSet<VariableSetDescription> VarSet;
    typedef VariableSet<ExtensionVariableSetDescription> ExtVarSet;
    typedef typename VariableSetDescription::template CoefficientVectorRepresentation<> CVL;
    typedef typename VariableSetDescription::template CoefficientVectorRepresentation<yId,yId+1> CVLY;
    typedef typename VariableSetDescription::template CoefficientVectorRepresentation<uId,uId+1> CVLU;
    typedef typename VariableSetDescription::template CoefficientVectorRepresentation<pId,pId+1> CVLP;
    typedef typename CVL::type CoefficientVectorL;
    typedef typename CVLY::type CoefficientVectorLY;
    typedef typename CVLU::type CoefficientVectorLU;
    typedef typename CVLP::type CoefficientVectorLP;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<0,1>::type CoefficientVectorH01;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<0,2>::type CoefficientVectorH02;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<2,3> CVH23;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<2,3>::type CoefficientVectorH23;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<> CVH;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<yId,yId+1> CVHY;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<uId,uId+1> CVHU;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<pId,pId+1> CVHP;
    typedef typename CVH::type CoefficientVectorH;
    typedef typename CVHY::type CoefficientVectorHY;
    typedef typename CVHU::type CoefficientVectorHU;
    typedef typename CVHP::type CoefficientVectorHP;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;

    typedef MatrixAsTriplet<Scalar> Matrix;

    //    HierarchicalBasisErrorEstimator(Functional& f_, LF& lf_, LVars& lvars, ExtensionVariableSetDescription& extensionVariableSetDescription_, ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    //    : f(f_), lf(lf_), normFunctional(f_), extensionVariableSetDescription(extensionVariableSetDescription_), extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
    //    squaredFraction(fraction*fraction), verbose(verbose_)
    //    {}

    HierarchicalBasisErrorEstimator2(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
    ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    : /*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
      squaredFraction(fraction*fraction), verbose(verbose_)
    {
      initAssemblers();
    }

    virtual ~HierarchicalBasisErrorEstimator2(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<VarSet> const& xl = dynamic_cast<const Bridge::Vector<VarSet>&>(x_);
      Bridge::Vector<ExtVarSet> xe(extensionVariableSetDescription);
      assembleAt(xl,xe);
      Bridge::Vector<VarSet> const& dxl = dynamic_cast<const Bridge::Vector<VarSet>&>(dx_);

      // estimate error in dual variable
      typename ExtensionVariableSetDescription::VariableSet tmpRep(extensionVariableSetDescription), mySol(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet tmpRepL(variableSetDescription);

      All all(*assll,false);
      Ahh ahh(*asshh,false);
      Alh alh(*asslh,false);

      Matrix mAll = all.template get<Matrix>();
      mAll.transpose();
      Matrix mAhh = ahh.template get<Matrix>();
      mAhh.transpose();
      Matrix mAlh = alh.template get<Matrix>();
      mAlh.transpose();

      tmpRep *= 0;
      at_c<yId>(tmpRep.data) = at_c<0>(rhs0.data);
//      //MyH1SemiNorm myNorm;
//      //std::cout << "||rhs||=" << myNorm(tmpRep) << std::endl;
//      LinearSolverHA(stateOp,DirectType::MUMPS,MatrixProperties::GENERAL).apply(sol0,rhs0);
//      //directInverseOperator(stateOp,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhs0,sol0);
//      tmpRep *= 0;
      at_c<pId>(tmpRep.data) = at_c<0>(sol0.data);
      at_c<pId>(mySol.data) = at_c<0>(sol0.data);
      //std::cout << "||sol||=" << myNorm(tmpRep) << std::endl;

      // estimate error in control variable
      // rhs vectors
      CoefficientVectorLU rhs1L(asshl->template rhs<uId,uId+1>());
      CoefficientVectorHU rhs1H(asshh->template rhs<uId,uId+1>());
      Bhl bhl(*asshl,false);
      Blh blh(*asslh,false);
      Bll bll(*assll,false);
      Bhh bhh(*asshh,false);

      MatrixAsTriplet<Scalar> tmpmat = blh.template get<MatrixAsTriplet<Scalar> >();
      std::vector<Scalar> tmpx, tmpy;
      IstlInterfaceDetail::toVector(sol0,tmpy);
      IstlInterfaceDetail::toVector(rhs1L,tmpx);
      tmpmat.usmtv(-1.0,tmpy,tmpx);
      IstlInterfaceDetail::fromVector(tmpx,rhs1L);
      //std::cout << "rhs1L=" << at_c<0>(rhs1L.data)[0] << std::endl;

      tmpx.clear(); tmpy.clear();
      IstlInterfaceDetail::toVector(sol0,tmpy);
      IstlInterfaceDetail::toVector(rhs1H,tmpx);
      tmpmat = bhh.template get<MatrixAsTriplet<Scalar> >();
      tmpmat.usmtv(-1.0,tmpy,tmpx);
      IstlInterfaceDetail::fromVector(tmpx,rhs1H);
      //        bophh.applyscaleadd(-1.0,sol0,rhs1H);
      //std::cout << "rhs1H=" << at_c<0>(rhs1H.data)[0] << std::endl;

      Muuhh muuhh(*asshh,true);
      Muull muull(*assll,true);
      Muulh muulh(*asslh,true);
      CoefficientVectorHU solhu(CVHU::init(extensionVariableSetDescription));
      CoefficientVectorLU sollu(CVLU::init(variableSetDescription));
      LinearSolverLM(muull,DirectType::MUMPS,MatrixProperties::GENERAL).apply(sollu,rhs1L);
      //directInverseOperator(muull,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhs1L,sollu);
      muulh.applyscaleadd(-1.0,sollu,rhs1H);
      LinearSolverHM(muuhh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(solhu,rhs1H);
      //      directInverseOperator(muuhh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhs1H,solhu);

      //std::cout << "sollu=" << at_c<0>(sollu.data)[0] << std::endl;
      //std::cout << "solhu=" << at_c<0>(solhu.data)[0] << std::endl;
      at_c<uId>(mySol.data) = at_c<0>(solhu.data);

      // estimate error in state variable
      CoefficientVectorLY solly(CVLP::init(variableSetDescription));
      CoefficientVectorHY solhy(CVHP::init(extensionVariableSetDescription));
      CoefficientVectorLP rhslp(assll->template rhs<pId,pId+1>());
      CoefficientVectorHP rhshp(asshh->template rhs<pId,pId+1>());

      bll.applyscaleadd(-1.0,sollu,rhslp);
      blh.applyscaleadd(-1.0,sollu,rhshp);
      bhl.applyscaleadd(-1.0,solhu,rhslp);
      bhh.applyscaleadd(-1.0,solhu,rhshp);

      LinearSolverLA(all,DirectType::MUMPS,MatrixProperties::GENERAL).apply(solly,rhslp);
      //directInverseOperator(all,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhslp,solly);
      alh.applyscaleadd(-1.0,solly,rhshp);
      LinearSolverHA(ahh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(solhy,rhshp);
      //      directInverseOperator(ahh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhshp,solhy);
      at_c<yId>(mySol.data) = at_c<0>(solhy.data);

      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      EnergyError energyError(normFunctional,xl.get(),mySol);
      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());

      eeAssembler.assemble(linearization(energyError,xl.get()), EnergyErrorAssembler::RHS);
      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
      typename EnergyError::AnsatzVars::VariableSet mde( energyError.getVariableSetDescription() );
      mde = distError;
      //std::cout << "||mde||=" << myNorm(mde) << std::endl;

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

    template <typename... Args>
    void initFunctionals(const Args&... args)
    {
      Fll.reset(new FLL(args...));
      Flh.reset(new FLH(args...));
      Fhl.reset(new FHL(args...));
      Fhh.reset(new FHH(args...));
    }

  private:
    void initAssemblers()
    {
      assll.reset(new Assll(extensionVariableSetDescription.spaces));
      asslh.reset(new Asslh(extensionVariableSetDescription.spaces));
      asshl.reset(new Asshl(extensionVariableSetDescription.spaces));
      asshh.reset(new Asshh(extensionVariableSetDescription.spaces));
    }

    void assembleAt(const Bridge::Vector<VarSet>& xl, const Bridge::Vector<ExtVarSet>& xe)
    {
      assll->assemble(linearization(*Fll,xl.get()));
      asslh->assemble(linearization(*Flh,xl.get()));
      asshl->assemble(linearization(*Fhl,xe.get()));
      asshh->assemble(linearization(*Fhh,xe.get()));
    }

    std::unique_ptr<FLL> Fll;
    std::unique_ptr<FLH> Flh;
    std::unique_ptr<FHL> Fhl;
    std::unique_ptr<FHH> Fhh;
    std::unique_ptr<Assll> assll;
    std::unique_ptr<Asslh> asslh;
    std::unique_ptr<Asshl> asshl;
    std::unique_ptr<Asshh> asshh;
//    LF& lf;
//    LVars& lvars;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
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
