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

#ifndef ALGORITHM_3_ERROR_ESTIMATOR_EE
#define ALGORITHM_3_ERROR_ESTIMATOR_EE

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/timer/timer.hpp>
#include <boost/fusion/include/at_c.hpp>

#include "algorithm/newton_bridge.hh"
#include "algorithm/adaptationStrategy.hh"
#include "algorithm/errorDistribution.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "fem/forEach.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/variables.hh"
#include "utilities/enums.hh"
#include "errorEstimationTraits.hh"
#include "lagrangeLinearization.hh"

// forward declarations
namespace Dune
{
  struct InverseOperatorResult;
  template <class,int> class FieldVector;
}

namespace Kaskade
{
  struct Dummy{};

  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration, bool lump=false, int components=1, class ReferenceSolution = Dummy, class ReferenceOperator=Dummy>
  class YetAnotherHBErrorEstimator : public ErrorEstimatorBase<VariableSetDescription,typename ErrorDistribution<NormFunctional,ExtensionVariableSetDescription>::AnsatzVars::VariableSet,ExtensionSpace,RefinementStrategy>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    template <class AnsatzVars, class TestVars, class OriginVars> using LumpedFunctional = Functional<AnsatzVars,TestVars,OriginVars,lump>;
    typedef ErrorEstimationTraits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
    typedef ErrorEstimationTraits<LumpedFunctional,VariableSetDescription,ExtensionVariableSetDescription> LumpedTraits;
    typedef ErrorEstimatorBase<VariableSetDescription,typename ErrorDistribution<NormFunctional,ExtensionVariableSetDescription>::AnsatzVars::VariableSet,ExtensionSpace,RefinementStrategy> Base;
    using Base::squaredError;
    using Base::mde;
    using Base::errorDistribution;
    using Base::extensionSpace;
  public:
    // assemblers

    typedef typename Traits::Scalar Scalar;
    typedef typename VariableSetDescription::Grid Grid;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;
    typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
    typedef typename EnergyError::AnsatzVars::VariableSet ErrorRepresentation;


    YetAnotherHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    : Base(extensionSpace_,fraction,verbose), normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      squaredFraction(fraction*fraction), verbose(verbose_)
      {
        ass_EE.reset(new typename LumpedTraits::Assembler_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new typename Traits::Assembler_HE(variableSetDescription.spaces));
        ass_EH.reset(new typename Traits::Assembler_EH(extensionVariableSetDescription.spaces));
        ass_HH.reset(new typename Traits::Assembler_HH(variableSetDescription.spaces));
      }

    virtual ~YetAnotherHBErrorEstimator(){}

    void operator()(AbstractLinearization const& lin, AbstractFunctionSpaceElement const& x_, AbstractFunctionSpaceElement const& dx_, int step, AbstractFunctionSpaceElement const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      computeError(xl.get());
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed())  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      // state equation -> DLY
      typename Traits::A_HH a_HH(*ass_HH,false);
      typename Traits::A_HE a_HE(*ass_HE,false);
      typename LumpedTraits::A_EE a_EE(*ass_EE,false);
      typename Traits::B_HH b_HH(*ass_HH,false);
      typename Traits::B_HE b_HE(*ass_HE,false);
      typename Traits::B_EH b_EH(*ass_EH,false);
      typename LumpedTraits::B_EE b_EE(*ass_EE,false);

      typename Traits::ExtensionVectorP rhsPE(ass_EE->template rhs<Traits::adjointId,Traits::adjointId+1>());
      rhsPE *= -1.0;
      typename Traits::VectorY dy(at_c<Traits::stateId>(dxl.get().data).coefficients()), solYH = dy; solYH *= 0;
      typename Traits::VectorU du(at_c<Traits::controlId>(dxl.get().data).coefficients());
      typename Traits::ExtensionVectorY solYE(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));

      b_HE.applyscaleadd(-1.0,du,rhsPE);
      a_HE.applyscaleadd(-1.0,dy,rhsPE);

      JacobiPreconditionerForTriplets<Scalar,typename LumpedTraits::ExtensionVectorY,typename LumpedTraits::ExtensionVectorP> jacobi(a_EE);
      jacobi.apply(solYE,rhsPE);


      // error in adjoint equation
      typename Traits::Lyy_HH lyy_HH(*ass_HH,false);
      typename Traits::Lyy_HE lyy_HE(*ass_HE,false);
      typename Traits::Lyy_EH lyy_EH(*ass_EH,false);
      typename LumpedTraits::Lyy_EE lyy_EE(*ass_EE,false);
      TransposedOperator<typename Traits::A_HE> at_EH(a_HE);//(*ass_EH,false);

      typename Traits::VectorY rhsYH(ass_HH->template rhs<Traits::stateId,Traits::stateId+1>());
      rhsYH *= -1.0;
      typename Traits::ExtensionVectorY rhsYE(ass_EE->template rhs<Traits::stateId, Traits::stateId+1>());
      rhsYE *= -1.0;

      lyy_HH.applyscaleadd(-1.0,dy,rhsYH);
      lyy_HE.applyscaleadd(-1.0,dy,rhsYE);
      lyy_EH.applyscaleadd(-1.0,solYE,rhsYH);
      lyy_EE.applyscaleadd(-1.0,solYE,rhsYE);

      typename Traits::ExtensionVectorP solPE(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorP solPH(Traits::VectorP_Initializer::init(variableSetDescription));

      jacobi.apply(solPE,rhsYE);
      at_EH.applyscaleadd(-1.0,solPE,rhsYH);

      boost::timer::cpu_timer mgTimer;
      MultiGridSolver<Grid,components>( a_HH, boost::fusion::at_c<0>(ass_HH->spaces())->gridManager().grid(),
                                        typename MultiGridSolver<Grid,components>::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy), true ).apply(at_c<0>(solPH.data),
                                                                                                                                                at_c<0>(rhsYH.data));
      std::cout << "mg: " << boost::timer::format(mgTimer.elapsed()) << std::endl;

      // variational equality
      TransposedOperator<typename Traits::B_HH> bt_HH(b_HH);
      TransposedOperator<typename Traits::B_EH> bt_HE(b_EH);
      TransposedOperator<typename Traits::B_HE> bt_EH(b_HE);
      TransposedOperator<typename LumpedTraits::B_EE> bt_EE(b_EE);
      typename Traits::Luu_HH luu_HH(*ass_HH,false);
      typename Traits::Luu_EH luu_EH(*ass_EH,false);
      typename Traits::Luu_HE luu_HE(*ass_HE,false);
      typename LumpedTraits::Luu_EE luu_EE(*ass_EE,false);

      typename Traits::VectorU rhsUH(ass_HH->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUH *= -1;
      typename Traits::ExtensionVectorU rhsUE(ass_EE->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUE *= -1;

      luu_HH.applyscaleadd(-1.0,du,rhsUH);
      luu_HE.applyscaleadd(-1.0,du,rhsUE);

      bt_HH.applyscaleadd(-1.0,solPH,rhsUH);
      bt_HE.applyscaleadd(-1.0,solPH,rhsUE);
      bt_EH.applyscaleadd(-1.0,solPE,rhsUH);
      bt_EE.applyscaleadd(-1.0,solPE,rhsUE);

      typename Traits::VectorU solUH(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorU solUE(Traits::ExtensionVectorU_Initializer::init(extensionVariableSetDescription));

      JacobiPreconditionerForTriplets<Scalar,typename LumpedTraits::ExtensionVectorU,typename LumpedTraits::ExtensionVectorU>(luu_EE).apply(solUE,rhsUE);
      luu_EH.applyscaleadd(-1.0,solUE,rhsUH);
      boost::timer::cpu_timer chebTimer;
      ChebyshevPreconditioner<typename Traits::Luu_HH> cheb(luu_HH,chebySteps);
      cheb.initForMassMatrix_TetrahedralQ1Elements();
      cheb.apply(solUH,rhsUH);
      std::cout << "cheb: " << boost::timer::format(chebTimer.elapsed()) << std::endl;

      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);
      solYH *= 0.0; // DLY in state equation -> no error transport into space of linear finite elements
      at_c<Traits::stateId>(errorEstimate_H.data) = at_c<0>(solYH.data);
      at_c<Traits::controlId>(errorEstimate_H.data) = at_c<0>(solUH.data);
      at_c<Traits::stateId>(errorEstimate_E.data) = at_c<0>(solYE.data);
      at_c<Traits::controlId>(errorEstimate_E.data) = at_c<0>(solUE.data);
//      std::string savefilename = createFileName("estimate_E",".vtu",false);
//      std::cout << "Error estimator: writing estimates" << std::endl;
//      writeVTKFile(extensionVariableSetDescription.gridView,errorEstimate_E,savefilename,IoOptions(),2);
//      std::cout << "Error estimator: writing lower order estimates" << std::endl;
//      writeVTKFile(variableSetDescription.gridView,errorEstimate_H,savefilename,IoOptions(),1);
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));

      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
      energyError.considerStateVariable(true);
      energyError.considerControlVariable(true);
      energyError.considerAdjointVariable(false);
      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());

      eeAssembler.assemble(linearization(energyError,xl.get()), EnergyErrorAssembler::RHS);
      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
      mde.reset(new ErrorRepresentation(energyError.getVariableSetDescription()) );
      *mde = distError;

//      if(verbose)
//      {
//        std::string name = "errorDistribution_";
//        name += std::to_string(step);
//        writeVTKFile(mde->descriptions.gridView,*mde,name);
//      }

      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde->data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      squaredError = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, ErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    template <typename... Args>
    void initFunctionals(const Args&... args)
    {
      F_HH.reset(new typename Traits::Functional_HH(args...));
      F_HE.reset(new typename Traits::Functional_HE(args...));
    }
    
    template <typename... Args>
    void initExtensionFunctionals(const Args&... args)
    {
      F_EH.reset(new typename Traits::Functional_EH(args...));
      F_EE.reset(new typename LumpedTraits::Functional_EE(args...));
    }

    void setReference(ReferenceOperator const& Aref, ReferenceSolution const& ref, typename Traits::Vector& refcv)
    {
      referenceSolution = &ref;
      referenceOperator = &Aref;
      referenceCoeffVec = &refcv;
    }

    void computeError(typename Traits::VarSet const& x)
    {
      if(referenceSolution == nullptr) return;
      ReferenceSolution tmp(*referenceSolution);
      interpolateGloballyFromFunctor<PlainAverage>(boost::fusion::at_c<0>(tmp.data),[&x](Cell const& cell, Dune::FieldVector<Scalar,dim> const& xLocal)
      {
        return boost::fusion::at_c<0>(x.data).value(cell.geometry().global(xLocal));
      });

      interpolateGloballyFromFunctor<PlainAverage>(boost::fusion::at_c<1>(tmp.data),[&x](Cell const& cell, Dune::FieldVector<Scalar,dim> const& xLocal)
      {
        return boost::fusion::at_c<1>(x.data).value(cell.geometry().global(xLocal));
      });
//      interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<2>(tmp.data),boost::fusion::at_c<2>(x.data));
      tmp -= *referenceSolution;

      boost::fusion::at_c<0>(referenceCoeffVec->data) = boost::fusion::at_c<0>(tmp.data).coefficients();
      boost::fusion::at_c<1>(referenceCoeffVec->data) = boost::fusion::at_c<1>(tmp.data).coefficients();
      // boost::fusion::at_c<2>(referenceCoeffVec->data) = boost::fusion::at_c<2>(tmp.data).coefficients();

      auto tmp2 = *referenceCoeffVec;
      referenceOperator->apply(*referenceCoeffVec,tmp2);
      std::cout << "REAL ERROR: " << (*referenceCoeffVec)*tmp2 << std::endl;
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HH->assemble(linearization(*F_HH,xl.get()));
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EH->assemble(linearization(*F_EH,xe.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<typename Traits::Functional_HH> F_HH;
    std::unique_ptr<typename Traits::Functional_HE> F_HE;
    std::unique_ptr<typename Traits::Functional_EH> F_EH;
    std::unique_ptr<typename LumpedTraits::Functional_EE> F_EE;
    std::unique_ptr<typename Traits::Assembler_HH> ass_HH;
    //typename Traits::Assembler_HH* ass_HH;
    std::unique_ptr<typename Traits::Assembler_HE> ass_HE;
    std::unique_ptr<typename Traits::Assembler_EH> ass_EH;
    std::unique_ptr<typename LumpedTraits::Assembler_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    Scalar squaredFraction;
    bool verbose;
    size_t chebySteps = 30, mgSteps = 25, mgSmoothingSteps = 20;
    Scalar relativeAccuracy = 1e-3;
    ReferenceSolution const* referenceSolution = nullptr;
    ReferenceOperator const* referenceOperator = nullptr;
    typename Traits::Vector* referenceCoeffVec = nullptr;
    //std::function<Scalar(typename Traits::VarSet)> compareError;
  };

  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration, bool lump=false, int components=1>
  class YetAnotherHBErrorEstimator_Elasticity : public ErrorEstimatorBase<VariableSetDescription,typename ErrorDistribution<NormFunctional,ExtensionVariableSetDescription>::AnsatzVars::VariableSet,ExtensionSpace,RefinementStrategy>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    template <class AnsatzVars, class TestVars, class OriginVars> using LumpedFunctional = Functional<AnsatzVars,TestVars,OriginVars,lump>;
    typedef ErrorEstimationTraits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
    typedef ErrorEstimationTraits<LumpedFunctional,VariableSetDescription,ExtensionVariableSetDescription> LumpedTraits;
    typedef ErrorEstimatorBase<VariableSetDescription,typename ErrorDistribution<NormFunctional,ExtensionVariableSetDescription>::AnsatzVars::VariableSet,ExtensionSpace,RefinementStrategy> Base;
    using Base::squaredError;
    using Base::mde;
    using Base::errorDistribution;
    using Base::extensionSpace;
  public:
    // assemblers

    typedef typename Traits::Scalar Scalar;
    typedef typename VariableSetDescription::Grid Grid;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;
    typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
    typedef typename EnergyError::AnsatzVars::VariableSet ErrorRepresentation;


    YetAnotherHBErrorEstimator_Elasticity(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false)
    : Base(extensionSpace_,fraction,verbose), normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      squaredFraction(fraction*fraction), verbose(verbose_)
      {
        ass_EE.reset(new typename LumpedTraits::Assembler_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new typename Traits::Assembler_HE(variableSetDescription.spaces));
        ass_EH.reset(new typename Traits::Assembler_EH(extensionVariableSetDescription.spaces));
        ass_HH.reset(new typename Traits::Assembler_HH(variableSetDescription.spaces));
      }

    virtual ~YetAnotherHBErrorEstimator_Elasticity(){}

    void operator()(AbstractLinearization const& lin, AbstractFunctionSpaceElement const& x_, AbstractFunctionSpaceElement const& dx_, int step, AbstractFunctionSpaceElement const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed())  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      // state equation -> DLY
      typename Traits::A_HH a_HH(*ass_HH,false);
      typename Traits::A_HE a_HE(*ass_HE,false);
      typename LumpedTraits::A_EE a_EE(*ass_EE,false);
      typename Traits::B_HH b_HH(*ass_HH,false);
      typename Traits::B_HE b_HE(*ass_HE,false);
      typename Traits::B_EH b_EH(*ass_EH,false);
      typename LumpedTraits::B_EE b_EE(*ass_EE,false);

      typename Traits::ExtensionVectorP rhsPE(ass_EE->template rhs<Traits::adjointId,Traits::adjointId+1>());
      rhsPE *= -1.0;
      typename Traits::VectorY dy(at_c<Traits::stateId>(dxl.get().data).coefficients()), solYH = dy; solYH *= 0;
      typename Traits::VectorU du(at_c<Traits::controlId>(dxl.get().data).coefficients());
      typename Traits::ExtensionVectorY solYE(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));

      b_HE.applyscaleadd(-1.0,du,rhsPE);
      a_HE.applyscaleadd(-1.0,dy,rhsPE);

      JacobiPreconditionerForTriplets<Scalar,typename LumpedTraits::ExtensionVectorY,typename LumpedTraits::ExtensionVectorP> jacobi(a_EE);
      jacobi.apply(solYE,rhsPE);


      // error in adjoint equation
      typename Traits::Lyy_HH lyy_HH(*ass_HH,false);
      typename Traits::Lyy_HE lyy_HE(*ass_HE,false);
      typename Traits::Lyy_EH lyy_EH(*ass_EH,false);
      typename LumpedTraits::Lyy_EE lyy_EE(*ass_EE,false);
      typename Traits::Luy_HH luy_HH(*ass_HH,false);
      typename Traits::Luy_HE luy_HE(*ass_HE,false);
      typename Traits::Luy_EH luy_EH(*ass_EH,false);
      typename LumpedTraits::Luy_EE luy_EE(*ass_EE,false);
      TransposedOperator<typename Traits::Luy_HH> lyu_HH(luy_HH);
      TransposedOperator<typename Traits::Luy_EH> lyu_HE(luy_EH);
//      typename Traits::Lyu_HH lyu_HH(*ass_HH,false);
//      typename Traits::Lyu_HE lyu_HE(*ass_HE,false);
      TransposedOperator<typename Traits::A_HE> at_EH(a_HE);//(*ass_EH,false);

      typename Traits::VectorY rhsYH(ass_HH->template rhs<Traits::stateId,Traits::stateId+1>());
      rhsYH *= -1.0;
      typename Traits::ExtensionVectorY rhsYE(ass_EE->template rhs<Traits::stateId, Traits::stateId+1>());
      rhsYE *= -1.0;

      lyy_HH.applyscaleadd(-1.0,dy,rhsYH);
      lyy_HE.applyscaleadd(-1.0,dy,rhsYE);
      lyy_EH.applyscaleadd(-1.0,solYE,rhsYH);
      lyy_EE.applyscaleadd(-1.0,solYE,rhsYE);
      lyu_HH.applyscaleadd(-1.0,du,rhsYH);
      lyu_HE.applyscaleadd(-1.0,du,rhsYE);

      typename Traits::ExtensionVectorP solPE(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorP solPH(Traits::VectorP_Initializer::init(variableSetDescription));

      jacobi.apply(solPE,rhsYE);
      at_EH.applyscaleadd(-1.0,solPE,rhsYH);

      boost::timer::cpu_timer mgTimer;
      MultiGridSolver<Grid,components>( a_HH, boost::fusion::at_c<0>(ass_HH->spaces())->gridManager().grid(),
                                        typename MultiGridSolver<Grid,components>::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy), true ).apply(at_c<0>(solPH.data),
                                                                                                                                                at_c<0>(rhsYH.data));
      std::cout << "mg: " << boost::timer::format(mgTimer.elapsed()) << std::endl;


      // variational equality
      TransposedOperator<typename Traits::B_HH> bt_HH(b_HH);
      TransposedOperator<typename Traits::B_EH> bt_HE(b_EH);
      TransposedOperator<typename Traits::B_HE> bt_EH(b_HE);
      TransposedOperator<typename LumpedTraits::B_EE> bt_EE(b_EE);
      typename Traits::Luu_HH luu_HH(*ass_HH,false);
      typename Traits::Luu_EH luu_EH(*ass_EH,false);
      typename Traits::Luu_HE luu_HE(*ass_HE,false);
      typename LumpedTraits::Luu_EE luu_EE(*ass_EE,false);

      typename Traits::VectorU rhsUH(ass_HH->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUH *= -1;
      typename Traits::ExtensionVectorU rhsUE(ass_EE->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUE *= -1;

      luu_HH.applyscaleadd(-1.0,du,rhsUH);
      luu_HE.applyscaleadd(-1.0,du,rhsUE);

      bt_HH.applyscaleadd(-1.0,solPH,rhsUH);
      bt_HE.applyscaleadd(-1.0,solPH,rhsUE);
      bt_EH.applyscaleadd(-1.0,solPE,rhsUH);
      bt_EE.applyscaleadd(-1.0,solPE,rhsUE);

      luy_HH.applyscaleadd(-1.0,dy,rhsUH);
      luy_HE.applyscaleadd(-1.0,dy,rhsUE);
      luy_EH.applyscaleadd(-1.0,solYE,rhsUH);
      luy_EE.applyscaleadd(-1.0,solYE,rhsUE);

      typename Traits::VectorU solUH(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorU solUE(Traits::ExtensionVectorU_Initializer::init(extensionVariableSetDescription));

      JacobiPreconditionerForTriplets<Scalar,typename LumpedTraits::ExtensionVectorU,typename LumpedTraits::ExtensionVectorU>(luu_EE).apply(solUE,rhsUE);
      luu_EH.applyscaleadd(-1.0,solUE,rhsUH);
      boost::timer::cpu_timer chebTimer;
      ChebyshevPreconditioner<typename Traits::Luu_HH> cheb(luu_HH,chebySteps);
      cheb.initForMassMatrix_TetrahedralQ1Elements();
      cheb.apply(solUH,rhsUH);
      std::cout << "cheb: " << boost::timer::format(chebTimer.elapsed()) << std::endl;

      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);
      solYH *= 0.0; // DLY in state equation -> no error transport into space of linear finite elements
      at_c<Traits::stateId>(errorEstimate_H.data) = at_c<0>(solYH.data);
      at_c<Traits::controlId>(errorEstimate_H.data) = at_c<0>(solUH.data);
      at_c<Traits::stateId>(errorEstimate_E.data) = at_c<0>(solYE.data);
      at_c<Traits::controlId>(errorEstimate_E.data) = at_c<0>(solUE.data);
//      std::string savefilename = createFileName("estimate_E",".vtu",false);
//      std::cout << "Error estimator: writing estimates" << std::endl;
//      writeVTKFile(extensionVariableSetDescription.gridView,errorEstimate_E,savefilename,IoOptions(),2);
//      std::cout << "Error estimator: writing lower order estimates" << std::endl;
//      writeVTKFile(variableSetDescription.gridView,errorEstimate_H,savefilename,IoOptions(),1);
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));

      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
      energyError.considerStateVariable(true);
      energyError.considerControlVariable(true);
      energyError.considerAdjointVariable(false);
      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());

      eeAssembler.assemble(linearization(energyError,xl.get()), EnergyErrorAssembler::RHS);
      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
      mde.reset(new ErrorRepresentation(energyError.getVariableSetDescription()) );
      *mde = distError;

//      if(verbose)
//      {
//        std::string name = "errorDistribution_";
//        name += std::to_string(step);
//        writeVTKFile(mde->descriptions.gridView,*mde,name);
//      }

      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde->data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      squaredError = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, ErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    template <typename... Args>
    void initFunctionals(const Args&... args)
    {
      F_HH.reset(new typename Traits::Functional_HH(args...));
      F_HE.reset(new typename Traits::Functional_HE(args...));
    }

        template <typename... Args>
    void initExtensionFunctionals(const Args&... args)
    {
      F_EH.reset(new typename Traits::Functional_EH(args...));
      F_EE.reset(new typename LumpedTraits::Functional_EE(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HH->assemble(linearization(*F_HH,xl.get()));
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EH->assemble(linearization(*F_EH,xe.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<typename Traits::Functional_HH> F_HH;
    std::unique_ptr<typename Traits::Functional_HE> F_HE;
    std::unique_ptr<typename Traits::Functional_EH> F_EH;
    std::unique_ptr<typename LumpedTraits::Functional_EE> F_EE;
    std::unique_ptr<typename Traits::Assembler_HH> ass_HH;
    //typename Traits::Assembler_HH* ass_HH;
    std::unique_ptr<typename Traits::Assembler_HE> ass_HE;
    std::unique_ptr<typename Traits::Assembler_EH> ass_EH;
    std::unique_ptr<typename LumpedTraits::Assembler_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    Scalar squaredFraction;
    bool verbose;
    size_t chebySteps = 30, mgSteps = 25, mgSmoothingSteps = 20;
    Scalar relativeAccuracy = 1e-3;
  };
}

#endif
