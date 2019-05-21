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
//#include "algorithm/errorDistribution.hh"
//#include "algorithm/mynorm.hh"  // file does not exist
//#include "algorithm/errorDistribution3.hh"
#include "linalg/jacobiPreconditioner.hh"
// #include "linalg/minresSolver.hh"  // file does not exist
#include "fem/forEach.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/variables.hh"
#include "linalg/direct.hh"
//#include "linalg/blockDiagonalSchurPreconditioner.hh"
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
  namespace ErrorEstimator_Detail
  {

    template <class Functional>
    std::integral_constant<int,Functional::stateId> getStateIdImpl(decltype(Functional::stateId)*);

    template <class Functional>
    std::integral_constant<int,Functional::yIdx> getStateIdImpl(decltype(Functional::yIdx)*);

    template <class Functional>
    std::integral_constant<int,Functional::controlId> getControlIdImpl(decltype(Functional::controlId)*);

    template <class Functional>
    std::integral_constant<int,Functional::uIdx> getControlIdImpl(decltype(Functional::uIdx)*);

    template <class Functional>
    std::integral_constant<int,Functional::adjointId> getAdjointIdImpl(decltype(Functional::adjointId)*);

    template <class Functional>
    std::integral_constant<int,Functional::pIdx> getAdjointIdImpl(decltype(Functional::pIdx)*);

    template <class Functional>
    std::integral_constant<int,Functional::lIdx> getAdjointIdImpl(decltype(Functional::lIdx)*);

    template <class Functional>
    constexpr int getStateId()
    {
      typedef decltype(getStateIdImpl<Functional>(nullptr)) type;
      return type::value;
    }

    template <class Functional>
    constexpr int getControlId()
    {
      typedef decltype(getControlIdImpl<Functional>(nullptr)) type;
      return type::value;
    }

    template <class Functional>
    constexpr int getAdjointId()
    {
      typedef decltype(getAdjointIdImpl<Functional>(nullptr)) type;
      return type::value;
    }

    template < template <class,class,class> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription>
    struct Traits
    {
      typedef VariableSet<VariableSetDescription> VarSet;
      typedef VariableSet<ExtensionVariableSetDescription> ExtensionVarSet;

      typedef Functional<VariableSetDescription,VariableSetDescription,VariableSetDescription>                               Functional_HH;
      typedef Functional<VariableSetDescription,ExtensionVariableSetDescription,VariableSetDescription>                     Functional_HE;
      typedef Functional<ExtensionVariableSetDescription,VariableSetDescription,ExtensionVariableSetDescription>           Functional_EH;
      typedef Functional<ExtensionVariableSetDescription,ExtensionVariableSetDescription,ExtensionVariableSetDescription> Functional_EE;

      typedef typename Functional_HH::Scalar                                             Scalar;

      typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HH> >            Assembler_HH;
      typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HE> >            Assembler_HE;
      typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EH> >            Assembler_EH;
      typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EE> >            Assembler_EE;

      static constexpr int stateId = getStateId<Functional_HH>();
      static constexpr int controlId = getControlId<Functional_HH>();
      static constexpr int adjointId = getAdjointId<Functional_HH>();

      static constexpr int dim = VariableSetDescription::Grid::dimension;

      /*
       * operators:
       * the state operator is called A
       * the control enters the state equation via operator B
       * the derivatives of the lagrange functional are called Lyy, Lyu, Luy, Luu
       *
       * the endings indicate which ansatz and test spaces are used (L=ansatz space H=extension space), i.e. Luy_HE maps from the
       * ansatz space to the extension space
       */
      typedef AssembledGalerkinOperator<Assembler_HH> Operator_HH;
      typedef AssembledGalerkinOperator<Assembler_EH> Operator_EH;
      typedef AssembledGalerkinOperator<Assembler_HE> Operator_HE;
      typedef AssembledGalerkinOperator<Assembler_EE> Operator_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,stateId,stateId+1,stateId,stateId+1>         Lyy_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,stateId,stateId+1,stateId,stateId+1>         Lyy_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,stateId,stateId+1,stateId,stateId+1>         Lyy_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,stateId,stateId+1,stateId,stateId+1>         Lyy_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,stateId,stateId+1,controlId,controlId+1>     Lyu_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,stateId,stateId+1,controlId,controlId+1>     Lyu_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,stateId,stateId+1,controlId,controlId+1>     Lyu_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,stateId,stateId+1,controlId,controlId+1>     Lyu_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,controlId,controlId+1,stateId,stateId+1>     Luy_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,controlId,controlId+1,stateId,stateId+1>     Luy_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,controlId,controlId+1,stateId,stateId+1>     Luy_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,controlId,controlId+1,stateId,stateId+1>     Luy_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,controlId,controlId+1,controlId,controlId+1> Luu_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,controlId,controlId+1,controlId,controlId+1> Luu_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,controlId,controlId+1,controlId,controlId+1> Luu_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,controlId,controlId+1,controlId,controlId+1> Luu_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,adjointId,adjointId+1,controlId,controlId+1> B_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,adjointId,adjointId+1,controlId,controlId+1> B_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,adjointId,adjointId+1,controlId,controlId+1> B_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,adjointId,adjointId+1,controlId,controlId+1> B_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,adjointId,adjointId+1,stateId,stateId+1>     A_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,adjointId,adjointId+1,stateId,stateId+1>     A_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,adjointId,adjointId+1,stateId,stateId+1>     A_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,adjointId,adjointId+1,stateId,stateId+1>     A_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,controlId,controlId+1,adjointId,adjointId+1> BT_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,controlId,controlId+1,adjointId,adjointId+1> BT_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,controlId,controlId+1,adjointId,adjointId+1> BT_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,controlId,controlId+1,adjointId,adjointId+1> BT_EE;

      typedef AssembledGalerkinOperator<Assembler_HH,stateId,stateId+1,adjointId,adjointId+1>     AT_HH;
      typedef AssembledGalerkinOperator<Assembler_EH,stateId,stateId+1,adjointId,adjointId+1>     AT_EH;
      typedef AssembledGalerkinOperator<Assembler_HE,stateId,stateId+1,adjointId,adjointId+1>     AT_HE;
      typedef AssembledGalerkinOperator<Assembler_EE,stateId,stateId+1,adjointId,adjointId+1>     AT_EE;

      /**
       * Coefficient vectors
       */
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<>::type                               Vector;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<stateId,stateId+1>::type              VectorY;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<controlId,controlId+1>::type          VectorU;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<adjointId,adjointId+1>::type          VectorP;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<>::type                      ExtensionVector;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<stateId,stateId+1>::type     ExtensionVectorY;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<controlId,controlId+1>::type ExtensionVectorU;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<adjointId,adjointId+1>::type ExtensionVectorP;

      /*
       * Coefficient vector initializers
       */
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<>                                     Vector_Initializer;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<stateId,stateId+1>                    VectorY_Initializer;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<controlId,controlId+1>                VectorU_Initializer;
      typedef typename VariableSetDescription::template CoefficientVectorRepresentation<adjointId,adjointId+1>                VectorP_Initializer;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<>                            ExtensionVector_Initializer;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<stateId,stateId+1>           ExtensionVectorY_Initializer;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<controlId,controlId+1>       ExtensionVectorU_Initializer;
      typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<adjointId,adjointId+1>       ExtensionVectorP_Initializer;
    };
  }

  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverLA, class LinearSolverHA, class LinearSolverLM=LinearSolverLA, class LinearSolverHM=LinearSolverHA,
  bool lumpM=false, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class HierarchicalBasisErrorEstimator2 : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    static constexpr int noOfVariables = ExtensionVariableSetDescription::noOfVariables;
    template <class AnsatzVars, class TestVars, class OriginVars> using ErrorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<ErrorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
    //    typedef typename LF::AnsatzVars LVars;
    static constexpr int yId = 1, uId = 0, pId = 2;
  public:
    typedef Functional<VariableSetDescription,VariableSetDescription,VariableSetDescription,lumpM>                             Functional_HH;
    typedef Functional<ExtensionVariableSetDescription,VariableSetDescription,ExtensionVariableSetDescription,lumpM>         Functional_EH;
    typedef Functional<VariableSetDescription,ExtensionVariableSetDescription,VariableSetDescription,lumpM>                    Functional_HE;
    typedef Functional<ExtensionVariableSetDescription,ExtensionVariableSetDescription,ExtensionVariableSetDescription,lumpM>  Functional_EE;
    typedef HierarchicErrorEstimator<LinearizationAt<Functional_HH>,ExtensionVariableSetDescription,ExtensionVariableSetDescription> ErrorEstimator;
    //    typedef HierarchicErrorEstimator<LinearizationAt<Functional_HH>,ExtensionVariableSetDescription,ExtensionVariableSetDescription,HierarchicErrorEstimatorDetail::TakeAllD2<LinearizationAt<Functional_HH> > > ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> Assembler;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HH> >                                 Ass_HH;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HE> >                                 Ass_HE;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EH> >                                 Ass_EH;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EE> >                                 Ass_EE;

    typedef typename Traits::AT_HH      AT_HH;
    typedef typename Traits::AT_EE      AT_EE;
    typedef typename Traits::AT_EH      AT_EH;

    typedef AssembledGalerkinOperator<Ass_EE,pId,pId+1,yId,yId+1>                                              Ahh;
    typedef AssembledGalerkinOperator<Ass_HH,pId,pId+1,yId,yId+1>                                              All;
    typedef AssembledGalerkinOperator<Ass_HE,pId,pId+1,yId,yId+1>                                              Alh;
    typedef AssembledGalerkinOperator<Ass_EH,pId,pId+1,yId,yId+1>                                              Ahl;
    typedef AssembledGalerkinOperator<Ass_EE,yId,yId+1,yId,yId+1>                                              Myyhh;

    typedef AssembledGalerkinOperator<Ass_HH,pId,pId+1,uId,uId+1>                                              Bll;
    typedef AssembledGalerkinOperator<Ass_HE,pId,pId+1,uId,uId+1>                                              Blh;
    typedef AssembledGalerkinOperator<Ass_EH,pId,pId+1,uId,uId+1>                                              Bhl;
    typedef AssembledGalerkinOperator<Ass_EE,pId,pId+1,uId,uId+1>                                              Bhh;

    typedef AssembledGalerkinOperator<Ass_HH,uId,uId+1,pId,pId+1>                                              BTll;
    typedef AssembledGalerkinOperator<Ass_HE,uId,uId+1,pId,pId+1>                                              BTlh;
    typedef AssembledGalerkinOperator<Ass_EH,uId,uId+1,pId,pId+1>                                              BThl;
    typedef AssembledGalerkinOperator<Ass_EE,uId,uId+1,pId,pId+1>                                              BThh;

    typedef AssembledGalerkinOperator<Ass_HH,uId,uId+1,uId,uId+1>                                             Muull;
    typedef AssembledGalerkinOperator<Ass_HE,uId,uId+1,uId,uId+1>                                             Muulh;
    typedef AssembledGalerkinOperator<Ass_EH,uId,uId+1,uId,uId+1>                                             Muuhl;
    typedef AssembledGalerkinOperator<Ass_EE,uId,uId+1,uId,uId+1>                                             Muuhh;

    typedef typename Functional_HH::Scalar Scalar;
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
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;

    typedef MatrixAsTriplet<Scalar> Matrix;


    HierarchicalBasisErrorEstimator2(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
      initAssemblers();
      }

    virtual ~HierarchicalBasisErrorEstimator2(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      double ilutTol = 0.01;
      int ilutLfil = 100;
      double minresTol = 1e-3;

      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<VarSet> const& xl = dynamic_cast<const Bridge::Vector<VarSet>&>(x_);
      Bridge::Vector<ExtVarSet> xe(extensionVariableSetDescription);
      assembleAt(xl,xe);
      Bridge::Vector<VarSet> const& dxl = dynamic_cast<const Bridge::Vector<VarSet>&>(dx_);

      /***************************************************************************************/
      // estimate error in dual variable
      typename ExtensionVariableSetDescription::VariableSet tmpRep(extensionVariableSetDescription), errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet tmpRepL(variableSetDescription), errorEstimate_H(variableSetDescription);

      All all(*ass_HH,false);
      Ahh ahh(*ass_EE,false);
      Alh alh(*ass_HE,false);
      Ahl ahl(*ass_EH,false);

      Matrix mAll = all.template get<Matrix>();
      mAll.transpose();
      Matrix mAhh = ahh.template get<Matrix>();
      mAhh.transpose();
      Matrix mAlh = alh.template get<Matrix>();
      mAlh.transpose();
      Matrix mAhl = ahl.template get<Matrix>();
      mAhl.transpose();

      boost::timer::cpu_timer atimer;
      assembler.assemble(ErrorEstimator(LinearizationAt<Functional_HH>(*F_HH,xl.get()),dxl.get()));
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << "s" << std::endl;
      typedef AssembledGalerkinOperator<Assembler,yId,yId+1,pId,pId+1> StateOperator;
      StateOperator stateOp(assembler,false);
      CoefficientVectorHY rhs0(assembler.template rhs<yId,yId+1>());
      CoefficientVectorHP sol0(CVHP::init(extensionVariableSetDescription));


      tmpRep *= 0;
      at_c<yId>(tmpRep.data) = at_c<0>(rhs0.data);
      CoefficientVectorLP xlp(CVLP::init(variableSetDescription));
      at_c<0>(xlp.data) = at_c<pId>(dxl.get().data).coefficients();
      std::vector<Scalar> tmpRepVec, xlpVec;
      IstlInterfaceDetail::toVector(rhs0,tmpRepVec);
      if(verbose){
        std::cout << "initial rhs entries:" << std::endl;
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
      }
      IstlInterfaceDetail::toVector(xlp,xlpVec);
      mAhl.usmv(-1.0,xlpVec,tmpRepVec);

      if(verbose)
      {
        std::cout << "corrected rhs entries:" << std::endl;
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
      }
      IstlInterfaceDetail::fromVector(tmpRepVec,rhs0);
      Dune::InverseOperatorResult res;
      boost::timer::cpu_timer timer;
      if(verbose) std::cout << "first solve: ";
      //LinearSolverHA(mAhh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(sol0,rhs0);
      if(fast)
      {
        MatrixRepresentedOperator<Matrix,CoefficientVectorH23,CoefficientVectorH23> mAhhOp(mAhh);
        ILUTPreconditioner<MatrixRepresentedOperator<Matrix,CoefficientVectorH23,CoefficientVectorH23> > prec(mAhhOp,ilutLfil,ilutTol,1);

        MINRESSolverK<CoefficientVectorH23> minres(mAhhOp, prec, minresTol, 100, 1);
        minres.apply(sol0,rhs0,res);
      }
      else InverseLinearOperator<DirectSolver<CoefficientVectorHP,CoefficientVectorHP> >(DirectSolver<CoefficientVectorHP,CoefficientVectorHP>(mAhh,DirectType::UMFPACK3264,MatrixProperties::GENERAL)).apply(rhs0,sol0);
      if(verbose) std::cout << boost::timer::format(timer.elapsed()) << std::endl;
      //      InverseLinearOperator<DirectSolver<CoefficientVectorHP,CoefficientVectorHP> >(DirectSolver<CoefficientVectorHP,CoefficientVectorHP>(mAhh,DirectType::MUMPS,MatrixProperties::GENERAL)).apply(tmpRepVec);
      //      IstlInterfaceDetail::fromVector(tmpRepVec,sol0);

      //      //MyH1SemiNorm myNorm;
      //      //std::cout << "||rhs||=" << myNorm(tmpRep) << std::endl;
      //      LinearSolverHA(stateOp,DirectType::MUMPS,MatrixProperties::GENERAL).apply(sol0,rhs0);
      //      //directInverseOperator(stateOp,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhs0,sol0);
      //      tmpRep *= 0;
      at_c<pId>(tmpRep.data) = at_c<0>(sol0.data);
      at_c<pId>(errorEstimate_E.data) = at_c<0>(sol0.data);
      //std::cout << "||sol||=" << myNorm(tmpRep) << std::endl;
      Scalar maxVal = std::numeric_limits<Scalar>::min(),
          minVal = std::numeric_limits<Scalar>::max(),
          sum = 0;
      if(verbose){
        for(int i=0; i<at_c<pId>(errorEstimate_E.data).size();  ++i)
        {
          for(int j=0; j<at_c<pId>(errorEstimate_E.data).coefficients()[i].size(); ++j)
          {
            Scalar val = at_c<pId>(errorEstimate_E.data).coefficients()[i][j];
            sum += val*val;
            if(val > maxVal) maxVal = val;
            if(val < minVal) minVal = val;
          }
        }
        std::cout << "p: minimal value: " << minVal << std::endl;
        std::cout << "p: maximal value: " << maxVal << std::endl;
        std::cout << "p: l2-sum: " << sum << std::endl;
      }
      /***************************************************************************************/
      // estimate error in control variable
      // rhs vectors
      CoefficientVectorLU rhs1L(ass_EH->template rhs<uId,uId+1>());
      CoefficientVectorHU rhs1H(ass_EE->template rhs<uId,uId+1>());
      Bhl bhl(*ass_EH,false);
      Blh blh(*ass_HE,false);
      Bll bll(*ass_HH,false);
      Bhh bhh(*ass_EE,false);

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

      Muuhh muuhh(*ass_EE,true);
      Muull muull(*ass_HH,true);
      Muulh muulh(*ass_HE,true);
      CoefficientVectorHU solhu(CVHU::init(extensionVariableSetDescription));
      CoefficientVectorLU sollu(CVLU::init(variableSetDescription));
      //LinearSolverLM(muull,DirectType::MUMPS,MatrixProperties::GENERAL).apply(sollu,rhs1L);
      timer.start();
      if(verbose) std::cout << "second solve: ";
      JacobiPreconditioner<Muull>(muull).apply(sollu,rhs1L);
      //      directInverseOperator(muull,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(rhs1L,sollu);
      if(verbose) std::cout << boost::timer::format(timer.elapsed()) << std::endl;
      muulh.applyscaleadd(-1.0,sollu,rhs1H);
      JacobiPreconditioner<Muuhh>(muuhh).apply(solhu,rhs1H);
      //      LinearSolverHM(muuhh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solhu,rhs1H);
      //      directInverseOperator(muuhh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(rhs1H,solhu);

      //std::cout << "sollu=" << at_c<0>(sollu.data)[0] << std::endl;
      //std::cout << "solhu=" << at_c<0>(solhu.data)[0] << std::endl;
      at_c<uId>(errorEstimate_E.data) = at_c<0>(solhu.data);
      at_c<uId>(errorEstimate_H.data) = at_c<0>(sollu.data);
      /***************************************************************************************/
      // estimate error in state variable
      CoefficientVectorLY solly(CVLP::init(variableSetDescription));
      CoefficientVectorHY solhy(CVHP::init(extensionVariableSetDescription));
      CoefficientVectorLP rhslp(ass_HH->template rhs<pId,pId+1>());
      CoefficientVectorHP rhshp(ass_EE->template rhs<pId,pId+1>());

      std::vector<Scalar> tmprhslp, tmprhshp;
      IstlInterfaceDetail::toVector(rhslp,tmprhslp);
      IstlInterfaceDetail::toVector(rhshp,tmprhshp);
      if(verbose){
        std::cout << "initial rhsl entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhslp.size(); ++i){
          if(std::isnan(tmprhslp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhslp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }
      if(verbose){
        std::cout << "initial rhsh entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhshp.size(); ++i){
          if(std::isnan(tmprhshp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhshp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }

      bll.applyscaleadd(-1.0,sollu,rhslp);
      if(verbose){
        std::cout << "1: rhsl entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhslp.size(); ++i){
          if(std::isnan(tmprhslp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhslp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }
      bhl.applyscaleadd(-1.0,solhu,rhslp);
      if(verbose){
        std::cout << "2: rhsl entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhslp.size(); ++i){
          if(std::isnan(tmprhslp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhslp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }

      blh.applyscaleadd(-1.0,sollu,rhshp);
      if(verbose){
        std::cout << "1: rhsh entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhshp.size(); ++i){
          if(std::isnan(tmprhshp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhshp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }
      bhh.applyscaleadd(-1.0,solhu,rhshp);
      if(verbose) {
        std::cout << "2: rhsh entries:" << std::endl;
        Scalar maxRhs=std::numeric_limits<Scalar>::min(), minRhs=std::numeric_limits<Scalar>::max(), rhsSum=0;
        for(size_t i=0; i<tmprhshp.size(); ++i){
          if(std::isnan(tmprhshp[i])) std::cout << "nan at " << i << std::endl;

          Scalar val = tmprhshp[i];
          if(val < minRhs) minRhs = val;
          if(val > maxRhs) maxRhs = val;
          rhsSum += val*val;
        }
        std::cout << "max rhs: " << maxRhs << std::endl;
        std::cout << "min rhs: " << minRhs << std::endl;
        std::cout << "sum: " << rhsSum << std::endl;
      }

      //      LinearSolverLA(ahh,DirectType::MUMPS,MatrixProperties::GENERAL).apply(solhy,rhshp);
      //      ahl.applyscaleadd(-1.0,solhy,rhslp);
      //      LinearSolverLA(all,DirectType::MUMPS,MatrixProperties::GENERAL).apply(solly,rhslp);
      timer.start();
      if(verbose) std::cout << "third solve: ";
      LinearSolverLA(all,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solly,rhslp);
      if(verbose) std::cout << boost::timer::format(timer.elapsed()) << std::endl;
      if(verbose)std::cout << "first solve for last variable finished." << std::endl;
      alh.applyscaleadd(-1.0,solly,rhshp);
      timer.start();
      if(verbose) std::cout << "fourth solve: ";
      if(fast)
      {
        ILUTPreconditioner<Ahh> prec2(ahh,ilutLfil,ilutTol,1);

        MINRESSolverK<CoefficientVectorH23> minres2(ahh, prec2, minresTol, 100, 1);
        minres2.apply(solhy,rhshp,res);
      }
      else LinearSolverHA(ahh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solhy,rhshp);
      if(verbose) std::cout << boost::timer::format(timer.elapsed()) << std::endl;
      at_c<yId>(errorEstimate_E.data) = at_c<0>(solhy.data);
      at_c<yId>(errorEstimate_H.data) = at_c<0>(solly.data);

      if(verbose){
        Scalar tmpMax = maxVal, tmpMin = minVal, tmpSum = sum;
        maxVal = std::numeric_limits<Scalar>::min();
        minVal = std::numeric_limits<Scalar>::max();
        sum = 0;

        for(int i=0; i<at_c<uId>(errorEstimate_E.data).size();  ++i)
        {
          for(int j=0; j<at_c<uId>(errorEstimate_E.data).coefficients()[i].size(); ++j)
          {
            Scalar val = at_c<uId>(errorEstimate_E.data).coefficients()[i][j];
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

        for(int i=0; i<at_c<yId>(errorEstimate_E.data).size();  ++i)
        {
          for(int j=0; j<at_c<yId>(errorEstimate_E.data).coefficients()[i].size(); ++j)
          {
            Scalar val = at_c<yId>(errorEstimate_E.data).coefficients()[i][j];
            sum += val*val;
            if(val > maxVal) maxVal = val;
            if(val < minVal) minVal = val;
          }
        }
        std::cout << "y: minimal value: " << minVal << std::endl;
        std::cout << "y: maximal value: " << maxVal << std::endl;
        std::cout << "y: l2-sum: " << sum << std::endl;

        tmpSum += sum;
        std::cout << "minimal value: " << std::min(tmpMin,minVal) << std::endl;
        std::cout << "maximal value: " << std::max(tmpMax,maxVal) << std::endl;
        std::cout << "l2-sum: " << tmpSum << std::endl;
      }
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
//      Scalar tol = squaredFraction*totalErrorSquared;
//      Scalar sum = 0;
//      std::sort(errorDistribution.begin(), errorDistribution.end(), [](std::pair<double,int> const& a, std::pair<double,int> const& b) -> bool
//          {
//        return a.first > b.first;
//          });
//      std::vector<size_t> cellIds;
//      for(size_t i=0; i<errorDistribution.size(); ++i)
//      {
//        sum += errorDistribution[i].first;
//        cellIds.push_back(errorDistribution[i].second);
//        if(sum > tol) break;
//      }
//
//      forEachCell(variableSetDescription.gridView, [&](Cell const& cell)
//          {
//        if(std::find(cellIds.begin(),cellIds.end(),variableSetDescription.gridView.indexSet().index(cell))!=cellIds.end()) extensionSpace.gridManager().mark(1,cell);
//          });
//      extensionSpace.gridManager().adaptAtOnce();
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
      F_HH.reset(new Functional_HH(args...));
      F_HE.reset(new Functional_HE(args...));
      F_EH.reset(new Functional_EH(args...));
      F_EE.reset(new Functional_EE(args...));
    }

  private:
    void initAssemblers()
    {
      ass_HH.reset(new Ass_HH(extensionVariableSetDescription.spaces));
      ass_HE.reset(new Ass_HE(extensionVariableSetDescription.spaces));
      ass_EH.reset(new Ass_EH(extensionVariableSetDescription.spaces));
      ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
    }

    void assembleAt(const Bridge::Vector<VarSet>& xl, const Bridge::Vector<ExtVarSet>& xe)
    {
      ass_HH->assemble(linearization(*F_HH,xl.get()));
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EH->assemble(linearization(*F_EH,xe.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<Functional_HH> F_HH;
    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EH> F_EH;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Ass_HH> ass_HH;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EH> ass_EH;
    std::unique_ptr<Ass_EE> ass_EE;
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
    bool verbose, fast;
    // AdjustRHS<CoefficientVector> adjustRHS;
  };


  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverHA, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class StateEquationHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // Functionals
    typedef typename Traits::Functional_HH Functional_HH;
    typedef typename Traits::Functional_HE Functional_HE;
    typedef typename Traits::Functional_EE Functional_EE;

    // assemblers
    typedef typename Traits::Assembler_HE Ass_HE;
    typedef typename Traits::Assembler_EE Ass_EE;

    // operators
    typedef typename Traits::B_HE Blh;
    typedef typename Traits::A_HE Alh;
    typedef typename Traits::A_EE Ahh;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    StateEquationHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new Ass_HE(extensionVariableSetDescription.spaces));
      }

    virtual ~StateEquationHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      Blh blh(*ass_HE,false);
      Alh alh(*ass_HE,false);
      Ahh ahh(*ass_EE,false);
      typename Traits::ExtensionVectorP rhs0(ass_EE->template rhs<Traits::adjointId,Traits::adjointId+1>());
      rhs0 *= -1.0;
      typename Traits::ExtensionVectorY sol0(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorU sollu(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::VectorY solly(Traits::VectorY_Initializer::init(variableSetDescription));
      at_c<0>(solly.data) = at_c<Traits::stateId>(dxl.get().data).coefficients();
      at_c<0>(sollu.data) = at_c<Traits::controlId>(dxl.get().data).coefficients();
      blh.applyscaleadd(-1.0,sollu,rhs0);
      alh.applyscaleadd(-1.0,solly,rhs0);
      boost::timer::cpu_timer ctimer;
      LinearSolverHA(ahh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol0,rhs0);
      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << std::endl;
      at_c<Traits::stateId>(errorEstimate_E.data) = at_c<0>(sol0.data);

      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HE.reset(new Functional_HE(args...));
      F_EE.reset(new Functional_EE(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };


  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverHA, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class VariationalEquationHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // Functionals
    typedef typename Traits::Functional_HE Functional_HE;
    typedef typename Traits::Functional_EE Functional_EE;

    // assemblers
    typedef typename Traits::Assembler_HE Ass_HE;
    typedef typename Traits::Assembler_EE Ass_EE;

    // operators
    typedef typename Traits::Luu_HE MLH;
    typedef typename Traits::Luu_EE MHH;
    typedef typename Traits::BT_HE BTLH;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;

    VariationalEquationHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new Ass_HE(extensionVariableSetDescription.spaces));
      }

    virtual ~VariationalEquationHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable

      MLH mlh(*ass_HE,false);
      MHH mhh(*ass_EE,false);
      BTLH btlh(*ass_HE,false);

      typename Traits::ExtensionVectorU sol(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorU sollu(at_c<Traits::controlId>(dxl.get().data.coefficients()));
      typename Traits::VectorP sollp(at_c<Traits::adjointId>(dxl.get().data.coefficients()));
      typename Traits::ExtensionVectorU rhs(ass_EE->template rhs<Traits::controlId,Traits::controlId+1>());
      rhs *= -1.0;
      mlh.applyscaleadd(-1.0,sollu,rhs);
      btlh.applyscaleadd(-1.0,sollp,rhs);

      LinearSolverHA(mhh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol,rhs);
      at_c<Traits::control>(errorEstimate_E.data) = at_c<0>(sol.data);

      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HE.reset(new Functional_HE(args...));
      F_EE.reset(new Functional_EE(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };


  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverHA, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class AdjointEquationHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // assemblers
    typedef typename Traits::Assembler_HE Ass_HE;
    typedef typename Traits::Assembler_EE Ass_EE;

    // operators
    typedef typename Traits::Lyy_HE Lyy_HE;
    typedef typename Traits::AT_HE AT_HE;
    typedef typename Traits::AT_EE AT_EE;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    AdjointEquationHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new Ass_HE(extensionVariableSetDescription.spaces));
      }

    virtual ~AdjointEquationHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      Lyy_HE Lyy_HE(*ass_HE,false);
      AT_HE at_HE(*ass_HE,false);
      AT_EE at_EE(*ass_EE,false);
      typename Traits::ExtensionVectorY rhs0(ass_EE->template rhs<Traits::stateId,Traits::stateId+1>());
      rhs0 *= -1.0;
      typename Traits::ExtensionVectorP sol0(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorP sollp(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::VectorY solly(Traits::VectorY_Initializer::init(variableSetDescription));
      at_c<0>(solly.data) = at_c<Traits::stateId>(dxl.get().data).coefficients();
      at_c<0>(sollp.data) = at_c<Traits::adjointId>(dxl.get().data).coefficients();
      Lyy_HE.applyscaleadd(-1.0,solly,rhs0);
      at_HE.applyscaleadd(-1.0,sollp,rhs0);
      boost::timer::cpu_timer ctimer;
      LinearSolverHA(at_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol0,rhs0);
      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << std::endl;
      at_c<Traits::stateId>(errorEstimate_E.data) = at_c<0>(sol0.data);

      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HE.reset(new typename Traits::Functional_HE(args...));
      F_EE.reset(new typename Traits::Functional_EE(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<typename Traits::Functional_HE> F_HE;
    std::unique_ptr<typename Traits::Functional_EE> F_EE;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };


  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverLA, class LinearSolverHA, class LinearSolverLU, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class AdjointEquationLinearPropagationHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // assemblers

    // operators
    typedef typename Traits::Lyy_HE Lyy_HE;
    typedef typename Traits::AT_HE AT_HE;
    typedef typename Traits::AT_EE AT_EE;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    AdjointEquationLinearPropagationHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new typename Traits::Assembler_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new typename Traits::Assembler_HE(extensionVariableSetDescription.spaces));
        ass_EH.reset(new typename Traits::Assembler_EH(extensionVariableSetDescription.spaces));
        ass_HH.reset(new typename Traits::Assembler_HH(extensionVariableSetDescription.spaces));
      }

    virtual ~AdjointEquationLinearPropagationHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      /**
       * error in adjoint equation
       */
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      Lyy_HE Lyy_HE(*ass_HE,false);
      AT_HE at_HE(*ass_HE,false);
      AT_EE at_EE(*ass_EE,false);
      typename Traits::ExtensionVectorY rhs0(ass_EE->template rhs<Traits::stateId,Traits::stateId+1>());
      rhs0 *= -1.0;
      typename Traits::ExtensionVectorP sol0(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorP sollp(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::VectorY solly(Traits::VectorY_Initializer::init(variableSetDescription));
      at_c<0>(solly.data) = at_c<Traits::stateId>(dxl.get().data).coefficients();
      at_c<0>(sollp.data) = at_c<Traits::adjointId>(dxl.get().data).coefficients();
      Lyy_HE.applyscaleadd(-1.0,solly,rhs0);
      at_HE.applyscaleadd(-1.0,sollp,rhs0);
      boost::timer::cpu_timer ctimer;
      LinearSolverHA(at_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol0,rhs0);
      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << std::endl;

      std::cout << "ERR4 sol0 " << (sol0*sol0) << std::endl;
      /**
       * propagation through variational equality and state equation
       */
      // variational equality
      typename Traits::VectorU rhsU(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::BT_EH bt_EH(*ass_EH,false);
      std::cout << "ERR4 rhsU0 " << (rhsU*rhsU) << std::endl;
      bt_EH.applyscaleadd(-1.0,sol0,rhsU);
      std::cout << "ERR4 rhsU1 " << (rhsU*rhsU) << std::endl;

      typename Traits::Luu_HH luu(*ass_HH,false);
      typename Traits::VectorU solU(Traits::VectorU_Initializer::init(variableSetDescription));
      std::cout << "ERR4 solU0 " << (solU*solU) << std::endl;
      LinearSolverLU(luu,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solU,rhsU);

      std::cout << "ERR4 solU " << (solU*solU) << std::endl;
      // state equation
      typename Traits::VectorP rhsP(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::B_HH b_HH(*ass_HH,false);
      b_HH.applyscaleadd(-1.0,solU,rhsP);
      typename Traits::A_HH a_HH(*ass_HH,false);
      typename Traits::VectorY solY(Traits::VectorY_Initializer::init(variableSetDescription));
      LinearSolverLA(a_HH,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solY,rhsP);
      std::cout << "ERR4 solY " << (solY*solY) << std::endl;
      at_c<Traits::stateId>(errorEstimate_H.data) = at_c<0>(solY.data);
      at_c<Traits::controlId>(errorEstimate_H.data) = at_c<0>(solU.data);
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HH.reset(new typename Traits::Functional_HH(args...));
      F_HE.reset(new typename Traits::Functional_HE(args...));
      F_EH.reset(new typename Traits::Functional_EH(args...));
      F_EE.reset(new typename Traits::Functional_EE(args...));
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
    std::unique_ptr<typename Traits::Functional_EE> F_EE;
    std::unique_ptr<typename Traits::Assembler_HH> ass_HH;
    std::unique_ptr<typename Traits::Assembler_HE> ass_HE;
    std::unique_ptr<typename Traits::Assembler_EH> ass_EH;
    std::unique_ptr<typename Traits::Assembler_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };

  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverLA, class LinearSolverHA, class LinearSolverHU, class LinearSolverLU, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration, bool lump=false>
  class AnotherHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,lump>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // assemblers

    // operators
    typedef typename Traits::Lyy_HE Lyy_HE;
    typedef typename Traits::AT_HE AT_HE;
    typedef typename Traits::AT_EE AT_EE;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    AnotherHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new typename Traits::Assembler_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new typename Traits::Assembler_HE(variableSetDescription.spaces));
        ass_EH.reset(new typename Traits::Assembler_EH(extensionVariableSetDescription.spaces));
        ass_HH.reset(new typename Traits::Assembler_HH(variableSetDescription.spaces));
      }

    virtual ~AnotherHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      std::cout << "ids: " << Traits::stateId << " " << Traits::controlId << " " << Traits::adjointId << std::endl;
      /***************************************************************************************/
      /**
       * error in adjoint equation
       */
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      Lyy_HE Lyy_HE(*ass_HE,false);
      AT_HE at_HE(*ass_HE,false);
      AT_EE at_EE(*ass_EE,false);
      typename Traits::ExtensionVectorY rhs0(ass_EE->template rhs<Traits::stateId,Traits::stateId+1>());
      rhs0 *= -1.0;
      typename Traits::ExtensionVectorP sol0(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
      typename Traits::VectorP sollp(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::VectorY solly(Traits::VectorY_Initializer::init(variableSetDescription));
      at_c<0>(solly.data) = at_c<Traits::stateId>(dxl.get().data).coefficients();
      at_c<0>(sollp.data) = at_c<Traits::adjointId>(dxl.get().data).coefficients();
      Lyy_HE.applyscaleadd(-1.0,solly,rhs0);
      at_HE.applyscaleadd(-1.0,sollp,rhs0);
      boost::timer::cpu_timer ctimer;
      //JacobiPreconditioner<AT_EE>(at_EE).apply(sol0,rhs0);
      LinearSolverHA(at_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol0,rhs0);
      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << std::endl;

      std::cout << "|solPH|^2=" << (sol0*sol0) << std::endl;

      /**
       * propagation through variational equality and state equation
       */
      // variational equality
      typename Traits::VectorU rhsUL(ass_HH->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUL *= -1.0;
      typename Traits::ExtensionVectorU rhsUH(ass_EE->template rhs<Traits::controlId,Traits::controlId+1>());
      rhsUH *= -1.0;
      std::cout << "|rhsUH|^2=" << (rhsUH*rhsUH) << std::endl;
      std::cout << "|rhsUL|^2=" << (rhsUL*rhsUL) << std::endl;
      typename Traits::BT_EH bt_EH(*ass_EH,false);
      typename Traits::BT_EE bt_EE(*ass_EE,false);
      bt_EH.applyscaleadd(-1.0,sol0,rhsUL);
      bt_EE.applyscaleadd(-1.0,sol0,rhsUH);
      std::cout << "|rhsUH|^2=" << (rhsUH*rhsUH) << std::endl;
      std::cout << "|rhsUL|^2=" << (rhsUL*rhsUL) << std::endl;
      typename Traits::Luu_EE luu_EE(*ass_EE,true);
      typename Traits::Luu_EH luu_EH(*ass_EH,false);
      typename Traits::Luu_HH luu_HH(*ass_HH,false);
      typename Traits::VectorU solUL(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorU solUH(Traits::ExtensionVectorU_Initializer::init(extensionVariableSetDescription));
      std::cout << "computing u_h" << std::endl;
      JacobiPreconditioner<typename Traits::Luu_EE>(luu_EE).apply(solUH,rhsUH);//LinearSolverHU(luu_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solUH,rhsUH);
      std::cout << "|solUH|^2=" << (solUH*solUH) << std::endl;
      //luu_EH.applyscaleadd(-1.0,solUH,rhsUL);
      std::cout << "|rhsUL|^2=" << (rhsUL*rhsUL) << std::endl;
      std::cout << "computing u_l" << std::endl;
      JacobiPreconditioner<typename Traits::Luu_HH>(luu_HH).apply(solUL,rhsUL);//LinearSolverLU(luu_HH,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solUL,rhsUL);
      std::cout << "|solUH|^2=" << (solUH*solUH) << std::endl;
      std::cout << "|solUL|^2=" << (solUL*solUL) << std::endl;

      // state equation
      typename Traits::VectorP rhsPL(ass_HH->template rhs<Traits::adjointId,Traits::adjointId+1>());
      rhsPL *= -1.0;
      typename Traits::ExtensionVectorP rhsPH(ass_EE->template rhs<Traits::adjointId,Traits::adjointId+1>());
      rhsPH *= -1.0;
      typename Traits::B_HH b_HH(*ass_HH,false);
      typename Traits::B_HE b_HE(*ass_HE,false);
      typename Traits::B_EH b_EH(*ass_EH,false);
      typename Traits::B_EE b_EE(*ass_EE,false);
      b_HH.applyscaleadd(-1.0,solUL,rhsPL);
      b_HE.applyscaleadd(-1.0,solUL,rhsPH);
      b_EH.applyscaleadd(-1.0,solUH,rhsPL);
      b_HH.applyscaleadd(-1.0,solUH,rhsPH);

      typename Traits::A_HH a_HH(*ass_HH,false);
      typename Traits::A_EH a_EH(*ass_EH,false);
      typename Traits::A_EE a_EE(*ass_EE,false);

      typename Traits::VectorY solYL(Traits::VectorY_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorY solYH(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));
      JacobiPreconditioner<typename Traits::A_EE>(a_EE).apply(solYH,rhsPH);//LinearSolverHA(a_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solYH,rhsPH);

      a_EH.applyscaleadd(-1.0,solYH,rhsPL);
      JacobiPreconditioner<typename Traits::A_HH>(a_HH).apply(solYL,rhsPL);//LinearSolverLA(a_HH,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(solYL,rhsPL);
      std::cout << "|solYH|^2=" << (solYH*solYH) << std::endl;
      std::cout << "|solYL|^2=" << (solYL*solYL) << std::endl;
      at_c<Traits::stateId>(errorEstimate_H.data) = at_c<0>(solYL.data);
      at_c<Traits::controlId>(errorEstimate_H.data) = at_c<0>(solUL.data);
      at_c<Traits::stateId>(errorEstimate_E.data) = at_c<0>(solYH.data);
      at_c<Traits::controlId>(errorEstimate_E.data) = at_c<0>(solUH.data);
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,squaredFraction*totalErrorSquared);
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
      F_HH.reset(new typename Traits::Functional_HH(args...));
      F_HE.reset(new typename Traits::Functional_HE(args...));
      F_EH.reset(new typename Traits::Functional_EH(args...));
      F_EE.reset(new typename Traits::Functional_EE(args...));
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
    std::unique_ptr<typename Traits::Functional_EE> F_EE;
    std::unique_ptr<typename Traits::Assembler_HH> ass_HH;
    std::unique_ptr<typename Traits::Assembler_HE> ass_HE;
    std::unique_ptr<typename Traits::Assembler_EH> ass_EH;
    std::unique_ptr<typename Traits::Assembler_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };


  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverHA, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class StupidHBErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // Functionals
    typedef typename Traits::Functional_HH Functional_HH;
    typedef typename Traits::Functional_HE Functional_HE;
    typedef typename Traits::Functional_EE Functional_EE;

    // assemblers
    typedef typename Traits::Assembler_HE Ass_HE;
    typedef typename Traits::Assembler_EE Ass_EE;

    // operators
    typedef typename Traits::Operator_HE A_HE;
    typedef typename Traits::Operator_EE A_EE;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    StupidHBErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new Ass_HE(variableSetDescription.spaces));
      }

    virtual ~StupidHBErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      A_HE alh(*ass_HE,false);
      A_EE ahh(*ass_EE,false);
      typename Traits::ExtensionVector rhs_H(ass_EE->rhs());
      rhs_H *= -1.0;
      typename Traits::ExtensionVector sol_H(Traits::ExtensionVector_Initializer::init(extensionVariableSetDescription));
      typename Traits::Vector sol_L(Traits::Vector_Initializer::init(variableSetDescription));
      sol_L = dxl.get();
      alh.applyscaleadd(-1.0,sol_L,rhs_H);
      boost::timer::cpu_timer ctimer;
      LinearSolverHA(ahh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol_H,rhs_H);
      if(verbose) std::cout << "ERROR ESTIMATOR: computation time: " << boost::timer::format(atimer.elapsed(), 2, "%t") << std::endl;
      errorEstimate_E = sol_H;

      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HE.reset(new Functional_HE(args...));
      F_EE.reset(new Functional_EE(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EE> ass_EE;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };



  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSystemSolver_H, class LinearSystemSolver_L, class LinearSolverA_H, class LinearSolverA_L, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class MartinsErrorEstimator : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    template <class AnsatzVars, class TestVars, class OriginVars> using EstimatorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;
    template <class AnsatzVars, class TestVars, class OriginVars> using LumpedFunctional = Functional<AnsatzVars,TestVars,OriginVars,true>;
    typedef ErrorEstimator_Detail::Traits<EstimatorFunctional,VariableSetDescription,ExtensionVariableSetDescription> FullTraits;
    typedef ErrorEstimator_Detail::Traits<LumpedFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
  public:
    // Functionals
    typedef typename FullTraits::Functional_HH Functional_HH;
    typedef typename Traits::Functional_HE Functional_HE;
    typedef typename Traits::Functional_EH Functional_EH;
    typedef typename Traits::Functional_EE Functional_EE;

    // assemblers
    typedef typename Traits::Assembler_EE Ass_EE;
    typedef typename Traits::Assembler_HE Ass_HE;
    typedef typename Traits::Assembler_EH Ass_EH;
    typedef typename FullTraits::Assembler_HH Ass_HH;

    // operators
    typedef typename FullTraits::Operator_HH H_HH;
    typedef typename Traits::Operator_HE H_HE;
    typedef typename Traits::Operator_EH H_EH;
    typedef typename Traits::Operator_EE H_EE;

    typedef typename FullTraits::B_HH B_HH;
    typedef typename Traits::B_HE B_HE;
    typedef typename Traits::B_EH B_EH;
    typedef typename Traits::B_EE B_EE;

    typedef typename FullTraits::Lyy_HH Hyy_HH;
    typedef typename Traits::Lyy_EH Hyy_EH;
    typedef typename Traits::Lyy_HE Hyy_HE;
    typedef typename Traits::Lyy_EE Hyy_EE;

    typedef typename FullTraits::A_HH A_HH;
    typedef typename Traits::A_HE A_HE;
    typedef typename Traits::A_EE A_EE;

    typedef typename FullTraits::AT_HH AT_HH;
    typedef typename Traits::AT_HE AT_HE;
    typedef typename Traits::AT_EE AT_EE;

    typedef typename Traits::Scalar Scalar;
    static constexpr int dim = VariableSetDescription::Grid::dimension;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;
    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Entity Cell;



    MartinsErrorEstimator(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
        ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
        ass_HE.reset(new Ass_HE(variableSetDescription.spaces));
        ass_EH.reset(new Ass_EH(extensionVariableSetDescription.spaces));
        ass_HH.reset(new Ass_HH(variableSetDescription.spaces));
      }

    virtual ~MartinsErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename Traits::VarSet> const& xl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(x_);
      Bridge::Vector<typename Traits::ExtensionVarSet> xe(extensionVariableSetDescription);
      boost::timer::cpu_timer atimer;
      assembleAt(xl,xe);
      if(verbose) std::cout << "ERROR ESTIMATOR: assembly time: " << boost::timer::format(atimer.elapsed(), 2, "%t")  << std::endl;
//      Bridge::Vector<typename Traits::VarSet> const& dxl = dynamic_cast<const Bridge::Vector<typename Traits::VarSet>&>(dx_);

      /***************************************************************************************/
      H_EE h_EE(*ass_EE,false);
      H_EH h_EH(*ass_EH,false);
      H_HH h_HH(*ass_HH,false);
      H_HE h_HE(*ass_HE,false);

      typename Traits::ExtensionVector rhs_h(ass_EE->rhs());
      rhs_h *= -1.0;
      typename Traits::ExtensionVector sol_h(Traits::ExtensionVector_Initializer::init(extensionVariableSetDescription));

      LinearSystemSolver_H(h_EE,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol_h,rhs_h);

      typename Traits::Vector rhs_l(Traits::Vector_Initializer::init(variableSetDescription)),
          sol_l(Traits::Vector_Initializer::init(variableSetDescription));
      h_EH.applyscaleadd(-1.0,sol_h,rhs_l);

      LinearSystemSolver_L(h_HH,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(sol_l,rhs_l);

      // compute w_y
      B_HH bll(*ass_HH,false);
      B_HE blh(*ass_HE,false);
      B_EH bhl(*ass_EH,false);
      B_EE bhh(*ass_EE,false);
      A_HH all(*ass_HH,false);
      A_HE alh(*ass_HE,false);
      A_EE ahh(*ass_EE,false);

      typename Traits::VectorP rhsP_l(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorP rhsP_h(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));

      typename Traits::VectorU wu_l(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorU wu_h(Traits::ExtensionVectorU_Initializer::init(extensionVariableSetDescription));

      at_c<0>(wu_l.data) = at_c<Traits::controlId>(sol_l.data);
      at_c<0>(wu_h.data) = at_c<Traits::controlId>(sol_h.data);
      bhh.applyscaleadd(-1.0,wu_h,rhsP_h);
      blh.applyscaleadd(-1.0,wu_l,rhsP_h);
      bhl.applyscaleadd(-1.0,wu_h,rhsP_l);
      bll.applyscaleadd(-1.0,wu_l,rhsP_l);

      typename Traits::VectorY wy_l(Traits::VectorY_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorY wy_h(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));

      LinearSolverA_L(all,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(wy_l,rhsP_l);
      alh.applyscaleadd(-1.0,wy_l,rhsP_h);
      LinearSolverA_H(ahh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(wy_h,rhsP_h);

      // compute w_y
      typename Traits::VectorY rhsY_l(Traits::VectorY_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorY rhsY_h(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));

      typename Traits::VectorP wp_l(Traits::VectorP_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorP wp_h(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));

      Hyy_HH hyyll(*ass_HH,false);
      Hyy_EH hyyhl(*ass_EH,false);
      Hyy_HE hyylh(*ass_HE,false);
      Hyy_EE hyyhh(*ass_EE,false);

      AT_HH atll(*ass_HH,false);
      AT_HE atlh(*ass_HE,false);
      AT_EE athh(*ass_EE,false);

      hyyhh.applyscaleadd(-1.0,wy_h,rhsY_h);
      hyylh.applyscaleadd(-1.0,wy_l,rhsY_h);
      hyyhl.applyscaleadd(-1.0,wy_h,rhsY_l);
      hyyll.applyscaleadd(-1.0,wy_l,rhsY_l);

      LinearSolverA_L(atll,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(wp_l,rhsY_l);
      atlh.applyscaleadd(-1.0,wp_l,rhsY_h);
      LinearSolverA_H(athh,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(wp_h,rhsY_h);

      // compute error indicator
      typename Traits::ExtensionVector weights(Traits::ExtensionVector_Initializer::init(extensionVariableSetDescription));
      at_c<Traits::controlId>(weights.data) = at_c<0>(wu_h.data);
      at_c<Traits::stateId>(weights.data) = at_c<0>(wy_h.data);
      at_c<Traits::adjointId>(weights.data) = at_c<0>(wp_h.data);

      typename Traits::ExtensionVector errorIndicator(Traits::ExtensionVector_Initializer::init(extensionVariableSetDescription));
      typename Traits::Vector tmp(Traits::Vector_Initializer::init(variableSetDescription)),
          tmp2(Traits::Vector_Initializer::init(variableSetDescription));

      h_EE.apply(weights,errorIndicator);
      h_EH.apply(weights,tmp);
      LinearSystemSolver_L(h_HH,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(tmp2,tmp);
      h_HE.applyscaleadd(-1.0,tmp2,errorIndicator);

      typename ExtensionVariableSetDescription::VariableSet errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet errorEstimate_H(variableSetDescription);      // estimate error in state variable
      errorEstimate_E = errorIndicator;

      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());

      eeAssembler.assemble(linearization(energyError,xl.get()), EnergyErrorAssembler::RHS);
      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
      typename EnergyError::AnsatzVars::VariableSet mde( energyError.getVariableSetDescription() );
      mde = distError;
      std::cout << "nError: " << distError.dim() << std::endl;
      //std::cout << "||mde||=" << myNorm(mde) << std::endl;

      if(verbose)
      {
        std::string name = "errorDistribution_";
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      CellIterator cend = variableSetDescription.gridView.template end<0>();
      for (CellIterator ci=variableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
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
      F_HE.reset(new Functional_HE(args...));
      F_EH.reset(new Functional_EH(args...));
      F_EE.reset(new Functional_EE(args...));
      F_HH.reset(new Functional_HH(args...));
    }

  private:
    void assembleAt(const Bridge::Vector<typename Traits::VarSet>& xl, const Bridge::Vector<typename Traits::ExtensionVarSet>& xe)
    {
      ass_HE->assemble(linearization(*F_HE,xl.get()));
      ass_EH->assemble(linearization(*F_EH,xe.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
      ass_HH->assemble(linearization(*F_HH,xl.get()));
    }

    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EH> F_EH;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Functional_HH> F_HH;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EH> ass_EH;
    std::unique_ptr<Ass_EE> ass_EE;
    std::unique_ptr<Ass_HH> ass_HH;
    NormFunctional& normFunctional;
    VariableSetDescription& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    Scalar squaredFraction;
    Scalar totalErrorSquared;
    std::vector<std::pair<double,int> > errorDistribution;
    bool verbose, fast;
  };



  template <template <class,class,class,bool> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace,
  class NormFunctional, class LinearSolverLA, class LinearSolverHA, class LinearSolverLM=LinearSolverLA, class LinearSolverHM=LinearSolverHA,
  bool lumpM=false, template <class> class RefinementStrategy = Adaptivity::ErrorEquilibration>
  class HierarchicalBasisErrorEstimator3 : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
    static constexpr int noOfVariables = ExtensionVariableSetDescription::noOfVariables;
    template <class AnsatzVars, class TestVars, class OriginVars> using ErrorFunctional = Functional<AnsatzVars,TestVars,OriginVars,false>;

    typedef ErrorEstimator_Detail::Traits<ErrorFunctional,VariableSetDescription,ExtensionVariableSetDescription> Traits;
    static constexpr int controlId = Traits::controlId;
    static constexpr int stateId = Traits::stateId;
    static constexpr int adjointId = Traits::adjointId;

  public:
    typedef typename Traits::Scalar Scalar;
    typedef Functional<VariableSetDescription,VariableSetDescription,VariableSetDescription,lumpM>                             Functional_HH;
    typedef Functional<ExtensionVariableSetDescription,VariableSetDescription,ExtensionVariableSetDescription,lumpM>         Functional_EH;
    typedef Functional<VariableSetDescription,ExtensionVariableSetDescription,VariableSetDescription,lumpM>                    Functional_HE;
    typedef Functional<ExtensionVariableSetDescription,ExtensionVariableSetDescription,ExtensionVariableSetDescription,lumpM>  Functional_EE;
    typedef HierarchicErrorEstimator<LinearizationAt<Functional_HH>,ExtensionVariableSetDescription,ExtensionVariableSetDescription> ErrorEstimator;
    //    typedef HierarchicErrorEstimator<LinearizationAt<Functional_HH>,ExtensionVariableSetDescription,ExtensionVariableSetDescription,HierarchicErrorEstimatorDetail::TakeAllD2<LinearizationAt<Functional_HH> > > ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> Assembler;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HH> >                                 Ass_HH;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_HE> >                                 Ass_HE;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EH> >                                 Ass_EH;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional_EE> >                                 Ass_EE;

    typedef typename Traits::AT_EE      AT_EE;
    typedef typename Traits::AT_HE      AT_HE;

    typedef typename Traits::A_HH       A_HH;
    typedef typename Traits::A_EE       A_EE;
    typedef typename Traits::A_EH       A_EH;

    typedef typename Traits::Luu_HH     M_HH;
    typedef typename Traits::Luu_EE     M_EE;

    typedef typename Traits::BT_EH      BT_EH;
    typedef typename Traits::BT_EE      BT_EE;

    typedef typename Traits::B_HH       B_HH;
    typedef typename Traits::B_HE       B_HE;
    typedef typename Traits::B_EH       B_EH;
    typedef typename Traits::B_EE       B_EE;

    static constexpr int dim = Traits::dim;

    HierarchicalBasisErrorEstimator3(NormFunctional& normFunctional_, VariableSetDescription& variableSetDescription_, ExtensionVariableSetDescription& extensionVariableSetDescription_,
        ExtensionSpace& extensionSpace_, Scalar fraction=0.7, bool verbose_=false, bool fast_ = false)
    : RefinementStrategy<typename VariableSetDescription::Grid>(extensionSpace_.gridManager(),fraction,verbose),/*lf(lf_), lvars(lvars_),*/ normFunctional(normFunctional_), variableSetDescription(variableSetDescription_), extensionVariableSetDescription(extensionVariableSetDescription_),
      extensionSpace(extensionSpace_), assembler(extensionVariableSetDescription.spaces),
      squaredFraction(fraction*fraction), verbose(verbose_), fast(fast_)
      {
      initAssemblers();
      }

    virtual ~HierarchicalBasisErrorEstimator3(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step, AbstractVector const&)
    {
      boost::timer::cpu_timer overallTimer;
      using boost::fusion::at_c;
      if(verbose) std::cout << "ERROR ESTIMATOR: Start." << std::endl;
      Bridge::Vector<typename VariableSetDescription::VariableSet> const& xl = dynamic_cast<const Bridge::Vector<typename VariableSetDescription::VariableSet>&>(x_);
      Bridge::Vector<typename ExtensionVariableSetDescription::VariableSet> xe(extensionVariableSetDescription);
      assembleAt(xl,xe);
      Bridge::Vector<typename VariableSetDescription::VariableSet> const& dxl = dynamic_cast<const Bridge::Vector<typename VariableSetDescription::VariableSet>&>(dx_);

      /***************************************************************************************/
      // estimate error in dual variable
      typename ExtensionVariableSetDescription::VariableSet tmpRep(extensionVariableSetDescription), errorEstimate_E(extensionVariableSetDescription);
      typename VariableSetDescription::VariableSet tmpRepL(variableSetDescription), errorEstimate_H(variableSetDescription);

      AT_EE atee(*ass_EE,false);
      AT_HE athe(*ass_HE,false);

      typename Traits::ExtensionVectorP dpe(Traits::ExtensionVectorP_Initializer::init(extensionVariableSetDescription));
//      typename Traits::VectorP dph(at_c<adjointId>(dxl.get().data));
      typename Traits::VectorP dph(Traits::VectorP_Initializer::init(variableSetDescription));
      at_c<0>(dph.data) = at_c<adjointId>(dxl.get().data).coefficients();
      typename Traits::ExtensionVectorY rye(ass_EE->template rhs<stateId,stateId+1>());
      rye *= -1.;

      athe.applyscaleadd(-1.0,dph,rye);
      JacobiPreconditioner<AT_EE>(atee).apply(dpe,rye);

      at_c<adjointId>(errorEstimate_E.data) = at_c<0>(dpe.data);

      /***************************************************************************************/
      // estimate error in control variable

      M_HH mhh(*ass_HH,false);
      M_EE mee(*ass_EE,false);

      BT_EH bteh(*ass_EH,false);
      BT_EE btee(*ass_EE,false);

      typename Traits::VectorU ruh(ass_HH->template rhs<controlId,controlId+1>()),
                                duh(Traits::VectorU_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorU rue(ass_EE->template rhs<controlId,controlId+1>()),
                                         due(Traits::ExtensionVectorU_Initializer::init(extensionVariableSetDescription));
      ruh *= -1.;
      rue *= -1.;

      bteh.applyscaleadd(-1.0,dpe,ruh);
      btee.applyscaleadd(-1.0,dpe,rue);

      JacobiPreconditioner<M_HH>(mhh).apply(duh,ruh);
      JacobiPreconditioner<M_EE>(mee).apply(due,rue);

      at_c<controlId>(errorEstimate_E.data) = at_c<0>(due.data);
      at_c<controlId>(errorEstimate_H.data) = at_c<0>(duh.data);
      /***************************************************************************************/
      // estimate error in state variable

      B_HH bhh(*ass_HH,false);
      B_EH beh(*ass_EH,false);
      B_HE bhe(*ass_HE,false);
      B_EE bee(*ass_EE,false);

      A_HH ahh(*ass_HH,false);
      A_EE aee(*ass_EE,false);
      A_EH aeh(*ass_EH,false);

      typename Traits::VectorP rph(ass_HH->template rhs<adjointId,adjointId+1>());
      typename Traits::ExtensionVectorP rpe(ass_EE->template rhs<adjointId,adjointId+1>());
      rph *= -1.;
      rpe *= -1.;

      bhh.applyscaleadd(-1.0,duh,rph);
      bhe.applyscaleadd(-1.0,duh,rpe);
      beh.applyscaleadd(-1.0,due,rph);
      bee.applyscaleadd(-1.0,due,rpe);

      typename Traits::VectorY dyh(Traits::VectorY_Initializer::init(variableSetDescription));
      typename Traits::ExtensionVectorY dye(Traits::ExtensionVectorY_Initializer::init(extensionVariableSetDescription));

      JacobiPreconditioner<A_EE>(aee).apply(dye,rpe);
      aeh.applyscaleadd(-1.0,dye,rph);
      JacobiPreconditioner<A_HH>(ahh).apply(dyh,rph);

      at_c<stateId>(errorEstimate_E.data) = at_c<0>(dye.data);
      at_c<stateId>(errorEstimate_H.data) = at_c<0>(dyh.data);
      /***************************************************************************************/
      // Transfer error indicators to cells.
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
      errorDistribution.clear();
      errorDistribution.resize(is.size(0),std::make_pair(0.0,0));
      Scalar errLevel(0.0), minRefine(0.2);

      //      typedef ErrorDistribution<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      typedef ErrorDistribution3<NormFunctional,ExtensionVariableSetDescription> EnergyError;
      //      EnergyError energyError(normFunctional,xl.get(),errorEstimate_E);
      EnergyError energyError(normFunctional,xl.get(),errorEstimate_H,errorEstimate_E);
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
        name += std::to_string(step);
        writeVTKFile(mde.descriptions.gridView,mde.descriptions,mde,name);
      }
      // Fehler bzgl L_xx : ( sum_i estSol_i * (L_xx * estSol)_i )^{1/2}  (nur die x - Komponente, nicht p)
      // Lokalisierung: einfachste Idee: v_i = estSol_i * (L_xx * estSol)_i
      // Aufteilen von v_i auf die einzelnen Zellen
      auto cend = extensionVariableSetDescription.gridView.template end<0>();
      for (auto ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));

      totalErrorSquared = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, HierarchicErrorEstimator_Detail::add);

      if(verbose) std::cout << "overall error estimation time: " << boost::timer::format(overallTimer.elapsed()) << std::endl;
    }

    void refineGrid()
    {
      if(totalErrorSquared > 0) this->refineGrid_impl(errorDistribution,sqrt(totalErrorSquared));
//      Scalar tol = squaredFraction*totalErrorSquared;
//      Scalar sum = 0;
//      std::sort(errorDistribution.begin(), errorDistribution.end(), [](std::pair<double,int> const& a, std::pair<double,int> const& b) -> bool
//          {
//        return a.first > b.first;
//          });
//      std::vector<size_t> cellIds;
//      for(size_t i=0; i<errorDistribution.size(); ++i)
//      {
//        sum += errorDistribution[i].first;
//        cellIds.push_back(errorDistribution[i].second);
//        if(sum > tol) break;
//      }
//
//      forEachCell(variableSetDescription.gridView, [&](Cell const& cell)
//          {
//        if(std::find(cellIds.begin(),cellIds.end(),variableSetDescription.gridView.indexSet().index(cell))!=cellIds.end()) extensionSpace.gridManager().mark(1,cell);
//          });
//      extensionSpace.gridManager().adaptAtOnce();
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
      F_HH.reset(new Functional_HH(args...));
      F_HE.reset(new Functional_HE(args...));
      F_EH.reset(new Functional_EH(args...));
      F_EE.reset(new Functional_EE(args...));
    }

  private:
    void initAssemblers()
    {
      ass_HH.reset(new Ass_HH(extensionVariableSetDescription.spaces));
      ass_HE.reset(new Ass_HE(extensionVariableSetDescription.spaces));
      ass_EH.reset(new Ass_EH(extensionVariableSetDescription.spaces));
      ass_EE.reset(new Ass_EE(extensionVariableSetDescription.spaces));
    }

    void assembleAt(const Bridge::Vector<typename VariableSetDescription::VariableSet>& xh, const Bridge::Vector<typename ExtensionVariableSetDescription::VariableSet>& xe)
    {
      ass_HH->assemble(linearization(*F_HH,xh.get()));
      ass_HE->assemble(linearization(*F_HE,xh.get()));
      ass_EH->assemble(linearization(*F_EH,xe.get()));
      ass_EE->assemble(linearization(*F_EE,xe.get()));
    }

    std::unique_ptr<Functional_HH> F_HH;
    std::unique_ptr<Functional_HE> F_HE;
    std::unique_ptr<Functional_EH> F_EH;
    std::unique_ptr<Functional_EE> F_EE;
    std::unique_ptr<Ass_HH> ass_HH;
    std::unique_ptr<Ass_HE> ass_HE;
    std::unique_ptr<Ass_EH> ass_EH;
    std::unique_ptr<Ass_EE> ass_EE;
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
    bool verbose, fast;
    // AdjustRHS<CoefficientVector> adjustRHS;
  };
}

#endif
