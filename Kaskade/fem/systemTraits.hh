#ifndef SYSTEMTRAITS_HH
#define SYSTEMTRAITS_HH

#include "fem/istlinterface.hh"
#include "mg/multiGridSolver.hh"

/// \cond undocumented

namespace Dune
{
  template <class,class> class Preconditioner;
}

namespace Kaskade
{
  template <class,class,class,class> class PreconditionerFactory;

  namespace SystemTraits_Detail
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
  }


  template <class Functional>
  constexpr int getStateId()
  {
    typedef decltype(SystemTraits_Detail::getStateIdImpl<Functional>(nullptr)) type;
    return type::value;
  }

  template <class Functional>
  constexpr int getControlId()
  {
    typedef decltype(SystemTraits_Detail::getControlIdImpl<Functional>(nullptr)) type;
    return type::value;
  }

  template <class Functional>
  constexpr int getAdjointId()
  {
    typedef decltype(SystemTraits_Detail::getAdjointIdImpl<Functional>(nullptr)) type;
    return type::value;
  }


  template <class Functional, class Assembler>
  struct OptimalControlTraits
  {
    // variable ids
    static constexpr int yIdx = getStateId<Functional>();
    static constexpr int uIdx = getControlId<Functional>();
    static constexpr int pIdx = getAdjointId<Functional>();

    static constexpr int stateId = yIdx;
    static constexpr int adjointId = pIdx;
    static constexpr int controlId = uIdx;

    // grid
    typedef typename Assembler::Grid Grid;

    // spaces
    typedef typename Assembler::Spaces Spaces;

    // variable set descriptions
    typedef typename Assembler::AnsatzVariableSet AnsatzVariableSetDescription;
    typedef typename Assembler::TestVariableSet   TestVariableSetDescription;

    // spaces associated with variables
    static constexpr int stateSpaceId   = boost::fusion::result_of::value_at_c<typename AnsatzVariableSetDescription::Variables,stateId>::type::spaceIndex;
    static constexpr int controlSpaceId = boost::fusion::result_of::value_at_c<typename AnsatzVariableSetDescription::Variables,controlId>::type::spaceIndex;
    static constexpr int adjointSpaceId = boost::fusion::result_of::value_at_c<typename AnsatzVariableSetDescription::Variables,adjointId>::type::spaceIndex;

    // components of each variable
    static constexpr int stateComponents   = AnsatzVariableSetDescription::template Components<stateId>::m;
    static constexpr int controlComponents = AnsatzVariableSetDescription::template Components<controlId>::m;
    static constexpr int adjointComponents = AnsatzVariableSetDescription::template Components<adjointId>::m;
    static_assert(stateComponents == adjointComponents, "state and adjoint variables should have the same number of components");

    // operators
    typedef AssembledGalerkinOperator<Assembler> SystemOperator;
    typedef AssembledGalerkinOperator<Assembler,pIdx,pIdx+1,yIdx,yIdx+1> StateOperator;
    typedef AssembledGalerkinOperator<Assembler,pIdx,pIdx+1,uIdx,uIdx+1> ControlOperator;
    typedef AssembledGalerkinOperator<Assembler,yIdx,yIdx+1,pIdx,pIdx+1> AdjointOperator;
    typedef AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,uIdx,uIdx+1> NormUOperator;
    typedef AssembledGalerkinOperator<Assembler,yIdx,yIdx+1,yIdx,yIdx+1> NormYOperator;
    typedef AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,yIdx,yIdx+1> NormUYOperator;
    typedef AssembledGalerkinOperator<Assembler,yIdx,yIdx+1,uIdx,uIdx+1> NormYUOperator;
    typedef AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,pIdx,pIdx+1> ControlOperatorT;

    // domain and range
    typedef typename SystemOperator::Domain Domain;
    typedef typename SystemOperator::Range Range;
    typedef typename Functional::AnsatzVars::Scalar Scalar;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<yIdx,yIdx+1> CreateVectorY;
    typedef typename CreateVectorY::type VectorY;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<uIdx,uIdx+1> CreateVectorU;
    typedef typename CreateVectorU::type VectorU;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<pIdx,pIdx+1> CreateVectorP;
    typedef typename CreateVectorP::type VectorP;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<> CreateCoefficientVector;
    typedef typename CreateCoefficientVector::type CoefficientVector;
    typedef VectorY StateCoefficientVector;
    typedef VectorP AdjointCoefficientVector;
    typedef VectorU ControlCoefficientVector;
    typedef CreateVectorY StateCoefficientVectorCreator;
    typedef CreateVectorP AdjointCoefficientVectorCreator;
    typedef CreateVectorU ControlCoefficientVectorCreator;

    // preconditioners
    typedef Dune::Preconditioner<VectorY,VectorP> StatePreconditioner;
    typedef Dune::Preconditioner<VectorP,VectorY> AdjointPreconditioner;
    typedef Dune::Preconditioner<VectorU,VectorU> NormUPreconditioner;

    // preconditioner factories
    typedef PreconditionerFactory<Assembler,VectorY,VectorY,Scalar> NormYPreconditionerFactory;
    typedef PreconditionerFactory<Assembler,VectorU,VectorU,Scalar> NormUPreconditionerFactory;
    typedef PreconditionerFactory<Assembler,VectorP,VectorY,Scalar> AdjointPreconditionerFactory;
    typedef PreconditionerFactory<Assembler,VectorY,VectorP,Scalar> StatePreconditionerFactory;

    // multi grid solvers
    typedef MultiGridSolver<Grid,stateComponents> StateMGSolver;
    typedef MultiGridSolver<Grid,controlComponents> ControlMGSolver;

    // multi grid preconditioners
    typedef MultiGridPreconditioner<Grid,stateComponents> StateMGPreconditioner;
    typedef MultiGridPreconditioner<Grid,controlComponents> ControlMGPreconditioner;
  };

}


/// \endcond

#endif // SYSTEMTRAITS_HH
