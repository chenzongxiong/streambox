#ifndef ERRORESTIMATIONTRAITS_HH
#define ERRORESTIMATIONTRAITS_HH

#include <cmath>
#include <memory>

#include "fem/assemble.hh"
#include "fem/variables.hh"

#include "abstract_interface.hh"
#include "errorDistribution.hh"

namespace Kaskade
{
  namespace ErrorEstimator_Detail
  {
    double add(const double& x1, const std::pair<double,int>&  x2)
    {
      return std::fabs(x1) + std::fabs(x2.first);
    }
  }

  template < template <class,class,class> class Functional, class VariableSetDescription, class ExtensionVariableSetDescription>
  struct ErrorEstimationTraits
  {    
    typedef VariableSet<VariableSetDescription> VarSet;
    typedef VariableSet<ExtensionVariableSetDescription> ExtensionVarSet;

    typedef Functional<VariableSetDescription,VariableSetDescription,VariableSetDescription>                            Functional_HH;
    typedef Functional<VariableSetDescription,ExtensionVariableSetDescription,VariableSetDescription>                   Functional_HE;
    typedef Functional<ExtensionVariableSetDescription,VariableSetDescription,ExtensionVariableSetDescription>          Functional_EH;
    typedef Functional<ExtensionVariableSetDescription,ExtensionVariableSetDescription,ExtensionVariableSetDescription> Functional_EE;


    typedef typename Functional_HH::AnsatzVars::VariableSet DomainElement;
    typedef typename Functional_EE::AnsatzVars::VariableSet DomainExtensionElement;
    typedef typename Functional_HH::TestVars::VariableSet ImageElement;
    typedef typename Functional_EE::TestVars::VariableSet ImageExtensionElement;
    typedef typename Functional_HH::OriginVars::VariableSet OriginElement;
    typedef typename Functional_EE::OriginVars::VariableSet OriginExtensionElement;
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

  template <class VariableSetDescription, class ErrorRepresentation, class ExtensionSpace, template <class> class RefinementStrategy>
  class ErrorEstimatorBase : public AbstractHierarchicalErrorEstimator, public RefinementStrategy<typename VariableSetDescription::Grid>
  {
  public:
    ErrorEstimatorBase(ExtensionSpace& extensionSpace_, double fraction, bool verbose = false)
      : RefinementStrategy<typename VariableSetDescription::Grid>( extensionSpace_.gridManager(), fraction, verbose ),
        extensionSpace(extensionSpace_)
    {}

    double estimatedAbsoluteError() const final { return sqrt(squaredError); }

    size_t gridSize() const final { return extensionSpace.gridManager().grid().size(0); }

    void refineGrid() final { if( squaredError > 0 ) this->refineGrid_impl(*mde,errorDistribution,squaredError);  }

  protected:
    ExtensionSpace& extensionSpace;
    std::vector<std::pair<double,int> > errorDistribution = {}; // store error and corresponding cell index
    std::unique_ptr<ErrorRepresentation> mde;
    double squaredError = 0;
  };
}

#endif // ERRORESTIMATIONTRAITS_HH
