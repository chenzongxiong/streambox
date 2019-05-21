#ifndef COMP_STEP
#define COMP_STEP

#include <memory> // std::unique_ptr
#include <vector>

#include <boost/timer/timer.hpp>
#include <boost/fusion/include/at_c.hpp>

#include "algorithm/opt_interface.hh"
#include "algorithm/dune_bridge.hh"
#include "algorithm/newton_bridge.hh"
#include "linalg/cg.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/icc0precond.hh"
#include "linalg/directPreconditioner.hh"
#include "linalg/blockDiagonalSchurPreconditioner.hh"
//#include "linalg/preconditionerFactory.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/chebyshev.hh"
#include "fem/systemTraits.hh"
#include "mg/multiGridSolver.hh"
#include "utilities/save_file.hh"
#include "lagrangeLinearization.hh"

namespace Kaskade
{
  template <class,class,class,class> class PreconditionerFactory;

  template <class Functional, class Assembler> class NormalStepPreconditioner;
  template <class Functional, class Assembler, int components=1> class TangentSpacePreconditioner;
  template <class Functional, class Assembler, int components=1> class TangentSpacePreconditioner2;
template <class Functional, class Assembler, int components=1> class NormalStepPreconditioner2;
template <class Functional, class Assembler, int components=1> class NormalStepPreconditioner3;
  template <class Functional, class Assembler, int components=1,bool exactConstraint=false> class InexactTangentSpacePreconditioner;


  template <class X, class Xstar, int yIdx, int uIdx, int pIdx>
  class Scaledl2ScalarProduct : public DualPairing<X,Xstar>
  {
    typedef typename DualPairing<X,Xstar>::field_type Scalar;
  public:
    Scaledl2ScalarProduct(Scalar sy_, Scalar su_, Scalar sp_) : sy(sy_), su(su_), sp(sp_) {}

    virtual ~Scaledl2ScalarProduct(){}

    virtual Scalar operator()(X const& x, Xstar const& y) const { return sy * prod<yIdx>(x,y) + su * prod<uIdx>(x,y) + sp * prod<pIdx>(x,y); }

  private:
    template <int i> Scalar prod(X const& x, Xstar const& y) const { return boost::fusion::at_c<i>(x.data) * boost::fusion::at_c<i>(y.data); }

    Scalar sy = 1.0, su = 1.0, sp = 1.0;
  };


//  template <class Functional, class NormalStepAssembler, class TangentialStepAssembler>
//  class StollReesPreconditioner : public Dune::Preconditioner<typename OptimalControlTraits<Functional,NormalStepAssembler>::Domain, typename OptimalControlTraits<Functional,NormalStepAssembler>::Range>
//  {
//    typedef OptimalControlTraits<Functional,NormalStepAssembler> Traits;
//    typedef OptimalControlTraits<Functional,TangentialStepAssembler> TTraits;
//    typedef typename Functional::Scalar Scalar;
//    typedef MatrixAsTriplet<Scalar> Matrix;
//    typedef SchurComplement_v1<NormalStepAssembler, TangentialStepAssembler,
//                               IstlInterfaceDetail::BlockInfo<Traits::adjointId,Traits::adjointId+1,Traits::stateId,Traits::stateId+1>,
//                               IstlInterfaceDetail::BlockInfo<Traits::stateId,Traits::stateId+1,Traits::stateId,Traits::stateId+1>
//                              > SchurComplement;
//  public:
//    typedef typename Traits::Domain Domain;
//    typedef typename Traits::Range Range;
//    typedef NormalStepAssembler Assembler;

//    StollReesPreconditioner()
//      : Myy(nullptr), Muu(nullptr), A(nullptr), B(nullptr), S(nullptr), chebyy(nullptr), chebuu(nullptr)
//    {}

//    StollReesPreconditioner(NormalStepAssembler const& normalStepAssembler, TangentialStepAssembler const& tangentialStepAssembler, size_t maxSteps = 50) : StollReesPreconditioner()
//    {
//      update(normalStepAssembler,tangentialStepAssembler,maxSteps);
//    }

//    void update(NormalStepAssembler const& normalStepAssembler, TangentialStepAssembler const& tangentialStepAssembler, size_t maxSteps = 50)
//    {
//      Myy.reset( new typename TTraits::NormYOperator(tangentialStepAssembler) );
//      Muu.reset( new typename TTraits::NormUOperator(tangentialStepAssembler) );
//      A.reset( new typename Traits::StateOperator(normalStepAssembler) );
//      B.reset( new typename Traits::ControlOperator(normalStepAssembler) );
//      S.reset( new SchurComplement(normalStepAssembler,tangentialStepAssembler) );
//      chebyy.reset( new ChebyshevPreconditioner<typename TTraits::NormYOperator>(*Myy,maxSteps) );
//      chebuu.reset( new ChebyshevPreconditioner<typename TTraits::NormUOperator>(*Muu,maxSteps) );
//      chebyy->initForMassMatrix_TetrahedralQ1Elements();
//      chebuu->initForMassMatrix_TetrahedralQ1Elements();
//    }

//    virtual void pre(Domain&, Range&){}
//    virtual void post(Domain&){}

//    virtual void apply(Domain& x, Range const& y)
//    {
//      using namespace boost::fusion;

//      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data)), solY(rhsY); solY *= 0;
//      directInverseOperator(*Myy,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(rhsY,solY); solY *= 0.5; // scaling for positive definiteness of non-standard inner product
//      rhsY *= -1;
//      Myy->applyscaleadd(1.0,solY,rhsY);
//      at_c<Traits::yIdx>(x.data) = at_c<0>(rhsY.data);

//      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data)), solU(rhsU); solU *= 0;
//      chebuu->apply(solU,rhsU); solU *= 0.5;
//      rhsU *= -1;
//      Muu->applyscaleadd(1.0,solU,rhsU);
//      at_c<Traits::uIdx>(x.data) = at_c<0>(rhsU.data);

//      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(x.data)), solP(rhsP); solP *= 0;
//      rhsP *= -1;
//      A->applyscaleadd(1.0,solY,rhsP);
//      B->applyscaleadd(1.0,solU,rhsP);
//      at_c<Traits::pIdx>(x.data) = at_c<0>(rhsP.data);
//    }


//    virtual void apply2(Domain& x, Range const& y)
//    {
//      using namespace boost::fusion;

//      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data)), solY(rhsY); solY *= 0;
//      directInverseOperator(*Myy,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(rhsY,solY);
//      //      chebyy->apply(solY,rhsY);

//      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data)), solU(rhsU); solU *= 0;
//      chebuu->apply(solU,rhsU);

//      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data)), solP(rhsP); solP *= 0;
//      A->applyscaleadd(-1.0,solY,rhsP);
//      B->applyscaleadd(-1.0,solU,rhsP);

//      S->solve(solP,rhsP);

//      at_c<Traits::stateId>(x.data) = at_c<0>(solY.data);
//      at_c<Traits::controlId>(x.data) = at_c<0>(solU.data);
//      at_c<Traits::adjointId>(x.data) = at_c<0>(solP.data);
//    }

//  private:
//    std::unique_ptr<typename TTraits::NormYOperator> Myy;
//    std::unique_ptr<typename TTraits::NormUOperator> Muu;
//    std::unique_ptr<typename Traits::StateOperator> A;
//    std::unique_ptr<typename Traits::ControlOperator> B;
//    std::unique_ptr<SchurComplement> S;
//    std::unique_ptr<ChebyshevPreconditioner<typename TTraits::NormYOperator> > chebyy;
//    std::unique_ptr<ChebyshevPreconditioner<typename TTraits::NormUOperator> > chebuu;
//  };


//  template <class Functional,class Assembler>
//  class SchurPreconditioner : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
//  {
//    typedef OptimalControlTraits<Functional,Assembler> Traits;
//    typedef typename Functional::Scalar Scalar;
//    typedef MatrixAsTriplet<Scalar> Matrix;
//    typedef SchurPreconditionerDetail::InvertLumpedMatrix<Scalar> Solver;
//  public:
//    SchurPreconditioner(Assembler const& assembler)
//    : Myy(assembler,true), Muu(assembler,true)
//    {
//      std::cout << "Schurpreconditioner: init schur complement...";
//      // compute approximate schur complement
//      typename Traits::StateOperator A(assembler,false);
//      typename Traits::AdjointOperator At(assembler,false);
//      typename Traits::ControlOperator B(assembler,false);
//      typename Traits::ControlOperatorT Bt(assembler,false);

//      Matrix s0, s1;
//      Matrix MA = A.template get<Matrix>(), MAT = At.template get<Matrix>(), MMyy = Myy.template get<Matrix>();
//      for(size_t i=0; i<MAT.ridx.size(); ++i)
//      {
//        size_t row = MAT.ridx[i];
//        size_t newIndex = std::abs(std::distance(MMyy.ridx.begin(), std::find(MMyy.ridx.begin(), MMyy.ridx.end(), row)));
//        assert(newIndex < MMyy.ridx.size());

//        MAT.data[i] /= MMyy.data[newIndex];
//      }

//      std::vector<std::vector<Scalar> > At_cols;
//      MAT.toColumns(At_cols);
//      auto At_cols2(At_cols);
//      for(size_t i=0; i<At_cols.size(); ++i)
//      {
//        MA.ax(At_cols2[i],At_cols[i]);
//        s0.addColumn(At_cols2[i],i);
//      }

//      Matrix MB = B.template get<Matrix>(), MBT = Bt.template get<Matrix>(), MMuu = Muu.template get<Matrix>();
//      for(size_t i=0; i<MBT.ridx.size(); ++i)
//      {
//        size_t row = MBT.ridx[i];
//        size_t newIndex = std::abs(std::distance(MMuu.ridx.begin(), std::find(MMuu.ridx.begin(), MMuu.ridx.end(), row)));
//        assert(newIndex < MMuu.ridx.size());

//        MBT.data[i] /= MMuu.data[newIndex];
//      }


//      std::vector<std::vector<Scalar> > Bt_cols;
//      MBT.toColumns(Bt_cols);
//      auto Bt_cols2(Bt_cols);
//      for(size_t i=0; i<Bt_cols.size(); ++i)
//      {
//        MB.ax(Bt_cols2[i],Bt_cols[i]);
//        s1.addColumn(Bt_cols2[i],i);
//      }

//      s0 += s1;

//      S.reset(new MatrixRepresentedOperator<Matrix,typename Traits::VectorP, typename Traits::VectorP>(std::move(s0)));
//      std::cout << "done" << std::endl;
//    }

//    ~SchurPreconditioner(){}

//    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
//    virtual void post(typename Traits::Domain&){}

//    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
//    {
//      using namespace boost::fusion;

//      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data)), solY(rhsY);
//      Solver solverY(Myy);
//      solverY.apply(solY,rhsY);

//      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data)), solU(rhsU);
//      Solver solverU(Muu);
//      solverU.apply(solU,rhsU);

//      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data)), solP(rhsP);
//      directInverseOperator(*S,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(rhsP,solP);

//      at_c<Traits::yIdx>(x.data) = at_c<0>(solY.data);
//      at_c<Traits::uIdx>(x.data) = at_c<0>(solU.data);
//      at_c<Traits::pIdx>(x.data) = at_c<0>(solP.data);
//    }

//  private:
//    typename Traits::NormYOperator Myy;
//    typename Traits::NormUOperator Muu;
//    std::unique_ptr<MatrixRepresentedOperator<Matrix,typename Traits::VectorP, typename Traits::VectorP> > S;
//  };


//  template <class PreconditionerFactory,class VariableSet>
//  class PreconditionerAsPDESolver : public AbstractNormalDirection, public Dune::Preconditioner<typename PreconditionerFactory::Operator::Domain, typename PreconditionerFactory::Operator::Range>
//  {
//  public:
//    typedef typename PreconditionerFactory::Operator::Assembler Assembler;
//    typedef typename PreconditionerFactory::Operator::Domain Domain;
//    typedef typename PreconditionerFactory::Operator::Range Range;
//    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> BridgeLinearization;

//    PreconditionerAsPDESolver(PreconditionerFactory& prec_) : prec(prec_)
//    {}

//    virtual ~PreconditionerAsPDESolver(){}

//    virtual void pre(Domain& x, Range& b){}
//    virtual void post(Domain& x){}
//    virtual void apply(Domain& x, Range const& y)
//    {
//      preconditioner->apply(x,y);
//      //      P->apply(x,y);
//    }

//  private:
//    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& linearization)
//    {
//      A.reset(new AssembledGalerkinOperator<Assembler>(dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler()));
//      preconditioner.reset(new NormalStepPreconditioner<typename Assembler::Functional::Functional,Assembler>(dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler()));
//      P.reset(prec.create(*A).release());
//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
//      linearization.evald(*combinedrhs);
//      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*combinedrhs);
//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));

//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
//      y = conrhs;

//      P->apply(x,y);

//      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
//      cor = x;
//      cor *= -1.0;
//    }

//    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction,
//        AbstractLinearization const& lin) const
//    {
//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
//      lin.evald(*combinedrhs);
//      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*combinedrhs);

//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
//      x *= 0.0;
//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));

//      y=conrhs;

//      P->apply(x,y);

//      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);

//      cor = x;

//      cor *= -1.0;
//    }

//    PreconditionerFactory& prec;
//    std::unique_ptr<Dune::Preconditioner<Domain,Range> > P;
//    std::unique_ptr<NormalStepPreconditioner<typename Assembler::Functional::Functional,Assembler> > preconditioner;
//    std::unique_ptr<AssembledGalerkinOperator<Assembler> > A;
//  };

  template <class Operator, class PrecondAssembler, class PreconditionerFactory, class VariableSet>
  class PreconditionerAsNormalSolver : public AbstractNormalDirection, public Dune::Preconditioner<typename Operator::Domain, typename Operator::Range>
  {
  public:
    typedef typename Operator::Assembler Assembler;
    static constexpr int components = Assembler::Grid::dimension;
    typedef typename Assembler::Functional::Functional Functional;
    typedef typename PrecondAssembler::Functional::Functional PreconditionerFunctional;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Operator::Domain Domain;
    typedef typename Operator::Range Range;
    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> BridgeLinearization;

    virtual ~PreconditionerAsNormalSolver() {}

    PreconditionerAsNormalSolver(PreconditionerFunctional const& pf_, PreconditionerFactory& prec_) : pf(pf_), prec(prec_) {}

    virtual void pre(Domain &x, Range &b) {}
    virtual void apply (Domain &v, const Range &d) {
      //normalStepPreconditioner->apply(v,d);
      tangentialStepPreconditioner->apply(v,d);
//      P->apply(v,d);
    }

    virtual void post (Domain &x) {}

  private:
    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& lin, AbstractFunctionSpaceElement*,AbstractFunctionSpaceElement*)
    {
      A.reset(new Operator(dynamic_cast<BridgeLinearization&>(lin).getValidAssembler()));
      std::cout << "init P" << std::endl;
      //if(std::is_same<Operator,NormalStepAssembledGalerkinOperator<Assembler>>::value) std::cout << "normal step operator" << std::endl;
      // else std::cout << "standard operator" << std::endl;
      P.reset(new DirectPreconditioner<Operator>(*A));
      //P.reset(prec.create(*A).release());

      PrecondAssembler assembler(A->getAssembler().spaces());
      assembler.assemble(linearization(pf,Bridge::getImpl<VariableSet>(lin.getOrigin())));
      std::cout << "init nsp" << std::endl;
      normalStepPreconditioner.reset(new NormalStepPreconditioner<PreconditionerFunctional,PrecondAssembler>(assembler));
      tangentialStepPreconditioner.reset(new InexactTangentSpacePreconditioner<PreconditionerFunctional,PrecondAssembler,components>(assembler));
      std::cout << "done init p" << std::endl;
      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> derivativerhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
      *derivativerhs *= 0.0;
      *constraintrhs *= 0.0;
      lin.evald(*combinedrhs);
      derivativerhs->axpy(1.0,*combinedrhs,"primal");
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      VariableSet& derrhs=Bridge::getImpl<VariableSet>(*derivativerhs);
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);
      using namespace boost::fusion;

      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));

      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
      y = derrhs;

      P->apply(x,y);

      VariableSet& adj=Bridge::getImpl<VariableSet>(adjointCorrection);
      adj = x;
      adj *= -1.0;

      std::cout << "adjoint correction: " << x*x << std::endl;

      x *= 0;
      y = conrhs;
      P->apply(x,y);


      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;
    }


    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction, 
        AbstractLinearization const& lin, AbstractFunctionSpaceElement* residual) const
    {

      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
      *constraintrhs *= 0.0;
      lin.evald(*combinedrhs);
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      if(residual != nullptr) *constraintrhs += *residual;
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);
      std::cout << "RHS 2=" << boost::fusion::at_c<2>(conrhs.data).coefficients()*boost::fusion::at_c<2>(conrhs.data).coefficients() << std::endl;

      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
      x *= 0.0;
      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));

      y=conrhs;

      P->apply(x,y);

      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;
    }

    PreconditionerFunctional const& pf;
    PreconditionerFactory& prec;
    std::unique_ptr<Dune::Preconditioner<Domain,Range> > P;
    std::unique_ptr<Operator> A;
    std::unique_ptr<NormalStepPreconditioner<PreconditionerFunctional,PrecondAssembler> > normalStepPreconditioner;
    std::unique_ptr<InexactTangentSpacePreconditioner<PreconditionerFunctional,PrecondAssembler,components> > tangentialStepPreconditioner;
};

  template <class Assembler_, class PrecondAssembler, class Domain_, class Range_, class VariableSet, int components=Assembler_::Grid::dimension>
  class PPCGAsNormalSolver : public AbstractNormalDirection, public Dune::Preconditioner<Domain_, Range_>
  {
  public:
    typedef Domain_ Domain;
    typedef Range_ Range;
    typedef Assembler_ Assembler;
    typedef typename Assembler::Functional::Functional Functional;
    typedef typename PrecondAssembler::Functional::Functional PreconditionerFunctional;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Traits::Scalar Scalar;
    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> BridgeLinearization;
    typedef IterateType::CG<Domain,Range> Solver;

    PPCGAsNormalSolver(PreconditionerFunctional const& pf_, double relativeAccuracy_=1e-9, size_t maxSteps_=2000, int verbosity_=1, double eps_ = 1e-15)
    : pf(pf_), initialRelativeAccuracy(relativeAccuracy_), relativeAccuracy(initialRelativeAccuracy), eps(eps_), maxSteps(maxSteps_), verbosity(verbosity_)
    {}

    virtual ~PPCGAsNormalSolver(){}
    virtual void pre(Domain &x, Range &b) {}
    virtual void post (Domain &x) {}

    virtual void apply (Domain &v, const Range &d)
    {
      normalStepPreconditioner->apply(v,d);
      //tangentialStepPreconditioner->apply(v,d);
    }

    virtual void setRelativeAccuracy(double relativeAccuracy_)
    {
      relativeAccuracy = std::min(relativeAccuracy_,initialRelativeAccuracy);
    }

    virtual void setEps(double eps_) { eps = eps_; }

    void setMultiGridSteps(size_t mgSteps_) { mgSteps = mgSteps_; }
    void setMultiGridSmoothingSteps(size_t mgSmoothingSteps_) { mgSmoothingSteps = mgSmoothingSteps_; }
    void setChebyshevSteps(size_t chebySteps_) { chebySteps = chebySteps_; }

  private:
    void setNormalStepParameters() const
    {
        normalStepPreconditioner->setParameter(mgSteps, mgSmoothingSteps, chebySteps, 1e-9);
    }

    void setTangentialStepParameters() const
    {
      normalStepPreconditioner->setParameter(mgSteps, mgSmoothingSteps, chebySteps, 1e-9);
//        normalStepPreconditioner->setParameter(mgStepsT, mgSmoothingStepsT, chebyStepsT, 0);
    }

    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& lin, AbstractFunctionSpaceElement* correctionResidual, AbstractFunctionSpaceElement* adjointResidual)
    {
      A.reset(new AssembledGalerkinOperator<Assembler>(dynamic_cast<BridgeLinearization&>(lin).getValidAssembler()));

      PrecondAssembler assembler(A->getAssembler().spaces());
      assembler.assemble(linearization(pf,Bridge::getImpl<VariableSet>(lin.getOrigin())),PrecondAssembler::MATRIX);
      //normalStepPreconditioner.reset(new NormalStepPreconditioner<Functional,Assembler>(dynamic_cast<BridgeLinearization&>(lin).getValidAssembler(), false, sqrt(relativeAccuracy)));
      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> derivativerhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
      *derivativerhs *= 0.0;
      *constraintrhs *= 0.0;
      lin.evald(*combinedrhs);
      // rhs for computation of adjoint correction
      derivativerhs->axpy(1.0,*combinedrhs,"primal");
      VariableSet& derrhs=Bridge::getImpl<VariableSet>(*derivativerhs);
      // rhs for computation of correction
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);

      // compute adjoint correction
      if(verbosity > 0) std::cout << "computing adjoint correction" << std::endl;
      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
      y = derrhs;

      StrakosTichyEnergyErrorTerminationCriterion<double> terminationCriterion(relativeAccuracy,maxSteps,eps);
      terminationCriterion.lookahead(10);
      normalStepPreconditioner.reset(new NormalStepPreconditioner3<PreconditionerFunctional,PrecondAssembler,components>(assembler,conrhs.descriptions,derrhs.descriptions));
      setNormalStepParameters();
      //tangentialStepPreconditioner.reset(new TangentSpacePreconditioner2<PreconditionerFunctional,PrecondAssembler,components>(assembler,conrhs.descriptions,derrhs.descriptions, mgSteps, mgSmoothingSteps, chebySteps));
      Solver solver(*A, *normalStepPreconditioner, dualPairing, terminationCriterion, verbosity, eps);
      solver.apply(x,y);
      VariableSet& adj=Bridge::getImpl<VariableSet>(adjointCorrection);
      adj = x;
      adj *= -1.0;

      if( adjointResidual != nullptr )
      {
        y = derrhs;
        A->applyscaleadd(-1.0,x,y);
        Bridge::getImpl<VariableSet>(*adjointResidual) = y;
      }

      x = 0;
      y = conrhs;

      // compute initial iterate
      Domain x1(x); x1 = 0;
      typename Traits::VectorP rhs(boost::fusion::at_c<Traits::adjointId>(y.data));
      typename Traits::VectorY sol(Traits::CreateVectorY::init(derrhs.descriptions));
      normalStepPreconditioner->applyStatePreconditioner(sol,rhs);
      boost::fusion::at_c<Traits::stateId>(x1.data) = boost::fusion::at_c<0>(sol.data);
      A->applyscaleadd(-1.0,x1,y);
      // compute correction
      if(verbosity > 0) std::cout << "computing correction" << std::endl;
      solver.apply(x,y);

      x += x1;

      // compute residual
      if(correctionResidual != nullptr)
      {
        y = conrhs;
        A->applyscaleadd(-1.0,x,y);
        Bridge::getImpl<VariableSet>(*correctionResidual) = y;
      }

      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;

      setTangentialStepParameters();
    }


    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction,
        AbstractLinearization const& lin, AbstractFunctionSpaceElement* correctionResidual) const
    {
      setNormalStepParameters();
      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
      *constraintrhs *= 0.0;
      lin.evald(*combinedrhs);
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      using namespace boost::fusion;
      if(correctionResidual != nullptr) constraintrhs->axpy(-1.0,*correctionResidual,"dual");
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);

      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
      
      y=conrhs;

      // compute  initial iterate
      Domain x1(x); x1 = 0;
      typename Traits::VectorP rhs(boost::fusion::at_c<Traits::adjointId>(y.data));
      typename Traits::VectorY sol(Traits::CreateVectorY::init(conrhs.descriptions));
      normalStepPreconditioner->applyStatePreconditioner(sol,rhs);
      boost::fusion::at_c<Traits::stateId>(x1.data) = boost::fusion::at_c<0>(sol.data);
      A->applyscaleadd(-1.0,x1,y);

      // compute correction
      StrakosTichyEnergyErrorTerminationCriterion<double> terminationCriterion(relativeAccuracy,maxSteps,eps);
      terminationCriterion.lookahead(15);
      Solver solver(*A, *normalStepPreconditioner, dualPairing, terminationCriterion, verbosity, eps);
      solver.apply(x,y);

      x += x1;

      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;
    }

    PreconditionerFunctional const& pf;
    mutable std::unique_ptr<NormalStepPreconditioner3<PreconditionerFunctional,PrecondAssembler,components> > normalStepPreconditioner;
    std::unique_ptr<AssembledGalerkinOperator<Assembler> > A;
    DefaultDualPairing<Domain,Range> dualPairing;
    Scalar initialRelativeAccuracy, relativeAccuracy, eps;
    size_t maxSteps;
    int verbosity;
    size_t mgSteps = 500, mgSmoothingSteps = 25, chebySteps = 50;
    size_t mgStepsT = 20, mgSmoothingStepsT = 20, chebyStepsT = 50;
  };

  template <class Assembler_, class PrecondAssembler, class Domain_, class Range_, class VariableSet, int components=Assembler_::Grid::dimension>
  class DirectNormalSolver : public AbstractNormalDirection, public Dune::Preconditioner<Domain_, Range_>
  {
  public:
    typedef Domain_ Domain;
    typedef Range_ Range;
    typedef Assembler_ Assembler;
    typedef typename Assembler::Functional::Functional Functional;
    typedef typename PrecondAssembler::Functional::Functional PreconditionerFunctional;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Traits::Scalar Scalar;
    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> BridgeLinearization;

    explicit DirectNormalSolver(PreconditionerFunctional const& pf_, DirectType directType_=DirectType::UMFPACK3264, MatrixProperties properties_=MatrixProperties::GENERAL, int verbosity_=1)
    : directType(directType_), properties(properties_), verbosity(verbosity_)
    {}

    virtual ~DirectNormalSolver(){}
    virtual void pre(Domain &x, Range &b) {}
    virtual void post (Domain &x) {}

    virtual void apply (Domain &v, const Range &d) { directSolver.apply(d,v); }

    virtual void setRelativeAccuracy(double){}

  private:
    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& lin,
                                                        AbstractFunctionSpaceElement* normalStepResidual, AbstractFunctionSpaceElement* adjointResidual)
    {
      A.reset(new AssembledGalerkinOperator<Assembler>(dynamic_cast<BridgeLinearization&>(lin).getValidAssembler()));

      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.initZeroVector());
      std::unique_ptr<AbstractFunctionSpaceElement> derivativerhs(correction.initZeroVector());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.initZeroVector());
      lin.evald(*combinedrhs);
      // rhs for computation of adjoint correction
      derivativerhs->axpy(1.0,*combinedrhs,"primal");
      VariableSet& derrhs=Bridge::getImpl<VariableSet>(*derivativerhs);
      // rhs for computation of correction
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);

      // compute adjoint correction
      if(verbosity > 0) std::cout << "PrecondType::DIRECT NORMAL SOLVER: Computing adjoint correction." << std::endl;
      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
      y = derrhs;

      directSolver = directInverseOperator(*A,DirectType::UMFPACK3264,MatrixProperties::GENERAL);
      directSolver.apply(y,x);
      VariableSet& adj=Bridge::getImpl<VariableSet>(adjointCorrection);
      adj = x;
      adj *= -1.0;

      if( adjointResidual != nullptr )
      {
        y = derrhs;
        A->applyscaleadd(-1.0,x,y);
        Bridge::getImpl<VariableSet>(*adjointResidual) = y;
      }

      x = 0;
      y = conrhs;
      // compute correction
      if(verbosity > 0) std::cout << "PrecondType::DIRECT NORMAL SOLVER: Computing correction." << std::endl;
      directSolver.apply(y,x);

      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;

      if( normalStepResidual != nullptr )
      {
        y = conrhs;
        A->applyscaleadd(-1.0,x,y);
        Bridge::getImpl<VariableSet>(*normalStepResidual) = y;
      }
    }


    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction,
        AbstractLinearization const& lin, AbstractFunctionSpaceElement* residual) const
    {
      if( verbosity > 0 ) std::cout << "PrecondType::DIRECT NORMAL SOLVER: Computing simplified correction." << std::endl;
      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.initZeroVector());
      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.initZeroVector());
      lin.evald(*combinedrhs);
      constraintrhs->axpy(1.0,*combinedrhs,"dual");
      using namespace boost::fusion;
      if(residual!=nullptr) constraintrhs->axpy(1.0,*residual,"dual");
      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);

      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));

      y=conrhs;
      directSolver.apply(y,x);
      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
      cor = x;
      cor *= -1.0;
    }

//    std::unique_ptr<TangentSpacePreconditioner<PreconditionerFunctional,PrecondAssembler,components> > tangentialStepPreconditioner;
    std::unique_ptr<AssembledGalerkinOperator<Assembler> > A;
    DirectType directType;
    MatrixProperties properties;
    int verbosity;
    InverseLinearOperator<DirectSolver<typename AssembledGalerkinOperator<Assembler>::Domain,typename AssembledGalerkinOperator<Assembler>::Range> > directSolver;
  };



  template <class Functional, class Assembler>
  class NormalStepPreconditioner : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Traits::Scalar Scalar;
    typedef typename Assembler::Grid Grid;
 public:
    NormalStepPreconditioner(Assembler const& assembler, bool onlyLowerTriangle=false, Scalar tolerance=1e-6)
    : Mu(assembler,onlyLowerTriangle),
      B(assembler,false), A(assembler,false),
      PA(new DirectPreconditioner<typename Traits::StateOperator>(A)),
      cheb(Mu,30),
      tol(tolerance)//,
    {
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    ~NormalStepPreconditioner(){}

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&){}

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      using namespace boost::fusion;
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));
      PA->apply(dp,rhsY);

      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      typename Traits::VectorU du(rhsU);
      B.applyscaleaddTransposed(-1.0,dp,rhsU);
      //Lu.apply(du,rhsU);
      //PMu->apply(du,rhsU);
      cheb.apply(du,rhsU);

      //cg.apply(du,rhsU);

      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      typename Traits::VectorY dy(at_c<Traits::yIdx>(y.data));
      B.applyscaleadd(-1.0,du,rhsP);
      PA->apply(dy,rhsP);

      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
//      directInverseOperator(A,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(y,x);
      PA->apply(x,y);
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
      PA->apply(x,y);
    }

  private:
    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperator B;
    typename Traits::StateOperator A;

    std::unique_ptr<typename Traits::StatePreconditioner> PA;
//    std::unique_ptr<typename Traits::AdjointPreconditioner> PAt;
    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    Scalar tol;
 };

  template <class Functional, class Assembler, int components>
  class TangentSpacePreconditioner : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef typename Assembler::Scalar Scalar;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Assembler::Grid Grid;
    static constexpr int dim = Functional::dim;

  public:
    TangentSpacePreconditioner(Assembler const& assembler, bool onlyLowerTriangle=false)
    : grid(boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid()),
      At(assembler,onlyLowerTriangle),
      Mu(assembler,onlyLowerTriangle), Bt(assembler,onlyLowerTriangle),
      B(assembler,onlyLowerTriangle), A(assembler,onlyLowerTriangle),
      PA(new DirectPreconditioner<typename Traits::StateOperator>(A)),
      PAt(new DirectPreconditioner<typename Traits::AdjointOperator>(At)),
      cheb(Mu,50)
    {
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    virtual ~TangentSpacePreconditioner(){}

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&){}

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      std::cout << "in preconditioner" << std::endl;
      std::cout << "adjoint equation" << std::endl;
      using namespace boost::fusion;
      std::cout << "read rhs" << std::endl;
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      std::cout << "init dp" << std::endl;
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));
      std::cout << "solve" << std::endl;
      PAt->apply(dp,rhsY);

      std::cout << "control equation" << std::endl;
      std::cout << "read rhs" << std::endl;
      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      std::cout << "init du" << std::endl;
      typename Traits::VectorU du(rhsU);
      std::cout << "apply Bt" << std::endl;
      Bt.applyscaleadd(-1.0,dp,rhsU);
      std::cout << "invert Mu" << std::endl;
//      directInverseOperator(Mu,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(rhsU,du);
      cheb.apply(du,rhsU);

      std::cout << "state equation" << std::endl;
      typename Traits::Range tmpY(y);
      std::cout << "read rhs " << std::endl;
      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      std::cout << "init dy" << std::endl;
      typename Traits::VectorY dy(at_c<Traits::yIdx>(y.data));
      std::cout << "apply B" << std::endl;
      B.applyscaleadd(-1.0,du,rhsP);
      std::cout << "solve" << std::endl;
      PA->apply(dy,rhsP);

      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
      PA->apply(x,y);
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
      PAt->apply(x,y);
    }

  private:
    Grid const& grid;

    typename Traits::AdjointOperator At;
    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperatorT Bt;
    typename Traits::ControlOperator B;
    typename Traits::StateOperator A;

    std::unique_ptr<typename Traits::StatePreconditioner> PA;
    std::unique_ptr<typename Traits::AdjointPreconditioner> PAt;
    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    //std::unique_ptr<typename Traits::NormUPreconditioner> PMu;
  };

  template <class Functional, class Assembler, int components>
  class TangentSpacePreconditioner2 : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef typename Assembler::Scalar Scalar;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Assembler::Grid Grid;
    typedef MultiGridPreconditioner<Grid,components> MGPreconditioner;
    typedef typename Functional::AnsatzVars AnsatzVariableSetDesc;
    typedef typename Functional::TestVars TestVariableSetDesc;
    static constexpr int dim = Functional::dim;

  public:
    TangentSpacePreconditioner2(Assembler const& assembler, AnsatzVariableSetDesc const& ansatzVars, TestVariableSetDesc const& testVars,
                                size_t mgSteps = 10, size_t mgSmoothingSteps = 10, size_t chebySteps = 10, bool onlyLowerTriangle=false)
    : At(assembler,onlyLowerTriangle),
      Mu(assembler,onlyLowerTriangle), Bt(assembler,onlyLowerTriangle),
      B(assembler,onlyLowerTriangle), A(assembler,onlyLowerTriangle),
      cheb(Mu,chebySteps),
      mgA( A, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      mgAt( At, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      avar(ansatzVars), tvar(testVars)
    {
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    virtual ~TangentSpacePreconditioner2(){}

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&){}

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      using namespace boost::fusion;
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));

      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(rhsY.data);
      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
      at_c<0>(dp.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();

      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      typename Traits::VectorU du(rhsU);
      Bt.applyscaleadd(-1.0,dp,rhsU);
      cheb.apply(du,rhsU);

      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      typename Traits::VectorY dy(at_c<Traits::yIdx>(y.data));
      B.applyscaleadd(-1.0,du,rhsP);

      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(rhsP.data);
      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
      at_c<0>(dy.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();

      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
      using namespace boost::fusion;
      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
      using namespace boost::fusion;
      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();
    }

  private:
    typename Traits::AdjointOperator At;
    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperatorT Bt;
    typename Traits::ControlOperator B;
    typename Traits::StateOperator A;

    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    MGPreconditioner mgA, mgAt;
    AnsatzVariableSetDesc avar;
    TestVariableSetDesc tvar;
  };

//  template <class Functional, class Assembler, int components>
//  class NormalStepPreconditioner2 : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
//  {
//    typedef typename Assembler::Scalar Scalar;
//    typedef OptimalControlTraits<Functional,Assembler> Traits;
//    typedef typename Assembler::Grid Grid;
//    typedef MultiGridSolver<Grid,components> MGSolver;
//    typedef typename Functional::AnsatzVars AnsatzVariableSetDesc;
//    typedef typename Functional::TestVars TestVariableSetDesc;
//    static constexpr int dim = Functional::dim;

//  public:
//    NormalStepPreconditioner2(Assembler const& assembler, AnsatzVariableSetDesc const& ansatzVars, TestVariableSetDesc const& testVars,
//                              size_t mgSteps = 500, size_t mgSmoothingSteps = 10, size_t chebySteps = 20, double relativeAccuracy = 1e-6, bool onlyLowerTriangle=false)
//    : At(assembler,onlyLowerTriangle),
//      Mu(assembler,onlyLowerTriangle), Bt(assembler,onlyLowerTriangle),
//      B(assembler,onlyLowerTriangle), A(assembler,onlyLowerTriangle),
//      cheb(Mu,chebySteps),
//      mgA( A, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy) ),
//      mgAt( At, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy) ),
//      avar(ansatzVars), tvar(testVars)
//    {
//      cheb.initForMassMatrix_TetrahedralQ1Elements();
//    }

//    virtual ~NormalStepPreconditioner2(){}

//    void setParameter(size_t mgSteps, size_t mgSmoothingSteps, size_t chebySteps, double relativeAccuracy)
//    {
//        mgA.setParameter(typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy));
//        mgAt.setParameter(typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy));
//        cheb.setSteps(chebySteps);
//    }

//    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
//    virtual void post(typename Traits::Domain&){}

//    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
//    {
//      using namespace boost::fusion;
//      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
//      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));

//      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
//      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
//      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(rhsY.data);
//      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
//      at_c<0>(dp.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();

//      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
//      typename Traits::VectorU du(rhsU);
//      Bt.applyscaleadd(-1.0,dp,rhsU);
//      cheb.apply(du,rhsU);

//      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
//      typename Traits::VectorY dy(at_c<Traits::yIdx>(y.data));
//      B.applyscaleadd(-1.0,du,rhsP);

//      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(rhsP.data);
//      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
//      at_c<0>(dy.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();

//      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
//      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
//      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
//    }

//    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
//    {
//      using namespace boost::fusion;
//      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
//      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
//      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
//      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
//      at_c<0>(x.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();
//    }

//    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
//    {
//      using namespace boost::fusion;
//      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
//      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
//      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
//      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
//      at_c<0>(x.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();
//    }

//  private:
//    typename Traits::AdjointOperator At;
//    typename Traits::NormUOperator Mu;
//    typename Traits::ControlOperatorT Bt;
//    typename Traits::ControlOperator B;
//    typename Traits::StateOperator A;

//    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
//    MGSolver mgA, mgAt;
//    AnsatzVariableSetDesc avar;
//    TestVariableSetDesc tvar;
//  };


  template <class Functional, class Assembler, int components>
  class NormalStepPreconditioner3 : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef typename Assembler::Scalar Scalar;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Assembler::Grid Grid;
    typedef MultiGridSolver<Grid,components> MGSolver;
    typedef typename Functional::AnsatzVars AnsatzVariableSetDesc;
    typedef typename Functional::TestVars TestVariableSetDesc;
    static constexpr int dim = Functional::dim;

  public:
    NormalStepPreconditioner3(Assembler const& assembler, AnsatzVariableSetDesc const& ansatzVars, TestVariableSetDesc const& testVars,
                              size_t mgSteps = 500, size_t mgSmoothingSteps = 10, size_t chebySteps = 20, double relativeAccuracy = 1e-6, bool onlyLowerTriangle=false)
    :
      Mu(assembler,onlyLowerTriangle), B(assembler,onlyLowerTriangle), Bt(B), A(assembler,onlyLowerTriangle),
      cheb(Mu,chebySteps),
      mgA( A, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy) ),
      mgAt( A, boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid(), typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy), true ),
      avar(ansatzVars), tvar(testVars)
    {
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    virtual ~NormalStepPreconditioner3(){}

    void setParameter(size_t mgSteps, size_t mgSmoothingSteps, size_t chebySteps, double relativeAccuracy)
    {
        mgA.setParameter(typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy));
        mgAt.setParameter(typename MGSolver::Parameter(mgSteps,mgSmoothingSteps,relativeAccuracy));
        cheb.setSteps(chebySteps);
    }

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&){}

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      using namespace boost::fusion;
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));

      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(rhsY.data);
      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
      at_c<0>(dp.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();

      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      typename Traits::VectorU du(rhsU);
      Bt.applyscaleadd(-1.0,dp,rhsU);
      cheb.apply(du,rhsU);

      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      typename Traits::VectorY dy(at_c<Traits::yIdx>(y.data));
      B.applyscaleadd(-1.0,du,rhsP);

      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(rhsP.data);
      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
      at_c<0>(dy.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();

      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
      using namespace boost::fusion;
      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::pIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
      mgA.apply(at_c<Traits::yIdx>(tmpx.data).coefficients(),at_c<Traits::pIdx>(tmpy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::yIdx>(tmpx.data).coefficients();
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
      using namespace boost::fusion;
      typename Traits::AnsatzVariableSetDescription::VariableSet tmpx(avar);
      typename Traits::TestVariableSetDescription::VariableSet tmpy(tvar);
      at_c<Traits::yIdx>(tmpy.data).coefficients() = at_c<0>(y.data);
      mgAt.apply(at_c<Traits::pIdx>(tmpx.data).coefficients(),at_c<Traits::yIdx>(tmpy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::pIdx>(tmpx.data).coefficients();
    }

  private:
    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperator B;
    TransposedOperator<typename Traits::ControlOperator> Bt;
    typename Traits::StateOperator A;

    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    MGSolver mgA, mgAt;
    AnsatzVariableSetDesc avar;
    TestVariableSetDesc tvar;
  };

  template <class Functional, class Assembler, int components, bool exactConstraint>
  class InexactTangentSpacePreconditioner : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef typename Assembler::Scalar Scalar;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Assembler::Grid Grid;
    static constexpr int dim = Functional::dim;
    typedef MultiGridPreconditioner<Grid,components> MGPreconditioner;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,components,components> > BCRSMat;
    typedef typename Functional::AnsatzVars AnsatzVariableSetDesc;
    typedef typename Functional::TestVars TestVariableSetDesc;
    typedef MatrixRepresentedOperator< MatrixAsTriplet<double>, typename Traits::VectorP, typename Traits::VectorY > AdjointOperator;

  public:
    InexactTangentSpacePreconditioner(Assembler const& assembler, AnsatzVariableSetDesc aDesc, TestVariableSetDesc tdesc, size_t mgSteps = 10, size_t mgSmoothingSteps = 10, size_t chebySteps = 10, bool onlyLowerTriangle=false)
    : grid(boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid()),
      Mu(assembler,onlyLowerTriangle),
      B(assembler,onlyLowerTriangle),
      A(assembler,onlyLowerTriangle),
      At(assembler,onlyLowerTriangle),
      cheb(Mu,chebySteps),
      mgA( A, grid, typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      mgAt( At, grid, typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      mgAExact(A,grid, typename MultiGridSolver<Grid,components>::Parameter(500,25,1e-9)),
      mgAtExact(At, grid, typename MultiGridSolver<Grid,components>::Parameter(500,25,1e-9)),
      ansatzDescription(aDesc), testDescription(tdesc)
    {
      timerA.stop(), timerAt.stop(), timerMu.stop(), timerB.stop(), timerCopy.stop();
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    virtual ~InexactTangentSpacePreconditioner(){}

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&)
    {
      std::cout << "timer A : " << boost::timer::format(timerA.elapsed()) << std::endl;
      std::cout << "timer At: " << boost::timer::format(timerAt.elapsed()) << std::endl;
      std::cout << "timer Mu: " << boost::timer::format(timerMu.elapsed()) << std::endl;
      std::cout << "timer B : " << boost::timer::format(timerB.elapsed()) << std::endl;
      std::cout << "timer cp: " << boost::timer::format(timerCopy.elapsed()) << std::endl;
    }

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      using namespace boost::fusion;
      timerCopy.resume();
      std::cout << "inexact adjoint equation" << std::endl;
      std::cout << "read rhs1 " << std::endl;
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      std::cout << "init dp1" << std::endl;
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));
//      std::cout << "read dp2" << std::endl;
//      typename Functional::AnsatzVars::VariableSet mgAx(ansatzDescription);
//      std::cout << "read rhs2" << std::endl;
//      typename Functional::TestVars::VariableSet mgAy(testDescription);
//      std::cout << "copy" << std::endl;
//      at_c<Traits::yIdx>(mgAy.data).coefficients() = at_c<Traits::yIdx>(y.data);
      timerCopy.stop();
      timerAt.resume();
      std::cout << "solve" << std::endl;
      directInverseOperator(At,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(dp,rhsY);
//      if( exactConstraint ) mgAtExact.apply( at_c<Traits::pIdx>(mgAx.data).coefficients(), at_c<Traits::yIdx>(mgAy.data).coefficients() );
//      else mgAt.apply(at_c<Traits::pIdx>(mgAx.data).coefficients(), at_c<Traits::yIdx>(mgAy.data).coefficients());
      timerAt.stop();
      timerCopy.resume();
//      std::cout << "copy" << std::endl;
//      at_c<0>(dp.data) = at_c<Traits::pIdx>(mgAx.data).coefficients();
      //      mgAt.apply(at_c<Traits::pIdx>(x.data),at_c<0>(rhsY.data));

      std::cout << "control equation" << std::endl;
      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      typename Traits::VectorU du(rhsU);
      timerCopy.stop();
      timerB.resume();
      B.applyscaleaddTransposed(-1.0,dp,rhsU);
      timerB.stop();
      timerMu.resume();
      cheb.apply(du,rhsU);
      timerMu.stop();

      std::cout << "adjoint equation" << std::endl;
      timerCopy.resume();
      typename Traits::Range tmpY(y);
      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      timerCopy.stop();
      timerB.resume();
      std::cout << "apply B" << std::endl;
      B.applyscaleadd(-1.0,du,rhsP);
      timerB.stop();
      timerCopy.resume();
            std::cout << "read dp2" << std::endl;
            typename Functional::AnsatzVars::VariableSet mgAx(ansatzDescription);
            std::cout << "read rhs2" << std::endl;
            typename Functional::TestVars::VariableSet mgAy(testDescription);
      at_c<Traits::pIdx>(tmpY.data) = at_c<0>(rhsP.data);
      at_c<Traits::pIdx>(mgAy.data) = at_c<0>(rhsP.data);
      timerCopy.stop();
      timerA.resume();
      std::cout << "invert A" << std::endl;
      if( exactConstraint ) mgAExact.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      else mgA.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      timerA.stop();
      timerCopy.resume();
      at_c<Traits::yIdx>(x.data) = at_c<Traits::yIdx>(mgAx.data).coefficients();
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
      timerCopy.stop();
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
      std::cout << "apply state preconditioner" << std::endl;
      using namespace boost::fusion;
      typename Functional::AnsatzVars::VariableSet mgAx(ansatzDescription);
      typename Functional::TestVars::VariableSet mgAy(testDescription);
      at_c<Traits::pIdx>(mgAy.data) = at_c<0>(y.data);
      mgA.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::yIdx>(mgAx.data).coefficients();
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
      std::cout << "apply adjoint preconditioner" << std::endl;
      abort();
     // PAt->apply(x,y);
    }

  private:
    Grid const& grid;

    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperator B;
    typename Traits::StateOperator A;
    typename Traits::AdjointOperator At;
    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    MGPreconditioner mgA, mgAt;
    MultiGridSolver<Grid,components> mgAExact, mgAtExact;
    AnsatzVariableSetDesc ansatzDescription;
    TestVariableSetDesc testDescription;
    boost::timer::cpu_timer timerA, timerAt, timerMu, timerB, timerCopy;
  };

  template <class Functional, class Assembler, int components, bool exactConstraint>
  class InexactTangentSpacePreconditionerILU : public Dune::Preconditioner<typename OptimalControlTraits<Functional,Assembler>::Domain, typename OptimalControlTraits<Functional,Assembler>::Range>
  {
    typedef typename Assembler::Scalar Scalar;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
    typedef typename Assembler::Grid Grid;
    static constexpr int dim = Functional::dim;
    typedef MultiGridPreconditioner<Grid,components> MGPreconditioner;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,components,components> > BCRSMat;
    typedef typename Functional::AnsatzVars AnsatzVariableSetDesc;
    typedef typename Functional::TestVars TestVariableSetDesc;
    typedef MatrixRepresentedOperator< MatrixAsTriplet<double>, typename Traits::VectorP, typename Traits::VectorY > AdjointOperator;

  public:
    InexactTangentSpacePreconditionerILU(Assembler const& assembler, AnsatzVariableSetDesc const& aDesc, TestVariableSetDesc const& tdesc, size_t mgSteps = 10, size_t mgSmoothingSteps = 10, size_t chebySteps = 10, bool onlyLowerTriangle=false)
    : grid(boost::fusion::at_c<0>(assembler.spaces())->gridManager().grid()),
      Mu(assembler,onlyLowerTriangle),
      B(assembler,onlyLowerTriangle),
      A(assembler,onlyLowerTriangle),
      At(assembler,onlyLowerTriangle),
      cheb(Mu,chebySteps),
      mgA( A, grid, typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      mgAt( At, grid, typename MGPreconditioner::Parameter(mgSteps,mgSmoothingSteps) ),
      mgAExact(A,grid, typename MultiGridSolver<Grid,components>::Parameter(500,25,1e-9)),
      mgAtExact(At, grid, typename MultiGridSolver<Grid,components>::Parameter(500,25,1e-9)),
      ansatzDescription(aDesc), testDescription(tdesc)
    {
      timerA.stop(), timerAt.stop(), timerMu.stop(), timerB.stop(), timerCopy.stop();
      cheb.initForMassMatrix_TetrahedralQ1Elements();
    }

    virtual ~InexactTangentSpacePreconditionerILU(){}

    virtual void pre(typename Traits::Domain&, typename Traits::Range&){}
    virtual void post(typename Traits::Domain&)
    {
      std::cout << "timer A : " << boost::timer::format(timerA.elapsed()) << std::endl;
      std::cout << "timer At: " << boost::timer::format(timerAt.elapsed()) << std::endl;
      std::cout << "timer Mu: " << boost::timer::format(timerMu.elapsed()) << std::endl;
      std::cout << "timer B : " << boost::timer::format(timerB.elapsed()) << std::endl;
      std::cout << "timer cp: " << boost::timer::format(timerCopy.elapsed()) << std::endl;
    }

    virtual void apply(typename Traits::Domain& x, typename Traits::Range const& y)
    {
      using namespace boost::fusion;
      timerCopy.resume();
      typename Traits::VectorY rhsY(at_c<Traits::yIdx>(y.data));
      typename Traits::VectorP dp(at_c<Traits::pIdx>(y.data));
      typename Functional::AnsatzVars::VariableSet mgAx(ansatzDescription);
      typename Functional::TestVars::VariableSet mgAy(testDescription);
      at_c<Traits::yIdx>(mgAy.data) = at_c<Traits::yIdx>(y.data);
      timerCopy.stop();
      timerAt.resume();
      if( exactConstraint ) mgAtExact.apply( at_c<Traits::pIds>(mgAx.data).coefficients(), at_c<Traits::yIdx>(mgAy.data).coefficients() );
      else mgAt.apply(at_c<Traits::pIdx>(mgAx.data).coefficients(), at_c<Traits::yIdx>(mgAy.data).coefficients());
      timerAt.stop();
      timerCopy.resume();
      at_c<0>(dp.data) = at_c<Traits::pIdx>(mgAx.data).coefficients();
      //      mgAt.apply(at_c<Traits::pIdx>(x.data),at_c<0>(rhsY.data));

      typename Traits::VectorU rhsU(at_c<Traits::uIdx>(y.data));
      typename Traits::VectorU du(rhsU);
      timerCopy.stop();
      timerB.resume();
      B.applyscaleaddTransposed(-1.0,dp,rhsU);
      timerB.stop();
      timerMu.resume();
      cheb.apply(du,rhsU);
      timerMu.stop();

      timerCopy.resume();
      typename Traits::Range tmpY(y);
      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(y.data));
      timerCopy.stop();
      timerB.resume();
      B.applyscaleadd(-1.0,du,rhsP);
      timerB.stop();
      timerCopy.resume();
      at_c<Traits::pIdx>(tmpY.data) = at_c<0>(rhsP.data);
      at_c<Traits::pIdx>(mgAy.data) = at_c<0>(rhsP.data);
      timerCopy.stop();
      timerA.resume();
      if( exactConstraint ) mgAExact.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      else mgA.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      timerA.stop();
      timerCopy.resume();
      at_c<Traits::yIdx>(x.data) = at_c<Traits::yIdx>(mgAx.data).coefficients();
      at_c<Traits::uIdx>(x.data) = at_c<0>(du.data);
      at_c<Traits::pIdx>(x.data) = at_c<0>(dp.data);
      timerCopy.stop();
    }

    void applyStatePreconditioner(typename Traits::VectorY& x, typename Traits::VectorP const& y)
    {
      using namespace boost::fusion;
      typename Functional::AnsatzVars::VariableSet mgAx(ansatzDescription);
      typename Functional::TestVars::VariableSet mgAy(testDescription);
      at_c<Traits::pIdx>(mgAy.data) = at_c<0>(y.data);
      mgA.apply(at_c<Traits::yIdx>(mgAx.data).coefficients(),at_c<Traits::pIdx>(mgAy.data).coefficients());
      at_c<0>(x.data) = at_c<Traits::yIdx>(mgAx.data).coefficients();
    }

    void applyAdjointPreconditioner(typename Traits::VectorP& x, typename Traits::VectorY const& y)
    {
     // PAt->apply(x,y);
    }

  private:
    Grid const& grid;

    typename Traits::NormUOperator Mu;
    typename Traits::ControlOperator B;
    typename Traits::StateOperator A;
    typename Traits::AdjointOperator At;
    ChebyshevPreconditioner<typename Traits::NormUOperator> cheb;
    MGPreconditioner mgA, mgAt;
    MultiGridSolver<Grid,components> mgAExact, mgAtExact;
    AnsatzVariableSetDesc const& ansatzDescription;
    TestVariableSetDesc const& testDescription;
    boost::timer::cpu_timer timerA, timerAt, timerMu, timerB, timerCopy;
  };
//  template <class PreconditionerFactory_, class PDEPreconditionerFactory, class VariableSet>
//  class OptimalControlNormalSolver : public AbstractNormalDirection, public Dune::Preconditioner<typename PreconditionerFactory_::Operator::Domain,typename  PreconditionerFactory_::Operator::Range>
//  {
//  public:

//    typedef typename PreconditionerFactory_::Operator::Assembler Assembler;
//    typedef typename PreconditionerFactory_::Operator::Domain Domain;
//    typedef typename PreconditionerFactory_::Operator::Range Range;
//    //    typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<>::type Domain;
//    //    typedef Domain Range;
//    typedef typename Assembler::Functional::Functional Functional;
//    typedef Bridge::ConnectedKaskadeLinearization<Functional> BridgeLinearization;

//    static constexpr int yIdx = Functional::yIdx;
//    static constexpr int uIdx = Functional::uIdx;
//    static constexpr int pIdx = Functional::lIdx;
//    typedef AssembledGalerkinOperator<Assembler,pIdx,pIdx+1,yIdx,yIdx+1> StateOperator;
//    typedef AssembledGalerkinOperator<Assembler,pIdx,pIdx+1,uIdx,uIdx+1> ControlOperator;
//    typedef AssembledGalerkinOperator<Assembler,yIdx,yIdx+1,pIdx,pIdx+1> AdjointOperator;
//    typedef AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,uIdx,uIdx+1> NormUOperator;
//    typedef AssembledGalerkinOperator<Assembler,yIdx,yIdx+1,yIdx,yIdx+1> NormYOperator;
//    typedef AssembledGalerkinOperator<Assembler,uIdx,uIdx+1,pIdx,pIdx+1> ControlOperatorT;
//    typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<yIdx,yIdx+1> CreateVectorY;
//    typedef typename CreateVectorY::type VectorY;
//    typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<uIdx,uIdx+1> CreateVectorU;
//    typedef typename CreateVectorU::type VectorU;
//    typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<pIdx,pIdx+1> CreateVectorP;
//    typedef typename CreateVectorP::type VectorP;

//    typedef PreconditionerFactory<Assembler,VectorU,VectorU,double> UPrecondFactory;
//    typedef PreconditionerFactory<Assembler,VectorP,VectorY,double> AdjointPrecondFactory;
//    typedef PreconditionerFactory<Assembler,VectorY,VectorP,double> PDEPrecondFactory;

//    virtual ~OptimalControlNormalSolver() {}

//    OptimalControlNormalSolver(PreconditionerFactory_& prec_, PDEPreconditionerFactory& pdePrec_)
//    : A(nullptr), preconditioner(nullptr), tolerance(1e-6), maxSteps(100)
//    {}

//    virtual void pre(Domain &x, Range &b) {}
//    virtual void post (Domain &x) {}

//    virtual void apply(Domain &x, Range const& y)
//    {
//      preconditioner->apply(x,y);
//    }

//  private:
//    void solveAdjointEquation(AbstractFunctionSpaceElement& adjointCorrection, VariableSet const& variableSet) const
//    {
//      using namespace boost::fusion;
//      VectorY adjointRhs(CreateVectorY::init(variableSet.descriptions));
//      at_c<0>(adjointRhs.data) = at_c<yIdx>(variableSet.data).coefficients();
//      VectorP adjointSolution(CreateVectorP::init(variableSet.descriptions));

//      preconditioner->applyAdjointPreconditioner(adjointSolution,adjointRhs);
//      adjointSolution *= -1.0;

//      VariableSet& adj=Bridge::getImpl<VariableSet>(adjointCorrection);
//      at_c<pIdx>(adj.data).coefficients() = at_c<0>(adjointSolution.data);
//    }

//    void solveSystem(AbstractFunctionSpaceElement& correction, VariableSet const& variableSet) const
//    {
//      using namespace boost::fusion;

//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(variableSet.descriptions));
//      x *= 0.0;
//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(variableSet.descriptions));
//      y = variableSet;

//      // compute initial iterate
//      VectorY tmpY(CreateVectorY::init(variableSet.descriptions));
//      VectorP tmpRhs(CreateVectorP::init(variableSet.descriptions));
//      at_c<0>(tmpRhs.data) = at_c<pIdx>(y.data);
//      preconditioner->applyStatePreconditioner(tmpY,tmpRhs);
//      at_c<yIdx>(x.data) = at_c<0>(tmpY.data);

//      // compute dn
//      DefaultDualPairing<Domain,Range> dualpairing;
//      Dune::InverseOperatorResult res;
//      TPCG<Domain,Range> tpcg(*A, *preconditioner, dualpairing, tolerance, maxSteps, 1);
//      tpcg.apply(x,y,res);
//      //      directInverseOperator(*A,DirectType::UMFPACK3264,MatrixProperties::GENERAL).apply(y,x);


//      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
//      cor = x;
//      cor *= -1.0;
//    }

//    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& linearization)
//    {
//      using namespace boost::fusion;

//      Assembler& assembler = dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler();
//      boost::timer::cpu_timer timer1;
//      preconditioner.reset(new NormalStepPreconditioner<Functional,Assembler>(assembler));
//      std::cout << "preconditioner reset: " << boost::timer::format(timer1.elapsed());
//      boost::timer::cpu_timer timer2;
//      A.reset(new AssembledGalerkinOperator<Assembler>(assembler));
//      std::cout << "A reset: " << boost::timer::format(timer2.elapsed());

//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone()), derivativerhs(correction.clone()), constraintrhs(correction.clone());
//      *derivativerhs *= 0.0;
//      *constraintrhs *= 0.0;
//      linearization.evald(*combinedrhs);
//      derivativerhs->axpy(1.0,*combinedrhs,"primal");
//      constraintrhs->axpy(1.0,*combinedrhs,"dual");

//      // compute dy_n
//#ifdef TESTOUTPUT
//      boost::timer::cpu_timer systemTimer;
//#endif
//      solveSystem(correction, Bridge::getImpl<VariableSet>(*constraintrhs));
//#ifdef TESTOUTPUT
//      std::cout << "computation of normal step: " << boost::timer::format(systemTimer.elapsed());
//      boost::timer::cpu_timer adjointTimer;
//#endif
//      // compute dp
//      solveAdjointEquation(adjointCorrection, Bridge::getImpl<VariableSet>(*derivativerhs));
//#ifdef TESTOUTPUT
//      std::cout << "computation of adjoint variable: " << boost::timer::format(adjointTimer.elapsed());
//#endif
//    }


//    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction,
//        AbstractLinearization const& lin) const
//    {
//#ifdef TESTOUTPUT
//      boost::timer::cpu_timer timer;
//#endif
//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone()), constraintrhs(correction.clone());
//      *constraintrhs *= 0.0;
//      lin.evald(*combinedrhs);
//      constraintrhs->axpy(1.0,*combinedrhs,"dual");

//      solveSystem(correction,Bridge::getImpl<VariableSet>(*constraintrhs));
//#ifdef TESTOUTPUT
//      std::cout << "computation of simplified normal step: " << boost::timer::format(timer.elapsed());
//#endif
//    }

//    std::unique_ptr<AssembledGalerkinOperator<Assembler> > A;
//    std::unique_ptr<NormalStepPreconditioner<Functional,Assembler> > preconditioner;
//    double tolerance;
//    size_t maxSteps;
//  };


//  template <class Operator, class VariableSet>
//  class CGasNormalSolver : public AbstractNormalDirection, public Dune::Preconditioner<typename Operator::Domain,typename Operator::Range>
//  {
//  public:

//    typedef typename Operator::Assembler Assembler;
//    typedef typename Assembler::Functional::Functional Functional;
//    typedef typename Operator::Domain Domain;
//    typedef typename Operator::Range Range;
//    typedef typename Operator::Scalar Scalar;
//    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> BridgeLinearization;

//    virtual ~CGasNormalSolver() {}

//    CGasNormalSolver(Scalar accuracy_, int maxSteps_, int verbose_=0)
//    : P(nullptr), P2(nullptr), accuracy(accuracy_), maxSteps(maxSteps_), verbose(verbose_)
//    {}

//    //    CGasNormalSolver(Dune::Preconditioner<Domain,Range>& preconditioner, Scalar accuracy_, int maxSteps_, int verbose_=0)
//    //    : P(&preconditioner), accuracy(accuracy_), maxSteps(maxSteps_), verbose(verbose_)
//    //    {}
//    //
//    //    void setPreconditioner(Dune::Preconditioner<Domain,Range>& preconditioner)
//    //    {
//    //      P = &preconditioner;
//    //    }

//    virtual void pre(Domain &x, Range &b) {}
//    virtual void apply (Domain &v, const Range &d)
//    {
//      P2->apply(v,d);
//    }
//    virtual void post (Domain &x) {}

//  private:

//    virtual void  computeCorrectionAndAdjointCorrection(AbstractFunctionSpaceElement& correction, AbstractFunctionSpaceElement& adjointCorrection, AbstractLinearization& linearization)
//    {
//      /// ToDo: ???????????
//      A.reset(new AssembledGalerkinOperator<Assembler>(dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler()));

//      //P.reset(prec.create(*A).release());
//      P.reset(new NormalStepPreconditioner<Functional,Assembler>(dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler()));
//      P2.reset(new InexactTangentSpacePreconditioner<Functional,Assembler>(dynamic_cast<BridgeLinearization&>(linearization).getValidAssembler()));
//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
//      std::unique_ptr<AbstractFunctionSpaceElement> derivativerhs(correction.clone());
//      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
//      *derivativerhs *= 0.0;
//      *constraintrhs *= 0.0;
//      linearization.evald(*combinedrhs);
//      derivativerhs->axpy(1.0,*combinedrhs,"primal");
//      constraintrhs->axpy(1.0,*combinedrhs,"dual");
//      VariableSet& derrhs=Bridge::getImpl<VariableSet>(*derivativerhs);
//      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);
//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
//      x *= 0.0;

//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
//      y = derrhs;

//      DefaultDualPairing<Domain,Range> dualpairing;
//      TPCG<Domain,Range> tpcg(*A, *P, dualpairing, accuracy, maxSteps, verbose);
//      Dune::InverseOperatorResult res;

//      tpcg.apply(x,y,res);

//      VariableSet& adj=Bridge::getImpl<VariableSet>(adjointCorrection);
//      adj = x;
//      adj *= -1.0;

//      y = conrhs;

//      tpcg.apply(x,y,res);


//      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
//      cor = x;
//      cor *= -1.0;
//    }


//    virtual void computeSimplifiedCorrection(AbstractFunctionSpaceElement& correction,
//        AbstractLinearization const& lin) const
//    {
//      /// ToDo:
//      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(correction.clone());
//      std::unique_ptr<AbstractFunctionSpaceElement> constraintrhs(correction.clone());
//      *constraintrhs *= 0.0;
//      lin.evald(*combinedrhs);
//      constraintrhs->axpy(1.0,*combinedrhs,"dual");
//      VariableSet& conrhs=Bridge::getImpl<VariableSet>(*constraintrhs);

//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));
//      x *= 0.0;
//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(conrhs.descriptions));

//      y=conrhs;

//      DefaultDualPairing<Domain,Range> dualpairing;
//      TPCG<Domain,Range> tpcg(*A, *P, dualpairing, accuracy, maxSteps, verbose);
//      Dune::InverseOperatorResult res;

//      tpcg.apply(x,y,res);

//      VariableSet& cor=Bridge::getImpl<VariableSet>(correction);
//      cor = x;
//      cor *= -1.0;
//    }

//    Scalar accuracy;
//    int maxSteps, verbose;
//    std::unique_ptr<Dune::Preconditioner<Domain,Range> > P;
//    std::unique_ptr<InexactTangentSpacePreconditioner<Functional,Assembler> > P2;
//    //    std::unique_ptr<Dune::Preconditioner<Domain,Range> > P;
//    std::unique_ptr<AssembledGalerkinOperator<Assembler> > A;
//  };

  template <class> class ProjectionOnTangentSpace;
  template <class,int,int,int,int> class MGProjectionOnTangentSpace;

  template <class FSE, class CoeffVector, int rbegin, int rend>
  struct CopyFromFunctionSpaceElement
  {
    static void apply(FSE const& fse, CoeffVector& v)
    {
      boost::fusion::at_c<rbegin>(v.data) = boost::fusion::at_c<rbegin>(fse.data).coefficients();
      CopyFromFunctionSpaceElement<FSE,CoeffVector,rbegin+1,rend>::apply(fse,v);
    }
  };

  template <class FSE, class CoeffVector, int rbegin>
  struct CopyFromFunctionSpaceElement<FSE,CoeffVector,rbegin,rbegin>
  {
    static void apply(FSE const&, CoeffVector&) {}
  };

  template <class FSE, class CoeffVector, int rbegin, int rend>
  struct AddToFunctionSpaceElement
  {
    static void apply(CoeffVector const& v, FSE& fse)
    {
      boost::fusion::at_c<rbegin>(fse.data).coefficients() += boost::fusion::at_c<rbegin>(v.data);
      AddToFunctionSpaceElement<FSE,CoeffVector,rbegin+1,rend>::apply(v,fse);
    }
  };

  template <class FSE, class CoeffVector, int rbegin>
  struct AddToFunctionSpaceElement<FSE,CoeffVector,rbegin,rbegin>
  {
    static void apply(CoeffVector const&, FSE&) {}
  };

  template<int rbegin, int rend, class FSE, class CoeffVector>
  void copyFromFunctionSpaceElement(FSE const& fse, CoeffVector& coeffVector)
  {
    CopyFromFunctionSpaceElement<FSE,CoeffVector,rbegin,rend>::apply(fse,coeffVector);
  }

  template<int rbegin, int rend, class FSE, class CoeffVector>
  void addToFunctionSpaceElement(CoeffVector const& coeffVector,FSE& fse)
  {
    AddToFunctionSpaceElement<FSE,CoeffVector,rbegin,rend>::apply(coeffVector, fse);
  }

  template <class Assembler, class Lin, int rbegin, int rend, int cbegin, int cend>
  void ddxpy_template(Lin& lin, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x)
  {
    static_assert(cbegin<=cend,"incompatible template arguments");
    static_assert(rbegin<=rend,"incompatible template arguments");
    typedef typename Assembler::Functional::Functional Functional;
    boost::timer::cpu_timer t1;
    AssembledGalerkinOperator<Assembler,rbegin,rend,cbegin,cend> gop(dynamic_cast<Bridge::ConnectedKaskadeLinearization<Functional>&>(lin.getTangentialLinearization()).getValidAssembler());
    std::cout << "gop init: " << boost::timer::format(t1.elapsed()) << std::endl;
    boost::timer::cpu_timer t2;
    typename Functional::AnsatzVars::template CoefficientVectorRepresentation<rbegin,rend>::type tmpy
        (Functional::AnsatzVars::template CoefficientVectorRepresentation<rbegin,rend>::init(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet> const& >(x).get().descriptions));
    typename Functional::TestVars::template CoefficientVectorRepresentation<cbegin,cend>::type tmpx
        (Functional::TestVars::template CoefficientVectorRepresentation<cbegin,cend>::init(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet> const& >(x).get().descriptions));
    std::cout << "init vecs: " << boost::timer::format(t2.elapsed()) << std::endl;
    boost::timer::cpu_timer t3;
    auto const& tx = dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet> const& >(x).get();
    copyFromFunctionSpaceElement<cbegin,cend>(tx,tmpx);
    std::cout << "copy1: " << boost::timer::format(t3.elapsed()) << std::endl;
    //    boost::fusion::at_c<0>(tmpx.data) = boost::fusion::at_c<0>(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet> const& >(x).get().data).coefficients();
//    boost::fusion::at_c<1>(tmpx.data) = boost::fusion::at_c<1>(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet> const& >(x).get().data).coefficients();
    boost::timer::cpu_timer t4;
    gop.apply(tmpx,tmpy);
    std::cout << "gop apply: " << boost::timer::format(t4.elapsed()) << std::endl;

    boost::timer::cpu_timer t5;
    auto& ty = dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet>& >(y).get();
    addToFunctionSpaceElement<rbegin,rend>(tmpy,ty);
    std::cout << "copy2: " << boost::timer::format(t5.elapsed()) << std::endl;
    //    boost::fusion::at_c<0>(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet>& >(y).get().data).coefficients() += boost::fusion::at_c<0>(tmpy.data);
//    boost::fusion::at_c<1>(dynamic_cast<Bridge::Vector<typename Functional::AnsatzVars::VariableSet>& >(y).get().data).coefficients() += boost::fusion::at_c<1>(tmpy.data);
  }

  template<class Assembler, class Preconditioner, class VariableSet, int components = 1, CGImplementationType cgImpl = CGImplementationType::HYBRID>
  class ProjectedAPCGSolver : public AbstractTangentialSpace
  {
  public:
    typedef typename Preconditioner::Domain Domain;
    typedef typename Preconditioner::Range Range;
    typedef typename Preconditioner::Assembler NormalStepAssembler;
    typedef OptimalControlTraits<typename Assembler::Functional::Functional,Assembler> Traits;
    typedef double Scalar;
    typedef MGProjectionOnTangentSpace<NormalStepAssembler,components,Traits::stateId,Traits::controlId,Traits::adjointId> Projection;
    typedef typename Assembler::Functional::Functional Functional;
    static constexpr int yIdx = Traits::yIdx;
    static constexpr int uIdx = Traits::uIdx;
    static constexpr int pIdx = Traits::pIdx;
    typedef CGBase<Domain,Range,cgImpl> Solver;
    typedef Bridge::ConnectedKaskadeLinearization<typename NormalStepAssembler::Functional::Functional> NormalBridgeLinearization;
    typedef Bridge::ConnectedKaskadeLinearization<typename Assembler::Functional::Functional> TangentialBridgeLinearization;
    typedef typename Assembler::Functional::Functional TF;

    ProjectedAPCGSolver(Preconditioner& P_, Scalar accuracy_, int verbose_ = 0, size_t maxSteps_ = 500, double eps_ = 1e-12)
    : P(P_), accuracy(accuracy_), eps(eps_), verbose(verbose_), maxSteps(maxSteps_), terminationCriterion(new StrakosTichyEnergyErrorTerminationCriterion<double>(accuracy,maxSteps,eps))
    {}

    virtual ~ProjectedAPCGSolver() {}

    virtual void setRelativeAccuracy(double accuracy_) final
    {
      if(verbose > 0) std::cout << "TANGENTIAL SPACE: updating relative accuracy from " << accuracy << " to " << accuracy_ << std::endl;
      accuracy=accuracy_;
    }

    void usePreconditionerForNorm() { terminationCriterion.reset(new StrakosTichyPTerminationCriterion<double>(accuracy,maxSteps,eps)); }

    void setLookAHead(size_t d) { terminationCriterion->lookahead(d); }

    virtual int nSolutionVectors() const final { return 2; }
    virtual bool localConvergenceLikely() final
    {
      return !encounteredNonConvexity;
    }

    virtual AbstractFunctionSpaceElement& getCorrectRhs() final
        {
      return *localtmpvec;
        }

    virtual void setEps(double eps_) { eps = eps_; }

    virtual void setLipschitzConstant(double omega) { omegaL = omega; }

  private:
    virtual int computeBasis(std::vector<std::shared_ptr<AbstractFunctionSpaceElement> >& correction,
        LagrangeLinearization& lin,
        AbstractFunctionSpaceElement const& normalStep,
        double nu0, AbstractFunctionSpaceElement* residual)
        {
          using namespace boost::fusion;

      auto& normalAssembler = dynamic_cast<NormalBridgeLinearization&>(lin.getNormalLinearization()).getValidAssembler();
      auto& tangentialAssembler = dynamic_cast<TangentialBridgeLinearization&>(lin.getTangentialLinearization()).getValidAssembler();

      std::cout << "assemblers: normal: " << normalAssembler.nrows(0,1) << ", " << normalAssembler.nrows(1,2) << ", " << normalAssembler.nrows(2,3) << std::endl;
      std::cout << "assemblers: tangential: " << tangentialAssembler.nrows(0,1) << ", " << tangentialAssembler.nrows(1,2) << ", " << tangentialAssembler.nrows(2,3) << std::endl;

      LagrangeOperator<NormalStepAssembler,Assembler,yIdx,uIdx,pIdx> M(normalAssembler,tangentialAssembler);


      Range rhs( dynamic_cast<NormalBridgeLinearization&>(lin.getNormalLinearization()).getValidAssembler().rhs() );
      at_c<pIdx>(rhs.data) *= 0.;
      VariableSet xv = Bridge::getImpl<VariableSet>(normalStep);
      at_c<pIdx>(xv.data) *= 0.;
      Domain x(xv);
//      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(xv.descriptions));
      M.applyscaleadd(nu0,x,rhs);
      x *= 0;
      at_c<pIdx>(rhs.data) *= 0.;

/*      std::unique_ptr<AbstractFunctionSpaceElement> combinedrhs(normalStep.clone());
      std::unique_ptr<AbstractFunctionSpaceElement> derivativerhs(normalStep.clone());
      *derivativerhs *= 0.0;
      std::cout << "begin computation of basis" << std::endl;
      boost::timer::cpu_timer evalTimer;
      // derivativerhs = L_x(x_0,p_0)+nu_0 L_xx(x_0,p_0)dn
      lin.evald(*combinedrhs);
      std::cout << "evalTimer: " << boost::timer::format(evalTimer.elapsed()) << std::endl;
      boost::timer::cpu_timer evalTimer2;
      derivativerhs->axpy(1.0,*combinedrhs,"primal");
      std::cout << "evalTimer2: " << boost::timer::format(evalTimer2.elapsed()) << std::endl;
      *combinedrhs *= 0.0;
      boost::timer::cpu_timer evalTimer3;
      ddxpy_template<Assembler,LagrangeLinearization,0,2,0,2>(lin,*combinedrhs,normalStep);
//      lin.ddxpy(*combinedrhs,normalStep,0,2,0,2);
      std::cout << "evalTimer3: " << boost::timer::format(evalTimer3.elapsed()) << std::endl;
      boost::timer::cpu_timer evalTimer4;
      derivativerhs->axpy(nu0,*combinedrhs,"primal");
      std::cout << "evalTimer4: " << boost::timer::format(evalTimer4.elapsed()) << std::endl;
      VariableSet& derrhs=Bridge::getImpl<VariableSet>(*derivativerhs);
      Domain x(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));*/

      DefaultDualPairing<Domain,Range> dualpairing;

/*      MatrixAsTriplet<double> mat;
      lin.getNormalLinearization().getMatrixBlocks(mat);
      typedef typename VariableSet::Descriptions::template CoefficientVectorRepresentation<>::type Vec;
      MatrixRepresentedOperator<MatrixAsTriplet<double>, Vec, Vec> PP(mat);
      AlwaysRegularizeWithThetaGuess<MatrixRepresentedOperator<MatrixAsTriplet<double>, Vec, Vec>,Vec> f(PP,getImpl<Vec>(normalStep),nu0,omegaL);*/
      Solver solver(M, P, dualpairing, *terminationCriterion,  verbose, eps);
      Dune::InverseOperatorResult res;
//      Range y(VariableSet::Descriptions::template CoefficientVectorRepresentation<>::init(derrhs.descriptions));
//      y=derrhs;
      solver.apply(x, rhs, res);
      encounteredNonConvexity = solver.encounteredNonConvexity();

/*      if(residual != nullptr)
      {
        y = derrhs;
        M.applyscaleadd(-1.0,x,y);
        derrhs = y;
        Bridge::getImpl<VariableSet>(*residual) = derrhs;
      }*/


      VariableSet& cor=Bridge::getImpl<VariableSet>(*(correction[0]));
      cor = x;
      cor *= -1.0;

      // project on tangent space
//      if(solver.getSolutionEnergyNormEstimate() > eps)
//      {
//        projection.reset(new Projection (dynamic_cast<NormalBridgeLinearization&>(lin.getNormalLinearization()).getValidAssembler()));
//        projection->apply(cor);
//      }

      return 1;
    }

    Preconditioner& P;
    double accuracy, eps;
    int verbose;
    std::unique_ptr<AbstractFunctionSpaceElement> localtmpvec;
    std::unique_ptr<Projection> projection;
    size_t maxSteps;
    bool encounteredNonConvexity = false;
    std::unique_ptr<PCGTerminationCriterion<double> > terminationCriterion;
    double omegaL;
 };

  struct DefaultProjectionPolicy
  {
    template <class Arg>
    void projectOnConstraint(Arg&) const {}
  };

  template <class Preconditioner>
  struct ChooseProjectionPolicy
  {
    typedef DefaultProjectionPolicy type;
  };

  template <class Assembler> class ProjectionOnTangentSpace;

  template <class Functional, class Assembler, int components>
  struct ChooseProjectionPolicy<InexactTangentSpacePreconditioner<Functional,Assembler,components> >
  {
    typedef ProjectionOnTangentSpace<Assembler> type;
  };

  template <class Assembler>
  class ProjectionOnTangentSpace
  {
    typedef typename Assembler::Functional::Functional Functional;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
  public:
    ProjectionOnTangentSpace(Assembler const& assembler)
    : A(assembler,false), B(assembler,false),
      PA(new DirectPreconditioner<typename Traits::StateOperator>(A))
    {}

    template <class VarSetDesc>
    void apply(VarSetDesc& x)
    {
      using namespace boost::fusion;
      typename Traits::VectorP rhsP(at_c<Traits::pIdx>(x.data));
      rhsP *= 0.0;
      typename Traits::VectorU du(at_c<Traits::uIdx>(x.data));
      typename Traits::VectorY dy(at_c<Traits::yIdx>(x.data));
      B.applyscaleadd(-1.0,du,rhsP);
      PA->apply(dy,rhsP);

      at_c<Traits::yIdx>(x.data) = at_c<0>(dy.data);
    }

  private:
    typename Traits::StateOperator A;
    typename Traits::ControlOperator B;

    std::unique_ptr<typename Traits::StatePreconditioner> PA;
  };

  template <class Assembler, int nComponents, int stateId, int controlId, int adjointId>
  class MGProjectionOnTangentSpace
  {
    typedef typename Assembler::Functional::Functional Functional;
    typedef OptimalControlTraits<Functional,Assembler> Traits;
  public:
    MGProjectionOnTangentSpace(Assembler const& assembler, typename MultiGridSolver<typename Assembler::Grid,nComponents>::Parameter parameter = typename MultiGridSolver<typename Assembler::Grid,nComponents>::Parameter(500,25,1e-9))
      : A(assembler,false), B(assembler,false),
        mgSolver(A, assembler.getGridManager().grid(), parameter)
    {}

    template <class VarSet>
    void apply(VarSet& x)
    {
      using namespace boost::fusion;

      VarSet y(x.descriptions);

      typename Traits::VectorP rhsP(Traits::CreateVectorP::init(x.descriptions));
      typename Traits::VectorU du(Traits::CreateVectorU::init(x.descriptions));
      at_c<0>(du.data) = at_c<Traits::controlId>(x.data).coefficients();
      B.applyscaleadd(-1.0,du,rhsP);
      at_c<Traits::adjointId>(y.data) = at_c<Traits::adjointId>(x.data);
      mgSolver.apply(at_c<stateId>(x.data).coefficients(),at_c<Traits::adjointId>(y.data).coefficients());
    }

  private:
    typename Traits::StateOperator A;
    typename Traits::ControlOperator B;
    MultiGridSolver<typename Assembler::Grid,nComponents> mgSolver;
  };

} // namespace Kaskade
#endif
