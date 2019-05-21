#ifndef IP_AUX_HH
#define IP_AUX_HH

#include <memory>
#include "ipfunctional.hh"
#include "fem/variables.hh"
#include "fem/functional_aux.hh"
#include "dune_bridge.hh"
#include "newton_damped.hh"
#include "homotopy_base.hh"
#include "fem/functional_manipulation.hh"
#include <limits.h>
#include <boost/timer/timer.hpp>
#include "fem/special_quadrature.hh"

/**
 * @file 
 * @brief Much of the implementation for barrier methods
 * @author Anton Schiela
 *
 */
namespace Bridge{
/// Linearization for Constrained Optimization with Barrier methods
/** Combines the Functional to be minimized and appropriate barrier terms
 *  Barrier terms are assembled separately via a different type of quadrature rule.
 */
template<class BarrierFunctional, 
         class VectorImpl, 
         class ImageImpl=VectorImpl,
         class QuadRule =  TrapezoidalRule<typename BarrierFunctional::AnsatzVars::Grid::ctype, BarrierFunctional::AnsatzVars::Grid::dimension> >
class IPLinearization
{
public:
  typedef BarrierFunctional IPFunctional;
  typedef typename BarrierFunctional::AnsatzVars::Grid Grid;
  typedef typename Grid::LeafGridView GridView;
  typedef VectorImpl DomainElement;
  typedef ImageImpl ImageElement;
  typedef typename BarrierFunctional::Scalar Scalar;

  typedef KaskadeLinearization<typename BarrierFunctional::OptimalityFunctional, VectorImpl> Unconstrained;
  typedef KaskadeLinearization<BarrierFunctional, VectorImpl, QuadRule> Barrier;
  typedef Dune::LinearOperator<DomainElement, ImageElement> OperatorType;

  IPLinearization(BarrierFunctional const& fu, DomainElement const& x_) 
    : blin(fu,x_), vlin(fu.unconstrainedFunctional,x_), bfun(fu), prepared(false)
  {
    flush();
    bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
    prepared=true; 
  }

    void precompute() {
    if(!prepared) 
    { 
      bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
      prepared=true; 
    }
      vlin.precompute();
      blin.precompute();
    }


/// Return gap parameter
  double getMu() const {return bfun.mu; }

  int cols(int cbegin, int cend) const { return vlin.cols(cbegin, cend);}
  int rows(int cbegin, int cend) const { return vlin.rows(cbegin, cend);}

  void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int colbegin, int colend) const 
  { 
    matb.resize(0);
    if(!prepared) 
    { 
      bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
      prepared=true; 
    }
    vlin.getMatrixBlocks(mat, rbegin, rend, colbegin, colend); 
    blin.getMatrixBlocks(matb, rbegin, rend, colbegin, colend); 
    mat+=matb;
  }
  void getRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const 
  { 
    vlin.getRHSBlocks(rhs, rbegin, rend); 
    if(BarrierFunctional::parameterLin!=1) 
    {  
      if(!prepared) 
      { 
        bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
        prepared=true; 
      }
      blin.getRHSBlocks(rhsb, rbegin, rend); 
      int iend=std::min(rhs.size(),rhsb.size());
      if(BarrierFunctional::parameterLin==0) for(int i=0; i<iend;++i)  rhs[i]+=rhsb[i]; 
      if(BarrierFunctional::parameterLin==2) for(int i=0; i<iend;++i)  rhs[i]=rhsb[i]; 
    }
  }

  DomainElement const& getOriginImpl() const {return vlin.getOriginImpl(); }

  int nRowBlocks() const { return vlin.nRowBlocks(); }
  int nColBlocks() const { return vlin.nColBlocks(); }

  void flush() { vlin.flush(); blin.flush(); matb.resize(0); rhsb.resize(0); prepared=false; }

  void touch() { vlin.touch(); blin.touch(); }

  double getValue() const
  {
    if(!prepared) 
    { 
      bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
      prepared=true; 
    }
    return vlin.getValue()+blin.getValue();
  }

  double getL1norm() const
  {
    if(!prepared) 
    { 
      bfun.prepareConstraintsCache(vlin.getOriginImpl()); 
      prepared=true; 
    }
    return vlin.getL1norm()+blin.getL1norm();
  }

  Barrier const& getBarrier() const {return blin;}
  Unconstrained const& getUnconstrained() const {return vlin;}

private:
  mutable MatrixAsTriplet<Scalar> matb;
  mutable std::vector<Scalar> rhsb;
  Barrier blin;
  Unconstrained vlin;
  BarrierFunctional const& bfun;
  mutable bool prepared;
};



///Traits class to choose the right Bridge::Linearization class
template<class Descr, class Functional, class BarrierFu, int paralin>
  class LinearizationTraits<VariableSet<Descr>, IPFunctional<Functional,BarrierFu, paralin> >
  {
    typedef IPFunctional<Functional,BarrierFu, paralin> Func;
  public:
    typedef IPLinearization<Func,typename Func::AnsatzVars::VariableSet> Linearization;
  };

}

template<class BarrierFunction, int yIdx>
struct FirstOrderEquation
{
  double H,r,mu,upper;
  template<class Cache, class Arg>
  void computeRHS(Cache const& cache,Arg const& arg,double const& trial)
  {
    r= cache.template d1<yIdx>(arg)+H*trial;
 }

  double operator()(double y) const { return H*y-BarrierFunction::db(mu,mu/(upper-y))-r; }
};

template<class BarrierFunction,int yIdx>
struct SecondOrderEquation
{
  double H,r,mu,upper;
  template<class Cache, class Arg>
  void computeRHS(Cache const& cache,Arg  const& arg,double const& trial)
  {
    r=cache.d0()+H*trial;
  }
  double operator()(double y) const { return H*y-BarrierFunction::b(mu,mu/(upper-y))-r; }
};

#include "algorithm/opt_aux/src/include/zeroin.cpp"

/// Pointwise damping strategy
template<template <class BF, int yIdx> class PolyEquation>
class PointwiseCorrection
{
public:
template<class DomainElement, class QuadraticModelBarrier, class UnconstrainedLinearizationAt>
void correct(DomainElement& tx, const DomainElement& x, QuadraticModelBarrier const& blin, UnconstrainedLinearizationAt const& ulin, double mu)
{
  using namespace boost::fusion;
  boost::timer::cpu_timer refTimer;

  typedef typename UnconstrainedLinearizationAt::Functional                Unconstrained;
  typedef typename Unconstrained::AnsatzVars::Spaces                       Spaces;
  typedef typename Unconstrained::Scalar                                       Scalar;

  int const yIdx  = Unconstrained::yIdx;
  int const ySIdx = boost::fusion::result_of::value_at_c<typename Unconstrained::AnsatzVars::Variables,yIdx>::type::spaceIndex;

  typedef typename result_of::template value_at_c<typename DomainElement::Functions, yIdx>::type::Space YSpace;
  int const dim = YSpace::dim;
  YSpace const& yspace(at_c<yIdx>(tx.data).space());

  if(Unconstrained::ConstraintsCache::template bounds<yIdx>::upper==false) return;

  typedef ShapeFunctionCache<typename YSpace::Grid,Scalar> SfCache;
  SfCache sfCache;
  typedef typename result_of::as_vector<typename result_of::transform<Spaces, GetEvaluators<SfCache> >::type>::type SFEvaluators;
  SFEvaluators evaluators(transform(tx.descriptions.spaces,GetEvaluators<SfCache>(&sfCache)));
  
  typename QuadraticModelBarrier::DomainCache bdomc(blin.template createDomainCache());
  typename Unconstrained::ConstraintsCache                                          uconc(ulin.getFunctional().createConstraintsCache(tx));
  typename Unconstrained::DomainCache         udomc(ulin.template createDomainCache(8));

  typedef typename QuadraticModelBarrier::Linearization::Functional::BarrierFunction BarrierFunction;
    
  VariationalArg<Scalar,dim> unitarg;
  unitarg.value=1.0;
  unitarg.gradient=0.0;

  PolyEquation<BarrierFunction,yIdx> pe;
  pe.mu = mu;

  int nnodes(0),maxiter(0),sumiter(0);
  double ratio(2.0);
  std::vector<int> visited(yspace.degreesOfFreedom(),0);

  std::cout << "PwD " << std::flush;

  typedef typename YSpace::Grid::LeafGridView GridView;

  typedef typename GridView::template Codim<0>::Iterator CellIterator;

  for (CellIterator ci=yspace.gridView().template begin<0>(); ci!=yspace.gridView().template end<0>(); ++ci) 
  {

    moveEvaluatorsToCell(evaluators,*ci);
    bdomc.moveTo(*ci);
    uconc.moveTo(*ci);
    udomc.moveTo(*ci);

    std::vector<Dune::FieldVector<typename YSpace::Grid::ctype, dim> >const &  iNodes(at_c<ySIdx>(evaluators).shapeFunctions().interpolationNodes());

    for (int i=0; i<iNodes.size(); ++i) 
    {
      int globalIdx=at_c<ySIdx>(evaluators).globalIndices()[i];
      if(!visited[globalIdx])
      {
        visited[globalIdx]++;
        double trial   = (*at_c<yIdx>(tx.data))[globalIdx][0];
        double iterate = (*at_c<yIdx>( x.data))[globalIdx][0];
        if(trial > iterate)
        {
          moveEvaluatorsToIntegrationPoint(evaluators,iNodes[i]);
          uconc.evaluateAt(iNodes[i],evaluators);
          pe.upper = uconc.upperbound();
          assert(iterate <= pe.upper);

          if(pe.upper < 1e30)
          {
            nnodes++;
            udomc.evaluateAt(iNodes[i],evaluators);
            bdomc.evaluateAt(iNodes[i],evaluators);

// Initialize Polynomial Equation
// H y+b'(y)=r
// H = Lyy (Hessian of the Lagrangian)
// r = b'(y_)+b''(y_)delta y+Lyy y_+

            pe.H     = udomc.template d2<yIdx,yIdx,dim>(unitarg,unitarg)[0][0];
            pe.computeRHS(bdomc,unitarg,trial);
                
            int iterations;
            double upperB = std::min(trial,pe.upper);
            double lowerB = (upperB+iterate)/2.0;

// Here the solution of the pointwise damping equation takes place

            double solution = zeroin(lowerB,upperB,pe,(upperB-lowerB)*1e-12+4*std::numeric_limits<Scalar>::epsilon()*(1+upperB),iterations,100);

            (*at_c<yIdx>(tx.data))[globalIdx][0] = solution;

// Statistics and Warnings

            maxiter = std::max(maxiter,iterations);
            sumiter += iterations;
            ratio=std::min(ratio,std::fabs((solution-iterate)/(trial-iterate)));

//             if(iterations >= 100)
//             {
//               std::cout << "Warning: too many iterations: " << pe.upper-solution << std::endl;
//               std::cout << "Interval: [" << upperB << "," << lowerB << "] " << solution << std::endl;
//               std::cout << "mu:" << pe.mu << " r:" << pe.r << " H:" << pe.H << std::endl;
//             }
//             if(std::fabs((solution-iterate)/(trial-iterate)) < 1e-3 
//                || std::fabs((solution-iterate)/(trial-iterate) > 1.0 
//                             && std::fabs(trial-iterate)>1e-12 )  )
//             {
//               std::cout << "Small ratio:" << std::fabs((solution-iterate)/(trial-iterate)) << std::endl;
//               std::cout << "Interval: [" << upperB << "," << lowerB << "] " << trial << " " << solution << std::endl;
//               std::cout << trial - iterate << "="<< trial-solution << "+" << solution - iterate << std::endl;
//               std::cout << "mu:" << pe.mu << " r:" << pe.r << " H:" << pe.H << std::endl;
//             }
          }
        }
      }
    }      
  }
  std::cout << (double)(refTimer.elapsed().user)/1e9 << " sec.";
  std::cout << "MaxRatio: " << ratio << std::endl;
}
};

template<class Functional, class DomainElement,class Implementation=Bridge::IPLinearization<Functional, DomainElement> > 
class PWDampingChart: public AbstractChart
{
  typedef typename Functional::OptimalityFunctional Unconstrained;
  typedef Functional Combined;

public:

  PWDampingChart(bool includeDuals_=true) :includeDuals(includeDuals_){}

private:

  void addPerturbation(AbstractFunctionSpaceElement & trialIterate, 
                      AbstractFunctionSpaceElement const& correction,
                       AbstractLinearization const& lin,
                       std::vector<AbstractFunctionSpaceElement* > basis = std::vector<AbstractFunctionSpaceElement* >()) const
  {
    trialIterate = lin.getOrigin();

    if(includeDuals)
      trialIterate += correction;
    else
      for(int i=0; i<trialIterate.nComponents();++i)
        if(trialIterate.getRole(i)!="dual")
          trialIterate.axpy(1.0,correction,i);

    using namespace boost::fusion;
    DomainElement const& x=Bridge::getImpl<DomainElement>(lin.getOrigin());
    DomainElement& tx=Bridge::getImpl<DomainElement>(trialIterate);


    Bridge::Linearization<Implementation> const&  limpl(dynamic_cast<Bridge::Linearization<Implementation> const& >(lin));
    double mu = limpl.getLinImpl().getMu();
    
    LinearizationAt<Unconstrained> ulin(limpl.getLinImpl().getUnconstrained().getLinImpl().getFunctional(),tx);
    QuadraticModel<LinearizationAt<Combined>,false> blin(linearization(limpl.getLinImpl().getBarrier().getLinImpl().getFunctional(),x),
                                                         Bridge::getImpl<DomainElement>(correction));
    PointwiseCorrection<FirstOrderEquation> pc;

    pc.correct(tx, x, blin, ulin, mu);
  }
private:
  bool includeDuals;
};

/// Pointwise damping strategy
template<class DomainElement, class QBarrierModel, class ParameterModel, class LModel>
void corrTangent(DomainElement& tx,const DomainElement& x, QBarrierModel const& blin, ParameterModel const& plin, LModel const& ulin, double muLast, double mu)
  {
    using namespace boost::fusion;

    typedef typename LModel::Functional Unconstrained;
    typedef typename LModel::AnsatzVars::Spaces    Spaces;
    typedef typename LModel::AnsatzVars::IndexSet  IndexSet;
    typedef typename Unconstrained::Entity  Entity;
    typedef typename LModel::Scalar  Scalar;
    int const yIdx = Unconstrained::yIdx;
    typedef typename result_of::template value_at_c<typename DomainElement::Functions, yIdx>::type::Space YSpace;
    int const dim = YSpace::dim;

    double deltaMu = muLast-mu;

    if(Unconstrained::ConstraintsCache::template bounds<yIdx>::upper==0) return;

    YSpace const& yspace(at_c<yIdx>(tx.data).space());
    typename YSpace::Evaluator ysfs(yspace);

    typedef ShapeFunctionCache<typename YSpace::Grid,Scalar> SfCache;
    SfCache sfCache;
    
    typedef typename result_of::as_vector<typename result_of::transform<Spaces, GetEvaluators<SfCache> >::type>::type Evaluators;
    Evaluators evaluators(transform(tx.descriptions.spaces,GetEvaluators<SfCache>(&sfCache)));
  
    typename QBarrierModel::DomainCache bdc(blin.template createDomainCache());
    typename ParameterModel::DomainCache pdc(plin.template createDomainCache());
    typename Unconstrained::ConstraintsCache udc(ulin.getFunctional().createConstraintsCache(tx));
    
    VariationalArg<Scalar,dim> unitarg;
    unitarg.value[0]=1.0;
    for(int i=0; i < dim; ++i) unitarg.gradient[0][i]=0.0;

    typedef typename IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;
    for (CellIterator ci=yspace.indexSet().template begin<0,Dune::All_Partition>();
         ci!=yspace.indexSet().template end<0,Dune::All_Partition>(); ++ci) 
    {
      ysfs.moveTo(*ci);
      moveEvaluatorsToCell(evaluators,*ci);
      bdc.moveTo(*ci);
      pdc.moveTo(*ci);
      udc.moveTo(*ci);
      std::vector<Dune::FieldVector<typename YSpace::Grid::ctype, dim> >const &  iNodes(ysfs.shapeFunctions().interpolationNodes());
      FirstOrderEquation<typename QBarrierModel::Linearization::Functional::BarrierFunction,yIdx> pe;
      for (int i=0; i<iNodes.size(); ++i) 
      {
        ysfs.evaluateAt(iNodes[i]);
        moveEvaluatorsToIntegrationPoint(evaluators,iNodes[i]);
        udc.evaluateAt(iNodes[i],evaluators);
        bdc.evaluateAt(iNodes[i],evaluators);
        pdc.evaluateAt(iNodes[i],evaluators);
        pe.aqp = udc.quadTermA();
        pe.aq=(bdc.template d1<yIdx,dim>(unitarg))[0]+deltaMu*(pdc.template d1<yIdx,dim>(unitarg))[0]+udc.quadTermB();                 
        pe.mu= mu;
        double yc2;
        double trial= -(*at_c<yIdx>(tx.data))[ysfs.globalIndices()[i]][0]+udc.upperbound();
        double iterate= -(*at_c<yIdx>(x.data))[ysfs.globalIndices()[i]][0]+udc.upperbound();
        double upper = std::max(trial,iterate);
        double lower = std::min(trial,iterate);
        int iterations;
        yc2=bisection(std::max(0.0,lower),upper,pe,upper*1e-6+4*std::numeric_limits<Scalar>::epsilon()*(1+upper),iterations);
        (*at_c<yIdx>(tx.data))[ysfs.globalIndices()[i]][0]= udc.upperbound()-yc2;
      }      
    }
  }

struct StepEquation
{
  double k;   // order of predictor (k=0 : classical, k=1 : tangential)
  double err; // estimated iteration error of last step
  double eta; // order coefficent : (k=0 : slope, k=1 : curvature); 
  double mu;  // homotopy parameter
  double alpha; // assumed negative exponent for model of slope of path (usually alpha = 1/2)
  double Theta; // desired contraction
  double omega; // aff. cov. Lipschitz constant for Newton's method
  double beta;  // assumed exponent for model of omega (usually beta = 1/2) 
  double operator()(double sigma) const 
{ return -(err+eta*std::pow(mu,k+1)*std::pow(sigma,-alpha-k)*std::pow(1-sigma,k+1)-Theta/omega*std::pow(sigma,beta)); }
};

template<class IPF, class DomainVector>
class InteriorPointTangentialPredictor : public HomotopyBase
{
  typedef typename IPF::OptimalityFunctional Unconstrained;
  typedef IPF Barrier;
  typedef Bridge::IPLinearization<IPF, DomainVector> Implementation;

public:
  InteriorPointTangentialPredictor(NewtonsMethod& n_, AbstractNorm const& norm_, InteriorPointParameters& p_, 
                          AbstractLinearSolver& solver_, AbstractNorm const* normPlain_=0,
                          NewtonsMethod* finalsolver_=0) : 
    HomotopyBase(n_, p_),
    norm(norm_),
    solver(solver_),
    stp(0),
    normPlain(normPlain_),
    finalsolver(finalsolver_)
  {
    if(!normPlain_) normPlain=&norm;
    if(!finalsolver_) finalsolver=&n_;
  }

  virtual void initialize() { 
    iterate=trialIterate->clone(); 
    predictor=trialIterate->clone(); 
    tangent=trialIterate->clone(); 
    ptangent=trialIterate->clone(); 
    nNewton=0;
  };
  virtual void finalize() {};


  virtual double muOnSuccess(int step) {
    double b = pp.omega*pp.accuracyCorrector;
    double c = pp.eta*pp.omega*p.mu;
    double a = c+p.desiredContraction;
    double sigma = (b+sqrt(b*b+4*a*c))/(2*a);
    sigma = sigma*sigma;

    StepEquation f;
    f.k=1;
    f.err = pp.accuracyCorrector;
    f.eta = pp.eta;
    f.mu = p.mu;
    f.alpha = 0.5;
    f.Theta = p.desiredContraction;
    f.omega = pp.omega;
    f.beta = 0.5;
    int iter;
    double sig = bisection(1e-20,1.0-1e-15,f,1e-6,iter);
    sigma = sig; 
    sigma = std::min(sigma,p.reductionFactorWorst);
    sigma = std::max(sigma,p.reductionFactorBest);
    if(step==2) sigma=p.reductionStart;

    return sigma*p.mu;
  };

  virtual double muOnFailure() { 
    p.sigma=0.5+p.sigma/2.0;    
    return p.mu*p.sigma;
  };

  virtual double muFinal() { return p.muFinal;};

  virtual void updateModelOfHomotopy() 
  {
    *predictor -= *trialIterate;
    double n(norm(*predictor));
    *iterate *= 0.0;
     std::unique_ptr<AbstractFunctional> fuPtr(functional->getParameterLinFunctional(makePars(p.muTrial.value())));
     std::unique_ptr<AbstractLinearization> paraLinearization(fuPtr->getLinearization(*trialIterate));
     solver.solve(*tangent,*paraLinearization);
     double nD=norm(*tangent);
     double nDp=norm(*ptangent);

     if(stp>=2) 
     {
       *ptangent -= *tangent;
       std::cout << "Difference in Tangent vectors:" << norm(*ptangent)/std::max(nD,nDp) << std::endl; 
     }

     std::unique_ptr<AbstractFunctional> fuVPtr(functional->getLinFunctionValue(makePars(p.muTrial.value())));
     std::unique_ptr<AbstractLinearization> paraLinFu(fuVPtr->getLinearization(*trialIterate));
     paraLinFu->evald(*iterate);
     double slopeJ = iterate->applyAsDualTo(*tangent);
     pp.jpl=p.mu.value()*slopeJ;
     double nDP=(*normPlain)(*tangent);
     *ptangent -= *tangent;
     double nDiff=norm(*ptangent);

     pp.slope=nDP;
     pp.curvature=nDiff/p.deltaMu;
     pp.eta=pp.curvature;
     GuardedCovariantNewtonParameters const& cnp(dynamic_cast<GuardedCovariantNewtonParameters const&>(corrector.getParameters()));
     pp.omega=cnp.omega0.value();
     pp.accuracyCorrector=cnp.accuracyReached;
     std::cout << "Norm:" << n/p.deltaMu << " ||x'||_2:" << nDP << " Dmu:" << p.deltaMu << " ||x'||_x:"<< pp.eta << " omega:" << pp.omega << std::endl;
  };

  virtual void computePredictor(int step)
  {
    stp=step;
    if(p.mu==p.muTrial) return;
    if(step==2) 
    {
      std::unique_ptr<AbstractFunctional> fuPtr(functional->getParameterLinFunctional(makePars(p.muTrial.value())));
      std::unique_ptr<AbstractLinearization> paraLinearization(fuPtr->getLinearization(*trialIterate));
      solver.solve(*tangent,*paraLinearization);
    }
    *predictor=*iterate;
    predictor->axpy(p.mu-p.muTrial,*tangent);
    DomainVector& tx=Bridge::getImpl<DomainVector>(*iterate);
    DomainVector& px=Bridge::getImpl<DomainVector>(*predictor);
    DomainVector const& cx=Bridge::getImpl<DomainVector>(*tangent);
    std::unique_ptr<AbstractFunctional> fuPtr(functional->getFunctional(makePars(p.mu.value())));
    std::unique_ptr<AbstractLinearization> linPtr(fuPtr->getLinearization(*iterate));
    std::unique_ptr<AbstractLinearization> paraLinearization(fuPtr->getLinearization(*iterate));

    Bridge::Linearization<Implementation> const&  limpl(dynamic_cast<Bridge::Linearization<Implementation> const& >(*linPtr));
    Bridge::Linearization<Implementation> const&  plimpl(dynamic_cast<Bridge::Linearization<Implementation> const& >(*paraLinearization));
    
    LinearizationAt<Barrier> blin(limpl.getLinImpl().getBarrier().getLinImpl().getFunctional(),tx);
    LinearizationAt<Barrier> plin(plimpl.getLinImpl().getBarrier().getLinImpl().getFunctional(),tx);
    LinearizationAt<Unconstrained> ulin(limpl.getLinImpl().getUnconstrained().getLinImpl().getFunctional(),tx);
    Dune::BlockVector<Dune::FieldVector<typename Unconstrained::Scalar, 1> > temp(*boost::fusion::at_c<Unconstrained::yIdx>(px.data));
    corrTangent(px, tx, blin, plin, ulin, p.muTrial, p.mu);
    for (int i=0; i<temp.size(); ++i) 
      if((*boost::fusion::at_c<Barrier::yIdx>(cx.data))[i][0] < 0) (*boost::fusion::at_c<Barrier::yIdx>(px.data))[i][0] = temp[i];
    *trialIterate=*predictor;      
    *predictor-=*iterate;
    std::unique_ptr<AbstractFunctional> fuPtr2(functional->getParameterLinFunctional(makePars(p.muTrial.value())));
    std::unique_ptr<AbstractLinearization> paraLinearization2(fuPtr2->getLinearization(*trialIterate));
    solver.solve(*ptangent,*paraLinearization2);
  }


  virtual void initializeCorrector() 
  { 
    if(p.muTrial < p.muFinal) corrector.setDesiredAccuracy(p.accfinal);
    else
    {
      if(pp.eta.isValid()) corrector.setDesiredAccuracy(1000000*std::max(p.accfinal*0.25,lengthOfPath()*p.relDistanceToPath));
      else corrector.setDesiredAccuracy(sqrt(p.muTrial)*p.relDistanceToPath);
      corrector.setDesiredRelativeAccuracy(p.relDistanceToPath);//std::sqrt(p.muTrial/p.mu));
    }
  };
  virtual void finalizeCorrector() {    nNewton += corrector.stepsPerformed(); pp.newtonSum=nNewton;};

  virtual void finalizeHomotopy()
  {
    std::cout << "----------------------- Final Correction Step ---------------------" << std::endl;
    finalsolver->resetParameters();
    finalsolver->setDesiredAccuracy(p.accfinal); 
    finalsolver->setDesiredRelativeAccuracy(1.0); 
    std::unique_ptr<AbstractFunctional> fuPtr(functional->getFunctional(makePars(p.muTrial.value())));
    std::cout << "muTrial: " << p.muTrial << std::endl;
    finalsolver->reportOnIteration(true);
    finalsolver->solve(fuPtr.get(), *trialIterate);
    p.correctorTermination=corrector.getParameters().termination;
  }

  virtual void updateIterate() { 
    *iterate=*trialIterate;
    p.j=functional->getFunctional(makePars(p.muTrial.value()))->getLinearization(*iterate)->eval();
    jModelL.update(p.muTrial,p.j);
    jModelL.fixUpdate();
    p.jp=jModelL.getValue(p.muTrial);
    std::cout << "Est. Error in Functional:" << p.jp << std::endl;
    pp.timemeasured=overalltime.elapsed();
  };

  virtual void recoverIterate() { *trialIterate=*iterate; }


  virtual double lengthOfPath() 
  { 
    if(pp.eta.isValid())
      return 2*p.mu*pp.slope; 
    else
      return 1e300;
  };

  virtual void logQuantities()
  {
    HomotopyBase::logQuantities();
    pp.logStep();
    printDiagnosis();
  }

  void printDiagnosis()
  {
    std::ofstream eout("eta.log"); pp.eta.print(eout); 
    std::ofstream epout("slope.log"); pp.slope.print(epout); 
    std::ofstream cpout("curvature.log"); pp.curvature.print(cpout); 
    std::ofstream nout("newtonSum.log"); pp.newtonSum.print(nout); 
    std::ofstream oout("omega.log"); pp.omega.print(oout); 
    std::ofstream sout("sigma.log"); p.sigma.print(sout); 
    std::ofstream outmu("mu.log"); p.muTrial.print(outmu); 
    std::ofstream outj("j.log"); p.j.print(outj); 
    std::ofstream outjp("jp.log"); p.jp.print(outjp); 
    std::ofstream outjpl("jplin.log"); pp.jpl.print(outjpl); 
    std::ofstream outac("ac.log"); pp.accuracyCorrector.print(outac); 
    std::ofstream outtm("time.log"); pp.timemeasured.print(outtm); 
  }

private:
  std::unique_ptr<AbstractFunctionSpaceElement> iterate, predictor, tangent, ptangent;

  AbstractNorm const& norm;
  AbstractLinearSolver& solver;
  InteriorPointParametersSqrt pp;
  JModelLin jModelL;
  int stp;
  AbstractNorm const* normPlain;
  NewtonsMethod* finalsolver;
  int nNewton;
  boost::timer::cpu_timer overalltime;
};


/// Policy class that performs a pointwise damped Newton step
// template<class Functional, class DomainElement,class Implementation=Bridge::IPLinearization<Functional, DomainElement> > 
// class StepPolicyPointwiseDamping : public StepPolicyPlainNewton
// {
// public:
//   typedef typename Functional::OptimalityFunctional Unconstrained;
//   typedef Functional Combined;

//   virtual void setLinearization(AbstractLinearization& lin, AbstractLinearSolver& solver) 
//   {
//     linPtr= &lin;
//     linearSolver=&solver;
//     mu=dynamic_cast<Bridge::Linearization<Implementation>& >(lin).getLinImpl().getMu();
//   }

//   virtual void getTrialIterate(
//     AbstractFunctionSpaceElement& trialIterate, 
//     AbstractFunctionSpaceElement const& correction, 
//     AbstractFunctionSpaceElement const& iterate,
//     double damping)
//   {
//     trialIterate = iterate;
//     trialIterate.axpy(damping,correction);
//     dcorrection = correction.clone();
//     *dcorrection*=damping;
//     correctIterate(trialIterate, *dcorrection, *linPtr, iterate, mu); 
//   }
// private:

//   std::unique_ptr<AbstractFunctionSpaceElement> dcorrection;

// /* Makes a pointwise damping correction. The following parameters are used:
//  * - trialIterate: iterate after (simplified) Newton step -> iterate after pointwise damping
//  * - correction:   (simplified) Newton correction
//  * - lin       :   point of linearization/ iterate before Newton step
//  * - siterate  :   point, where rhs is evaluated, iterate after Newton step
//  *
//  */
//   void correctIterate(AbstractFunctionSpaceElement & trialIterate, 
//                       AbstractFunctionSpaceElement const& correction,
//                       AbstractLinearization const& lin, 
//                       AbstractFunctionSpaceElement const& siterate, 
//                       double mu)
//   {
//     using namespace boost::fusion;
//     DomainElement const& x=Bridge::getImpl<DomainElement>(lin.getOrigin());
//     DomainElement const& sx=Bridge::getImpl<DomainElement>(siterate);
//     DomainElement& tx=Bridge::getImpl<DomainElement>(trialIterate);

//     Bridge::Linearization<Implementation> const&  limpl(dynamic_cast<Bridge::Linearization<Implementation> const& >(lin));
    
//     LinearizationAt<Unconstrained> ulin(limpl.getLinImpl().getUnconstrained().getLinImpl().getFunctional(),tx);
//     QuadraticModel<LinearizationAt<Combined>,true> blin(linearization(limpl.getLinImpl().getBarrier().getLinImpl().getFunctional(),sx),
//                                                        linearization(limpl.getLinImpl().getBarrier().getLinImpl().getFunctional(),x),
//                                                        Bridge::getImpl<DomainElement>(correction));
//     PointwiseCorrection<FirstOrderEquation> pc;

//     pc.correct(tx, sx, blin, ulin, mu);
//   }
//   double mu;

// };


#endif
