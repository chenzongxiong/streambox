#ifndef DUNE_BRIDGE_HH
#define DUNE_BRIDGE_HH

#include <cmath>
#include <memory>
#include <utility>
#include <dune/istl/operators.hh>
#include "linalg/linearsystem.hh"
#include "linalg/triplet.hh"
#include "fem/assemble.hh"
#include "fem/linearspace.hh"
#include "fem/functional_aux.hh"
#include "algorithm/abstract_interface.hh"

namespace Kaskade
{
  // forward declaration
  namespace Bridge { template <class> class Linearization; }

template<class Functional>
struct DiscreteBlockTraits
{
  typedef typename Functional::AnsatzVars::VariableSet Vars;
  typedef typename Functional::Scalar Scalar;

  static constexpr bool anyPresent = false;
  static void getMatrixBlock(MatrixAsTriplet<Scalar>& mat,
      Functional const& fu,
      Vars const& origin,
      int row,
      int col)
  {}
  static void getRHSBlock(std::vector<Scalar>& rhs,
      Functional const& fu,
      Vars const& origin,
      int row)
  {}

  static Scalar getValue(Functional const& fu,
      Vars const& origin)
  {
    return 0.0;
  }
};

template <typename Scalar>
std::ostream& operator<<(std::ostream &s, std::vector<Scalar> const& vec)
{
  s << "[\n";
  for(size_t i=0; i<vec.size(); ++i)
    s << "  " << vec[i] << std::endl;
  s << "]\n";
  return s;
}

template<class Functional, class DomainElement>
struct InDomainTraits
{
  static bool inDomain(DomainElement const& x){ return true; }
};


namespace Bridge{
  template <class> class Vector;

  template<class Implementation> Implementation& getImpl(AbstractFunctionSpaceElement& v);
  template<class Implementation> Implementation const& getImpl(AbstractFunctionSpaceElement const& v);

//  /// Bridge::Linearization class that uses a VariationalFunctionalAssembler to create linear systems
//  /** Implements AbstractLinearization */
//  template<class Functional>
//  class KaskadeLinearization : public AbstractLinearization, public SparseLinearSystem
//  {
//  public:
//
//    static const int nThreads = 16;
//
//    typedef typename Functional::AnsatzVars::VariableSet DomainElement;
//    typedef typename Functional::TestVars::VariableSet ImageElement;
//    typedef typename Functional::Scalar Scalar;
//    typedef LinearizationAt<Functional> Implementation;
//    typedef VariationalFunctionalAssembler<Implementation> Assembler;
//    typedef typename DomainElement::Descriptions::template CoefficientVectorRepresentation<>::type CoefficientVector;
//     typedef Dune::LinearOperator<CoefficientVector, CoefficientVector> OperatorType;
//     typedef OperatorType Operator;
//
//    KaskadeLinearization() = delete;
//
//    /// Creation of a linearization for a functional fu at x_
//    KaskadeLinearization(Functional const& fu_, DomainElement const& x_)
//    : x(x_), fu(fu_), lin(fu_,x), ass(x.descriptions.spaces)
//    {
//      flush();
//    }
//
//    virtual ~KaskadeLinearization() { changed(); flushconn.disconnect(); }
//
//    /// Number of columns of components [cbegin, cend)
//    int cols(int cbegin, int cend) const
//    {
//      if(cbegin < cend)
//        return x.descriptions.degreesOfFreedom(cbegin,cend);
//      return 0;
//    }
//
//    /// Number of rows of components [rbegin, rend)
//    int rows(int rbegin, int rend) const
//    {
//      if(rbegin < rend)
//        return x.descriptions.degreesOfFreedom(rbegin,rend);
//      return 0;
//    }
//
//    void precompute() {
//      doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);
//    }
//
//    /// write blocks of the hessian matrix into mat
//    void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int colbegin, int colend) const
//    {
//      doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);
//      mat = ass.template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend,colbegin,colend);
//
//      if(DiscreteBlockTraits<Functional>::anyPresent)
//      {
//        MatrixAsTriplet<Scalar> matD;
//        getDiscreteMatrixBlocks(matD,rbegin,rend,colbegin,colend);
//        mat+=matD;
//      }
//    }
//
//    /// write components of the gradient into rhs
//    void getRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const
//    {
//      doAssemble(Assembler::VALUE | Assembler::RHS);
//      rhs.resize(rows(rbegin,rend),0.0);
//      ass.toSequence(rbegin,rend,rhs.begin());
//      if(DiscreteBlockTraits<Functional>::anyPresent)
//      {
//        addDiscreteRHSBlocks(rhs,rbegin,rend);
//      }
//    }
//
//    MatrixAsTriplet<Scalar> assembleGradientCheck(int const lastBlockId)
//        {
//      int const firstBlockId = 0;
//      std::cout << "assembling for gradient check" << std::flush;
//
//      MatrixAsTriplet<Scalar> result, reference;
//      flush();
//      getMatrixBlocks(reference, 0, lastBlockId, 0, lastBlockId);
//      std::vector<Scalar> ref_sol, tmp;
//      size_t numRows = rows(firstBlockId, lastBlockId), numCom = Functional::TestVars::template Components<0>::m;
//      std::cout << "1: numRows: " << numRows << ", " << numCom << std::endl;
//      Scalar const increment = 1.0e-9;
//
//      tmp.resize(numRows, 0.0);
//      ref_sol.resize(numRows, 0.0), tmp.resize(numRows,0.0);
//      flush();
//      doAssemble(Assembler::RHS);
//      ass.toSequence(firstBlockId,lastBlockId, ref_sol.begin());
//      for(size_t i=0; i<boost::fusion::at_c<0>(x.data).coefficients().N(); ++i)
//      {
//        for(int j=0; j<numCom; ++j){
//          tmp.clear();
//          tmp.resize(numRows,0.0);
//          boost::fusion::at_c<0>(x.data).coefficients()[i][j] += increment;
//          flush();
//          doAssemble(Assembler::RHS);
//          boost::fusion::at_c<0>(x.data).coefficients()[i][j] -= increment;
//
//          ass.toSequence(firstBlockId, lastBlockId, tmp.begin());
//
//          for(size_t k=0; k<tmp.size(); ++k)
//          {
//            result.addEntry(k,i*numCom+j, (tmp[k] - ref_sol[k])/increment );
//          }
//        }
//        std::cout << "." << std::flush;
//      }
//      std::cout << "done." << std::flush << std::endl << std::endl;
//
//      result *= -1.0;
//      result += reference;
//
//      std::cout << "Error:" << std::endl;
//      std::cout << result << std::endl;
//
//      return result;
//        }
//
//    std::vector<Scalar> checkRealGradient(int const lastBlockId)
//         {
//      std::vector<Scalar> result, reference;
//      int const firstBlockId = 0;
//      getRHSBlocks(reference, firstBlockId, lastBlockId);
//      size_t numRows = rows(firstBlockId, lastBlockId);
//      Scalar const increment = 1.0e-9;
//
//      result.resize(numRows,0.0);
//      flush();
//      Scalar const value = getValue();
//      Scalar tmp(0);
//
//      for(size_t i=0; i<boost::fusion::at_c<0>(x.data).coefficients().N(); ++i)
//      {
//        boost::fusion::at_c<0>(x.data).coefficients()[i][0] += increment;
//        flush();
//        tmp = getValue();
//        boost::fusion::at_c<0>(x.data).coefficients()[i][0] -= increment;
//        result[i] = (tmp - value)/increment;
//      }
//      std::cout << "Error\n";
//      for(size_t i=0; i<result.size(); ++i){
//        std::cout << result[i] << " : " << reference[i] << " -> " << (result[i] + reference[i]) << std::endl;
//        result[i] += reference[i];
//        if(fabs(result[i]) > 1.0e-8)
//          std::cout << "  " << result[i] << std::endl;
//      }
//      // std::cout << result << std::endl;
//
//      return result;
//         }
//
//
//    /// return number of columns
//    constexpr int nColBlocks() const {return Functional::AnsatzVars::noOfVariables;}
//
//    /// return number of rows
//    constexpr int nRowBlocks() const { return Functional::TestVars::noOfVariables; }
//
//    /// return point of linearization
//    DomainElement const& getOriginImpl() const {return x;}
//
//    // return point of linearization
//    DomainElement getX() const { return x; }
//
//    void setX(DomainElement const& x_){ x = x_; }
//
//    /// return the implementation
//    Implementation const& getLinImpl() const {return lin; }
//
//    /// flush all data, gathered so far
//    void flush() { changed(); ass.flush((Assembler::VALUE | Assembler::RHS | Assembler::MATRIX)); }
//
//    /// return whether x is in the domain of definition
//    bool inDomain(DomainElement const& x) { return InDomainTraits<Functional,DomainElement>::inDomain(x); }
//
//    /// return the current value Functional(Origin)
//    double getValue() const
//    {
//      doAssemble(Assembler::VALUE);
//      Scalar value =ass.functional();
//      if(DiscreteBlockTraits<Functional>::anyPresent)
//        value += DiscreteBlockTraits<Functional>::getValue(fu,x);
//      return value;
//    }
//
//    double eval() const { return getValue(); }
//    double evalL1norm() const { return getL1norm(); }
//    void evald(AbstractFunctionSpaceElement& v, int rbegin=0, int rend=nRowBlocks()) const
//    {
//      assert(rbegin<=rend);
//      std::vector<Scalar> rhs(rows(rbegin,rend),0.0);
//      getRHSBlocks(rhs,rbegin,rend);
//      dynamic_cast<Vector<ImageElement>& >(v).read(rhs,rbegin,rend);
//    }
//    void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=nRowBlocks(), int cbegin=0, int cend=nColBlocks()) const
//    {
//      assert(cbegin<=cend);
//      assert(rbegin<=rend);
//      std::vector<Scalar> result, argument;
//      dynamic_cast<Vector<ImageElement> const& >(x).write(argument,cbegin,cend);
//      dynamic_cast<Vector<ImageElement>& >(y).write(result,rbegin,rend);
//      MatrixAsTriplet<Scalar> mat;
//      getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
//      mat.axpy(result, argument);
//      dynamic_cast<Vector<ImageElement>& >(y).read(result,rbegin,rend);
//    }
//    AbstractFunctionSpaceElement const& getOrigin() const
//    {
//      xptr.reset(new Vector<DomainElement>(x));
//      return *xptr;
//    }
//    void connectToSignalForFlush(boost::signals2::signal0<void>& sig) {
//      if(flushconn.connected()) flushconn.disconnect();
//      flushconn=sig.connect(boost::bind(&KaskadeLinearization::flush, this));
//    }
//
//    double getL1norm() const
//    {
//      doAssemble(Assembler::VALUE);
//      Scalar value =0;//ass.fL1norm();
//      if(DiscreteBlockTraits<Functional>::anyPresent)
//        value += std::fabs(DiscreteBlockTraits<Functional>::getValue(fu,x));
//      return value;
//    }
//
//    Assembler& getValidAssembler() const
//    {
//      doAssemble((Assembler::VALUE | Assembler::RHS | Assembler::MATRIX));
//      return ass;
//    }
//
//    Functional const& getFunctional() const
//    {
//      return fu;
//    }
//
//  private:
//    void addDiscreteRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const
//    {
//      for(int i=rbegin; i<rend; ++i)
//      {
//        std::vector<Scalar> rhsD;
//        DiscreteBlockTraits<Functional>::getRHSBlock(rhsD,fu,x,i);
//        if(rhsD.size() > 0)
//        {
//          assert(rhsD.size() == rows(i,i+1));
//          for(int k=rows(rbegin,i), l=0;k<rows(rbegin,i+1);++k, ++l)
//            rhs[k]+=rhsD[l];
//        }
//      }
//    }
//
//    void getDiscreteMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int colbegin, int colend) const
//    {
//      for(int i=rbegin; i<rend; ++i)
//        for(int j=colbegin; j<colend; ++j)
//        {
//          MatrixAsTriplet<Scalar> matD;
//          DiscreteBlockTraits<Functional>::getMatrixBlock(matD, fu,x,i,j);
//          matD.shiftIndices(rows(rbegin,i),cols(colbegin,j));
//          mat+=matD;
//        }
//    }
//
//    void doAssemble(int flags) const
//    {
//      int toDoFlag= ((~ass.valid()) & flags);
//      if(toDoFlag!=0)
//      {
//        ass.assemble(lin,toDoFlag,nThreads);
//      }
//    }
//
//
//    DomainElement x;
//    Functional const& fu;
//    Implementation lin;
//    mutable Assembler ass;
//    mutable std::unique_ptr<Vector<DomainElement> > xptr;
//    boost::signals2::connection flushconn;
//};

//  template<class Functional>
//  class LinearizationTraits<typename Functional::AnsatzVars::VariableSet, Functional>
//  {
//  public:
//    typedef KaskadeLinearization<Functional> Linearization;
//  };


  /// Bridge class for an error estimator, implements AbstractErrorEstimator
  /** An implementation must support:
   *
   * - typedef DomainElement  (Vector type)
   *
   * - typedef Estimate (Structure where error estimate and indicators are stored)
   *
   * - std::unique_ptr<AbstractErrorEstimate> createEstimate(DomainElement const&, DomainElement const&)
   */
  template<typename Implementation>
  class ErrorEstimator : public AbstractErrorEstimator
  {
    typedef typename Implementation::DomainElement VectorImplementation;
  public:
    typedef typename Implementation::Estimate Estimate;

    ErrorEstimator(std::unique_ptr<Implementation>& imp) : myImplementation(std::move(imp)) {}

    virtual std::unique_ptr<AbstractErrorEstimate> createEstimate(AbstractFunctionSpaceElement const& correction,
        AbstractLinearization const& lin) const
        {
      VectorImplementation const& ci = getImpl<VectorImplementation>(correction);
      typedef typename LinearizationTraits<VectorImplementation,typename Implementation::Functional>::Linearization LinImpl;
      return myImplementation->createEstimate(ci,dynamic_cast<Linearization<LinImpl> const&>(lin));
        }
  private:
    std::unique_ptr<Implementation> myImplementation;
  };

  template<class ErrorEst, class GridMan>
  ErrorEst extendMarkings(ErrorEst const&, GridMan &);


  /// Implements AbstractAdaptiveGrid, uses a bulk criterion as a marking strategy
  /** Gets a signal from the GridManager, if Grid is going to change. Emits a signal to inform others
   */
  template <class GridManager, class Estimate>
  class AdaptiveGrid : public AbstractAdaptiveGrid
  {
  public:
    AdaptiveGrid(GridManager& gm) 
    : gridMan(gm), nMarked(0)
    {
      gridChange=gm.signals.informAboutRefinement.connect(
          int(GridSignals::freeResources),
          boost::bind(&AdaptiveGrid<GridManager,Estimate>::onGridChange,this,_1));
    }

    virtual int size() { return gridMan.grid().size(0); }

    virtual void adapt() { gridMan.preAdapt(); gridMan.adapt(); gridMan.postAdapt(); nMarked = 0; };

    virtual void mark(AbstractErrorEstimate& estimate, double portion) 
    {
      Estimate cdmarks(CellData<typename GridManager::Grid,int>(markByBulkCriterion(dynamic_cast<Estimate&>(estimate).getData(),portion)));

      nMarked = cdmarks.getData().sum();
      //      if(cdmarks.getData().sum() >= 00)
      gridMan.mark(cdmarks.getData());
      //       else
      //       {
      //         Estimate
      //           extMarks(extendMarkings<Estimate, GridManager>(cdmarks,gridMan));
      //         gridMan.mark(extMarks.getData());
      //       }
    }

    virtual int getNMarked() { return nMarked; }

    virtual void flushMarks() { gridMan.flushMarks(); nMarked = 0; }

    virtual ~AdaptiveGrid() { gridChange.disconnect(); }

  private:
    void onGridChange(GridSignals::Status stat) 
    { 
      if(stat==GridSignals::BeforeRefinement) gridWillChange();
    }

    GridManager& gridMan;
    int nMarked;
    boost::signals2::connection gridChange;
  };


  //--------------------------------------------------------------------------------

  ///Convenience routine: makes an ErrorEstimator of the right type
  template<class Est>
  std::unique_ptr<ErrorEstimator<Est> > makeErrorEstimator(Est* est){return makeUP(new ErrorEstimator<Est>(makeUP(est)));}

}

}// namespace Kaskade
#endif 
