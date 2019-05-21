/*
 * kaskadeBridge.hh
 *
 *  Created on: 06.12.2013
 *      Author: Lars Lubkoll, Anton Schiela
 */

#ifndef KASKADE_BRIDGE_HH_
#define KASKADE_BRIDGE_HH_

#include <boost/signals2.hpp>
#include <boost/bind.hpp>

#include "algorithm/abstract_interface.hh"
#include "algorithm/dune_bridge.hh"
#include "algorithm/newton_bridge.hh"

namespace Kaskade
{
  template <class Value>
  void assignIfNegative(Value& val, Value newVal)
  {
    if(val < 0) val = newVal;
  }

  namespace Bridge
  {
    template <class Linearization>
    class ConnectedLinearization : public Linearization, public AbstractFlushConnection
    {
    public:
      template <typename... Args>
      ConnectedLinearization(const Args&... args) : Linearization(args...)
      {}

      virtual ~ConnectedLinearization() { changed(); flushconn.disconnect(); }

      virtual void flush() { changed(); Linearization::flush(); }

      virtual void connectToSignalForFlush(boost::signals2::signal<void ()>& sig)
      {
        if(flushconn.connected()) flushconn.disconnect();
        flushconn=sig.connect(boost::bind(&ConnectedLinearization<Linearization>::flush, this));
      }

      boost::signals2::signal<void ()> changed;
      boost::signals2::connection flushconn;
    };



    /// Bridge::Linearization class that uses a VariationalFunctionalAssembler to create linear systems
    /** Implements AbstractLinearization */
    template<class Functional>
    class KaskadeLinearization : public AbstractLinearization, public SparseLinearSystem
    {
    public:

      static const int nThreads = 32;

      typedef typename Functional::AnsatzVars::VariableSet DomainElement;
      typedef typename Functional::TestVars::VariableSet ImageElement;
      typedef typename Functional::Scalar Scalar;
      typedef LinearizationAt<Functional> Implementation;
      typedef VariationalFunctionalAssembler<Implementation> Assembler;
      typedef typename DomainElement::Descriptions::template CoefficientVectorRepresentation<>::type CoefficientVector;
      typedef Dune::LinearOperator<CoefficientVector, CoefficientVector> OperatorType;
      typedef OperatorType Operator;

      /*      KaskadeLinearization()
      : x(0), fu(nullptr), lin(*fu,x), ass(nullptr), xptr(nullptr)
      {
        flush();
      }

      /// Creation of a linearization for a functional fu at x_
      KaskadeLinearization(Functional const& fu_)
      : x(0), fu(&fu_), lin(*fu,x), ass(nullptr), xptr(nullptr)
      {
        flush();
      }*/

      /// Creation of a linearization for a functional fu at x_
      KaskadeLinearization(Functional const& fu_, DomainElement const& x_)
        : x(x_), fu(&fu_), lin(*fu,x), ass(new Assembler(x.descriptions.spaces)), xptr(new Vector<DomainElement>(x))
      {
        flush();
      }

      /// Creation of a linearization for a functional fu at x_
      KaskadeLinearization(Functional const& fu_, DomainElement const& x_, std::shared_ptr<Assembler> const& ass_)
        : x(x_), fu(&fu_), lin(*fu,x), ass(ass_), xptr(new Vector<DomainElement>(x))
      {
        flush();
      }

      KaskadeLinearization(KaskadeLinearization const& other)
        : x(other.x), fu(other.fu), lin(*fu,x), ass(other.ass), xptr(new Vector<DomainElement>(x))
      {
      }

      virtual ~KaskadeLinearization() {}

      /// Number of columns of components [cbegin, cend)
      int cols(int cbegin=0, int cend=-1) const
      {
        assignIfNegative(cend,nColBlocks());
        assert(cbegin<cend);
        return x.descriptions.degreesOfFreedom(cbegin,cend);
      }

      /// Number of rows of components [rbegin, rend)
      int rows(int rbegin=0, int rend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assert(rbegin < rend);
        return x.descriptions.degreesOfFreedom(rbegin,rend);
      }

      void precompute() {
        doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);
      }

      /// write blocks of the hessian matrix into mat
      void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int cbegin, int cend) const
      {
        doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);
        mat = ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend,cbegin,cend);

        if(DiscreteBlockTraits<Functional>::anyPresent)
        {
          MatrixAsTriplet<Scalar> matD;
          getDiscreteMatrixBlocks(matD,rbegin,rend,cbegin,cend);
          mat+=matD;
        }
      }

      /// write components of the gradient into rhs
      void getRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const
      {
        doAssemble(Assembler::VALUE | Assembler::RHS);
        rhs.resize(rows(rbegin,rend),0.0);
        ass->toSequence(rbegin,rend,rhs.begin());
        if(DiscreteBlockTraits<Functional>::anyPresent)
        {
          addDiscreteRHSBlocks(rhs,rbegin,rend);
        }
      }


      /// return number of columns
      int nColBlocks() const { return Functional::AnsatzVars::noOfVariables; }

      /// return number of rows
      int nRowBlocks() const { return Functional::TestVars::noOfVariables; }

      /// return point of linearization
      AbstractFunctionSpaceElement const& getOrigin() const
      {
        return *xptr;
      }

      void setOrigin(AbstractFunctionSpaceElement const& x_)
      {
        x = Bridge::getImpl<DomainElement>(x_);
        xptr.reset(new Vector<DomainElement>(x));

        if(ass == nullptr)
        {
          ass.reset(new Assembler(x.descriptions.spaces));
        }
      }

      /// return the implementation
      Implementation const& getLinImpl() const {return lin; }

      /// flush all data, gathered so far
      void flush() { ass->flush( Assembler::VALUE | Assembler::RHS | Assembler::MATRIX ); }

      /// return whether x is in the domain of definition
      bool inDomain(DomainElement const& x) { return InDomainTraits<Functional,DomainElement>::inDomain(x); }

      /// return the current value Functional(Origin)
      double eval() const
      {
        doAssemble(Assembler::VALUE);
        Scalar value =ass->functional();
        if(DiscreteBlockTraits<Functional>::anyPresent)
          value += DiscreteBlockTraits<Functional>::getValue(*fu,x);
        return value;
      }

      double getValue() const { return eval(); }

      void evald(AbstractFunctionSpaceElement& v, int rbegin, int rend) const
      {
        assignIfNegative(rend,nRowBlocks());
        assert(rbegin<=rend);
        std::vector<Scalar> rhs(rows(rbegin,rend),0.0);
        getRHSBlocks(rhs,rbegin,rend);
        dynamic_cast<Vector<ImageElement>& >(v).read(rhs,rbegin,rend);
      }

      void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        d2axpy(1.0,y,x,rbegin,rend,cbegin,cend);
      }

      void d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        assert( cbegin <= cend );
        assert( rbegin <= rend );
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const&>(x).write(argument,cbegin,cend);
        dynamic_cast<Vector<ImageElement> const&>(y).write(result,rbegin,rend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.axpy(result, argument, a);
        dynamic_cast<Vector<ImageElement>&>(y).read(result,rbegin,rend);
      }

      void ddtxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        d2taxpy(1.0,y,x,rbegin,rend,cbegin,cend);
      }

      void d2taxpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());

        assert( cbegin <= cend );
        assert( rbegin <= rend );
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const&>(x).write(argument,rbegin,rend);
        dynamic_cast<Vector<ImageElement> const&>(y).write(result,cbegin,cend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.transpose().axpy(result, argument, a);
        dynamic_cast<Vector<ImageElement>&>(y).read(result,cbegin,cend);
      }

      double evalL1norm() const
      {
        doAssemble(Assembler::VALUE);
        Scalar value =0;//ass.fL1norm();
        if(DiscreteBlockTraits<Functional>::anyPresent)
          value += std::fabs(DiscreteBlockTraits<Functional>::getValue(*fu,x));
        return value;
      }

      Assembler const& getValidAssembler() const
      {
        doAssemble((Assembler::VALUE | Assembler::RHS | Assembler::MATRIX));
        return *ass;
      }

      Functional const& getFunctional() const
      {
        return *fu;
      }

    protected:
      void addDiscreteRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const
      {
        for(int i=rbegin; i<rend; ++i)
        {
          std::vector<Scalar> rhsD;
          DiscreteBlockTraits<Functional>::getRHSBlock(rhsD,*fu,x,i);
          if(rhsD.size() > 0)
          {
            assert(rhsD.size() == rows(i,i+1));
            for(int k=rows(rbegin,i), l=0;k<rows(rbegin,i+1);++k, ++l)
              rhs[k]+=rhsD[l];
          }
        }
      }

      void getDiscreteMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int cbegin, int cend) const
      {
        for(int i=rbegin; i<rend; ++i)
          for(int j=cbegin; j<cend; ++j)
          {
            MatrixAsTriplet<Scalar> matD;
            DiscreteBlockTraits<Functional>::getMatrixBlock(matD, *fu,x,i,j);
            matD.shiftIndices(rows(rbegin,i),cols(cbegin,j));
            mat+=matD;
          }
      }

      void doAssemble(int flags) const
      {
        int toDoFlag= ((~ass->valid()) & flags);
        if(toDoFlag!=0) ass->assemble(lin,toDoFlag,nThreads);
      }


      DomainElement x;
      Functional const* fu;
      Implementation lin;
      mutable std::shared_ptr<Assembler> ass;
      std::unique_ptr<Vector<DomainElement> > xptr;
    };

    template <class Functional> using ConnectedKaskadeLinearization = ConnectedLinearization<KaskadeLinearization<Functional> >;

    template<class Functional, int stateId=1, int adjointId=2>
    class NormalStepLinearization : public ConnectedKaskadeLinearization<Functional>
    {
      typedef ConnectedKaskadeLinearization<Functional> Base;
      typedef typename Base::Scalar Scalar;
      typedef typename Base::Assembler Assembler;
      typedef typename Base::DomainElement DomainElement;
      typedef typename Base::ImageElement ImageElement;
      using Base::nRowBlocks;
      using Base::nColBlocks;
    public:
      NormalStepLinearization() : Base() {}

      /// Creation of a linearization for a functional fu at x_
      NormalStepLinearization(Functional const& fu, DomainElement const& x) : Base(fu,x) {}

      /// Creation of a linearization for a functional fu at x_
      NormalStepLinearization(Functional const& fu, DomainElement const& x, std::shared_ptr<Assembler> const& assembler) : Base(fu,x,assembler) {}

      NormalStepLinearization(Base const& other) : Base(other) {}

      /// write blocks of the hessian matrix into mat
      void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if(rbegin==adjointId && rend==adjointId+1 && cbegin==adjointId && cend==adjointId+1) return;
        this->doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);
        if(rbegin==stateId && rend==stateId+1 && cbegin==stateId && cend==stateId+1)
        {
          mat = this->ass->template get<MatrixAsTriplet<Scalar> >(false,adjointId,adjointId+1,adjointId,adjointId+1);
          return;
        }

        if(rbegin<=stateId && rend>stateId && cbegin<=stateId && cend > stateId)
        {
          size_t rows0 = 0, cols0 = 0, rows1 = 0, cols1 = 0;
          MatrixAsTriplet<Scalar> tmp;
          if(rbegin < stateId)
          {
            tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,stateId,cbegin,cend);
            rows0 = tmp.nrows();
            mat += tmp;
          }
          if(cbegin < stateId)
          {
            tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,stateId,stateId+1,cbegin,stateId);
            cols0 = tmp.ncols();
            tmp.shiftIndices(rows0,0);
            mat += tmp;
          }
          tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,adjointId,adjointId+1,adjointId,adjointId+1);
          rows1 = tmp.nrows(), cols1 = tmp.ncols();
          tmp.shiftIndices(rows0,cols0);
          mat += tmp;
          if(cend > stateId+1)
          {
            tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,stateId,stateId+1,stateId+1,cend);
            tmp.shiftIndices(rows0,cols0+cols1);
            mat += tmp;
          }
          if(rend > stateId+1)
          {
            tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,stateId+1,rend,cbegin,cend);
            tmp.shiftIndices(rows0+rows1,0);
            mat += tmp;
          }

        }
        else
          mat = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend,cbegin,cend);


        if(DiscreteBlockTraits<Functional>::anyPresent)
        {
          MatrixAsTriplet<Scalar> matD;
          this->getDiscreteMatrixBlocks(matD,rbegin,rend,cbegin,cend);
          mat+=matD;
        }
      }

      void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if(rbegin==adjointId && rend==adjointId+1 && cbegin==adjointId && cend==adjointId+1) return;
        assert(cbegin<=cend);
        assert(rbegin<=rend);
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const& >(x).write(argument,cbegin,cend);
        dynamic_cast<Vector<ImageElement>& >(y).write(result,rbegin,rend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.axpy(result, argument);
        dynamic_cast<Vector<ImageElement>& >(y).read(result,rbegin,rend);
      }

      void d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if( rbegin == adjointId && rend == adjointId+1 && cbegin == adjointId && cend == adjointId+1 ) return;
        assert( cbegin <= cend );
        assert( rbegin <= rend );
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const&>(x).write(argument,cbegin,cend);
        dynamic_cast<Vector<ImageElement> const&>(y).write(result,rbegin,rend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.axpy(result, argument, a);
        dynamic_cast<Vector<ImageElement>&>(y).read(result, rbegin, rend);
      }
    };

    template<class Functional, int stateId=1, int adjointId=2>
    class TangentialStepLinearization : public ConnectedKaskadeLinearization<Functional>
    {
      typedef ConnectedKaskadeLinearization<Functional> Base;
      typedef typename Base::Scalar Scalar;
      typedef typename Base::Assembler Assembler;
      typedef typename Base::DomainElement DomainElement;
      typedef typename Base::ImageElement ImageElement;
      using Base::nRowBlocks;
      using Base::nColBlocks;
    public:
      TangentialStepLinearization() : Base() {}

      /// Creation of a linearization for a functional fu at x_
      TangentialStepLinearization(Functional const& fu, DomainElement const& x) : Base(fu,x) {}

      /// Creation of a linearization for a functional fu at x_
      TangentialStepLinearization(Functional const& fu, DomainElement const& x, std::shared_ptr<Assembler> const& assembler) : Base(fu,x,assembler) {}


      TangentialStepLinearization(Base const& other) : Base(other) {}

      /// write blocks of the hessian matrix into mat
      void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if(rbegin==adjointId && rend==adjointId+1 && cbegin==adjointId && cend==adjointId+1) return;
        this->doAssemble(Assembler::VALUE | Assembler::RHS | Assembler::MATRIX);

        if(rend > adjointId && cend > adjointId)
        {
          mat = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend-1,cbegin,cend-1);
          size_t rows = mat.nrows(), cols = mat.ncols();
          auto tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rend-1,rend,cbegin,cend-1);
          tmp.shiftIndices(rows,0);
          mat += tmp;
          tmp = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend-1,cend-1,cend);
          tmp.shiftIndices(0,cols);
          mat += tmp;
        }
        else
          mat = this->ass->template get<MatrixAsTriplet<Scalar> >(false,rbegin,rend,cbegin,cend);

        if(DiscreteBlockTraits<Functional>::anyPresent)
        {
          MatrixAsTriplet<Scalar> matD;
          this->getDiscreteMatrixBlocks(matD,rbegin,rend,cbegin,cend);
          mat+=matD;
        }
      }

      void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if(rbegin==adjointId && rend==adjointId+1 && cbegin==adjointId && cend==adjointId+1) return;
        assert(cbegin<=cend);
        assert(rbegin<=rend);
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const& >(x).write(argument,cbegin,cend);
        dynamic_cast<Vector<ImageElement>& >(y).write(result,rbegin,rend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.axpy(result, argument);
        dynamic_cast<Vector<ImageElement>& >(y).read(result,rbegin,rend);
      }

      void d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const
      {
        assignIfNegative(rend,nRowBlocks());
        assignIfNegative(cend,nColBlocks());
        if( rbegin == adjointId && rend == adjointId+1 && cbegin == adjointId && cend == adjointId+1 ) return;
        assert( cbegin <= cend );
        assert( rbegin <= rend );
        std::vector<Scalar> result, argument;
        dynamic_cast<Vector<ImageElement> const&>(x).write(argument,cbegin,cend);
        dynamic_cast<Vector<ImageElement> const&>(y).write(result,rbegin,rend);
        MatrixAsTriplet<Scalar> mat;
        getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
        mat.axpy(result, argument, a);
        dynamic_cast<Vector<ImageElement>&>(y).read(result,rbegin,rend);
      }
    };

    //template <class Functional, int stateId=1, int adjointId=2> using NormalStepLinearization = ConnectedLinearization<NormalStepKaskadeLinearization<Functional,stateId,adjointId> >;
    //template <class Functional, int stateId=1, int adjointId=2> using TangentialStepLinearization = ConnectedLinearization<TangentialStepKaskadeLinearization<Functional,stateId,adjointId> >;

  }
}

#endif /* BRIDGE_HH_ */
