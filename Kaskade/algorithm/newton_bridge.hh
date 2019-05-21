#ifndef NEWTON_BRIDGE_HH
#define NEWTON_BRIDGE_HH

/**
 * @file 
 * @brief Bridge classes that connect low level FE algorithms to higher level algorithms
 * @author Anton Schiela
 *
 * Usually these classes are implemented using the "pointer to
 * implementation (pimpl)" idiom. Newton's Method and similar
 * algorithms can be formulated in terms of abstract quantities by
 * including aInterface.hh. Moreover, since all base classes beginning
 * with "Abstract..." are not templates, the higher level algorithms
 * can be pre-compiled.
 */

#include <typeinfo>
#include <iostream>
#include <fstream>
#include <memory>
#include <numeric>
#include <utility>
#include <sys/wait.h>
#include <sys/types.h>

#include <boost/signals2.hpp>
#include <boost/bind.hpp>
#include <boost/fusion/include/transform.hpp>

#include "abstract_interface.hh"
#include "algorithm_base.hh"
#include "algorithm/dune_bridge.hh"
#include "linalg/linearsystem.hh"
#include "linalg/direct.hh"
//#include "linalg/mumps_solve.hh"
#include "linalg/umfpack_solve.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/istlinterface.hh"
#include "fem/norms.hh"
#include "io/iobase.hh"
#include "algorithm/kaskadeBridge.hh"

namespace Kaskade{

  /// Namespace of Classes that implement the abstract interface classes via the bridge-pattern
  namespace Bridge{

    template<class Implementation> Implementation& getImpl(AbstractFunctionSpaceElement& v);
    template<class Implementation> Implementation const& getImpl(AbstractFunctionSpaceElement const& v);

    /// Traits class to unify the assignment of std::vector to a Vector
    template<class VectorImpl>
    struct VectorTraits
    {
      static const int nComponents = 1;

      static void read(VectorImpl &out, std::vector<double> const& in, int vbegin, int vend)
      {
        out.resize(in.size());
        if(vend==-1)
          for(int i=0; i<in.size(); ++i) out[i]=in[i];
        else
          for(int i=vbegin; i<vend; ++i) out[i]=in[i];
      }

      static void write(VectorImpl const &in, std::vector<double> & out, int vbegin, int vend)
      {
        out.resize(in.size());
        if(vend==-1)
          for(int i=0; i<in.size(); ++i) out[i]=in[i];
        else
          for(int i=vbegin; i<vend; ++i) out[i]=in[i];
      }

      static void writeToFile(std::string const& name, bool append, VectorImpl const& x)
      {
//        std::_Ios_Openmode mode= std::ios::out;
//        if(append) mode= std::ios::app;
//        std::ofstream f(name.c_str(),mode);
//        f.setf(std::ios::scientific,std::ios::floatfield);
//        f.precision(16);
//        for(int i=0; i<x.size();++i)
//          f << x[i] << " ";
//        f << "\n";
      }

      static void print(VectorImpl const& out, std::string const& message)
      {
        std::cout << message << "[";
        for(int i=0; i<out.size();++i)
          std::cout << out[i] << " ";
        std::cout << "]" << std::endl;
      }

      static std::string getRole(const VectorImpl &out, int component) { return ""; }

      static void scale(VectorImpl &out, std::vector<double> const& lambda)
      {
        for(int i=0; i<out.size(); ++i) out[i]*=lambda[i];
      }
    };

    /// Specialization for variable sets.
    template<class Descr>
    struct VectorTraits<VariableSet<Descr> >
    {
      typedef VariableSet<Descr> VectorImpl;
      static const int nComponents = VectorImpl::Descriptions::noOfVariables;

      static void print(VectorImpl const& out, std::string const& message)
      {
        // Usually these vectors are too large for printing, do nothing!
      }


      static void read(VectorImpl &out, std::vector<double> const& in, int vbegin, int vend)
      {
        assert(in.size() >= out.descriptions.degreesOfFreedom(vbegin,vend));
        out.read(vbegin,vend,in.begin());
      }

      static void write(VectorImpl const& in, std::vector<double>& out, int vbegin, int vend)
      {
        out.resize(in.descriptions.degreesOfFreedom(vbegin,vend),0.0);
        in.write(vbegin,vend,out.begin());
      }

      static void writeToFile(std::string const& name, bool append, VectorImpl const& x, int order)
      {
//        std::cout << "write file of order " << order << std::endl;
//        boost::thread writeThread([&x,&name,&order]()
//        {
          writeVTKFile(x.descriptions.gridView,x,name,IoOptions(),order);
//        });
//        std::cout << "starting writeToFile" << std::endl;
//        //      std::cout << "Nothing written.." << std::endl;
//        pid_t pid1, pid2;
//        int status;
//
//        // Double fork, avoids zombie processes
//        if((pid1=fork())!=0)
//        {
//          std::cout << "waiting..." << std::flush;
//          waitpid(pid1, &status, 0);
//          std::cout << "done." << std::endl;
//        }
//        else
//        {
//          if((pid2=fork())!=0)
//          {
//            std::cout << "exit" << std::endl;
//            exit(0);
//          }
//          else
//          {
//            std::cout << "write before" << std::endl;
//            writeVTKFile(x.descriptions.gridView,x.descriptions,x,name,IoOptions(),order);
//            //writeAMIRAFile(x.descriptions.gridView,x.descriptions,x,name,IoOptions(),order);
//            std::cout << "exit" << std::endl;
//            exit(0);
//          }
//        }
      }


      struct ComponentWiseScalarMult {
        ComponentWiseScalarMult(std::vector<double>const& s_): s(s_), i(0) {}
        template <class Element> void operator()(Element& e) const { e *= s[i]; ++i; }
      private:
        std::vector<double>const& s;
        mutable int i;
      };

      static void scale(VectorImpl &out, std::vector<double> const& lambda)
      {
        boost::fusion::for_each(out.data, ComponentWiseScalarMult(lambda));
      }

      static std::string getRole(const VectorImpl &out, int component)
      {
        return out.descriptions.roles[component];
      }

    };


    /// Mathematical Vector that supports copy-on-write, implements AbstractFunctionSpaceElement
    /** An implementation must support the following three operations:
     *
     * - Implementation(Implementation const& )              (copy-constructor)
     * - operator*=(double lambda)                           (scaling operation: *this*=lambda, arbitrary return type)
     * - axpy(double alpha, Implementation const& t2)        (axpy:              *this+=alpha*t2, arbitrary return type)
     */
    template<typename Implementation>
    class Vector : public AbstractFunctionSpaceElement
    {
    public:
      template <class Creator>
      explicit Vector(Creator& creator) : implementation(creator) {}

      explicit Vector(Implementation const& gi) : implementation(gi) {}

      virtual std::unique_ptr<AbstractFunctionSpaceElement> clone() const
      {
        return std::unique_ptr<AbstractFunctionSpaceElement>(new Vector<Implementation>(*this));
      }

      virtual std::unique_ptr<AbstractFunctionSpaceElement> initZeroVector() const
      {
        auto zeroVector = clone();
        *zeroVector *= 0;
        return zeroVector;
      }

      /// Access to the data
      Implementation const & get() const { return implementation; }

      /// Access data
      Implementation& get() { return implementation; }

      virtual void read(std::vector<double> const& in, int vbegin, int vend)
      {
        VectorTraits<Implementation>::read(implementation,in,vbegin,vend);
      }

      virtual void read(std::vector<double> const& in)
      {
        VectorTraits<Implementation>::read(implementation,in,0,this->nComponents());
      }

      virtual void write(std::vector<double>& out, int vbegin, int vend) const
      {
        VectorTraits<Implementation>::write(implementation,out,vbegin,vend);
      }

      virtual void write(std::vector<double>& out) const
      {
        VectorTraits<Implementation>::write(implementation,out,0,this->nComponents());
      }

      virtual std::string getRole(int component) const
      {
        return VectorTraits<Implementation>::getRole(implementation,component);
      }

      int nComponents() const
      {
        return VectorTraits<Implementation>::nComponents;
      }

      virtual void print(std::string const& message="") const
      {
        VectorTraits<Implementation>::print(implementation,message);
      }

      virtual double doapplyAsDualTo(AbstractFunctionSpaceElement const& v, int vbegin, int vend) const
      {
        std::vector<double> v1,v2;
        VectorTraits<Implementation>::write(implementation,v1,vbegin,vend);
        VectorTraits<Implementation>::write(dynamic_cast<Vector<Implementation> const& >(v).get(),v2,vbegin,vend);
        //    if(v1.size() != v2.size()) std::cout << "Warning: different sizes: " << v1.size() << " " << v2.size() << std::endl;
        double sum=0.0;
        for(int i=0; i< std::min(v1.size(),v2.size()); ++i)
        {
          sum+=v1[i]*v2[i];
        }
        return sum;
      }


    private:
      ///Write access to the data. This causes copying
      Implementation& get_nonconst() { return implementation; }

      virtual Vector& doscale(std::vector<double>const& lambda)
      {
        VectorTraits<Implementation>::scale(implementation,lambda);
        return *this;
      }

      virtual Vector& doaxpy(double alpha, AbstractFunctionSpaceElement const& t2, int component) {
        std::vector<double> add1,add2;
        this->write(add1,component,component+1);
        VectorTraits<Implementation>::write(dynamic_cast<Vector<Implementation> const& >(t2).get(),add2,component,component+1);
        for(int i=0; i<add1.size();++i)
          add1[i]+=alpha*add2[i];
        this->read(add1,component,component+1);
        return *this;
      }

      virtual Vector& doaxpy_(double alpha, AbstractFunctionSpaceElement const& t2) {
        implementation.axpy(alpha,dynamic_cast<Vector<Implementation> const& >(t2).get());
        return *this;
      }

      virtual Vector& doassign(AbstractFunctionSpaceElement const& t2) {
        implementation = dynamic_cast<Vector<Implementation> const& >(t2).implementation;
        return *this;
      }

      virtual void doswap(AbstractFunctionSpaceElement& other)
      {
        std::swap(implementation,dynamic_cast<Vector<Implementation>& >(other).implementation);
      }

      //Hide assignment operator, rather use clone
      Vector<Implementation>& operator=(Vector<Implementation>const & v);

      // Private copy-constructor, which makes a shallow copy
      Vector(Vector<Implementation> const & v)
      : implementation(v.implementation) {}
      // Private copy-constructor, which makes a deep copy
      //  Vector(Vector<Implementation> const & v):
      //    implementation(new Implementation(*(v.implementation)))
      //  {}

      virtual void writeToFile(std::string const& file, bool append, int order=1) const { VectorTraits<Implementation>::writeToFile(file, append, implementation, order); };

    public:
      Implementation implementation;

      friend Implementation& getImpl<Implementation>(AbstractFunctionSpaceElement& v);
    };

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    /// Get the implementation of an AbstractFunctionSpaceElement
    template<class Implementation>
    Implementation& getImpl(AbstractFunctionSpaceElement& v)
    {
      return dynamic_cast<Vector<Implementation>& >(v).get_nonconst();

    }

    /// Get the implementation of an AbstractFunctionSpaceElement
    template<class Implementation>
    Implementation const& getImpl(AbstractFunctionSpaceElement const& v)
    {
      return dynamic_cast<Vector<Implementation>const& >(v).get();
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    /// Object that represents the linearization of a nonlinear functional, implements AbstractLineariztion
    /** Can be connected to an AdaptiveGrid to be able to know, when the grid changes
     * Emits a signal if it is invalidated. This can be used by a solver to flush data
     *
     * Currently there are two implementations available with an automatic switch
     * depending on the type of vector used. This switch is implemented via
     * LinearizationTraits:
     *
     * The first implementation is used in case of a Dune::BlockVector<T>. Then
     * the creation of a linear system is forwarded to the functional
     *
     * The second implementation is used in case of a VariableSet. Here a
     * VariationalFunctionalAssembler is used to create the linear system.
     *
     * Implementations have to provide the following routines
     * Example: DuneLinearization
     *
     * - RT getValue()
     *
     * - void getRHSBlocks(std::vector<RT>& rhs, int rbegin, int rend) const
     *
     * - void getMatrixBlocks(MatrixAsTriplet<RT>& mat, int rbegin, int rend, int colbegin, int colend) const
     *
     * - AbstractFunctionSpaceElement const& getOrigin()
     *
     * - int size()
     *
     * - void flush()
     *
     * - int nColBlocks();
     *
     * - int nRowBlocks();
     *
     * - int cols(int cbegin, int cend) const
     *
     * - int rows(int rbegin, int rend) const
     *
     * Linearizations can be defined to set the right hand side to 0 everywhere, except in [rhsstart,rhsfinal) this is useful on the optimization context
     */
//    template<typename Impl>
//    class Linearization : public AbstractLinearization, public SparseLinearSystem, public Impl::Operator
//    {
//      typedef typename Impl::Scalar Scalar;
//      typedef typename Impl::Operator Operator;
//    public:
//      typedef Impl Implementation;
//
//      Linearization() = delete;
//
//      Linearization(std::unique_ptr<Implementation>&& impl) : myImpl(std::move(impl)) {};
//
//      virtual int nRowBlocks() const {return myImpl->nRowBlocks(); }
//      virtual int nColBlocks() const {return myImpl->nColBlocks(); }
//
//      /// Implementation of AbstractLinearization
//      double eval() const { return myImpl->getValue();}
//
//      double evalL1norm() const { return myImpl->getL1norm(); }
//
//      double getValue() const { return eval();}
//
//      double getValue(AbstractFunctionSpaceElement const& dx) const { return eval(dx);}
//
//      void evald(AbstractFunctionSpaceElement& g, int rbegin=0, int rend=-1) const
//      {
//        if(rend==-1) rend=nRowBlocks();
//        assert(rbegin<=rend);
//        std::vector<Scalar> rhs(rows(rbegin,rend),0.0);
//        getRHSBlocks(rhs,rbegin,rend);
//        dynamic_cast<Vector<typename Implementation::ImageElement>& >(g).read(rhs,rbegin,rend);
//      }
//
//      void precompute()
//      {
//        myImpl->precompute();
//      }
//
//      void ddxpy(AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin=0, int rend=-1, int cbegin=0, int cend=-1) const {
//        if(rend==-1) rend=nRowBlocks();
//        if(cend==-1) cend=nColBlocks();
//        assert(cbegin<=cend);
//        assert(rbegin<=rend);
//        std::vector<Scalar> result, argument;
//        dynamic_cast<Vector<typename Implementation::ImageElement> const& >(x).write(argument,cbegin,cend);
//        dynamic_cast<Vector<typename Implementation::ImageElement>& >(y).write(result,rbegin,rend);
//        MatrixAsTriplet<Scalar> mat;
//        myImpl->getMatrixBlocks(mat,rbegin,rend,cbegin,cend);
//        mat.axpy(result, argument);
//        dynamic_cast<Vector<typename Implementation::ImageElement>& >(y).read(result,rbegin,rend);
//      }
//
//      AbstractFunctionSpaceElement const& getOrigin() const
//      {
//        optr.reset(new Vector<typename Implementation::DomainElement>(myImpl->getOriginImpl()));
//        return *optr;
//      }
//
//      // Implementation of Impl::OperatorType==Dune::LinearOperator
//      void applyscaleadd(typename Operator::domain_type::field_type scale, const typename Operator::domain_type& in,typename Operator::range_type& out) const
//      {
//        std::vector<Scalar> result(this->rows(0,this->nRowBlocks())), argument(this->cols(0,this->nColBlocks()));
//        in.write(argument.begin());
//        out.write(result.begin());
//        MatrixAsTriplet<Scalar> mat;
//        myImpl->getMatrixBlocks(mat,0,nRowBlocks(),0,nColBlocks());
//        mat.axpy(argument, result, scale);
//        out.read(result.begin());
//      }
//
//      void apply(const typename Operator::domain_type& in,typename Operator::range_type& out) const
//      {
//        std::vector<Scalar> result(this->rows(0,this->nRowBlocks())), argument(this->cols(0,this->nColBlocks()));
//        MatrixAsTriplet<Scalar> mat;
//        in.write(argument.begin());
//        myImpl->getMatrixBlocks(mat,0,nRowBlocks(),0,nColBlocks());
//        mat.ax(result, argument);
//        out.read(result.begin());
//      }
//
//      // Implementation of SparseLinearSystem
//      virtual void getMatrixBlocks(MatrixAsTriplet<Scalar>& mat, int rbegin, int rend, int colbegin, int colend) const
//      { myImpl->getMatrixBlocks(mat,rbegin,rend,colbegin,colend); }
//
//      virtual void getRHSBlocks(std::vector<Scalar>& rhs, int rbegin, int rend) const
//      {
//        myImpl->getRHSBlocks(rhs,rbegin,rend);
//      }
//
//      virtual int cols(int cbegin, int cend) const {return myImpl->cols(cbegin, cend);}
//      virtual int rows(int rbegin, int rend) const {return myImpl->rows(rbegin, rend);}
//      // Other..
//
//      typename Implementation::DomainElement const&
//      getOriginImpl() const { return myImpl->getOriginImpl(); }
//
//      void connectToSignalForFlush(boost::signals2::signal0<void>& sig) {
//        if(flushconn.connected()) flushconn.disconnect();
//        flushconn=sig.connect(boost::bind(&Linearization<Implementation>::onGridChange, this));
//      }
//
//      virtual ~Linearization() { changed(); flushconn.disconnect(); }
//
//
//      void onGridChange(){ changed(); myImpl->flush();}
//
//      void flush(){ changed(); myImpl->flush();}
//
//      void touch() { changed(); myImpl->touch(); }
//
//      Implementation const& getLinImpl() const {return *myImpl;}
//      Implementation* getLinImplPtr() {return myImpl.get();}
//
//      typename Implementation::Assembler const& getValidAssembler() const
//      {
//        return myImpl->getValidAssembler();
//      }
//
//    private:
//      std::unique_ptr<Implementation> myImpl;
//      mutable boost::signals2::connection flushconn;
//      mutable std::unique_ptr<Vector<typename Implementation::DomainElement> > optr;
//    };


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    /// Linear Solver. Can get a signal from Linearization in order to flush its data, implements AbstractLinearSolver
    /** Implementation must support the following member functions:
     *
     * - getRelativeAccuracy() : relative Accuracy of the Solution
     *
     * - solve(std::vector<std::vector<double> >& c, SparseLinearSystem& lin) : solve a linear system
     *
     * - resolve(std::vector<std::vector<double> >& c, SparseLinearSystem& lin) : solve a linear system, using an old factorization
     *
     * - onChangedLinearization() : do something, if the Linearization has changed, for example flush memory
     *
     * There is one pure virtual function to be overloaded
     *
     */
//    class NoParameters {};
//
//    template<typename Implementation, typename VectorImpl, class Parameters=NoParameters>
//    class FixedGridNewtonDirection : public AbstractNewtonDirection
//    {
//      //  typedef FixedGridNewtonDirection<Implementation,VectorImpl,Parameters>  Self;
//    public:
//      FixedGridNewtonDirection(std::shared_ptr<Implementation> imp)
//      : myImplementation(imp), p(0)
//      {
//      }
//
//      FixedGridNewtonDirection(std::shared_ptr<Implementation> imp, Parameters const& p_)
//      : myImplementation(imp), p(&p_) {
//      }
//
//      ~FixedGridNewtonDirection()
//      {
//        if(flushconn.connected()) flushconn.disconnect();
//      }
//
//      /// Solves always exactly
//      virtual void setRelativeAccuracy(double accuracy) {myImplementation->setRelativeAccuracy(accuracy);}
//
//      /// Always exact solution
//      virtual double getRelativeAccuracy() {return myImplementation->getRelativeAccuracy();}
//
//      /// Always exact solution
//      virtual double getAbsoluteAccuracy() {return myImplementation->getAbsoluteAccuracy();}
//
//      virtual bool improvementPossible() { return myImplementation->improvementPossible(); }
//
//    private:
//
//      template<class Impl, class Paras>
//      struct CallTraits
//      {
//        static void solve(Impl& impl, std::vector<double>& res, SparseLinearSystem& lin, Paras const* p)
//        {
//          impl.resetParameters(*p);
//          impl.solve(res,lin);
//        }
//
//        static void resolve(Impl& impl, std::vector<double> &res, SparseLinearSystem const& lin, Paras const* p)
//        {
//          impl.resetParameters(*p);
//          impl.resolve(res,lin);
//        }
//      };
//
//      template<class Impl>
//      struct CallTraits<Impl, NoParameters>
//      {
//        static void solve(Impl& impl, std::vector<double>& res, SparseLinearSystem& lin, NoParameters const* )
//        {
//          impl.solve(res,lin);
//        }
//
//        static void resolve(Impl& impl, std::vector<double>& res, SparseLinearSystem const& lin, NoParameters const* )
//        {
//          impl.resolve(res,lin);
//        }
//      };
//
//      void connectToSignalForFlush(boost::signals2::signal0<void>& sig)
//      {
//        if(flushconn.connected()) flushconn.disconnect();
//        flushconn=sig.connect(boost::bind(&FixedGridNewtonDirection::onGridChange, this));
//      }
//
//      virtual void doSolve(AbstractFunctionSpaceElement& correction, AbstractLinearization& lin)
//      {
//        if(Implementation::needMatrix)
//        {
//          connectToSignalForFlush(lin.changed);
//          std::vector<double> result;
//          CallTraits<Implementation,Parameters>::solve(*myImplementation,result,dynamic_cast<SparseLinearSystem&>(lin),p);
//          dynamic_cast<Vector<VectorImpl>& >(correction).read(result);
//        } else
//        {
//          std::cout << "Not implemented!" << std::endl;
//          //      myImplementation->solve(dynamic_cast<VectorImpl& >correction,dynamic_cast<VectorImpl& >lin,dynamic_cast<VectorImpl& >lin)
//        }
//      }
//
//      virtual void doResolve(AbstractFunctionSpaceElement& correction,
//          AbstractLinearization const& lin,
//          AbstractLinearization const& olin) const
//      {
//        if(Implementation::needMatrix)
//        {
//          std::vector<double>  result;
//          CallTraits<Implementation,Parameters>::resolve(*myImplementation,result,dynamic_cast<SparseLinearSystem const &>(lin),p);
//          dynamic_cast<Vector<VectorImpl>& >(correction).read(result);
//        } else
//        {
//          std::cout << "Not implemented!" << std::endl;
//          //      myImplementation->solve(correction_,lin_,olin_)
//        }
//      }
//
//      virtual void doResolve(AbstractFunctionSpaceElement& correction,
//          AbstractLinearization const& lin) const
//      {
//        if(Implementation::needMatrix)
//        {
//          std::vector<double>  result;
//          CallTraits<Implementation,Parameters>::resolve(*myImplementation,result,dynamic_cast<SparseLinearSystem const&>(lin),p);
//          dynamic_cast<Vector<VectorImpl>& >(correction).read(result);
//        } else
//        {
//          std::cout << "Not implemented!" << std::endl;
//          //      myImplementation->solve(correction_,lin_,lin_)
//        }
//      }
//
//      void onGridChange(){ myImplementation->onChangedLinearization(); }
//
//      std::shared_ptr<Implementation> myImplementation;
//      boost::signals2::connection flushconn;
//      Parameters const* p;
//    };

    /// Bridge class for Functionals. Its most prominent task is to create linearizations, implements AbstractFunctional
    template<typename FunctionalImpl, typename DomainImpl=typename FunctionalImpl::AnsatzVars::VariableSet, typename ImageImpl=DomainImpl>
    class Functional : public AbstractFunctional
    {
    public:
      typedef typename FunctionalImpl::Scalar Scalar;


      explicit Functional(FunctionalImpl* impl) : myImplementation(impl){}

      explicit Functional(std::unique_ptr<FunctionalImpl>&& imp) : myImplementation(imp.release()) {}

      template <typename... ConstructorArguments>
      explicit Functional(const ConstructorArguments&... args) : myImplementation( new FunctionalImpl(args...) ) {}

      std::unique_ptr<AbstractLinearization> getLinearization(AbstractFunctionSpaceElement const& x) const
                  {
                    return std::unique_ptr<AbstractLinearization>(new ConnectedKaskadeLinearization<FunctionalImpl>(*myImplementation,getImpl<DomainImpl>(x)));
        //return std::unique_ptr<AbstractLinearization>(new Linearization<LinImpl>(std::unique_ptr<LinImpl>(new LinImpl(*myImplementation,getImpl<DomainImpl>(x)))));
                  }

      virtual std::unique_ptr<AbstractFunctionSpaceElement> getImageVector(AbstractFunctionSpaceElement const& x) const
                  {
        return std::unique_ptr<AbstractFunctionSpaceElement>(new Vector<ImageImpl>(getImpl<ImageImpl>(x)));
        //return std::unique_ptr<AbstractFunctionSpaceElement>(new Vector<ImageImpl>(std::unique_ptr<ImageImpl>(new ImageImpl(getImpl<DomainImpl>(x)))));
                  }

      virtual bool inDomain(AbstractFunctionSpaceElement const& x) const
      {
        return myImplementation->inDomain(getImpl<DomainImpl>(x));
      }
    protected:
      std::unique_ptr<FunctionalImpl> myImplementation;
    };

    //template <class,int,int> class NormalStepLinearization;
    //template <class,int,int> class TangentialStepLinearization;

     /// Bridge class for Functionals. Its most prominent task is to create linearizations, implements AbstractFunctional
    template<typename FunctionalImpl, typename DomainImpl=typename FunctionalImpl::AnsatzVars::VariableSet, typename ImageImpl=DomainImpl>
    class KaskadeNormalStepFunctional : public Functional<FunctionalImpl,DomainImpl,ImageImpl>
    {
      typedef typename FunctionalImpl::Scalar Scalar;
      typedef Functional<FunctionalImpl,DomainImpl,ImageImpl> Base;
      typedef VariationalFunctionalAssembler<LinearizationAt<FunctionalImpl> > Assembler;
    public:
      explicit KaskadeNormalStepFunctional(std::unique_ptr<FunctionalImpl>&& imp) : Base(std::move(imp)) {}

      explicit KaskadeNormalStepFunctional(FunctionalImpl* impl) : Base(impl){}

      template <typename... ConstructorArguments>
      explicit KaskadeNormalStepFunctional(const ConstructorArguments&... args) : Base(args...) {}

      std::unique_ptr<AbstractLinearization> getLinearization(AbstractFunctionSpaceElement const& x) const
      {
        return std::unique_ptr<AbstractLinearization>(new NormalStepLinearization<FunctionalImpl>(*(this->myImplementation),getImpl<DomainImpl>(x),assembler));
      }

      void setAssembler(std::shared_ptr<Assembler>& assembler_)
      {
        assembler = assembler_;
      }

    private:
      std::shared_ptr<Assembler> assembler;
    };

    /// Bridge class for Functionals. Its most prominent task is to create linearizations, implements AbstractFunctional
    template<typename FunctionalImpl, typename DomainImpl=typename FunctionalImpl::AnsatzVars::VariableSet, typename ImageImpl=DomainImpl>
    class KaskadeTangentialStepFunctional : public Functional<FunctionalImpl,DomainImpl,ImageImpl>
    {
      typedef Functional<FunctionalImpl,DomainImpl,ImageImpl> Base;
      typedef VariationalFunctionalAssembler<LinearizationAt<FunctionalImpl> > Assembler;
    public:
      explicit KaskadeTangentialStepFunctional(std::unique_ptr<FunctionalImpl>&& imp) : Base(std::move(imp)) {}

      explicit KaskadeTangentialStepFunctional(FunctionalImpl* impl) : Base(impl){}

      template <typename... ConstructorArguments>
      explicit KaskadeTangentialStepFunctional(const ConstructorArguments&... args) : Base(args...) {}

      std::unique_ptr<AbstractLinearization> getLinearization(AbstractFunctionSpaceElement const& x) const
      {
        return std::unique_ptr<AbstractLinearization>(new TangentialStepLinearization<FunctionalImpl>(*(this->myImplementation),getImpl<DomainImpl>(x),assembler));
      }

      void setAssembler(std::shared_ptr<Assembler> const& assembler_)
      {
        assembler = assembler_;
      }

    private:
      std::shared_ptr<Assembler> assembler;
    };

    template<class Variables, class Functional>
    std::unique_ptr<AbstractFunctional> getFunctional(std::unique_ptr<Functional>&& F)
    {
      return std::unique_ptr<AbstractFunctional>(new Bridge::Functional<Functional,Variables>(std::move(F)));
    }

    template<class Variables, class Functional>
    std::unique_ptr<AbstractFunctional> getFunctional(Functional*&& F)
    {
      return std::unique_ptr<AbstractFunctional>(new Bridge::Functional<Functional,Variables>(std::unique_ptr<Functional>(F)));
    }


    /// A functional that may depend on parameters, implements AbstractParameterFunctional
    template<typename Func, typename VectorImpl>
    class ParameterFunctional : public AbstractParameterFunctional
    {
      typedef typename Func::OptimalityFunctional VarFu;
      typedef typename Func::Parameter Parameter;
      typedef typename Func::ParameterLinearization FuncPLin;
      typedef typename Func::ParameterValueLinearization FuncPVLin;

      template<class VariFu>
      struct ParaTraits
      {
        static void setParameter(Parameter , VariFu&) { std::cout << "Pns"; };
      };

      template<template<class Scalar, class T, class Obj, int nEFields, class S> class V, class T, class Obj, int nEFields>
      struct ParaTraits<V<typename VarFu::Scalar, T, Obj, nEFields, Parameter> >
      {
        static void setParameter(Parameter mu, VarFu& f) { f.setParameter(mu); };
      };
    public:
      ParameterFunctional(VarFu& F_) : F(F_) {};

      virtual std::unique_ptr<AbstractFunctional> getFunctional(AbstractParameters const& p) const
                  {
        ParaTraits<VarFu>::setParameter(dynamic_cast<Parameters<Parameter> const&>(p).getPars(),F);
        return std::unique_ptr<AbstractFunctional> (new Functional<Func,VectorImpl>(makeUP(new Func(dynamic_cast<Parameters<typename Func::Parameter> const&>(p).getPars(),F))));
                  }
    private:
      VarFu& F;
    };

    /// A functional that may depend on parameters, implements AbstractC1ParameterFunctional
    template<typename Func, typename VectorImpl>
    class C1ParameterFunctional : public AbstractParameterFunctional
    {
      typedef typename Func::OptimalityFunctional VarFu;
      typedef typename Func::Parameter Parameter;
      typedef typename Func::ParameterLinearization FuncPLin;
      typedef typename Func::ParameterValueLinearization FuncPVLin;

      template<class VariFu>
      struct ParaTraits
      {
        static void setParameter(Parameter , VariFu&) { std::cout << "Pns"; };
      };

      template<template<class Scalar, class T, class Obj, int nEFields, class S> class V, class T, class Obj, int nEFields>
      struct ParaTraits<V<typename VarFu::Scalar, T, Obj, nEFields, Parameter> >
      {
        static void setParameter(Parameter mu, VarFu& f) { f.setParameter(mu); };
      };
    public:
      C1ParameterFunctional(VarFu& F_) : F(F_) {};

      virtual std::unique_ptr<AbstractFunctional> getFunctional(AbstractParameters const& p) const
                  {
        ParaTraits<VarFu>::setParameter(dynamic_cast<Parameters<Parameter> const&>(p).getPars(),F);
        return std::unique_ptr<AbstractFunctional>(new Functional<Func,VectorImpl>(makeAP(new Func(dynamic_cast<Parameters<typename Func::Parameter> const&>(p).getPars(),F))));
                  }

      virtual std::unique_ptr<AbstractFunctional> getParameterLinFunctional(AbstractParameters const& p) const
                  {
        ParaTraits<VarFu>::setParameter(dynamic_cast<Parameters<Parameter> const&>(p).getPars(),F);
        return std::unique_ptr<AbstractFunctional>(new Functional<FuncPLin,VectorImpl>(makeAP(new FuncPLin(dynamic_cast<Parameters<typename FuncPLin::Parameter> const&>(p).getPars(),F))));
                  }

      virtual std::unique_ptr<AbstractFunctional> getLinFunctionValue(AbstractParameters const& p) const
                  {
        //     typename VarFu::FunctionalDiagonal FD;
        //     ParaTraits<typename VarFu::FunctionalDiagonal>::setParameter(dynamic_cast<Parameters<Parameter> const&>(p).getPars(),FD);
        //     std::unique_ptr<AbstractFunctional> ptr
        //       (new Functional<FuncPVLin,VectorImpl>
        //        (makeAP(new FuncPVLin(dynamic_cast<Parameters<typename FuncPVLin::Parameter> const&>(p).getPars(),FD))));
        ParaTraits<typename VarFu::Functional>::setParameter(dynamic_cast<Parameters<Parameter> const&>(p).getPars(),F);
        return std::unique_ptr<AbstractFunctional>(new Functional<FuncPVLin,VectorImpl>(makeUP(new FuncPVLin(dynamic_cast<Parameters<typename FuncPVLin::Parameter> const&>(p).getPars(),F))));
                  }
    private:
      VarFu& F;
    };


    ///Convenience routine: makes an unique_ptr of the right type
    template<class T>
    inline std::unique_ptr<T> makeUP(T* t) { return std::unique_ptr<T>(t); };


    struct SpaceTransfer
    {
      template <class Fu1, class Fu2>
      void operator()(Fu1& f1, Fu2 const& f2)
      {
        spaceTransfer(f1,f2);
      }
    };
  }
} 
#endif
