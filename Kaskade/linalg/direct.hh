/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DIRECT_HH
#define DIRECT_HH

#include <cmath>
#include <memory>
#include <type_traits>

#include <boost/timer/timer.hpp>

#include "dune/istl/operators.hh"
#include <dune/istl/solvers.hh>

#include "fem/istlinterface.hh"
#include "linalg/factorization.hh"

namespace Kaskade
{
  /// \internal
  namespace DirectSolver_Detail
  {
    class NoScalar{};
    class NoFieldType{};

    template <class T>
    typename T::Scalar hasScalar(typename T::Scalar*);
    template <class T>
    NoScalar hasScalar(...);

    template <class T>
    typename T::field_type hasFieldType(typename T::field_type*);
    template <class T>
    NoFieldType hasFieldType(...);

    template <class T>
    struct HasScalar
    {
      typedef decltype(hasScalar<T>(nullptr)) type;
      static constexpr bool value = !std::is_same<type,NoScalar>::value;
    };

    template <class T>
    struct HasFieldType
    {
      typedef decltype(hasFieldType<T>(nullptr)) type;
      static constexpr bool value = !std::is_same<type,NoFieldType>::value;
    };



    template <class T>
    struct ScalarType
    {
      typedef typename std::conditional< HasScalar<T>::value, typename HasScalar<T>::type,
                                           typename std::conditional< HasFieldType<T>::value, typename HasFieldType<T>::type, void>::type >::type type;

    };
  }
  /// \endinternal

  /**
   * \ingroup direct
   * \brief Dune::InverseOperator and Dune::Preconditioner interface for direct solvers.
   *
   * This keeps a factorization during the lifetime of the object, however, due to
   * shared data, efficient copying is possible.
   */
  template <class Domain_, class Range_>
  class DirectSolver: public Dune::InverseOperator<Domain_,Range_>,
		      public Dune::Preconditioner<Domain_,Range_>
  {
  public:
    typedef typename DirectSolver_Detail::ScalarType<Domain_>::type Scalar;
    typedef Domain_ Domain;
    typedef Range_ Range;

    /**
     * \brief Default constructor. 
     * 
     * A default constructed DirectSolver implements a (pretty useless) zero operator.
     */
    DirectSolver() {}

    /**
     * \brief Constructs a direct solver from an assembled linear operator. 
     * 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     */
    template <class AssembledGOP>
    explicit DirectSolver(AssembledGOP const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL):
        solver(getFactorization(directType,properties,A.template get<MatrixAsTriplet<typename AssembledGOP::Scalar> >()))
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Constructs a direct solver from a triplet matrix. 
     */
    template <class FieldType>
    explicit DirectSolver(MatrixAsTriplet<FieldType> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL):
        solver(getFactorization(directType,properties,A))
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Constructs a direct solver from a BCRS matrix. 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     */
    template <class FieldType>
    explicit DirectSolver(Dune::BCRSMatrix<Dune::FieldMatrix<FieldType,1,1>> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL):
        solver(getFactorization<FieldType,int>(directType,properties,A))
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Constructs a direct solver from a NumaBCRS matrix. 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     */
    template <class FieldType, int n, int m, class Index>
    explicit DirectSolver(NumaBCRSMatrix<Dune::FieldMatrix<FieldType,n,m>,Index> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL):
        solver(getFactorization<FieldType,n,m,Index,int>(directType,properties,A))
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * The provided right hand side \arg b is not modified (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) { apply(x,b,1e-10,res); }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an extension to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods.
     */
    void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) const { apply(x,b,1e-10,res); }

    /**
     * Solves the system for the given right hand side \arg b, which is
     * guaranteed not to be overwritten (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, double reduction, Dune::InverseOperatorResult& res)
    {
      const_cast<DirectSolver<Domain,Range> const*>(this)->apply(x,b,res);
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an addition to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods, since the factorization is not modified during the solution phase.
     */
    void apply(Domain& x, Range const& b, double reduction, Dune::InverseOperatorResult& res) const
    {
      boost::timer::cpu_timer timer;

      std::vector<Scalar> rhs(b.dim());
      vectorToSequence(b,begin(rhs));
      apply(rhs,reduction,res);
      vectorFromSequence(x,begin(rhs));
    }

    /**
     * \brief The input vector v contains the rhs. It will be overwritten by the solution.
     *
     * Factorization is not touched by solving. As an addition to the
     * Dune::InverseOperator interface, we provide const versions of
     * the apply methods.
     */
    void apply(std::vector<Scalar> &v, double reduction, Dune::InverseOperatorResult& res) const
    {
      boost::timer::cpu_timer timer;

      // Check for solver availability. Otherwise return zero solution.
      if (solver) 
      {
        #ifndef NDEBUG
        size_t rhsNan = 0, solNan = 0;
        size_t rhsInf = 0;
        for (auto const& x: v)
        {
          if (std::isnan(x)) ++rhsNan;
          if (std::isinf(x)) ++rhsInf;
        }
        if (rhsNan>0) std::cerr << __FILE__ << ':' << __LINE__ << ": " << rhsNan << " rhs entries are nan\n";
        if (rhsInf>0) std::cerr << __FILE__ << ':' << __LINE__ << ": " << rhsInf << " rhs entries are inf\n";
        #endif
        solver->solve(v);
        #ifndef NDEBUG
        for (int i=0; i<v.size(); ++i)
          if (std::isnan(v[i])) ++solNan;
          if (solNan>0) std::cerr << __FILE__ << ':' << __LINE__ << ": " << solNan << " solution entries are nan\n";
          #endif
      } else
        std::fill(v.begin(),v.end(),0.0);

      // Write solver statistics. Currently dummy.
      res.clear();
      res.iterations = 1;
      res.reduction = 1e-10; // dummy!
      res.converged = true;
      res.conv_rate = 1e-10;
      res.elapsed = (double)(timer.elapsed().user)/1e9;
    }

    virtual void  apply (std::vector<Scalar> &v)
    {
      Dune::InverseOperatorResult dummy_res;
      apply(v,1e-10,dummy_res);
    }

    virtual void pre (Domain &x, Range &b) {}

    virtual void apply (Domain &v, const Range &d)
    {
      Dune::InverseOperatorResult dummy_res;
      apply(v,d, 1e-10,dummy_res);
    }

    virtual void post (Domain &x) {}

  private:
    std::shared_ptr<Factorization<Scalar>> solver;
  };

  //---------------------------------------------------------------------
  
  // partial specialization for tiny fixed-size matrices
  template <class S, int n>
  class DirectSolver<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>: public Dune::InverseOperator<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>,
                                                                     public Dune::Preconditioner<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>
  {
  public:
    typedef S Scalar;
    typedef Dune::FieldVector<S,n> Domain;
    typedef Dune::FieldVector<S,n> Range;

    /**
     * \brief Default constructor. 
     * 
     * A default constructed DirectSolver implements a (pretty useless) zero operator.
     */
    DirectSolver() {}

    /**
     * \brief Constructs a direct solver from a FieldMatrix matrix. 
     */
    template <class FieldType>
    explicit DirectSolver(Dune::FieldMatrix<FieldType,n,n> const& A):
        inverse(A)
    {
      inverse.invert();
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * The provided right hand side \arg b is not modified (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) { apply(x,b,1e-10,res); }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an extension to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods.
     */
    void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) const { apply(x,b,1e-10,res); }

    /**
     * Solves the system for the given right hand side \arg b, which is
     * guaranteed not to be overwritten (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, double reduction, Dune::InverseOperatorResult& res)
    {
      const_cast<DirectSolver<Domain,Range> const*>(this)->apply(x,b,res);
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an addition to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods, since the factorization is not modified during the solution phase.
     */
    void apply(Domain& x, Range const& b, double reduction, Dune::InverseOperatorResult& res) const
    {
      x = inverse*b;
      
      // Write solver statistics. Currently dummy.
      res.clear();
      res.iterations = 1;
      res.reduction = 1e-10; // dummy!
      res.converged = true;
      res.conv_rate = 1e-10;
      res.elapsed = 0;
    }


    virtual void        pre (Domain &x, Range &b) {}

    virtual void        apply (Domain &v, const Range &d)
    {
      v = inverse*d;
    }

    virtual void        post (Domain &x) {}

  private:
    Dune::FieldMatrix<Scalar,n,n> inverse;
  };

  //---------------------------------------------------------------------

  /**
   * \ingroup linalgsolution
   * \brief Dune::LinearOperator interface for inverse operators.
   */
  template <class InverseOperator>
  class InverseLinearOperator: public Dune::LinearOperator<typename InverseOperator::Range, typename InverseOperator::Domain>
  {
  public:
    typedef typename InverseOperator::Domain Domain;
    typedef typename InverseOperator::Range Range;
    typedef typename InverseOperator::Scalar Scalar;

    InverseLinearOperator() = default;

    InverseLinearOperator(InverseOperator const& op_) : op(op_)
    {}

    virtual void apply(Domain const& x, Range& y) const
    {
      Dune::InverseOperatorResult result;
      Domain rhs(x);
      op.apply(y,rhs,result);
    }

    virtual void applyscaleadd(Scalar alpha, Domain const& x, Range& y) const
    {
      Range ynew(y);
      apply(x,ynew);
      y.axpy(alpha,ynew);
    }

  private:
    InverseOperator op;
  };

  //---------------------------------------------------------------------

  /**
   * \ingroup direct
   * \brief convenience function for constructing a DirectInverseOperator
   */
  template <class GOP, int firstRow, int lastRow, int firstCol, int lastCol>
  InverseLinearOperator<DirectSolver<typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Domain,
                                     typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Range> >
  directInverseOperator(AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol> const& A,
                        DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
  {
    typedef typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Domain Domain;
    typedef typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Range Range;
    return InverseLinearOperator<DirectSolver<Domain,Range> >(DirectSolver<Domain,Range>(A,directType,properties));
  }

  template <class Matrix, class Domain, class Range>
  InverseLinearOperator<DirectSolver<Domain,Range> >
  directInverseOperator(MatrixRepresentedOperator<Matrix,Domain,Range> const& A,
                        DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
  {
    return InverseLinearOperator<DirectSolver<Domain,Range> >(DirectSolver<Domain,Range>(A,directType,properties));
  }

} // namespace Kaskade

//---------------------------------------------------------------------

#endif
