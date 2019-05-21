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
/*
 * functor.hh
 *
 *  Created on: 01.10.2012
 *      Author: lars lubkoll
 */

#ifndef FUNCTOR_HH_
#define FUNCTOR_HH_

#include <memory>

// forward declaration
namespace Dune
{
  template <class,class> class BCRSMatrix;
  template <class,class> class LinearOperator;
}

namespace Kaskade
{
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * Functor * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /// Abstract functor implementation base class.
  template <typename ReturnType_, typename... Arguments>
  class FunctorImpl
  {
  public:
    typedef ReturnType_ ReturnType;

    virtual ~FunctorImpl(){}

    virtual ReturnType operator()(Arguments...) = 0;

    virtual FunctorImpl* clone() const = 0;
  };

  /// Implementation of FunctorImpl. Allows the use of function pointers as well as functors.
  template <class ParentFunctor, typename Fun, typename... Arguments>
  class FunctorHandler : public FunctorImpl<typename ParentFunctor::ReturnType, Arguments...>
  {
  public:
    typedef typename ParentFunctor::ReturnType ReturnType;

    FunctorHandler(const Fun& fun_) : fun(fun_){}

    FunctorHandler* clone() const { return new FunctorHandler(*this); }

    ReturnType operator()(Arguments... args) { return fun(args...); }

  private:
    Fun const& fun;
  };

  /// Functor with first class semantics.
  /**
   * Defined by return type and arguments of the callable entity.
   * See "Modern C++ Design", chapter 5, by Andrei Alexandrescu for details.
   */
  template <typename ReturnType_, typename... Arguments>
  class Functor
  {
    typedef FunctorImpl<ReturnType_,Arguments...> Impl;
  public:
    typedef ReturnType_ ReturnType;

    Functor() : impl(nullptr){}

    Functor(Functor const& functor) : impl(functor.impl->clone()){}

    Functor(const std::unique_ptr<Impl>& impl_) : impl(impl_->clone()){}

    template <class Fun>
    Functor(const Fun& fun) : impl(new FunctorHandler<Functor,Fun,Arguments...>(fun)){}

    Functor& operator=(const Functor& functor) { impl.reset(functor.impl->clone()); return *this; }

    ReturnType operator()(Arguments... args) { return (*impl)(args...); }

    ReturnType operator()(Arguments... args) const { return (*impl)(args...); }

  private:
    std::unique_ptr<Impl> impl;
  };


  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * ApplyOperator * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /// Store operator, independent of operator type.
  template <class Scalar, typename Domain, typename Range, class Assembler>
  class ApplyOperatorImpl
  {
  public:
    ApplyOperatorImpl(Assembler const& assembler_, bool onlyLowerTriangle_, int rbegin_, int rend_, int cbegin_, int cend_)
     : assembler(assembler_), onlyLowerTriangle(onlyLowerTriangle_), rbegin(rbegin_), rend(rend_), cbegin(cbegin_), cend(cend_)
    {}

    virtual ~ApplyOperatorImpl(){}

    virtual void apply(Domain const&, Range&) const = 0;

    virtual void applyscaleadd(Scalar,Domain const&,Range&) const = 0;

    virtual ApplyOperatorImpl* clone() const = 0;

    template <class Matrix>
    Matrix get() const { return assembler.template get<Matrix>(onlyLowerTriangle,rbegin,rend,cbegin,cend); }

    template <class Matrix>
    std::unique_ptr<Matrix> getPointer() const { return assembler.template get<Matrix>(onlyLowerTriangle,rbegin,rend,cbegin,cend); }

  private:
    Assembler const& assembler;
    bool onlyLowerTriangle;
    int rbegin, rend, cbegin, cend;
  };

  /// Implementation of FunctorImpl. Allows the use of function pointers as well as functors.
  template <class ParentApplyOperator, typename Operator>
  class ApplyOperatorHandler : public ApplyOperatorImpl<typename ParentApplyOperator::Scalar, typename ParentApplyOperator::Domain, typename ParentApplyOperator::Range, typename Operator::Assembler>
  {
  public:
    typedef typename ParentApplyOperator::Scalar Scalar;
    typedef typename ParentApplyOperator::Domain Domain;
    typedef typename ParentApplyOperator::Range  Range;
    typedef typename Operator::Assembler         Assembler;
    typedef ApplyOperatorImpl<Scalar,Domain,Range,Assembler> Base;

    ApplyOperatorHandler(const Operator& op_)
    : Base(op_.getAssembler(), op_.getOnlyLowerTriangle(), Operator::BlockInfo::firstRow, Operator::BlockInfo::lastRow, Operator::BlockInfo::firstCol, Operator::BlockInfo::lastCol),
      op_(op_)
    {}

    ApplyOperatorHandler* clone() const { return new ApplyOperatorHandler(*this); }

    void apply(Domain const& x, Range& y) const { op_.apply(x,y); }

    void applyscaleadd(Scalar alpha, Domain const& x, Range& y) const { op_.applyscaleadd(alpha,x,y); }

  private:
    Operator op_;
  };

  /// Operator
  template <typename Scalar_, typename Domain_, typename Range_,typename Assembler_>
  class ApplyOperator : public Dune::LinearOperator<Domain_,Range_>
  {
  public:
    typedef Scalar_     Scalar;
    typedef Domain_     Domain;
    typedef Range_      Range;
    typedef Assembler_  Assembler;

  private:
    typedef ApplyOperatorImpl<Scalar,Domain,Range,Assembler> Impl;

  public:
    ApplyOperator() : impl(nullptr){}

    ApplyOperator(ApplyOperator const& other) : impl(other.impl->clone()){}

    ApplyOperator(const std::unique_ptr<Impl>& impl_) : impl(impl_->clone()){}

    template <class Operator>
    ApplyOperator(const Operator& op) : impl(new ApplyOperatorHandler<ApplyOperator,Operator>(op)){}

    ApplyOperator& operator=(const ApplyOperator& other) { impl.reset(other.impl->clone()); }

    void apply(Domain const& x, Range& y) const { impl->apply(x,y); }

    void applyscaleadd(Scalar alpha, Domain const& x, Range& y) const { impl->applyscaleadd(alpha,x,y); }

    template <class Matrix>
    Matrix get() const { return impl->template get<Matrix>(); }

    template <class Matrix>
    std::unique_ptr<Matrix> getPointer() const { return impl->template getPointer<Matrix>(); }

  private:
    std::unique_ptr<Impl> impl;
  };
}

#endif /* FUNCTOR_HH_ */
