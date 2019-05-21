/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SYMMETRICOPERATORS_HH
#define SYMMETRICOPERATORS_HH

namespace Kaskade
{
  /**
   * \ingroup linalgbasic
   * \brief Abstract base class for dual pairing of \f$ X \f$ and its dual space \f$ X^* \f$.
   */
  template <class X, class Xstar>
  class DualPairing {
  public:
    typedef X domain_type;
    typedef typename X::field_type field_type;
    
    virtual field_type operator()(X const& x, Xstar const& y) const = 0;
  };
  
  /**
   * \ingroup linalgbasic
   * \brief Default implementation of dual pairing (relies on scalar product operator * being defined)
   */
  template <class X, class Xstar>
  class DefaultDualPairing: public DualPairing<X,Xstar> {
    typedef typename DualPairing<X,Xstar>::field_type field_type;
  public:
    virtual field_type operator()(X const& x, Xstar const& y) const { /* std::cerr << "default dual pairing = " << x*y << "\n"; */ return x*y; }
  };
  
  /**
   * \ingroup linalgbasic
   * \brief Zero dual pairing implementation (mainly useful whenever a dual pairing is needed for formal language reasons, but never used)
   */
  template <class X, class Xstar>
  class ZeroDualPairing: public DualPairing<X,Xstar> {
    typedef typename DualPairing<X,Xstar>::field_type field_type;
  public:
    virtual field_type operator()(X const& , Xstar const& ) const { return 0; }
  };
  
  //--------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic
   * \brief Interface for symmetric linear operators
   * 
   * In addition to the Dune::LinearOperator interface, symmetric operators provide access to a
   * dual pairing \f$ \langle\cdot,\cdot\rangle \f$ (w.r.t. which they are symmetric) and 
   * offer the simultaneous evaluation of \f$ Ax \f$ and \f$ \langle x,Ax\rangle \f$, which can be 
   * implemented more efficiently than the sequential evaluation of the two quantities. This
   * efficiency gain is exploited, e.g., in the conjugate gradient solver \ref Pcg.
   */
  template <class X, class Xstar>
  class SymmetricLinearOperator: public Dune::LinearOperator<X,Xstar>
  {
  public:
    typedef typename X::field_type field_type;
    
    
    /**
     * \brief operator application
     * Computes \f$ y \leftarrow Ax \f$ and returns the dual pairing \f$ \langle x,y \rangle \f$.
     * 
     * The default implementation just calls dp(x,y) after the operator application. When implementing this interface
     * in derived classes, consider a more efficient version.
     */
    virtual field_type applyDp(X const& x, Xstar& y) const
    {
      this->apply(x,y);
      return dp(x,y);
    }
    
    /**
     * \brief returns the dual pairing \f$ \langle x, y \rangle \f$ with respect to which the operator is symmetric
     */
    virtual field_type dp(X const& x, Xstar const& y) const = 0;
  };
    
  /**
   * \ingroup linalgbasic
   * \brief Wrapper class to present (hopefully) symmetric linear operators in the SymmetricLinearOperator interface.
   * 
   * Note that the "simultaneous" evaluation of \f$ Ax \f$ and \f$ \langle x,Ax\rangle \f$ is just the sequential one,
   * i.e. there is no performance gain.
   */
  template <class X, class Xstar>
  class SymmetricLinearOperatorWrapper: public SymmetricLinearOperator<X,Xstar>
  {
  public:
    typedef typename X::field_type field_type;
    
    /**
     * \brief Constructor.
     *      
     * \param A the (hopefully symmetric) linear operator 
     * \param dualPairing the dual pairing with respect to which the operator is symmetric
     * 
     * Both arguments have to exist during the lifetime of the SymmetricLinearOperatorWrapper object.
     */
    SymmetricLinearOperatorWrapper(Dune::LinearOperator<X,Xstar> const& A_, DualPairing<X,Xstar> const& dualPairing_)
    : A(A_), dualPairing(dualPairing_) {}
    
    virtual void apply(X const& x, Xstar& y) const
    {
      A.apply(x,y);
    }
    
    virtual void applyscaleadd(field_type alpha, X const& x, Xstar& y) const
    {
      A.applyscaleadd(alpha,x,y);
    }
    
    /**
     * \brief returns the dual pairing \f$ \langle x, y \rangle \f$ with respect to which the operator is symmetric
     */
    virtual field_type dp(X const& x, Xstar const& y) const
    {
      return dualPairing(x,y);
    }
    
  private:
    Dune::LinearOperator<X,Xstar> const& A;
    DualPairing<X,Xstar> const& dualPairing;
  };
    
  //--------------------------------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \brief Interface for symmetric preconditioners.
   * 
   * In addition to the Dune::preconditioner interface, symmetric preconditioners provide the simultaneous
   * evaluation of \f$ x = By \f$ and the energy product \f$ \langle By,y \rangle \f$, which can be implemented
   * more efficiently than the sequential evaluation of both quantities.
   */
  template <class X, class Xstar>
  class SymmetricPreconditioner: public Dune::Preconditioner<X,Xstar> 
  {
  public:
    typedef typename X::field_type field_type;
    
    /**
     * \brief Preconditioner preparation.
     * 
     * The provided default implementation does nothing.
     */
    virtual void pre(X& , Xstar& ) {}
    
    /**
     * \brief Preconditioner cleanup.
     * 
     * The provided default implementation does nothing.
     */
    virtual void post(X& x) {}
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ and returns \f$ \langle By, y \rangle \f$.
     */
    virtual field_type applyDp(X& x, Xstar const& y) = 0;
    
    /** 
     * \brief Returns true if the target vector x has to be initialized to zero before calling apply or applyDp
     */
    virtual bool requiresInitializedInput() const = 0;
  };
  
  /**
   * \ingroup linalgbasic
   * \brief Wrapper class presenting a symmetric preconditioner interface for any preconditioner.
   * 
   * Use this only if you know that the preconditioner is in fact symmetric.
   */
  template <class X, class Xstar>
  class SymmetricPreconditionerWrapper: public SymmetricPreconditioner<X,Xstar>
  {
  public:
    typedef typename X::field_type field_type;
    
    /**
     * \brief Constructor.
     *      
     * \param B the (hopefully symmetric) preconditioner
     * \param dualPairing the dual pairing with respect to which the operator is symmetric
     * 
     * Both arguments have to exist during the lifetime of the SymmetricLinearOperatorWrapper object.
     */
    SymmetricPreconditionerWrapper(Dune::Preconditioner<X,Xstar>& B_, DualPairing<X,Xstar> const& dualPairing_)
    : B(B_), dualPairing(dualPairing_) {}
    
    virtual void pre(X& x, Xstar& y)
    {
      B.pre(x,y);
    }
    
    virtual void post(X& x)
    {
      B.post(x);
    }
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ .
     * 
     * The default implementation evaluates both quantities sequentially, without any performance gain.
     */
    virtual void apply(X& x, Xstar const& y)
    {
      B.apply(x,y);
    }
    
     /**
     * \brief Computes \f$ x \leftarrow By \f$ and returns \f$ \langle By, y \rangle \f$.
     * 
     * The default implementation evaluates both quantities sequentially, without any performance gain.
     */
    virtual field_type applyDp(X& x, Xstar const& y)
    {
      B.apply(x,y); 
      return dualPairing(x,y);
    }
    
    /** 
     * \brief Returns true if the target vector x has to be initialized to zero before calling apply or applyDp
     * 
     * Does not assume anything about the preconditioner and hence returns true in any case.
     */
    virtual bool requiresInitializedInput() const { return true; }

  private:
    Dune::Preconditioner<X,Xstar>& B;
    DualPairing<X,Xstar> const& dualPairing;
  };
  
} // end of namespace Kaskade

#endif
