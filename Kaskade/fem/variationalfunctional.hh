/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file
 * @brief  Documentation for a Variational Functional
 * @author Martin Weiser, Anton Schiela
 */
class unspecified;

/**
 * \brief Documentation of the concept of a quadratic variational functional
 * \ingroup concepts
 * The variational functional concept defines the interface that is accessed by the VariationalFunctionalAssembler. 
 * 
 * This is *not* a base class, only a documentation!
 *
 * The mathematical concept represented by this interface is that of a
 * variational functional
 * \f[ J(u) = \int_\Omega f(u_1(x),\nabla u_1(x),\dots,u_n(x),\nabla u_n(x)) \, dx + \int_{\partial\Omega} g(u_1(x),\dots,u_n(x)) \, dx \f]
 * with quadratic functions \f$f\f$ and \f$g\f$ to be evaluated at \f$0\f$.
 * The different variables \f$u_i\f$ may be scalar or vector-valued, and may live in different finite element
 * function spaces.
 *
 * Some local classes have to be written into VariationalFunctional,
 * where most of the functionality of the functional is defined. These
 * are:
 *
 * - DomainCache   (for evaluation of \f$f\f$ and its derivatives)
 * - BoundaryCache (for evaluation of \f$g\f$ and its derivatives)
 * - D1 (information on the block structure of \f$ f' \f$ and \f$g'\f$)
 * - D2 (information on the block structure of \f$f''\f$ and \f$g'\f$)
 *
 * Also read the documentation of LinearizationAt in
 * functional_aux.hh, an adapter class for the usage and
 * the definition of a variational functional. If only quadratic
 * functionals (linear problems) are to be defined, these classes are
 * not necessary.
 *
 * Const methods are required to be thread-safe.
 * Note that the space lists AnsatzVars::Spaces and TestVars::Spaces must coincide!
 */
class VariationalFunctional 
{
public:
  /**
   * \brief The scalar type to be used for this variational problem.
   */
  typedef unspecified Scalar;

  /**
   * \brief A description of the ansatz variables occuring in the
   * variational problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified AnsatzVars;


  /**
   * A description of the test variables occuring in the "variational"
   * problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified TestVars;

  /**
   * \brief A description of the variables defining the linearization point.
   * 
   * In most cases, this is simply AnsatzVars.
   * 
   * In some cases, however, the linear(ized) problem to be solved is posed in a different 
   * ansatz space than the current evaluation point, e.g. hierarchical error 
   * estimators computing corrections in an extension space. Then there is a 
   * difference between ansatz and test variables on one hand and the origin 
   * variables on the other. 
   * 
   * This should be an instance of VariableSetDescription<...>.
   */
  typedef unspecified OriginVars;

  /**
   * \brief The type of problem, either VariationalFunctional or
   * WeakFormulation.
   */
  static ProblemType const type;

  /**
   * The following typedef should be useful, when accessing Entities,
   * it is however not strictly required for a functional to work
   */
  typedef typename Vars::Grid::template Codim<0>::Entity Entity;

  /**
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent w.r.t. parallel
   * execution.
   */
  struct DomainCache : public EvalCacheBase    
  {

    /**
     * Construct all data that is constant in an entity. Default
     * implementation in EvalCacheBase: does nothing
     */
    template<class Entity>
    void moveTo(Entity const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default
     * implementation in EvalCacheBase: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$f(0)\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$f\f$ at \f$ 0 \f$ with respect to
     * the variable \f$ u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of \f$f\f$ at \$ 0 \f$ with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by \arg argT and \arg argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent
   * w.r.t. parallel execution. 
   */
  struct BoundaryCache : public EvalCacheBase    
  {
    /**
     *  \brief Switch for the assembler:
     *  - false: BoundaryCache is only evaluated at the domain boundary
     *  - true: BoundaryCache is evaluated also at interior faces, e.g. for computing face-gradient jumps for error estimation.
     */
    static const bool hasInteriorFaces=false;

/// Type of the Face Information
    typedef typename AnsatzVars::Grid::template Codim<0>::Entity::LeafIntersectionIterator FaceIterator;

    /**
     * Construct all data that is constant in an entity. Default implementation in EvalBaseCache: does nothing
     */
    template<class FaceIterator>
    void moveTo(FaceIterator const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default implementation in EvalBaseCache: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$g(0)\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$g\f$ at \f$ 0 \f$,
     * with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of 
     * \f$ g \f$ at \f$ 0 \f$, with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$ u_{\mathrm{col}}\f$ in direction of the test functions given
     * by argT and the ansatz functions, given by argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };


  /**
   * \brief Provides static information about the right hand side blocks.
   * 
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row block of
   * the first derivative.
   */
  template <int row>
  struct D1
  {
    /**
     * \brief Provides static information about the subvector blocks of the right hand side.
     * 
     * Is true if that block is statically present, i.e. if \f$ f \f$
     * depends on variable \f$ u_{\mathrm{row}} \f$.
     */
    static bool const present;
  };

  /**
   * \brief Provides static information about the submatrix blocks of the stiffness block matrix.
   * 
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the second derivative.
   */
  template <int row, int col>
  struct D2 
  {
    /**
     * \brief Specifies the presence of a submatrix block.
     * 
     * Is true if that block is statically present, i.e. if variable
     * \f$ u_{\mathrm{row}} \f$ and \f$ u_{\mathrm{col}}\f$ are
     * nonlinearly coupled.
     *
     * For second derivatives, which are assumed to be always
     * symmetric, only the lower triangular blocks need be present.
     * 
     * Announcing blocks as not present allows the assembler to save both memory and cput time
     * as such blocks are completely ignored.
     */
    static bool const present;

    /**
     * \brief Specifies whether the subblock is symmetric.
     * 
     * Should true if the block is conceptually symmetric, i.e. if
     * EvaluationCache::d2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     * 
     * Announcing a block to be symmetric allows the assembler to save both memory and cpu time 
     * by accessing only the lower half of the matrix block.
     */
    static bool const symmetric;

    /**
     * \brief Specifies whether only the diagonal of the subblock shall be assembled.
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
  };
  
  /**
   * \brief Creates a DomainCache.
   * 
   * \param[in] flags a bit field describing what values will be accessed. See \ref VariationalAssembler enums.
   */
  DomainCache createDomainCache(int flags=7) const {}
  /**
   * \brief Creates a BoundaryCache.
   */
  BoundaryCache createBoundaryCache(int flags=7) const {}


  /**
   * \brief This method defines a suitable quadrature order to be used on the
   * given cell. 
   * 
   * The maximum shape function order occuring in the
   * finite element spaces for this cell is provided for
   * convenience. If boundary is true, this refers to evaluation of
   * \f$g\f$, otherwise to the evaluation of \f$f\f$.
   */
  int integrationOrder(typename Variables::Grid::template Codim<0>::Entity const& cell,
                       int shapeFunctionOrder, bool boundary) const;
};




/// Documentation of the concept of a nonlinear variational functional
/**
 * \ingroup concepts
 * The nonlinear variational functional concept defines the interface that is accessed by the LinearizationAt
 * and linearizationAt adapters. It differs from the basic VariationalFunctional concept mainly in how the
 * DomainCache and the BoundaryCache are constructed.
 * 
 * This is *not* a base class, only a documentation!
 *
 * The mathematical concept represented by this interface is that of a
 * variational functional
 * \f[ J(u) = \int_\Omega f(u_1(x),\nabla u_1(x),\dots,u_n(x),\nabla u_n(x)) \, dx + \int_{\partial\Omega} g(u_1(x),\dots,u_n(x)) \, dx \f]
 * to be evaluated for a specific fixed value of \f$u\f$.
 * The different variables \f$u_i\f$ may be scalar or vector-valued, and may live in different finite element
 * function spaces.
 *
 * Some local classes have to be written into VariationalFunctional,
 * where most of the functionality of the functional is defined. These
 * are:
 *
 * - DomainCache   (for evaluation of \f$f\f$)
 * - BoundaryCache (for evaluation of \f$g\f$)
 * - D1 (information on the block structure of \f$ f' \f$)
 * - D2 (information on the block structure of \f$f''\f$)
 *
 * Also read the documentation of LinearizationAt in
 * functional_aux.hh, which provides helper classes for the usage and
 * the definition of a variational functional. If only quadratic
 * functionals (linear problems) are to be defined, these classes are
 * not necessary.
 *
 * Const methods are required to be thread-safe.
 * Note that the space lists AnsatzVars::Spaces and TestVars::Spaces must coincide!
 */
class NonlinearVariationalFunctional 
{
public:
  /**
   * \brief The scalar type to be used for this variational problem.
   */
  typedef unspecified Scalar;

  /**
   * \brief A description of the ansatz variables occuring in the
   * variational problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified AnsatzVars;


  /**
   * A description of the test variables occuring in the "variational"
   * problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified TestVars;
  
  /**
   * \brief A description of the type of the point of linearization
   */
  typedef unspecified OriginVars;

  /**
   * \brief The type of problem, either VariationalFunctional or
   * WeakFormulation.
   */
  static ProblemType const type;

  /**
   * The following typedef should be useful, when accessing Entities,
   * it is however not strictly required for a functional to work
   */
  typedef typename Vars::Grid::template Codim<0>::Entity Entity;

  /**
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent w.r.t. parallel
   * execution.
   */
  struct DomainCache : public EvalCacheBase    
  {

    /**
     * Constructs an evaluation cache from the associated functional
     * f_, a point of linearization u_, and flags from the assembler (see VariationalFunctionalAssembler)
     */
    DomainCache(Functional const& f, typename Vars::VariableSet const& u, int const flags=7)  {};

    /**
     * Construct all data that is constant in an entity. Default
     * implementation in EvalCacheBase: does nothing
     */
    template<class Entity>
    void moveTo(Entity const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default
     * implementation in EvalCacheBase: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$f(u(x))\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$f\f$  with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of \f$f\f$ with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by \arg argT and \arg argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent
   * w.r.t. parallel execution. 
   */
  struct BoundaryCache : public EvalCacheBase    
  {
    /**
     *  \brief Switch for the assembler:
     *  - false: BoundaryCache is only evaluated at the domain boundary
     *  - true: BoundaryCache is evaluated also at interior faces, e.g. for computing face-gradient jumps for error estimation.
     */
    static const bool hasInteriorFaces=false;

/// Type of the Face Information
    typedef typename AnsatzVars::Grid::template Codim<0>::Entity::LeafIntersectionIterator FaceIterator;

    /**
     * Constructs an evaluation cache from the associated functional f_, a point of linearization u_, and flags from the assembler 
     */
    BoundaryCache(Functional const& f, typename Vars::VariableSet const& u, int const flags=7)  {};

    /**
     * Construct all data that is constant in an entity. Default implementation in EvalBaseCache: does nothing
     */
    template<class FaceIterator>
    void moveTo(FaceIterator const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default implementation in EvalBaseCache: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$g(u(x))\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$g\f$,
     * with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of 
     * \f$g\f$, with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by argT and the ansatz functions, given by argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };


  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row block of
   * the first derivative.
   */
  template <int row>
  struct D1
  {
    /**
     * Is true if that block is statically present, i.e. if \f$ f \f$
     * depends on variable \f$u_{\mathrm{row}}\f$.
     */
    static bool const present;
  };

  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the second derivative.
   */
  template <int row, int col>
  struct D2 
  {
    /**
     * Is true if that block is statically present, i.e. if variable
     * \f$u_{\mathrm{row}}\f$ and \f$u_{\mathrm{col}}\f$ are
     * nonlinearly coupled.
     *
     * For second derivatives, which are assumed to be always
     * symmetric, only the lower triangular blocks need be present.
     */
    static bool const present;

    /**
     * Is true if the block is conceptually symmetric, i.e. whether
     * EvaluationCache::d2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     */
    static bool const symmetric;

    /**
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
  };

  /**
   * This method defines a suitable quadrature order to be used on the
   * given cell. The maximum shape function order occuring in the
   * finite element spaces for this cell is provided for
   * convenience. If boundary is true, this refers to evaluation of
   * \f$g\f$, otherwise to the evaluation of \f$f\f$.
   */
  int integrationOrder(typename Variables::Grid::template Codim<0>::Entity const& cell,
                       int shapeFunctionOrder, bool boundary) const;
};


/**
 * \ingroup functional
 *
 * \brief Concept for providing block information to hierarchical error estimator.
 * 
 * Concept that the fourth template parameter of the
 * HierarchicErrorEstimator must model. Note that this is no base
 * class, just an interface documentation.
 */
struct HierarchicErrorEstimatorD2Info 
{
  /**
   * Member template that models the VariationalFunctional::D2 concept.
   */
  template <int row, int col> class D2;
};
