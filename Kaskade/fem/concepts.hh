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
 * @brief  Some important concepts of Kaskade
 * @author Martin Weiser, Anton Schiela
 */
class unspecified;

/**
 * \ingroup concepts
 * \brief Management of degrees of freedom.
 *
 * LocalToGlobalMappers define how a global FE coefficient vector is
 * to be interpreted as a FE function in terms of shape functions on
 * each cell of the grid.
 *
 * Finite element functions are given by piecewise linear combinations
 * of shape functions on cells. On each cell \f$ T \f$, a finite
 * element function \f$ u \f$ is given as \f[ u|_T(x) = \sum_{j=1}^n
 * (K a_{I})_j \psi(x) \phi_{j}(x)  = \psi(x) (\phi^T(x) K a). \f]
 * Here, \f$ a \f$ denotes the global coefficient vector, \f$ I \in
 * {\bf N}^m, \, m\le n, \f$ is the set of \f$ m \f$ global degrees of freedom
 * that is associated to this cell, and \f$ K\in {\bf R}^{n\times m}
 * \f$ is a matrix (called combiner) describing which linear
 * combination of the \f$ n \f$ shape functions yields which global ansatz function
 * (FE basis function) restricted to the cell:
 * \f[ (\varphi|_T)_i = \sum_{j=1}^n K_{ji} (\psi\phi)_j \f]
 * 
 * 
 * Often, \f$ K \f$ is
 * simply the identity, or at least diagonal (e.g. for edge elements,
 * where the global orientation of edges has to be taken into
 * account), but it may be non-square and arbitrarily dense (e.g. for
 * hp-methods or hanging nodes). \f$ \psi(x) \f$ is a linear geometric
 * transformation of shape function values (and gradients). Often this
 * is the identity, but for vectorial shape functions it may depend on
 * the Jacobian at the evaluation point inside the cell (e.g. for edge
 * elements). Even in the most trivial case, \f$ \psi \f$ transforms
 * at least the shape function gradients due to the chain
 * rule. \f$\phi\f$ are the shape functions themselves. They may be
 * scalar or vectorial. Note that vectorial finite element functions
 * can perfectly well be defined in terms of scalar shape functions:
 * Each component is defined in terms of the scalar shape
 * functions. Often the set of shape functions is the same on each
 * cell, but this need not be the case (e.g. for hp-methods).
 *
 * Working with sets of evaluation points \f$ x_i \f$, the
 * representation formula may conveniently be written as \f[ u|_T =
 * \Psi \Phi K a_I, \f] with \f$ u \f$ being the column vector of
 * values \f$ u(x_i)\f$, \f$\Psi = \mathrm{diag}(\psi(x_i))\f$, and
 * \f$ \Phi_{ij} = \phi_j(x_i) \f$.
 *
 * LocalToGlobalMappers define the ingredients for the above
 * representation. Their interface is described here.
 */
class LocalToGlobalMapperConcept {
public:
  /// The type of scalars used.
  typedef unspecified RT;

  /// The grid type on which the finite element functions are defined.
  typedef unspecified Grid;

  /// The subset of the grid on which the finite element functions are defined.
  typedef unspecified IndexSet;

  /// The type of shape function sets living on the cells.
  typedef unspecified ShapeFunctionSet;
  
  /**
   * \brief A standard sequence with value type IndexSet::IndexType
   */
  typedef unspecified GlobalIndexRange;

  /**
   * \brief A standard sequence with value type std::pair<IndexSet::IndexType,int>.
   */
  typedef unspecified SortedIndexRange;

  /// The type of the geometrical transformation \f$ \psi \f$.
  typedef ConverterConcept Converter;

  /// The type of the algebraic combiner \f$ K \f$ associated with this cell.
  /** This type represents a matrix of small to moderate size. For
   * performance reasons, it is a class of its own with special
   * interface.
   */
  typedef CombinerConcept Combiner;

  /// The type of codimension 0 entities in the grid.
  typedef typename Grid::template Codim<0>::Entity Cell;

  /// Returns the shape function set living on the given cell.
  ShapeFunctionSet const& shapefunctions(Cell const& cell) const;

  /// Returns the shape function set living on the cell with given index.
  ShapeFunctionSet const& shapefunctions(size_t index) const;


  /// Returns the number of global degrees of freedom managed.
  size_t size() const;
  
  /**
   * \brief Returns the maximal polynomial order of shape functions encountered in any cell.
   */
  int maxOrder() const;

  /**
   * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
   */
  SortedIndexRange sortedIndices(Cell const& cell) const;

  /**
   * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
   */
  SortedIndexRange sortedIndices(size_t n) const;

  /** 
   * \brief The global indices \f$ I \f$ associated with this cell.
   * Returns a const sequence containing the global indices of the
   *  ansatz functions with nonvanishing support on this cell.
   */
  GlobalIndexRange globalIndices(Cell const& cell) const;

  /** 
   * \brief The global indices \f$ I \f$ associated with the cell with given index.
   * Returns a const sequence containing the global indices of the
   *  ansatz functions with nonvanishing support on this cell.
   */
  GlobalIndexRange globalIndices(size_t index) const;


  /**
   * \brief Returns the algebraic combinator.
   * \param cell the grid cell for which the combiner is requested
   * \param index the index of the cell
   */
  Combiner combiner(Cell const& cell, size_t index) const;

  /// Recomputes the node management, e.g. after grid modifications.
  void update();
};


/**
 * \ingroup concepts
 * \brief Geometrical transformation of shape function values.
 *
 * Converters \f$ \psi(x) \f$ transform shape function values defined
 * on the reference elements to their values on an actual cell. This
 * is nontrivial for vectorial shape functions, e.g. edge
 * elements. Converters do only depend on the geometry of the actual
 * cell, not on any combinatorial or algebraic information. They can
 * therefore instantiated without a valid LocalToGlobalMapperConcept
 * object.
 *
 * Here we assume that Grid is the grid type, RT the scalar type in
 * use, and k the number of components of the shape functions.
 *
 * Converters have to be assignable and copy constructible.
 */
class ConverterConcept {
public:
  /**
   * \brief Create a converter.
   *
   * Default constructed converters are good for nothing. Use moveTo()
   * to define the cell on which they do operate.
   */
  ConverterConcept();

  /// Convenience constructor. Equivalent to default constructor
  /// subsequented by moveTo.
  ConverterConcept(Grid::template Codim<0>::Entity const& cell);

  /**
   * \brief Defines the cell on which the converter lives.
   *
   * The cell that is given has to persist as long as the converter is
   * used on this cell.
   */
  void moveTo(Grid::template Codim<0>::Entity const& cell);

  /// Sets the argument \f$ x \f$ of the transformation \f$ \psi \f$.
  void setLocalPosition(Dune::FieldVector<typename Grid::ctype,Grid::dimension> const& x);

  /// Applies the transformation \f$ \psi(x) \f$ to shape function value.
  Dune::FieldMatrix<RT,k,1> global(Dune::FieldMatrix<RT,k,1> const& sf) const;

  /// Applies the transformation \f$ \psi(x) \f$ to shape function value and derivatives up to order deriv.
  VariationalArg<RT,dim,k> global(VariationalArg<RT,dim,k> const& sf, int deriv=1) const;

  /**
   * \brief Applies the transformation \f$ \psi \f$ to shape function value, gradient, and Hessian, returning a filled VariationalArg.
   */
  template <class Scalar>
  VariationalArg<Scalar,dim,1> global(std::tuple<Dune::FieldVector<Scalar,1>,
                                                 Dune::FieldMatrix<Scalar,1,dim>,
                                                 Tensor3<Scalar,1,dim,dim>> const& sf) const;
                                                 
  /// Transforms global values to local shape function values.
  Dune::FieldMatrix<RT,k,1> local(Dune::FieldMatrix<RT,k,1>) const;
};


/**
 * \ingroup concepts
 * \brief Algebraic combiner of global degrees of freedom to shape function coefficients.
 *
 * Combiners \f$ K \in \R^{n\times m} \f$ map the global degrees of freedom to linear
 * combinations of shape functions, which represent the FE ansatz
 * function on the given cell. To be precise, the global ansatz function \f$ \varphi \f$ 
 * restricted to the cell \f$ T \f$ is given as
 * \f[ (\varphi|_T)_i = \sum_{j=1}^n K_{ji} (\psi\phi_j) \f]
 * in terms of the \f$ n \f$ local shape functions \f$ \phi_j \f$ and the Converter \f$ \psi \f$. 
 * 
 * 
 * Often this is trivial (just the
 * identity), but may be diagonal containing 1 and -1 entries if
 * global orientation has to be taken into account (e.g. edge
 * elements), but may be a sparse or dense non-square matrix as well
 * (trace spaces, or hp methods, hanging nodes).
 *
 * For performance reasons, \f$ K \f$ is not simply represented as a
 * matrix, but as an object with special interface.
 *
 * Here we assume that RT the scalar type in use.
 *
 * Combiners have to be assignable and copy constructible.
 */
class CombinerConcept {
public:
  /// In-place computation of \f$ A \leftarrow A K \f$.
  template <int n, int m>
  void rightTransform(DynamicMatrix<Dune::FieldMatrix<RT,n,m>>& A) const;

  /**
   * In-place computation of row vectors \f$ v \leftarrow v K \f$.
   * Given a sequence of shape function values, this computes the
   * corresponding sequence of global ansatz function values.
   */
  template <int n, int m>
  void rightTransform(std::vector<VariationalArg<RT,n,m> >& v) const;

  /**
   * In-place computation of \f$ A \leftarrow K^+ A \f$. Given a set
   * of shape function coefficients, this computes the best fit of
   * global ansatz function coefficients.
   */
  template <int n, int m>
  void leftPseudoInverse(DynamicMatrix<Dune::FieldMatrix<RT,n,m> >& A) const;

  /// Implicit conversion to a sparse matrix.
  operator Dune::BCRSMatrix<Dune::FieldMatrix<RT,1,1>>() const;
};

//---------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Function is the interface for evaluatable functions \f$ f:\Omega\to\mathbb{R}^m \f$.
 *
 * Functions allow the evaluation of possibly vector valued
 * functions defined on a spatial domain given by a grid and provide
 * some metadata.
 *
 * Prominent models of the Function  concept are \ref FunctionSpaceElement s.
 */
class Function
{
public:
  /**
   * A scalar floating point type, which is the field over which the
   * return values are defined.
   */
  typedef unspecified Scalar;

  /**
   * The number of vectorial components of the result value.
   */
  static int const components = unspecified;

  /**
   * The type of the function values.
   */
  typedef Dune::FieldVector<Scalar,components> ValueType;

  /**
   * Returns the polynomial order of the function on the specified
   * cell (grid entity of codim 0). If the function is locally not
   * polynomial, std::numeric_limits<int>::max() is returned.
   */
  template <class Cell>
  int order(Cell const& cell) const;

  /**
   * Returns the value of the function at the specified local
   * coordinate of the given cell.
   */
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::ctype,
                  Cell::dimension> const& localCoordinate) const;
};

//---------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief A Weighing is a class that allows to compute weighted averages of values associated to cells of a grid.
 *
 * Weighings are used in interpolateGlobally and interpolateGloballyWeak to define
 * reasonable FE coefficients whenever the data to be interpolated cannot be represented
 * in the FE space, e.g. interpolating a discontinuous function in a space of continuous
 * functions. The conflicting interpolation values at one and the same interpolation node
 * are then averaged - the kind of averaging is specified by the Weighing.
 *
 * Predefined weighings are InverseVolume, Volume, PlainAverage.
 */
struct Weighing
{
  /**
   * \brief This defines a weight \f$ w(c) \f$ for each cell \f$ c \f$ of the grid.
   */
  template<class Cell>
  double operator()(Cell const& c) const;
  /**
   * \brief Defines the global scaling factor \f$ s \f$.
   */
  double scalingFactor() const;
  /**
   * \brief Defines the global offset \f$ o \f$.
   */
  double offset() const;
};

//---------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Allows to implement weighted (scaled) norms.
 *
 * This interface needs to be implemented for computing spatially
 * weighted norms. Plain, unscaled norms can be constructed using the
 * trivial IdentityScaling.
 */
class Scaling
{
public:
  /**
   * \tparam Cell a Dune grid entity type of codimension 0
   * \tparam Sequence a boost::fusion sequence containing the values to be scaled
   * \param cell the cell on which the scaling should be applied
   * \param localPos the position at which the scaling schould be applied, given in local coordinates
   * \param values the values that are to be scaled in place. This is a heterogeneous vector of
   *        Dune::FieldVector's that contain the values of a set of FE functions evaluated at the
   *        given point.
   */
  template <class Cell, class Sequence>
  void scale(Cell const& cell,
             Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localPos,
             Sequence& values) const;
};

//---------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Defines the block structure of defect systems to be solved by the hierarchical error estimator.
 *
 * The purpose of this class is to provide information to the HierarchicErrorEstimator about how the
 * blocks in the defect system are to be assembled. Blocks are indexed by (row,col), and individual
 * information is provided by the member template D2.
 *
 * Two bool flags need to be defined in D2<row,col>:
 * - present: whether the block is assembled at all
 * - lumped: whether only the diagonal of that block has to be assembled (true) or all entries (false)
 */
class HierarchicErrorEstimator_D2BlockInfo
{
  template <int row, int col>
  class D2
  {
    static bool const present;
    static bool const lumped;
  };
};

//---------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Defines a linear operator \f$A:\ X\rightarrow Y\f$ with matrix representation \f$M_A \in \mathbb{K}^{n\times m}\f$.
 *
 * This concept adds the following commutative diagram to the Dune::LinearOperator concept:
 *
 * \f[ \begin{CD}
 *       X @>A>> Y \\
 *     @AA{S^{-1}_X}A\  @VV{S_Y}V \\
 *     \mathbb{K}^n @>M_A>> \mathbb{K}^m \\
 *     \end{CD} \f]
 *
 * The mappings \f$S_X\f$ and \f$S_Y\f$ are the coordinate mappings of the domain and range space of \f$A\f$. \f$M_A\f$ is the matrix
 * representation corresponding to the canonical basis in \f$\mathbb{K}^n\f$ and \f$\mathbb{K}^m\ (\mathrm{dim}(X)=n, \mathrm{dim}(Y)=m) \f$.
 *
 * The matrix-representation \f$M_A\f$ can be accessed via MatrixGeneratedOperator<X,Y>::get<MatrixType>() thus leaving the choice of matrix type
 * to the user (i.e. Dune::BCRSMatrix, MatrixAsTriplet,...). Note that the set of matrix types that can be obtained depends on the 
 * actual model of this concept.
 *
 * For the application of \f$S^{-1}_X\f$ or \f$S_Y\f$ use the member functions:
 *  - template <class Vector> vectorToDomain(Vector const&, X&)
 *  - template <class Vector> rangeToVector(Y const&, Vector&)
 * The direct application of the inverse operations \f$S_X\f$ and \f$S^{-1}_Y\f$ is not supported.
 *
 */
template <class X, class Y>
class MatrixRepresentedOperator: public Dune::LinearOperator<X,Y>
{
public:
  /**  type of the operator's domain */
  typedef X domain_type;
  /** type of the operator's domain */
  typedef Y range_type;

  /// Get matrix-representation of the operator.
  template <class MatrixType>
  MatrixType get() const;

  /// Get pointer to matrix-represenation as there still exist matrix implementations that don't support move semantics.
  template <class MatrixType>
  std::shared_ptr<MatrixType> getPointer() const __attribute__((deprecated));

  /// Get coordinate vector \f$coord\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
  /**
   * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$coord=S_Y(y)\f$.
   *
   * The used vector type Vector must provide:
   * - its iterator type via typename Vector::iterator
   * - (possibly overloads of) the free functions:
   *    - typename Vector::iterator std::begin(Vector&);
   *    - typename Vector::iterator std::end(Vector&)
   */
  template <class Vector>
  void rangeToVector(Y const& y, Vector& coord) const;

  /// Get \f$x\in X\f$ from coordinate vector \f$coord\in\mathbb{K}^n\f$
  /**
   * Apply \f$S^{-1}_X\f$ to \f$coord\f$: \f$x=S^{-1}_X(x)\f$.
   *
   * The used vector type Vector must provide:
   * - its iterator type via typename Vector::const_iterator
   * - (possibly overloads of) the free functions:
   *    - typename Vector::const_iterator std::begin(Vector const&);
   *    - typename Vector::const_iterator std::end(Vector const&);
   */
  template <class Vector>
  void vectorToDomain(Vector const& coord, X& x) const;
};


/**
 * \ingroup concepts
 * \brief Defines an interface for integrating values computed over a grid.
 */
class Collector
{
public:
  /**
   * \brief Copy constructor.
   */
  Collector(Collector const&); 
  
    template <class Cell>
    int integrationOrder(Cell const& /* ci */, int shapeFunctionOrder) const;
    
  /**
   * \brief Collects values  computed on the grid.
   * 
   * \param cell a reference to a codim 0 entity in the grid
   * \param idx the index of the cell (within the used grid view)
   * \param weigth the integration weight of the quadrature rule on the cell
   * \param x the values computed at this integration point
   */
  template <class CellIterator, class Index, class Sequence>
  void operator()(CellIterator const& cell, Index idx, double weight, Sequence const& x);
  
  /**
   * \brief Merges the result of two concurrently filled collectors.
   * 
   * Leaves the other collector c in a valid, destructable but unspecified state.
   */
  void join(Collector& c);
};

/**
 * \ingroup concepts
 * \brief A function that supports efficient evaluation via an evaluator.
 * 
 * This interface is provided by finite element functions (see \ref FunctionSpaceElement), but
 * can also be provided by adaptors computing various values from finite element functions
 * on the fly.
 * 
 * \see interpolateGlobally
 * \see \ref WeakFunctionView
 */
class FunctionView
{
public:
  /**
   * \brief The type of finite element space the evaluators of which are required for evaluating this FunctionView.
   */
  using Space = unspecified;
  
  /**
   * \brief Provides a reference to the finite element space required for evaluating this FunctionView.
   */
  Space const& space() const;
  
  /**
   * \brief Evaluates the function view.
   */
  Dune::FieldVector<unspecified,unspecified> value(typename Space::Evaluator const& evaluator) const;
};

/**
 * \ingroup concepts
 * \brief A function that supports efficient evaluation via cell and local coordinate.
 * 
 * This interface is provided by finite element functions (see Kaskade::FunctionSpaceElement), but
 * can also be used to provide various analytically defined functions.
 * 
 * \see \ref Kaskade::interpolateGloballyWeak
 * \see \ref FunctionView
 */
class WeakFunctionView
{
public:
  /**
   * \brief Evaluates the function view.
   * \param cell a codim 0 entity in the grid
   * \param xloc the local coordinate at the evaluation point in the cell
   */
  Dune::FieldVector<unspecified,unspecified> value(Cell const& cell, Position const& xloc) const;
};



/**
 * \ingroup concepts
 * \brief A hyperelastic material law.
 * 
 * \see Kaskade::MaterialLawBase
 */
class HyperelasticMaterialLaw
{
public:
  /**
   * \brief The scalar field type (usually double).
   */
  using Scalar = unspecified;
  
  /**
   * \brief The spatial dimension (usually 2 or 3).
   */
  static int const dim = unspecified;
  
  /**
   * \brief The type of stress and strain tensors (usually Dune::FieldMatrix<Scalar,dim,dim>).
   */
  using Tensor = unspecified; 
      
  /**
   * \brief The type of stress and strain tensors (usually Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2>).
   */
  using VoigtTensor = unspecified; 
  
  /**
   * \brief Defines a new evaluation/linearization point.
   * \param e the strain tensor
   */
  void setLinearizationPoint(Tensor const& e);
      
  /**
   * \brief Evaluates the stored energy density \f$ W(\epsilon) \f$.
   */
  Scalar d0() const;
  
  /**
   * \brief Evaluates the first directional derivative \f$ W'(\epsilon)\epsilon_1 \f$.
   */
  Scalar d1(Tensor const& e1) const;
  
  /**
   * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
   */
  Scalar d2(Tensor const& e1, Tensor const& e2) const;
  
  /**
   * \brief Returns the 2nd Piola-Kirchhoff stress tensor \f$ \sigma = W'(E) \f$ corresponding to the current strain \f$ E \f$.
   * 
   * The following invariant holds for any strain tensor eps:
   * \code
   * d1(eps) == contraction(stress(),eps)
   * \endcode
   */
  Tensor stress() const;
  
  /**
   * \brief Returns the tangent stiffness tensor \f$ C(E) = W''(E)\f$ mapping strain tensor variations to stress tensor (2nd Piola-Kirchhoff) variations.
   * 
   * In this context, the symmetric matrices \f$ \sigma \f$ and \f$ \epsilon \f$ are interpreted as \f$d(d+1)/2\f$-vectors in Voigt notation as
   * returned by Kaskade::pack. Note that strain tensor off-diagonal entries appear with factor two in the vector.
   * 
   * Now, a matrix \f$ C\in\mathbb{R}^{(d+1)d/2\times(d+1)d/2} \f$ is returned such that for infinitesimal changes \f$ \delta\epsilon \f$ the 
   * corresponding change in \f$ \sigma(\epsilon+\delta\epsilon) = \sigma(\epsilon)+\delta\sigma\f$ is given as \f$ \delta\sigma = C \delta\epsilon \f$. 
   * Note that due to the symmetry of the stored energy Hessian \f$ W''(E) \f$, the Voigt representation of the stiffness tensor is a
   * symmetric matrix.
   * 
   * The following invariant holds for any strain tensors eps1 and eps2:
   * \code
   * d2(eps1,eps2) == pack(eps1)*(strainToStressMatrix()*pack(eps2))
   * \endcode
   * 
   * The Drucker/Hill stability criterion states that the returned matrix shall be positive semidefinite.
   */
   VoigtTensor strainToStressMatrix() const;
};

/**
 * \ingroup concepts
 * \brief A hyperelastic material law formulated in terms of the invariants \f$ I_1, I_2, I_3 \f$ of the Cauchy-Green strain tensor \f$ C \f$.
 * 
 * \see Kaskade::InvariantsMaterialLaw
 */
class InvariantsMaterialConcept
{
public:
  /**
   * \brief The scalar field type (usually double).
   */
  using Scalar = unspecified;
  
  /**
   * \brief The spatial dimension (usually 2 or 3).
   */
  static int const dim = unspecified;
  
  /**
   * \brief A dim-dimensional vector type.
   */
  using Invariants = Dune::FieldVector<Scalar,dim>;
  
  /**
   * \brief Constructs the material at shifted invariants \f$ i \f$.
   * 
   * The invariants of \f$ A \f$ are not given by their actual value \f$ I_k \f$, but in the shifted invariants \f$ i_k \f$ such that
   * \f$ i_k(I) = 0 \f$. This allows a numerically stable evaluation of the material laws in the vicinity of 
   * the reference configuration.
   * 
   * For \f$ d=2 \f$ we have \f$ i_1 = I_1-2 = \mathrm{tr}(A)-2 \f$ and \f$ i_2 = I_2-1 = |A| - 1 \f$. 
   * For \f$ d=3 \f$ we have \f$ i_1 = I_1-3 = \mathrm{tr}(A)-3 \f$, \f$ i_2 = I_2-3 = \frac{1}{2}(\mathrm{tr}(A)^2 - \mathrm{tr}(A^2)) - 3\f$, 
   * and \f$ i_3 = I_3-1 = |A| - 1 \f$.
   * 
   * \param i the shifted invariants
   */
  void InvariantsMaterialConcept(Invariants const& i, ... );
  
  /**
   * \brief Evaluates the stored energy density \f$ W(I) \f$.
   */
  Scalar d0() const;
  
  /**
   * \brief Evaluates the first directional derivative \f$ W'(I)\tilde I \f$.
   */
  Scalar d1(Invariants const& i1) const;
  
  /**
   * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
   */
  Scalar d2(Invariants const& i1, Invariants const& i2) const;  
};

