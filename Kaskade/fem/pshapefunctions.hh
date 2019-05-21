/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSHAPEFUNCTIONS_HH
#define PSHAPEFUNCTIONS_HH

#include <tuple>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

#include <boost/multi_array.hpp>


#include "fem/barycentric.hh"
#include "fem/fixdune.hh"
#include "fem/mllgeometry.hh"
#include "linalg/dynamicMatrix.hh"

namespace Kaskade
{
  /** \ingroup fem
   * \brief Shape functions
   *
   * This class provides a common interface to the access of shape functions
   *
   * \tparam ctype the scalar type used for coordinates (usually double)
   * \tparam dim the spatial dimension
   * \tparam T  Scalar type, the value field type of the shape functions
   * \tparam comp number of components
   *
   * A shape function lives on the reference element and may be scalar valued or vector valued:
   * \f[ \phi: R^d \to T^m \f]
   * with the spatial (grid) dimension \f$ d \f$, the value field type \f$ T \f$ (either a complex
   * or a real type), and the vectorial dimension \f$ m \f$ (comp) of the image space. Shape function
   * values are indeed defined not only on the reference element but on whole \f$ R^d \f$, even though
   * the values outside the reference elements may be of little interest.
   * It is possible to evaluate the shape function and its derivative.
   *
   */
  template <class ctype, int dim, class T=double, int comp=1> 
  class ShapeFunction
  {
  public:
    virtual ~ShapeFunction() {}

    /**
     * \brief Creates a copy of the actual derived shape function.
     */
    virtual std::unique_ptr<ShapeFunction> clone() const = 0;

    /**
     * \brief Evaluates the shape function (all components at once).
     */
    virtual Dune::FieldVector<T,comp> evaluateFunction(Dune::FieldVector<ctype,dim> const& xi) const = 0;

    /**
     * \brief Evaluates the derivative of the shape function (all components and all directions at once).
     */
    virtual Dune::FieldMatrix<T,comp,dim> evaluateDerivative(Dune::FieldVector<ctype,dim> const& xi) const = 0;

    /**
     * \brief Evaluates the second derivative of the shape function (all components and all directions at once).
     */
    virtual Tensor3<T,comp,dim,dim> evaluate2ndDerivative(Dune::FieldVector<ctype,dim> const& xi) const { 
      std::cerr << "NOT IMPLEMENTED: " << __FILE__ << ":" << __LINE__ << "\n"; 
      abort();
    }

    /**
     * \brief Returns a tuple (nominalOrder,codim,entity,index) giving detailed information about the location of the shape function. 
     * 
     * Each shape function is associated to a certain subentity of the element.
     *
     * nominalOrder is a nonnegative ordering parameter that is usually
     * the polynomial order of the shape function, but need not coincide.
     *
     * codim is the codimension of the subentity to which the shape
     * function is associated, entity is the number of the subentity,
     * and index is the number of the shape function among those that
     * are associated to the same subentity.
     */
    virtual std::tuple<int,int,int,int> location() const = 0;
  };

  //---------------------------------------------------------------------
  /** \ingroup fem
   *  \brief A set of shape functions
   *
   *
   * \tparam ctype scalar type for coordinates
   * \tparam dimension spatial dimension 
   * \tparam T  scalar type
   * \tparam comp number of components
   */
  template <class ctype, int dimension, class T, int comp=1>   
  class ShapeFunctionSet
  {
  public:
    /// Scalar field type
    typedef T Scalar;
    /// type of one shape function
    typedef ShapeFunction<ctype,dimension,T,comp> value_type;
    
    /**
     * \brief A matrix type mapping one-component coefficient vectors.
     */
    typedef DynamicMatrix<Dune::FieldMatrix<T,1,1> > Matrix;

    /// Constructor
    ShapeFunctionSet(Dune::GeometryType gt):
      order_(-1), size_(-1),
      gt_(gt),
      refElem(&Dune::ReferenceElements<ctype,dimension>::general(gt))
    {}

    ShapeFunctionSet(ShapeFunctionSet const& other)
    : iNodes(other.iNodes), projection(other.projection), order_(other.order_), size_(other.size_),
      gt_(other.gt_), refElem(&Dune::ReferenceElements<ctype,dimension>::general(gt_))
    {}

    /**
     * \brief Number of components of the shape function values.
     */
    static int const comps = comp;


    virtual ~ShapeFunctionSet() {}


    /// Random access to a shape function
    virtual value_type const& operator[](int i) const = 0;

    /**
     * \brief Number of shape functions in the set.
     */
    virtual int size() const { return size_; }

    /**
     * \brief Maximal polynomial order of shape functions.
     */
    int order() const { return order_; }

    /// Type of geometry on which the shape functions are defined.
    Dune::GeometryType type() const { return gt_; }

    /**
     * Returns a reference to the reference element on which this shape
     * function set is defined.
     */
    Dune::ReferenceElement<ctype,dimension> const& referenceElement() const { return *refElem; }

    /// A container type for holding interpolation points in the reference elements.
    typedef std::vector<Dune::FieldVector<ctype,dimension> > InterpolationNodes;

    /// A twodimensional array type for holding shape function values
    /// evaluated at a set of nodes.
    typedef DynamicMatrix<Dune::FieldMatrix<T,comp,1> > SfValueArray;

    /**
     * \brief Initialize the hierarchical projection matrix based on the given
     * lower order shape function set.
     */
    void initHierarchicalProjection(ShapeFunctionSet<ctype,dimension,T,comp> const* sfl)
    {
      // If no lower order shape function set is given, we assume an
      // empty set, i.e. the projection is just zero.
      if (!sfl || size()==0) {
        projection.setSize(size(),size());
//       projection = 0 // causes compilation error with dune-2.4.0 and clang++ on OS X (Darwin)
        projection.fill(0);
        return;
      }

      // Compute the hierarchic projection P = I_h E_l(x_h) I_l E_h(x_l),
      // where E_h(x_l) the evaluation of high order shape functions at
      // low order nodes, I_l the low order interpolation, E_l(x_h) the
      // evaluation of low order shape functions at high order nodes,
      // and I_h the high order shape function evaluation.

      // Compute the local restriction matrix A
      SfValueArray Ehxl;
      evaluate(sfl->interpolationNodes(),Ehxl);  // E_h(x_l)
      Matrix A;
      sfl->interpolate(Ehxl,A);                  // I_l  E_h(x_l)

      // Compute the local prolongation matrix B
      SfValueArray Elxh;
      sfl->evaluate(interpolationNodes(),Elxh);  // E_l(x_h)
      Matrix B;
      interpolate(Elxh,B);                       // I_h E_l(x_h)

      // Compute projection P=B*A
      MatMult(projection,B,A);
    }



    /**
     * @brief Interpolation points.
     *
     * Provides interpolation points such
     * that the shape function coefficients can be computed from
     * function values at interpolation nodes by multiplication by this
     * matrix.  The interpolation points are guaranteed to be inside the
     * reference element associated to this shape function set.
     */
    InterpolationNodes const& interpolationNodes() const { return iNodes; }


    /**
     * @brief Left-multiplies the provided matrix with the interpolation
     * matrix of the shape function set.
     *
     * Each column of A is interpreted as values of some function
     * evaluated at this shape function set's interpolation points (see
     * below). The columns of the output array IA then contain the shape
     * function coefficients such that the corresponding
     * linearcombination of shape functions "interpolates" that function
     * in the interpolation points. What "interpolation" means is up to
     * the actual implementation.
     *
     * IA is automatically resized if needed.
     *
     * Storage order: A[i][j] contains the value of function j at
     * interpolation node i.
     */
    virtual void interpolate(SfValueArray const& A, Matrix& IA) const = 0;

    /**
     * @brief Evaluate shape function set at a set of points. In
     * notation of the LocalToGlobalMapperConcept, this gives the matrix
     * \f$ \Phi \f$: the entry Phi[i][j] is the value of shape function
     * j evaluated at iNodes[i].
     *
     * @param iNodes the points at which the shape functions are to be evaluated.
     * \param phi    the array that is filled with shape function values. The array
     *               will be resized if needed.
     */
    void evaluate(InterpolationNodes const& iNodes, SfValueArray& phi) const
    {
      int s = size();
      phi.setSize(iNodes.size(),s);

      for (int j=0; j<s; ++j) {
        value_type const& sf = (*this)[j];
        for (int i=0; i<iNodes.size(); ++i) {
          Dune::FieldVector<T,comp> v = sf.evaluateFunction(iNodes[i]);
          for (int k=0; k<comp; ++k)
            phi[i][j][k] = v[k];  // todo: improve efficiency! cache directly!
        }
      }
    }

    /**
     * \brief Returns a square matrix that projects shape function
     * coefficients to a subspace spanned by shape functions of lower
     * order.  This is intended to be used to implement embedded error
     * estimators. The actual definition of this subspace depends on the
     * shape function set specified in the call to
     * initHierarchicalProjection().
     */
    Matrix const& hierarchicProjection() const
    {
      // check that the projection has been properly initialized.
      assert(projection.N()==size() && projection.M()==size());

      return projection;
    }

    virtual void removeShapeFunction(size_t index) {}

  protected:
    InterpolationNodes iNodes;
    Matrix             projection;
    int                order_;
    int                size_;

  private:
    Dune::GeometryType gt_;
    Dune::ReferenceElement<ctype,dimension> const* refElem;
   };

  //---------------------------------------------------------------------
  /**
   * \brief Restricted shape function set.
   * Introduces a new local ordering for the shape functions. To retrieve the original shape function id use getId(int newLocalId)
   * \todo docme: what is is good for? Which ordering? Why "restricted"?
   */
  template <class ShapeFunctionSet_>
  class RestrictedShapeFunctionSet: public ShapeFunctionSet_
  {
  public:
    /// Grid type
    //typedef typename ShapeFunctionSet_::Grid Grid; causes compiler error
    /// Scalar field type
    typedef typename ShapeFunctionSet_::Scalar Scalar;
    /// type of one shape function
    typedef typename ShapeFunctionSet_::value_type value_type;

    explicit RestrictedShapeFunctionSet(int order) : ShapeFunctionSet_(order) {}

    RestrictedShapeFunctionSet(ShapeFunctionSet_ const& other) : ShapeFunctionSet_(other) {}
    RestrictedShapeFunctionSet(RestrictedShapeFunctionSet const& other) : ShapeFunctionSet_(other) {}
    /// Copy constructor. Resets restriction ids
    RestrictedShapeFunctionSet(RestrictedShapeFunctionSet const& other, std::vector<int>* ids_) : ShapeFunctionSet_(other) {}

    virtual ~RestrictedShapeFunctionSet(){}

    /**
     * \param ids_ vector of used ids of the shape function set
     */
    void setRestriction(std::vector<int> const& ids_)
    {
      if(ids_.size()==0)
      {
        this->sf.clear();
        this->iNodes.clear();
        this->size_ = 0;
        return;
      }


      for(int i=this->sf.size()-1; i>-1; --i)
        if(std::find(ids_.begin(),ids_.end(),i)==ids_.end()) this->removeShapeFunction(i);

      this->size_ = ids_.size();
    }
  };

  //---------------------------------------------------------------------

  /**
   * \brief Base class for sets of shape function containers.
   *
   * \tparam ctype scalar type for coordinates
   * \tparam dimension spatial dimension
   * \tparam T the scalar shape function return value
   * \tparam comp the number of components of the shape functions' values
   */
  template <class ctype, int dimension, class T, int comp=1>
  class ShapeFunctionSetContainer
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,T,comp> value_type;

    virtual ~ShapeFunctionSetContainer() {}


    /// access a shape function via type and order
    virtual value_type const& operator() (Dune::GeometryType type, int order) const = 0;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * For a given shape function set \arg sfs which is defined on a
   * simplex and permits a simple coupling, and a given permutation \arg
   * vPermutation of the vertices of the simplex defining an affine
   * spatial transformation \f$ f \f$, this function returns a
   * permutation \arg sfPermutation of the shape functions and a sign
   * vector \arg sign, such that
   *
   * \f[ \mathrm{sign[k]}\phi_{\mathrm{sfPermutation[k]}}(f(x)) =
   * \phi_k(x). \f]
   */
  template <class ShapeFunctionSet>
  void computeSimplexSfPermutation(ShapeFunctionSet const& sfs,
      int const* vPermutation,
      int*       sfPermutation,
      int*       sign)
  {
    assert(vPermutation);
    assert(sfPermutation);
    assert(sign);

    // Ok, we rely on the evaluation & interpolation of the shape
    // function set. First we evaluate the shape functions on the
    // interpolation nodes transformed by f, then we interpolate these
    // values. This gives us a matrix which, if the shape function set
    // indeed allows simple coupling, is a signed permutation
    // matrix. Finally we extract the signs and the permutation from
    // this matrix.

    // Get the interpolation nodes. Make a copy because we need to
    // modify them.
    typedef typename ShapeFunctionSet::InterpolationNodes InterpolationNodes;
    InterpolationNodes in = sfs.interpolationNodes();


    // Transform the nodes. This is done by permuting their barycentric
    // coordinates. Remember that the interpolation nodes are given in
    // Cartesian, not in barycentric coordinates. This determines the
    // loop termination for j.
    for (int i=0; i<in.size(); ++i) {
      Dune::FieldVector<typename ShapeFunctionSet::Grid::ctype,ShapeFunctionSet::Grid::dimension+1> b = barycentric(in[i]);
      for (int j=0; j<ShapeFunctionSet::Grid::dimension; ++j)
        in[i][j] = b[vPermutation[j]];
    }


    // Evaluate the shape functions at the transformed nodes.
    typename ShapeFunctionSet::SfValueArray phi(in.size(),sfs.size());
    sfs.evaluate(in,phi);

    // Interpolate the shape function values.
    typename ShapeFunctionSet::Matrix P(sfs.size(),sfs.size());
    sfs.interpolate(phi,P);


#ifndef NDEBUG
    // Check that the permutation matrix is indeed a signed permutation
    bool ok = true;
    for (int i=0; i<sfs.size(); ++i) {
      double rowSum = 0;
      double colSum = 0;
      for (int j=0; j<sfs.size(); ++j) {
        if (std::abs(P[i][j])>1e-8 && std::abs(1-std::abs(P[i][j]))>1e-8) ok = false;
        rowSum += std::abs(P[i][j]);
        colSum += std::abs(P[j][i]);
      }
      if (std::abs(rowSum-1)>1e-8 || std::abs(colSum-1)>1e-8) ok = false;
    }
    if (!ok) {
      std::cout << "\nPermutation matrix:\n";
      for (int i=0; i<P.N(); ++i) {
        for (int j=0; j<P.M(); ++j)
          std::cout << std::setw(2) << P[i][j] << " ";
        std::cout << "\n";
      }
    }
#endif


    // Extract the permutation and the sign.
    for (int i=0; i<sfs.size(); ++i)
      for (int j=0; j<sfs.size(); ++j)
        if (std::abs(P[j][i])>0.5) {
          sfPermutation[i] = j; // XXX or the other way round?
          sign[i] = P[j][i]>0? 1: -1;
          break;
        }

  }
} // end of namespace Kaskade

#endif
