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

#ifndef HIERARCHICSHAPEFUNCTIONS_HH
#define HIERARCHICSHAPEFUNCTIONS_HH

/**
 * \file
 * \brief  define hierarchic shape functions for simplicial elements of
 *         arbitrary dimension and order
 * \author Martin Weiser
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include "fem/barycentric.hh"
#include "fem/combinatorics.hh"
#include "fem/gridcombinatorics.hh"
#include "fem/pshapefunctions.hh"
#include "linalg/simpleLAPmatrix.hh"

//---------------------------------------------------------------------
namespace Kaskade
{
  // forward declarations
  template <class,int,class> class HierarchicSimplexShapeFunctionSet;
  template <class,int,class> class HierarchicExtensionSimplexShapeFunctionSet;

  /**
   * \cond internals
   */
  namespace HierarchicDetail {

    /**
     * Returns the number of subentites of given codimension in the simplex of dimension \arg d.
     *
     * \param[in] d the spatial dimension (less or equal to 3)
     * \param[in] c the codimension of subentities
     */
    int nSubentities(int d, int c);


    /**
     * Returns the number of hierarchic shape functions of order p
     * associated to c-codimensional subentities of a d-dimensional
     * simplex.
     */
    int cdimSize(int d, int p, int c);


    /**
     * Returns the number of hierarchic shape functions of order \arg p on a
     * \arg d -dimensional simplex.
     */
    int size(int d, int p);


    /**
     * Returns a pair (c,e) where c is the codimension of the subentity to
     * which the k-th shape function is associated and e the number of the
     * subentity on which the shape function is located.
     */
    std::pair<int,int> codim(int d, int p, int k);


    /**
     * Computes the tuple index from the sequential index. For a given
     * spatial dimension and nominal order, any shape function can be
     * identified either by a sequential number (the sequential index) or
     * by the tuple (c,e,k), where c is the codimension of the subentity
     * to which the shape function is associated, e the subentity number
     * among all subentities with same codimension, and k a sequential
     * number of the shape function associated to the subentity.
     */
    std::tuple<int,int,int> tupleIndex(int d, int p, int k);


    /**
     * Computes the value of the Zumbusch shape function at given barycentric coordinates.
     *
     * \param[in] d is the dimension of the simplex
     * \param[in] p is the order of the shape function
     * \param[in] k is the number of the shape function (among those with same (d,p)
     * \param[in] b contains the barycentric coordinates
     */
    double shapeFunction(int d, int p, int k, double const b[]);


    /**
     * Computes the derivative of the Zumbusch shape function at given barycentric coordinates.
     *
     * \param[in] d is the dimension of the simplex
     * \param[in] p is the order of the shape function
     * \param[in] k is the number of the shape function (among those with same (d,p)
     * \param[in] b contains the barycentric coordinates
     * \param[out] grad is filled with the d+1 partial derivatives w.r.t. the coordinates
     */
    void dShapeFunction(int d, int p, int k, double const b[], double* grad);


    /**
     * Computes the actual polynomial order of the shape function. This
     * can be the nomial order or one more.
     */
    int actualOrder(int d, int p, int k);


    /**
     * Computes the number of the subentity of codimension \arg c in the
     * reference simplex that has the same vertex number as the normalized
     * vertex numbers \arg vGlobal of the subentity with number \arg e.
     *
     * \param[in] d the spatial dimension
     * \param[in] c the codimension of the subentity
     * \param[in] e the number of the subentity with codimension c in the reference simplex
     * \param[in] vGlobal global indices of the vertices of the reference simplex
     */
    int entityPermutation(int d, int c, int e, int const vGlobal[]);


    /**
     * Computes a globally unique index K of the ansatz function associated to the shape function
     * given by d,p,c,e,k on the reference simplex. The index K denotes the globally unique number
     * within the subentity and order p. Additional, the sign is returned, whether the shape function
     * or its negative is the ansatz function restricted to the simplex.
     *
     * \param[in] d spatial dimension
     * \param[in] p nominal order of shape function
     * \param[in] c codimension of subentity to which the shape function is associated
     * \param[in] e index of subentity in all subentities with codimension \arg c
     * \param[in] k index of shape function on its subentity
     * \param[in] vGlobal the global vertex indices of the reference simplex's vertices
     */
    std::tuple<int,int> sfPermutation(int d, int p, int c, int e, int k, int const vGlobal[]);


    template <class ctype, int dimension, class Scalar, bool restricted>
    struct ChooseShapeFunctionSet
    {
      typedef HierarchicSimplexShapeFunctionSet<ctype,dimension,Scalar> type;
    };

    template <class ctype, int dimension, class Scalar>
    struct ChooseShapeFunctionSet<ctype,dimension,Scalar,true>
    {
      typedef RestrictedShapeFunctionSet<HierarchicSimplexShapeFunctionSet<ctype,dimension,Scalar> > type;
    };

    template <class ctype, int dimension, class Scalar, bool restricted>
    struct ChooseExtensionShapeFunctionSet
    {
      typedef HierarchicExtensionSimplexShapeFunctionSet<ctype,dimension,Scalar> type;
    };

    template <class ctype, int dimension, class Scalar>
    struct ChooseExtensionShapeFunctionSet<ctype,dimension,Scalar,true>
    {
      typedef RestrictedShapeFunctionSet<HierarchicExtensionSimplexShapeFunctionSet<ctype,dimension,Scalar> > type;
    };

  } // End of namespace HierarchicDetail

  /**
   * \endcond
   */

  //---------------------------------------------------------------------

  /**
   * \brief Scalar hierarchic shape functions on the unit simplex.
   *
   * These are polynomial shape functions spanning hierarchical
   * polynomial spaces containing (but slightly extening) the full
   * polynomial spaces.
   *
   * The shape functions respect the geometrical symmetries of
   * simplices. In particular, they allow to construct global ansatz
   * spaces where every shape function restricted to a cell is given by
   * either one shape function or by the negative of one shape
   * function. No linear combinations of shape functions are required
   * (see Zumbusch [1]).
   *
   * On a simplex of dimension \f$ d \f$, the shape functions of \em
   * nomial order \f$ p \f$ spanning the set \f$ \Phi_p^d \f$ are
   * denoted by \f[ \phi_p^{d,j}, \quad j=0,\dots,\sharp_p^d-1 \f].
   *
   * The shape functions are constructed in terms of barycentric
   * coordinates \f$ \lambda \f$ as follows. For order \f$ p=0 \f$,
   * there are the \f$ d+1 \f$ shape functions \f[ \phi_0^{d,k} = \lambda_k,
   * \quad k=0,\dots, d. \f] For \f$ p=1 \f$, there are no shape
   * functions.
   *
   * For higher order, the shape functions are defined depending on the
   * simplex subentity to which they are associated. On vertices, no
   * additional shape functions are defined.
   *
   * On edges (codimension \f$ c = d-1 \f$) between the vertices \f$
   * j\f$ and \f$ k \f$, the only shape function of order \f$ p \ge 2
   * \f$ is given by \f[ \phi_p^{d,c,j,k} = \lambda_j \lambda_k
   * T_{p-2}(\lambda_j - \lambda_k), \f] where the \f$ T_p \f$ are
   * orthogonal polynomials of order \f$ p \f$, orthogonal w.r.t. the
   * scalar product \f[ (f,g) = \int_{-1}^1 ((1+x)(1-x)f)'
   * ((1+x)(1-x)g)'\,dx. \f] Note that this differs from [1], where the
   * shape functions are defined as \f$ (\lambda_j - \lambda_k)^p
   * \f$. The definition used here should lead to better conditioned
   * stiffness matrices (without the need to resort to redefining shape
   * functions by linear combinations as proposed in [1], see also [2]).
   *
   * On lower codimensional subentities (codimension \f$ c<d-1 \f$)
   * spanned by the vertices \f$ K = k_0, \dots, k_{d-c}\f$, the shape
   * functions are defined recursively by \f[ \phi_p^{d,c,K,j} =
   * \lambda_{k_0}\dots\lambda_{k_{d-c}} \phi, \quad
   * \phi\in \Phi_{p-(d-c+1)}^{d-c}. \f]
   *
   * Note that due to the definition of nominally zeroth order shape
   * functions as affine functions, the set of shape functions of orders
   * \f$ 0,\dots, p \f$ spans a space that can be slightly larger than
   * the space \f$ \mathbf{P}_p \f$ of polynomials of order up to \f$ p
   * \f$.
   *
   * \tparam ctype scalar type for coordinates
   * \tparam dimension spatial dimension
   * \tparam Scalar the value type (a scalar real number type)
   *
   * References:
   *
   * [1] Zumbusch, Grid.: Proc. of the Third Int. Conf. on Spectral and
   * High Order Methods, ICOSAHOM '95, Houston, Texas, June 5-9, 1995,
   * A. V. Ilin, L. R. Scott (eds.), Houston Journal of Mathematics
   * 1996, pp. 529-540
   *
   * [2] Solin, P. Segeth, K., Dolezel, I.: Higher-Order Finite Element
   * Methods. Chapman & Hall 2004.
   */
  template <class ctype, int dimension, class Scalar>
  class HierarchicSimplexShapeFunction: public ShapeFunction<ctype,dimension,Scalar>
  {
  public:
    static int const dim = dimension;
    static int const comps = 1;

    typedef ctype CoordType;
    typedef Scalar ResultType;

    /**
     * Default constructor. Constructs the constant shape function. Not
     * particularly useful, but it's often convenient to have a default
     * constructible class.
     */
    HierarchicSimplexShapeFunction() {}

    /**
     * \brief Creates the k-th shape function of given order.
     */
    explicit HierarchicSimplexShapeFunction(int order, int k):
                    order_(order), number_(k)
    {
      assert(order>=0);
      assert(0<=k && k<HierarchicDetail::size(dim,order));
      std::tie(codim_,entity_,entityIndex_) = HierarchicDetail::tupleIndex(dim,order,k);
      actualOrder = HierarchicDetail::actualOrder(dim,order_,k);
    }

    HierarchicSimplexShapeFunction(HierarchicSimplexShapeFunction const& other)
    : actualOrder(other.actualOrder), order_(other.order_), number_(other.number_),
      codim_(other.codim_), entity_(other.entity_), entityIndex_(other.entityIndex_)
    {}

    virtual  std::unique_ptr<ShapeFunction<ctype,dim,Scalar>> clone() const
    {
      return std::unique_ptr<ShapeFunction<ctype,dim,Scalar>>(new HierarchicSimplexShapeFunction(*this));
    }

    /**
     * Evaluates the shape function at point x. The point x is not
     * restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldVector<ResultType,comps> evaluateFunction(Dune::FieldVector<CoordType,dim> const& x) const
    {
      Dune::FieldVector<CoordType,dim+1> zeta = barycentric(x);
      return HierarchicDetail::shapeFunction(dim,order_,number_,&zeta[0]);
    }

    /**
     * \brief Evaluates the derivative of the shape function.
     *
     * The point \arg x is not restricted to be inside the unit simplex, but
     * the meaning of evaluating a shape function outside of the unit
     * simplex is questionable.
     */
    virtual Dune::FieldMatrix<Scalar,1,dim> evaluateDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      Dune::FieldVector<double,dim+1> zeta = barycentric(x);
      Dune::FieldVector<double,dim+1> bgrad;
      HierarchicDetail::dShapeFunction(dim,order_,number_,&zeta[0],&bgrad[0]);
      
      Dune::FieldMatrix<Scalar,1,dim> grad;
      for (int i=0; i<dim; ++i)
        grad[0][i] = bgrad[i]-bgrad[dim];
      
      
        // for (int i=0; i<dim; ++i) {
        //   double h = 1e-5;
        //   Dune::FieldVector<CoordType,dim> yl(x), yr(x);
        //   yl[i] -= h;
        //   yr[i] += h;
        //   double diff = (evaluateFunction(yr)-evaluateFunction(yl))/(2*h);
        //
        //   if (std::abs(diff-grad[0][i]) > 1e-4) {
          //     std::cout << "i=" << i << ": numdiff=" << diff << "  bdiff=" << grad[0][i] << '\n';
        //   }
        //   grad[0][i] = diff;
        // }
        
        
      return grad;
    }

    /**
     * \brief Evaluates the second derivative of the shape function.
     *
     * The point \arg x is not restricted to be inside the unit simplex, but
     * the meaning of evaluating a shape function outside of the unit
     * simplex is questionable.
     */
    virtual Tensor3<Scalar,1,dim,dim> evaluate2ndDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      static bool visited = false;
      if (!visited) 
      {
        visited = true;
        std::cerr << "Not implemented: Hierarchic shape function 2nd derivative at " << __FILE__ << ":" << __LINE__ << "\n";
      }
      return 0;
    }

    /**
     * \brief Information about the association of the shape function to
     * subentities of the simplex.
     *
     * Returns a tuple \f$ (p,c,e,k) \f$, where \f$ p \f$ is the nominal
     * order of the shape function, \f$ c \f$ is the codimension of the
     * associated subentity, \f$ e \f$ is the number of the subentity
     * (within the subentities of codimension \f$ c \f$), and \f$ k \f$
     * is the number of the shape function (within all shape functions
     * associated to that subentity).
     */
    virtual std::tuple<int,int,int,int> location() const
    {
      return std::make_tuple(order_,codim_,entity_,entityIndex_);
    }

    /**
     * \brief The actual polynomial order of the shape function.
     *
     * This can be the nominal order \f$ p \f$ or \f$ p+1 \f$.
     */
    int order() const { return actualOrder; }


  private:
    int actualOrder;

    int order_;
    int number_;
    int codim_;
    int entity_;
    int entityIndex_;
  };


  //---------------------------------------------------------------------

  /**
   * \brief A container of hierarchic shape functions (see \ref
   * HierarchicSimplexShapeFunction) up to a given nominal order on the
   * unit simplex of grid dimension.
   *
   * For Ord<0, this set is empty.
   */
  template <class ctype, int dimension, class Scalar>
  class HierarchicSimplexShapeFunctionSet : public ShapeFunctionSet<ctype,dimension,Scalar>
  {
  public:
    typedef HierarchicSimplexShapeFunction<ctype,dimension,Scalar> value_type;
    typedef typename ShapeFunctionSet<ctype,dimension,Scalar>::Matrix Matrix;

    /**
     * \brief Constructs a hierarchic shape function set.
     *
     * The shape functions span the set \f$ \bigcup_{p=0}^{\mathrm{ord}}
     * \Phi_p^d \f$ (see \ref * HierarchicSimplexShapeFunction).
     */
    HierarchicSimplexShapeFunctionSet(int ord)
     : ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dimension)),
       pinv(nullptr)
    {
      this->order_ = -1;

      // Empty dummy set.
      if (ord<0) {
        this->size_ = 0;
        return;
      }

      for (int o=0; o<=ord; ++o) {
        int n = HierarchicDetail::size(dimension,o);
        for (int i=0; i<n; ++i) {
          sf.push_back(value_type(o,i));
          this->order_ = std::max(this->order_,sf.back().order());
        }
      }

      // Use equidistant "interpolation" nodes. Note that we may have
      // more shape functions than nodes with the same order. To enable
      // a faithful reconstruction during grid transfer, we use
      // interpolation nodes for the maximal ACTUAL order of the shape
      // functions and a least squares "interpolation" approach.
      //
      // TODO: Evaluate the alternative of using interpolation nodes for
      // the maximal *NOMINAL* order of the shape functions and a least
      // coefficient l2-norm interpolation.  This does not allow a
      // faithful reconstruction (except for the lower dimensional
      // subspace of nominal order), but is faster due to the lower
      // number of interpolation points. Note that for some applications
      // the exact reproduction is necessary.
      int n = binomial(dimension+this->order_,dimension);
      this->iNodes.insert(this->iNodes.end(),n,Dune::FieldVector<ctype,dimension>());
      Dune::FieldVector<int,dimension> b(0);
      for (int i=0; i<n; ++i) {
        for (int j=0; j<dimension; ++j)
          this->iNodes[i][j] = static_cast<ctype>(b[j])/this->order_;
        next_multiindex(b.begin(),b.end(),this->order_);
      }

      // Compute the pseudoinverse for use in interpolation.
      computePseudoInverse();
//      pinv.reset(new SLAPMatrix<Scalar>(sf.size(),n));
//      SLAPMatrix<Scalar> A(n,sf.size());
//      for (int j=0; j<sf.size(); ++j) {
//        for (int i=0; i<n; ++i) {
//          A(i,j) = sf[j].evaluateFunction(this->iNodes[i]);
//        }
//      }
//
//      pseudoinverse(A,*pinv);

      this->size_ = sf.size();
    }

    HierarchicSimplexShapeFunctionSet(HierarchicSimplexShapeFunctionSet const& other) : ShapeFunctionSet<ctype,dimension,Scalar>(other),
      sf(other.sf), pinv(new SLAPMatrix<Scalar>(*other.pinv.get()))
    {}

    HierarchicSimplexShapeFunctionSet(HierarchicSimplexShapeFunctionSet&& other)
    : ShapeFunctionSet<ctype,dimension,Scalar>(other),
      sf(std::move(other.sf)), pinv(std::move(other.pinv))
    {}

    virtual ~HierarchicSimplexShapeFunctionSet(){}

    virtual value_type const& operator[](int i) const{ return sf[i]; }

    virtual Dune::GeometryType type() const { return Dune::GeometryType(Dune::GeometryType::simplex,dimension); }

    virtual void interpolate(typename ShapeFunctionSet<ctype,dimension,Scalar>::SfValueArray const& A,
        Matrix& IA) const
    {
      IA.setSize(sf.size(),A.M());
      if(pinv==nullptr) return;

      if (sf.size() > 0) {
//        std::cout << A.N() << ", " << A.M() << ", " << pinv->rows() << ", " << pinv->cols() << std::endl;
        assert(A.N()==pinv->cols());

        // Apply the pseudoinverse. TODO: Use BLAS level 3.
        for (int i=0; i<IA.N(); ++i)
          for (int j=0; j<IA.M(); ++j) {
            double r = 0;
            for (int k=0; k<A.N(); ++k)
              r += (*pinv)(i,k) * A[k][j];
            IA[i][j] = r;
          }
      }
    }

    void removeShapeFunction(size_t index)
    {
      assert(index>=0 && index < sf.size());
      sf.erase(sf.begin()+index);
      this->size_ = sf.size();
      computePseudoInverse();
    }

  protected:
    std::vector<value_type> sf;
    std::unique_ptr<SLAPMatrix<Scalar> > pinv;

  private:
    void computePseudoInverse()
    {
      pinv.reset(new SLAPMatrix<Scalar>(sf.size(),this->iNodes.size()));
      SLAPMatrix<Scalar> A(this->iNodes.size(),sf.size());
      for (int j=0; j<sf.size(); ++j) {
        for (int i=0; i<this->iNodes.size(); ++i) {
          A(i,j) = sf[j].evaluateFunction(this->iNodes[i]);
        }
      }

      pseudoinverse(A,*pinv);
    }
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  template <class ctype, int dimension, class Scalar, bool restricted = false>
  class HierarchicShapeFunctionSetContainer: public ShapeFunctionSetContainer<ctype,dimension,Scalar>
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,Scalar> value_type;

  private:
    typedef std::map<std::pair<Dune::GeometryType,int>,value_type const*> Container;
    typedef typename HierarchicDetail::ChooseShapeFunctionSet<ctype,dimension,Scalar,restricted>::type ShapeFunctionType;

  public:

    HierarchicShapeFunctionSetContainer(int maxOrder)
    {
      // Fill the container with shape function sets. This is not done
      // on demand in operator(), since then operator() is not easily
      // and efficiently made thread-safe.
      //
      // Currently only shape functions for simplices are defined.
      Dune::GeometryType simplex(Dune::GeometryType::simplex,dimension);

      value_type* sfsPrevious = 0;

      for (int i=-1; i<=maxOrder; ++i) {
        value_type* sfs = new ShapeFunctionType(i);
        sfs->initHierarchicalProjection(sfsPrevious);
        container.insert(typename Container::value_type(std::make_pair(simplex,i),sfs));
        sfsPrevious = sfs;
      }
    }

    virtual ~HierarchicShapeFunctionSetContainer()
    {
      for (typename Container::iterator i=container.begin(); i!=container.end(); ++i)
        delete i->second;
    }

    virtual value_type const& operator()(Dune::GeometryType type, int order) const
    {
      typename Container::const_iterator i = container.find(std::make_pair(type,order));
      assert(i!=container.end()); // should throw, not assert!
      return *i->second;
    }

  private:
    // A map storing Hierarchic shape function sets associated with
    // (reference element type, order). Due to runtime polymorphism,
    // pointers are stored in the map. This class owns the shape
    // function sets pointed to, and has to delete them in the
    // destructor.
    Container container;
  };


  /**
   * \brief Returns a Hierarchic shape function set for given reference element
   * type and given polynomial order. Singleton function.
   */
  template <class ctype, int dimension, class Scalar>
  typename HierarchicShapeFunctionSetContainer<ctype,dimension,Scalar>::value_type const&
  hierarchicShapeFunctionSet(Dune::GeometryType type, int order)
  {
    // This static variable will be destructed on program shutdown,
    // causing its contents to be deleted properly from the heap by
    // the destructor. No memory leaks should be induced.
    static HierarchicShapeFunctionSetContainer<ctype,dimension,Scalar> container(7);

    return container(type,order);
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  /**
   * \brief A container of hierarchic shape functions of nominal order Ord on
   * the unit simplex of grid dimension. For Ord<0, this set is empty.
   *
   * TODO: Common functionality of
   * HierarchicExtensionSimplexShapeFunctionSet and
   * HierarchicSimplexShapeFunctionSet should be factored out into a
   * base class.
   */
  template <class ctype, int dimension, class Scalar>
  class HierarchicExtensionSimplexShapeFunctionSet: public ShapeFunctionSet<ctype,dimension,Scalar>
  {
  public:
    typedef HierarchicSimplexShapeFunction<ctype,dimension,Scalar> value_type;
    typedef typename ShapeFunctionSet<ctype,dimension,Scalar>::Matrix Matrix;

    HierarchicExtensionSimplexShapeFunctionSet(int ord):
      ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dimension))
      {
      this->order_ = -1;

      // Empty dummy set.
      if (ord<0) {
        this->size_ = 0;
        return;
      }

      int n = HierarchicDetail::size(dimension,ord);
      if (n==0) {
        this->size_ = 0;
        return;
      }

      for (int i=0; i<n; ++i) {
        sf.push_back(value_type(ord,i));
        this->order_ = std::max(this->order_,sf.back().order());
      }

      // Use equidistant "interpolation" nodes. To enable a faithful
      // reconstruction during grid transfer, we use interpolation nodes
      // for the maximal ACTUAL order of the shape functions and a least
      // squares "interpolation" approach.
      //
      // TODO: Since this is only an *extension*, we usually have fewer
      // shape functions. Check whether this can be exploited by using
      // fewer interpolation nodes.
      int m = binomial(dimension+this->order_,dimension);
      this->iNodes.resize(m);
      Dune::FieldVector<int,dimension> b(0);
      for (int i=0; i<m; ++i) {
        for (int j=0; j<dimension; ++j)
          this->iNodes[i][j] = static_cast<ctype>(b[j])/this->order_;
        next_multiindex(b.begin(),b.end(),this->order_);
      }

      // Compute the pseudoinverse for use in interpolation.
      pinv.reset(new SLAPMatrix<Scalar>(sf.size(),m));
      SLAPMatrix<Scalar> A(m,sf.size());
      for (int j=0; j<sf.size(); ++j) {
        for (int i=0; i<m; ++i) {
          A(i,j) = sf[j].evaluateFunction(this->iNodes[i]);
        }
      }

      pseudoinverse(A,*pinv);

      this->size_ = sf.size();
    }

    HierarchicExtensionSimplexShapeFunctionSet(HierarchicExtensionSimplexShapeFunctionSet const& other) : ShapeFunctionSet<ctype,dimension,Scalar>(other),
      sf(other.sf), pinv(new SLAPMatrix<Scalar>(*other.pinv.get()))
    {}

    HierarchicExtensionSimplexShapeFunctionSet(HierarchicExtensionSimplexShapeFunctionSet&& other)
    : ShapeFunctionSet<ctype,dimension,Scalar>(other),
      sf(std::move(other.sf)), pinv(std::move(other.pinv))
    {}

    virtual value_type const& operator[](int i) const{ return sf[i]; }

    virtual Dune::GeometryType type() const { return Dune::GeometryType(Dune::GeometryType::simplex,dimension); }

    virtual void interpolate(typename ShapeFunctionSet<ctype,dimension,Scalar>::SfValueArray const& A,
        Matrix& IA) const
    {
      IA.setSize(sf.size(),A.M());

      if (sf.size() > 0) {
        assert(A.N()==pinv->cols());

        // Apply the pseudoinverse. TODO: Use BLAS level 3.
        for (int i=0; i<IA.N(); ++i)
          for (int j=0; j<IA.M(); ++j) {
            double r = 0;
            for (int k=0; k<A.N(); ++k)
              r += (*pinv)(i,k) * A[k][j];
            IA[i][j] = r;
          }
      }
    }

    void removeShapeFunction(size_t index)
    {
      assert(index>=0 && index < sf.size());
      sf.erase(sf.begin()+index);
      this->size_ = sf.size();
    }
  protected:
    std::vector<value_type> sf;
    std::unique_ptr<SLAPMatrix<Scalar> > pinv;
  };

  //---------------------------------------------------------------------

  template <class ctype, int dimension, class Scalar, bool restricted = false>
  class HierarchicExtensionShapeFunctionSetContainer: public ShapeFunctionSetContainer<ctype,dimension,Scalar>
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,Scalar> value_type;

  private:
    typedef std::map<std::pair<Dune::GeometryType,int>,value_type const*> Container;
    typedef typename HierarchicDetail::ChooseExtensionShapeFunctionSet<ctype,dimension,Scalar,restricted>::type ShapeFunctionType;

  public:

    HierarchicExtensionShapeFunctionSetContainer(int maxOrder)
    {
      // Fill the container with shape function sets. This is not done
      // on demand in operator(), since then operator() is not easily
      // and efficiently made thread-safe.
      //
      // Currently only shape functions for simplices are defined.
      Dune::GeometryType simplex(Dune::GeometryType::simplex,dimension);

      value_type* sfsPrevious = 0;

      for (int i=-1; i<=maxOrder; ++i) {
        value_type* sfs = new ShapeFunctionType(i);
        sfs->initHierarchicalProjection(sfsPrevious);
        container.insert(typename Container::value_type(std::make_pair(simplex,i),sfs));
        sfsPrevious = sfs;
      }
    }

    virtual ~HierarchicExtensionShapeFunctionSetContainer()
    {
      for (typename Container::iterator i=container.begin(); i!=container.end(); ++i)
        delete i->second;
    }

    virtual value_type const& operator()(Dune::GeometryType type, int order) const
    {
      typename Container::const_iterator i = container.find(std::make_pair(type,order));
      assert(i!=container.end()); // should throw, not assert!
      return *i->second;
    }

  private:
    // A map storing Hierarchic shape function sets associated with
    // (reference element type, order). Due to runtime polymorphism,
    // pointers are stored in the map. This class owns the shape
    // function sets pointed to, and has to delete them in the
    // destructor.
    Container container;
  };


  /**
   * Returns a Hierarchic shape function set for given reference element
   * type and given polynomial order. Singleton function.
   */
  template <class ctype, int dimension, class Scalar>
  typename HierarchicExtensionShapeFunctionSetContainer<ctype,dimension,Scalar>::value_type const&
  hierarchicExtensionShapeFunctionSet(Dune::GeometryType type, int order)
  {
    // This static variable will be destructed on program shutdown,
    // causing its contents to be deleted properly from the heap by
    // the destructor. No memory leaks should be induced.
    static HierarchicExtensionShapeFunctionSetContainer<ctype,dimension,Scalar> container(7);

    return container(type,order);
  }

} // End of namespace Kaskade

#endif
