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

#ifndef NEDELECSHAPEFUNCTIONS_HH
#define NEDELECSHAPEFUNCTIONS_HH

#include <map>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include "fem/barycentric.hh"
#include "fem/lagrangeshapefunctions.hh"
#include "fem/pshapefunctions.hh"

#include "linalg/dynamicMatrix.hh"

namespace Kaskade
{
  /**
   * \brief Vectorial Nedelec shape functions on the unit simplex. 
   * 
   * These are linear shape functions associated to edges. In barycentric coordinates \f$ \lambda \f$, 
   * with \f$ \lambda_i = x_i \f$ for \f$ i=0,\dots, \mathrm{dim}-1 \f$ and \f$ \lambda_\mathrm{dim} = 1-\sum_{i=0}^{\mathrm{dim}-1} x_i \f$ 
   * on the unit simplex, the shape function associated to the edge \f$ e = (p_i,p_j) \f$ is
   * \f[  \phi_e = \pm (\lambda_i \nabla \lambda_j - \lambda_j \nabla \lambda i). \f]
   * Here, \f$ p_k \f$ is the vertex with \f$ \lambda_k = 1 \f$. The sign of the shape function depends on
   * the numbering of vertices relative to the edge implemented by the grid manager.
   *
   * \tparam ctype the coordinate type (a real number type)
   * \tparam T the value type (a real number type)
   * \tparam D the dimension of the unit simplex over which the shape functions are defined
   */
  template <class ctype, int dimension, class T>
  class NedelecSimplexShapeFunction: public ShapeFunction<ctype,dimension,T,dimension>
  {
  public:
    static int const dim = dimension;
    static int const comps = dim;
    static int const order = 1;

    typedef ctype CoordType;
    typedef T ResultType;

    /**
     * Default constructor. Constructs the shape function associated
     * with the edge 0. Not particularly useful, but it's often
     * convenient to have a default constructible class.
     */
    NedelecSimplexShapeFunction()
    {
      *this = NedelecSimplexShapeFunction<ctype,dimension,T>(0);
    }

    /**
     * Creates a shape function associated to the edge e.
     *
     * @param e subentity index of the edge (codim=dim-1 subentity index)
     */
    explicit NedelecSimplexShapeFunction(int e)
    {
      // compute index local in element
      entity_ = e;

      // Compute barycentric dimensions of edge end points. Note that the origin, which 
      // has barycentric index dim, has the subentity number 0. Hence we need to shift and
      // wrap the subentity indices to obtain the barycentric indices.
      Dune::ReferenceElement<CoordType,dim> const &simplex = Dune::ReferenceElements<CoordType,dim>::simplex();
      int p0Id = simplex.subEntity(e,dim-1,0,dim);
      p0_ = p0Id==0? dim: p0Id-1;
      int p1Id = simplex.subEntity(e,dim-1,1,dim);
      p1_ = p1Id==0? dim: p1Id-1;

      // compute gradients of barycentric coordinate functions
      if (p0_==dim) dp0_ = -1;
      else {        dp0_ = 0; dp0_[p0_] = 1; }
      if (p1_==dim) dp1_ = -1;
      else {        dp1_ = 0; dp1_[p1_] = 1; }

      // compute position and direction of edge
      pos_ = simplex.position(e,dim-1);
      edge_ = simplex.position(p1Id,dim) - simplex.position(p0Id,dim);
    }
    
    virtual std::unique_ptr<ShapeFunction<ctype,dim,T,dim>> clone() const
    {
      return std::unique_ptr<ShapeFunction<ctype,dim,T,dim>>(new NedelecSimplexShapeFunction(*this));
    }

    /**
     * \brief Evaluates the shape function at point x. 
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldVector<T,dim> evaluateFunction(Dune::FieldVector<ctype,dim> const& x) const
    {
      Dune::FieldVector<CoordType,dim+1> bc = barycentric(x);
      return bc[p0_]*dp1_ - bc[p1_]*dp0_;
    }

    /**
     * \brief Evaluates the derivative of the shape function. 
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldMatrix<ResultType,dim,dim> evaluateDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      // derivative of bc[p0]*dp1 - bc[p1]*dp0 is dp1 * bc[p0]' - dp0 * bc[p1]', where the products are *outer products*
      // resulting in the proper matrices. Start with derivatives of the barycentric coordinates \lambda_i and \lambda_j
      Dune::FieldVector<T,dim> dbc0(0.0), dbc1(0.0);
      for (int dir=0; dir<dim; ++dir) 
      {
        if (p0_ == dir) dbc0 = -1.0;
        else            dbc0[p0_] = 1.0;
        if (p1_ == dir) dbc1 = -1.0;
        else            dbc1[p1_] = 1.0;
      }
      
      return outerProduct(dp1_,dbc0) - outerProduct(dp0_,dbc1);
      
      // old code, should be equivlaent...
//       // TODO: faster implementation needed
//       Dune::FieldMatrix<ResultType,dim,dim> d;
//       for (int c=0; c<dim; ++c)
//         for (int dir=0; dir<dim; ++dir) {
//           ResultType dbcp0_ddir=0, dbcp1_ddir=0;
//           if (p0_==dir) dbcp0_ddir = 1;
//           if (p0_==dim) dbcp0_ddir = -1;
//           if (p1_==dir) dbcp1_ddir = 1;
//           if (p1_==dim) dbcp1_ddir = -1;
// 
//           d[c][dir] = dbcp0_ddir*dp1_[c] - dbcp1_ddir*dp0_[c];
//         }
// 
//       return d;
    }

    /**
     * \brief Evaluates the second derivative of the shape function. 
     * 
     * The result is r[i][j][k] = \f$ \partial^2 \phi_i / (\partial x_j \partial x_k) = 0. \f$
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Tensor3<ResultType,1,dim,dim> evaluate2ndDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      // Nedelec shape functions are linear -- the second derivative vanishes.
      return Tensor3<ResultType,1,dim,dim>(0.0);
    }

    virtual std::tuple<int,int,int,int> location() const
    {
      return std::make_tuple(1,dim-1,entity_,0);
    }

    /**
     * \brief Returns the element-local position of the node associated to the shape function.
     */
    Dune::FieldVector<CoordType,dim> position() const { return pos_; }

    /**
     * \brief Returns the edge vector of the edge associated to this shape function.
     */
    Dune::FieldVector<CoordType,dim> edge() const { return edge_; }

  private:
    int entity_;
    // p0_ and p1_ contain the barycentric coordinate into which the
    // "from" and "to" points of the edge lie, respectively.
    int p0_, p1_;
    Dune::FieldVector<ResultType,dim> dp0_, dp1_;
    Dune::FieldVector<CoordType,dim> pos_, edge_; 
  };

  //--------------------------------------------------

  template <class ctype, int dimension, class T>
  class NedelecSimplexShapeFunctionSet: public ShapeFunctionSet<ctype,dimension,T,dimension>
  {
    static int const dim = dimension;
    static int const noEdges = dim*(dim+1)/2; // number of edges


  public:
    typedef NedelecSimplexShapeFunction<ctype,dim,T> value_type;

    NedelecSimplexShapeFunctionSet():
      ShapeFunctionSet<ctype,dim,T,dim>(Dune::GeometryType(Dune::GeometryType::simplex,dim))
      {
     auto const& simplex = Dune::ReferenceElements<ctype,dim>::simplex();
      assert(noEdges == simplex.size(dim-1));

      sf.reserve(noEdges);
      for (int i=0; i<noEdges; ++i)
        sf.push_back(value_type(i));

      // The interpolation uses the fact that for any edge vector e_i of
      // the reference simplex and any shape function phi_j evaluated on
      // this edge it holds that <e_i,phi_j>=0 if i!=j. Thus we first
      // store the tangential component <e_i,phi_i> for all shape
      // functions.
      this->iNodes.reserve(noEdges);
      tangentialComponents.reserve(noEdges);
      for (size_t i=0; i<noEdges; ++i) {
        this->iNodes.push_back(sf[i].position());
        tangentialComponents.push_back(sf[i].edge()*sf[i].evaluateFunction(sf[i].position()));
      }

      this->order_ = 1;
      this->size_ = sf.size();
      }

    virtual value_type const& operator[](int i) const{ return sf[i]; }
    virtual Dune::GeometryType type() const { return Dune::GeometryType(Dune::GeometryType::simplex,dim); }

    /**
     * @brief Left-multiplies the provided matrix with the interpolation
     * matrix of the shape function set.
     *
     * Each column of A is interpreted as values of some function
     * evaluated at this shape function set's interpolation points (edge midpoints). The
     * columns of the output matrix IA then contain the shape function
     * coefficients such that the corresponding linear combination of
     * shape functions "interpolates" that function in the interpolation
     * points. "Interpolation" means equal tangential components.
     */
    virtual void interpolate(typename ShapeFunctionSet<ctype,dim,T,dim>::SfValueArray const& A,
                             DynamicMatrix<Dune::FieldMatrix<T,1,1> >& IA) const
    {
      assert(A.N()==noEdges);

      IA.setSize(noEdges,A.M());
      for (int i=0; i<noEdges; ++i)
        for (int j=0; j<A.M(); ++j)
          IA[i][j] = sf[i].edge()*asVector(A[i][j])/tangentialComponents[i];
    }


  private:
    std::vector<value_type> sf;
    std::vector<T>          tangentialComponents;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  template <class ctype, int dimension, class T>
  class NedelecShapeFunctionSetContainer: public ShapeFunctionSetContainer<ctype,dimension,T,dimension>
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,T,dimension> value_type;

  private:
    typedef std::map<std::pair<Dune::GeometryType,int>,value_type const*> Container;

  public:

    NedelecShapeFunctionSetContainer()
    {
      // Fill the container with shape function sets. This is not done
      // on demand in operator(), since then operator() is not easily
      // and efficiently made thread-safe.
      //
      // Currently only shape functions for simplices are defined.
      Dune::GeometryType simplex(Dune::GeometryType::simplex,dimension);
      container.insert(typename Container::value_type(std::make_pair(simplex,1),new NedelecSimplexShapeFunctionSet<ctype,dimension,T>));
    }


    virtual ~NedelecShapeFunctionSetContainer()
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
    // A map storing Nedelec shape function sets associated with
    // (reference element type, order). Due to runtime polymorphism,
    // pointers are stored in the map. This class owns the shape
    // function sets pointed to, and has to delete them in the
    // destructor.
    Container container;
  };


  /**
   * Returns a Nedelec shape function set for given reference element
   * type and given polynomial order. Singleton function.
   */
  template <class ctype, int dimension, class T>
  typename NedelecShapeFunctionSetContainer<ctype,dimension,T>::value_type const& nedelecShapeFunctionSet(Dune::GeometryType type, int order)
  {
    // This static variable will be destructed on program shutdown,
    // causing its contents to be deleted properly from the heap by
    // the destructor. No memory leaks should be induced.
    static NedelecShapeFunctionSetContainer<ctype,dimension,T> container;
    return container(type,order);
  }
  
  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  
  /**
   * \brief Vectorial \f$ H(\text{div}) \f$ conforming shape functions on the unit simplex. 
   * 
   * These form a basis of \f$ \mathbb{P}_p^d \f$ of the space of vectorial polynomials of degree at
   * most \f$ p, p\ge 1 \f$. A shape function is represented as a linear combination of Lagrangian
   * shape functions.
   *
   * \tparam ctype the coordinate type (a real number type)
   * \tparam dimension the dimension of the unit simplex over which the shape functions are defined
   * \tparam T the value type (a real number type)
   */
  template <class ctype, int dimension, class T>
  class HdivSimplexShapeFunction: public ShapeFunction<ctype,dimension,T,dimension>
  {
  public:
    static int const dim = dimension;
    static int const comps = dim;
    static int const order = 1;

    typedef ctype CoordType;

    /**
     * \brief Default constructor. Constructs a vanishing "shape function". Not particularly useful, but it's often
     * convenient to have a default constructible class.
     * \todo makes only sense if it is assignable - any real need for default constructibility?
     */
    HdivSimplexShapeFunction() = default;

    /**
     * \brief Creates a shape function as a linear combination of Lagrangian shape functions.
     * \tparam Coefficients an STL container with value type Dune::FieldVector<T,comps>.
     * \param coeff contains coefficients of linear combinations of Lagrangian shape functions. Note that the coefficients are vectorial,
     *        as the Lagrangian shape functions are scalar, but we create vectorial shape functions.
     * \param loc the topological location information (nominalOrder,codim,entity,index)
     */
    template <class Coefficients>
    explicit HdivSimplexShapeFunction(Coefficients const& coeff_, std::tuple<int,int,int,int> loc_)
    : coeff(coeff_), loc(loc_), lsfs(lagrangeShapeFunctionSet<ctype,dimension>(Dune::GeometryType(Dune::GeometryType::simplex,dim),std::get<0>(loc)))
    {
      assert(coeff.size() == lsfs.size());
    }
    
    virtual std::unique_ptr<ShapeFunction<ctype,dim,T,dim>> clone() const
    {
      return std::unique_ptr<ShapeFunction<ctype,dim,T,dim>>(new HdivSimplexShapeFunction(*this));
    }

    /**
     * \brief Evaluates the shape function at point xi. 
     * 
     * The point xi is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldVector<T,dim> evaluateFunction(Dune::FieldVector<ctype,dim> const& xi) const
    {
      Dune::FieldVector<T,dim> phi(0);
      for (int i=0; i<coeff.size(); ++i)
          phi += lsfs[i].evaluateFunction(xi)[0] * coeff[i];
      return phi;
    }

    /**
     * \brief Evaluates the derivative of the shape function. 
     * 
     * The point xi is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is questionable.
     */
    virtual Dune::FieldMatrix<T,dim,dim> evaluateDerivative(Dune::FieldVector<CoordType,dim> const& xi) const
    {
      Dune::FieldMatrix<T,dim,dim> dphi(0);
      for (int i=0; i<coeff.size(); ++i)
      {
        auto lsf = lsfs[i].evaluateDerivative(xi)[0]; // has only one component
        for (int j=0; j<dim; ++j)
          dphi[j] +=  coeff[i][j] * lsf;
      }
      return dphi;
    }

    /**
     * \brief Evaluates the second derivative of the shape function. 
     * 
     * The result is r[i][j][k] = \f$ \partial^2 \phi_i / (\partial \xi_j \partial \xi_k) = 0. \f$
     * 
     * The point xi is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is questionable.
     */
    virtual Tensor3<T,dim,dim,dim> evaluate2ndDerivative(Dune::FieldVector<CoordType,dim> const& xi) const
    {
      Tensor3<T,dim,dim,dim> ddphi(0);
      for (int i=0; i<coeff.size(); ++i)
      {
        auto lsf = lsfs[i].evaluate2ndDerivative(xi)[0];// has only one component
        for (int j=0; j<dim; ++j)
          ddphi[j] +=  coeff[i][j] * lsf;
      }
      return ddphi;
    }

    virtual std::tuple<int,int,int,int> location() const
    {
      return loc;
    }

  private:
    std::vector<Dune::FieldVector<T,comps>>   coeff;  // coefficients of Lagrangian shape functions.
    std::tuple<int,int,int,int> loc;                  // location information
    ShapeFunctionSet<ctype,dim,T> const&  lsfs;       // Lagrangian shape function set
  };

  //--------------------------------------------------

  /**
   * \brief A shape function set of \f$ H(\text{div}) \f$ conforming functions.
   * 
   * For each Lagrangian interpolation point \f$ \xi_i \f$ of each face, there is one shape function with
   * non-vanishing normal component of size 1 exactly at this point. Those shape functions are associated
   * to the corresponding face. The remaining basis functions that complete the polynomial space all have
   * vanishing normal components and are associated to the cell.
   * 
   * Only the normal components at the face interpolation points are specified, the tangential components 
   * are well-defined but unspecified (implementation detail).
   */
  template <class ctype, int dimension, class T=double>
  class HdivSimplexShapeFunctionSet: public ShapeFunctionSet<ctype,dimension,T,dimension>
  {
    static int const dim = dimension; 


  public:
    typedef HdivSimplexShapeFunction<ctype,dim,T> value_type;

    /**
     * \brief Constructs a \f$ H(\text{div}) \f$ conforming shape function set of given polynomial order.
     * \param order the polynomial order of the shape functions (at least 1)
     */
    HdivSimplexShapeFunctionSet(int order);
    
    virtual value_type const& operator[](int i) const{ return sf[i]; }
    
    virtual Dune::GeometryType type() const { return Dune::GeometryType(Dune::GeometryType::simplex,dim); }

    /**
     * \brief Left-multiplies the provided matrix with the interpolation
     * matrix of the shape function set.
     *
     * Each column of A is interpreted as values of some function
     * evaluated at this shape function set's interpolation points (edge midpoints). The
     * columns of the output matrix IA then contain the shape function
     * coefficients such that the corresponding linear combination of
     * shape functions "interpolates" that function in the interpolation
     * points. "Interpolation" means equal tangential components.
     */
    virtual void interpolate(typename ShapeFunctionSet<ctype,dim,T,dim>::SfValueArray const& A,
                             DynamicMatrix<Dune::FieldMatrix<T,1,1> >& IA) const
    {
      std::cerr << "not implemented\n";
      abort(); // to be implemented
    }


  private:
    std::vector<value_type> sf;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  template <class ctype, int dimension, class T>
  class HdivShapeFunctionSetContainer: public ShapeFunctionSetContainer<ctype,dimension,T,dimension>
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,T,dimension> value_type;

  private:
    typedef std::map<std::pair<Dune::GeometryType,int>,value_type const*> Container;

  public:

    HdivShapeFunctionSetContainer()
    {
      // Fill the container with shape function sets. This is not done
      // on demand in operator(), since then operator() is not easily
      // and efficiently made thread-safe.
      //
      // Currently only shape functions for simplices are defined.
      Dune::GeometryType simplex(Dune::GeometryType::simplex,dimension);
      for (int order=1; order<6; ++order)
        container.insert(typename Container::value_type(std::make_pair(simplex,order),new HdivSimplexShapeFunctionSet<ctype,dimension,T>(order)));
    }


    virtual ~HdivShapeFunctionSetContainer()
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
    // A map storing Hdiv shape function sets associated with
    // (reference element type, order). Due to runtime polymorphism,
    // pointers are stored in the map. This class owns the shape
    // function sets pointed to, and has to delete them in the
    // destructor.
    Container container;
  };


  /**
   * \brief Returns a Hdiv shape function set for given reference element
   * type and given polynomial order. Singleton function.
   */
  template <class ctype, int dimension, class T=double>
  typename HdivShapeFunctionSetContainer<ctype,dimension,T>::value_type const& hdivShapeFunctionSet(Dune::GeometryType type, int order)
  {
    // This static variable will be destructed on program shutdown,
    // causing its contents to be deleted properly from the heap by
    // the destructor. No memory leaks should be induced.
    static HdivShapeFunctionSetContainer<ctype,dimension,T> container;
    return container(type,order);
  }
  
} // end of namespace Kaskade


#endif
