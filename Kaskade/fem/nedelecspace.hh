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

#ifndef NEDELECSPACE_HH
#define NEDELECSPACE_HH

#include <numeric>

#include "fem/combiner.hh"
#include "fem/functionspace.hh"
#include "fem/nedelecshapefunctions.hh"
#include "fem/views.hh"

//---------------------------------------------------------------------
namespace Kaskade
{
  
  
  /**
   * \ingroup fem
   * \brief A class mapping local vectorial shape function values and gradients to
   * global shape function values and gradients. 
   * 
   * This converter realizes the transform \f$ \psi(x) \f$ and gives the values and derivatives of global shape functions from local shape
   * functions: \f$ \phi(x) = C \hat \phi(\xi) \f$. Derived classes shall redefine the method update()
   * in order to provide \f$ C \f$.
   * 
   * \tparam Grid the FE grid class
   */
  template <class GridView>
  class VectorialConverterBase
  {
    typedef typename GridView::template Codim<0>::Entity Cell;
    static int const dim = GridView::dimension;
    
  public:
    VectorialConverterBase() = default;
    
    VectorialConverterBase(Cell const& cell): cell_(&cell) {}
    
    void moveTo(Cell const& cell) { cell_ = &cell; }
    
    void setLocalPosition(Dune::FieldVector<typename GridView::ctype,dim> const& x)
    {
      assert(cell_);
      assert(cell_->geometry().affine()); // works only with affine geometries (and hence constant C) - otherwise the spatial derivatives are wrong!
      Btinv = cell_->geometry().jacobianInverseTransposed(x);
      update();
    }
    
    /// Applies the transformation \f$ \psi(x) \f$ to shape function value.
    template <class Scalar>
    Dune::FieldMatrix<Scalar,dim,1> global(Dune::FieldMatrix<Scalar,dim,1> const& sf) const 
    { 
      return C * sf; 
    }
    
    /// Applies the transformation \f$ \psi(x) \f$ to shape function value and derivative.
    template <class Scalar>
    VariationalArg<Scalar,dim,dim> global(std::pair<Dune::FieldVector<Scalar,dim>,Dune::FieldMatrix<Scalar,dim,dim> > const& sf) const
    {
      VariationalArg<Scalar,dim,dim> phi;
      // global value: C*phi
      phi.value = Btinv * sf.first;
      
      // global derivative: C*dphi*inv(B)
      phi.gradient = C * sf.second * transpose(Btinv); 
      
      return phi;
    }
    
    /**
     * \brief Applies the transformation \f$ \psi \f$ to shape function value, derivative, and Hessian, returning a filled VariationalArg.
     */
    template <class Scalar>
    VariationalArg<Scalar,dim,dim> global(VariationalArg<Scalar,dim,dim> const& sf, int deriv) const
    {
      VariationalArg<Scalar,dim,dim> phi;
      // global value: C*phi
      phi.value = C * sf.value;
      
      // global derivative: C*dphi*inv(B)
      if (deriv >= 1) 
        phi.gradient = C * sf.gradient * transpose(Btinv); 
      
      assert(deriv<2 && "not yet implemented"); // not yet implemented
      
      return phi;
    }
    
    /**
     * \brief Applies the inverse transform \f$ \psi^{-1} \f$ to global shape function values, giving the local shape function value.
     */
    template <class Scalar>
    Dune::FieldMatrix<Scalar,dim,1> local(Dune::FieldMatrix<Scalar,dim,1> const& glob) const
    {
      Dune::FieldMatrix<Scalar,dim,1> loc;
      solve(C,loc,glob);
      return loc;
    }
    
    
  protected:
    /**
     * \brief Redefine this to set \f$ C \f$.
     */
    virtual void update() = 0;
    
    // Solves Ax=b. Correct :) implementation of LU decomposition. 
    // Should use Dune version, which is, however, buggy in some versions.
    void solve(Dune::FieldMatrix<typename GridView::ctype,dim,dim> A,
               Dune::FieldVector<typename GridView::ctype,dim>& x,
               Dune::FieldVector<typename GridView::ctype,dim> b) const
    {
      typedef typename GridView::ctype real;
      
      for (int i=0; i<dim-1; ++i) {
        // column pivoting
        real max = std::abs(A[i][i]);
        int midx = i;
        for (int k=i+1; k<dim; ++k)
          if (std::abs(A[k][i])>max) {
            max = std::abs(A[k][i]);
            midx = k;
          }
          // swap rows
          for (int j=i; j<dim; ++j)
            std::swap(A[i][j],A[midx][j]);
          std::swap(b[i],b[midx]);
        // perform elimination
        for (int k=i+1; k<dim; ++k) {
          real p = A[k][i]/A[i][i];
          for (int j=i; j<dim; ++j)
            A[k][j] -= p*A[i][j];
          b[k] -= p*b[i];
        }
      }
      // backsubstitution
      for (int i=dim-1; i>=0; --i) {
        for (int k=i+1; k<dim; ++k)
          b[i] -= A[i][k]*b[k];
        b[i] /= A[i][i];
      }
      x = b;
    }
    
    
    Cell const*                                     cell_;
    Dune::FieldMatrix<typename GridView::ctype,dim,dim> Btinv; // contains TRANSPOSED inverse jacobian
    Dune::FieldMatrix<typename GridView::ctype,dim,dim> C;
  };
  
  //--------------------------------------------------------------------------------------
  
  /**
   * \ingroup fem
   * \brief A class mapping local vectorial shape function values and gradients to
   * global shape function values and gradients in \f$ H(\text{curl}) \f$ conforming FE spaces. 
   * 
   * For vectorial \f$ H(\text{curl}) \f$ conforming finite elements, the values of the shape functions
   * transform as \f[ \phi(x) = B^{-T} \hat \phi(\xi), \f]
   * where \f$ B \f$ is the Jacobian of the affine transformation from reference element to actual element
   * and \f$ J \f$ its determinant. See, e.g., Solin, Segeth, Dolezel: Higher Order Finite Element Methods, Chapman & Hall 2004, p 174.
   * 
   * \tparam GridView the FE grid class
   */
  template <class GridView>
  class HcurlConverter: public VectorialConverterBase<GridView>
  {
  public:
    HcurlConverter() = default;
    
    HcurlConverter(typename GridView::template Codim<0>::Entity const& cell): VectorialConverterBase<GridView>(cell) {}
    
  private:
    virtual void update()
    {
      this->C = this->Btinv;
    }
  };
  
  //--------------------------------------------------------------------------------------
  
  /**
   * \ingroup fem
   * \brief A class mapping local vectorial shape function values and gradients to
   * global shape function values and gradients in \f$ H(\text{div}) \f$ conforming FE spaces. 
   * 
   * For vectorial \f$ H(\text{div}) \f$ conforming finite elements, the values of the shape functions
   * transform as \f[ \phi(x) = B^T J^{-1} \hat \phi(\xi), \f]
   * where \f$ B \f$ is the Jacobian of the affine transformation from reference element to actual element
   * and \f$ J \f$ its determinant. See, e.g., Solin, Segeth, Dolezel: Higher Order Finite Element Methods, Chapman & Hall 2004, p 176.
   * 
   * \tparam GridView the FE grid class
   */
  template <class GridView>
  class HdivConverter: public VectorialConverterBase<GridView>
  {
  public:
    HdivConverter() = default;
    
    HdivConverter(typename GridView::template Codim<0>::Entity const& cell): VectorialConverterBase<GridView>(cell) {}
    
  private:
    virtual void update()
    {
      auto& C = this->C;
      C = this->Btinv;
      C.invert();                      // C = B^T
      C *= this->Btinv.determinant();  // C = B^T * det(B^{-T}) = B^T / det(B^T)
    }
  };

  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------

  /**
   * \ingroup fem
   * \brief Degrees of freedom manager for H(div) conforming elements.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class HdivMapper
  {
  public:
    typedef ScalarType Scalar;

    typedef GV                          GridView;
    typedef typename GridView::Grid     Grid;
    typedef typename GridView::IndexSet IndexSet;

    static int const dim = Grid::dimension;

    /**
     * \brief Whether the ansatz functions have global support (i.e. lead to dense matrices).
     */
    static bool const globalSupport = false;

    typedef int ConstructorArgument;
    static int const continuity = -1;
    typedef RangeView<std::vector<size_t>::const_iterator> GlobalIndexRange;

    typedef std::pair<size_t,int>                  IndexPair;
    typedef std::vector<IndexPair>::const_iterator SortedIndexIterator;
    typedef RangeView<SortedIndexIterator>         SortedIndexRange;
  
    typedef typename Grid::template Codim<0>::Entity                            Cell;
    typedef typename HdivShapeFunctionSetContainer<typename Grid::ctype,dim,Scalar>::value_type ShapeFunctionSet;

    /**
     * \brief Constructs an HdivMapper for a given grid view. 
     * \param order ansatz function order (at least 1).
     */
    HdivMapper(GridView const& gridView_, int order_):
      gridView(gridView_), indexSet(gridView_.indexSet()), order(order_), 
      sfs(hdivShapeFunctionSet<typename Grid::ctype,dim,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex),order))
    {
      auto const& simplex = Dune::GenericReferenceElements<typename Grid::ctype,dim>::simplex();
      update();
    }

    virtual ~HdivMapper() {}

    /**
     * \brief Returns the set of shape functions defined on this cell.
     * \param cell the codim 0 entity of the grid for wich the shape functions are to be retrieved
     * \param contained if true, the method may assume that the cell is contained in the index set of the space.
     *                  (The other case occurs during interpolation between different grids).
     */
    ShapeFunctionSet const& shapefunctions(Cell const& cell, bool /* contained */) const
    {
      assert(cell.type().isSimplex());
      return sfs;
    }

    ShapeFunctionSet const& shapefunctions(size_t /* n */) const
    {
      return sfs;
    }

    /**
     * \brief Returns the maximal polynomial order of shape functions encountered in any cell.
     */
    int maxOrder() const
    {
      return order;
    }

    /**
     * \brief Returns an empty range just for initialization purposes, since RangeView is not default constructible.
     */
    GlobalIndexRange initGlobalIndexRange() const
    {
      return GlobalIndexRange(globIdx[0].begin(),globIdx[0].end());
    }

    /**
     * \brief Returns a const sequence containing the global indices of the shape functions associated to this cell. 
     * Global indices start at 0 and are consecutive.
     */
    GlobalIndexRange globalIndices(Cell const& cell) const
    {
      return globalIndices(indexSet.index(cell));
    }
    
    /**
     * \brief Returns a const sequence containing the global indices of the shape functions associated to the cell with index n. 
     * Global indices start at 0 and are consecutive.
     */
    GlobalIndexRange globalIndices(size_t n) const 
    {
      return GlobalIndexRange(globIdx[n].begin(),globIdx[n].end());
    }

    /**
     * \brief Returns an empty range just for initialization, since RangeView is not default constructible.
     */
    static SortedIndexRange initSortedIndexRange() 
    {
      static std::vector<IndexPair> dummy; // empty
      return SortedIndexRange(dummy.begin(),dummy.end());
    }
    
    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(Cell const& cell) const 
    {
      return sortedIndices(indexSet().index(cell));
    }
    
    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(size_t n) const 
    {
      return SortedIndexRange(sortedIdx[n].begin(),sortedIdx[n].end());
    }
    
    /**
     * \brief Returns the number of global degrees of freedom managed.
     */
    size_t size() const { return nDof; }

    /**
     * \brief A class mapping local shape function values and derivatives to global shape function values and derivatives.
     */
    typedef HdivConverter<GridView> Converter;

    /**
     * \brief A class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For Hdiv conforming elements, this is a diagonal matrix with entries 1 or -1
     * depending on the global vs. local orientation of the faces.
     */
    class Combiner: public DiagonalCombiner<Scalar> {
    public:
      /**
       * \param cell         the cell for which the combiner is to construct
       * \param gv           the grid view
       * \param is           the index set used for global orientation of edges by global vertex indices
       * \param index        the index of the current cell
       * \param sfs          the shape function set of the current cell
       */
      Combiner(Cell const& cell, GridView const& gv, IndexSet const& is, size_t index, ShapeFunctionSet const& sfs)
      : DiagonalCombiner<Scalar>(sfs.size())
      {
        // First compute the global indices of the neighboring cells for face orientation. If the 
        // neighbor has lower index, we flip our face functions.
        Scalar faceOrient[dim+1];                                  // we have dim+1 faces in a simplex
        for (auto fi=gv.ibegin(cell); fi!=gv.iend(cell); ++fi)     // and the same number of intersections in conforming grids
          if (fi->neighbor())
            faceOrient[fi->indexInInside()] = is.index(*(fi->outside()))>index? 1.0: -1.0;
          else
            faceOrient[fi->indexInInside()] = 1.0;

        // Assign an orientation to all shape functions depending on which face (if any) they are associated with.
        for (int i=0; i<orient.size(); ++i)
        {
          auto loc = sfs[i].location();
          int codim = std::get<1>(loc);
          if (codim>0)
            this->orient[i] = faceOrient[std::get<2>(loc)]; // face shape function: look up the face number
          else
            this->orient[i] = 1.0;                          // cell shape functions: orientation doesn't matter
        }
      }
    };

    /**
     * \brief Returns a combiner for the given cell.
     * \param cell the grid cell for which the combiner is requested
     * \param index the index of the cell
     */
    Combiner combiner(Cell const& cell, size_t index) const 
    { 
      assert(indexSet.index(cell)==index);
      return Combiner(cell,gridView,indexSet,index,shapefunctions(cell,true)); 
    }

    /**
     * \brief (Re)computes the internal data.
     *
     * This has to be called after grid modifications.
     * 
     * In Hdiv conforming spaces, there are two types of ansatz functions: (i) cell-interior functions (with 
     * vanishing normal components on the boundary) and (ii) face functions (with Lagrangian property on
     * normal component at the face/edge/vertex nodes). They are ordered 
     */
    void update() {
      // Precompute and cache all the global indices. First allocate the
      // memory needed to prevent frequent reallocation.
      globIdx.resize(indexSet.size(0));
      
      // Compute the number of cell-local and face functions.
      int nCellSf = 0;
      int nFaceSf = 0;
      for (int i=0; i<sfs.size(); ++i)
        if (std::get<1>(sf[i].location()) == 0) // codim 0 -> cell function
          ++nCellSf;
        else
          ++nFaceSf;
        
      nCellF = indexSet.size(0) * nCellSf;
      size t const nFaceF = indexSet.size(1) * nFaceSf;
      nDof = nCellF + nFaceF;
      
      // We sort the ansatz functions first cell local functions then face functions. This way,
      // the matrix gets an arrow shape. Within each cell and face, the shape functions are sorted 
      // according to their local number and their globally unique local number, respectively.
      // Step through all cells.
      auto end = gridView.template end<0>() ;
      for (auto ci=gridView.template begin<0>(); ci!=end; ++ci) 
      {
        size_t const cellIdx = indexSet.index(*ci);
        auto& gidx = globIdx[cellIdx];
        auto const& sf = shapefunctions(*ci);
        
        GlobalBarycentricPermutation<dimension> gbp(indexSet,*ci);
        
        // Allocate needed memory. The number of local ansatz functions the number of vectorial 
        // shape functions.
        gidx.resize(sf.size());
        sortedIdx[k].resize(gidx.size());
        for (int i=0; i<gidx.size(); ++i)
        {
          int nominalOrder, codim, entity, indexInEntity;
          std::tie(nominalOrder,codim,entity,indexInEntity) = sf[i].location();
          if (codim == 0) // codim 0 -> cell function
          {
            gidx[i] = cellIdx*nCellSf + indexInEntity;
          }
          else                                    // codim 1 -> face function
          {
            // obtain the permutation to globally unique numbering and the barycentric coordinates of 
            // the shape function's location to be permuted.
            std::array<int,dimension> pi = gbp.barycentricSubsetPermutation(entity);
            std::array<int,dimension+1> bc = barycentric(SimplexLagrangeDetail::tupleIndex(order,sf[i].location()),order);
            std::array<int,dimension> subBc;
            for (int i=0; i<dimension; ++i)
              subBc[i] = bc[pi[i]];
            
            // get local sequence number of (permuted) shape function's nodal point within the subentity.
            int k = SimplexLagrangeDetail::local(&subBc[0],dim-1,order);
            
            size_t faceIdx = indexSet.subIndex(*ci,entity,1);
            gidx[i] = nCellF + faceIdx*nFaceSf + k;
          }
          sortedIdx[cellIdx][i] = std::make_pair(gidx[i],i);
        }
        std::sort(sortedIdx[cellIdx].begin(),sortedIdx[cellIdx].end());
      }
    }

  private:
    GridView                   gridView;
    IndexSet const&            indexSet;
    int                        order;
    
    size_t                     nCellF;  // number of cell-local ansatz functions
    size_t                     nDof;

    // For each cell, this vector contains a collection with global
    // indices of each shape function on the cell. Note that up to their
    // sign, global shape functions and global ansatz functions
    // coincide.
    std::vector<std::vector<size_t> > globIdx;

    // In sortedIdx, for each cell there is a vector of (global index, local index) pairs, sorted ascendingly
    // by the global index. This could be computed outside on demand, but the sorting takes time and may be 
    // amortized over multiple assembly passes/matrix blocks if cached here.
    std::vector<std::vector<IndexPair> > sortedIdx;
    
    ShapeFunctionSet const& sfs;
  };
  
  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------
  
  /**
   * \ingroup fem
   * \brief Degrees of freedom manager for Nedelec edge elements.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class NedelecMapper
  {
  public:
    typedef ScalarType RT __attribute__((deprecated));
    typedef ScalarType Scalar;


    typedef GV    GridView;
    typedef typename GridView::Grid     Grid;
    typedef typename GridView::IndexSet IndexSet;

    static int const dim = Grid::dimension;

    /**
     * \brief Whether the ansatz functions have global support (i.e. lead to dense matrices).
     */
    static bool const globalSupport = false;

    typedef bool ConstructorArgument;
    static int const continuity = -1;
    typedef RangeView<std::vector<size_t>::const_iterator> GlobalIndexRange;

    typedef std::pair<size_t,int>                  IndexPair;
    typedef std::vector<IndexPair>::const_iterator SortedIndexIterator;
    typedef RangeView<SortedIndexIterator>         SortedIndexRange;
  
    typedef typename Grid::template Codim<0>::Entity                       Cell;
    typedef typename NedelecShapeFunctionSetContainer<typename Grid::ctype,dim,Scalar>::value_type ShapeFunctionSet;

    /**
     * \brief Constructs a NedelecMapper for a given grid view. The second
     * parameter is unused and only specified for interface consistency.
     */
    NedelecMapper(GridView const& gridView_, bool /* unused */):
      gridView(gridView_), indexSet(gridView_.indexSet())
    {
      Dune::GenericReferenceElement<typename Grid::ctype,dim> const &simplex = Dune::GenericReferenceElements<typename Grid::ctype,dim>::simplex();
      //    Dune::ReferenceSimplex<typename Grid::ctype,dim> simplex;

      vertexIndex_.reserve(simplex.size(dim-1));

      for (int e=0; e<simplex.size(dim-1); ++e)
        vertexIndex_.push_back(std::make_pair(simplex.subEntity(e,dim-1,0,dim),
            simplex.subEntity(e,dim-1,1,dim)));
      update();
    }

    virtual ~NedelecMapper() {}

    ShapeFunctionSet const& shapefunctions(Cell const& cell) const
    {
      assert(cell.type().isSimplex());
      return nedelecShapeFunctionSet<typename Grid::ctype,dim,Scalar>(cell.type(),1);
    }

  ShapeFunctionSet const& shapefunctions(size_t n) const
  {
    assert(cell.type().isSimplex());
    return nedelecShapeFunctionSet<typename Grid::ctype,dim,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dim),1);
  }

    /**
     * Returns an empty range just for initialization purposes, since
     * RangeView is not default constructible.
     */
    GlobalIndexRange initGlobalIndexRange() const
    {
      return GlobalIndexRange(globIdx[0].begin(),globIdx[0].end());
    }

    /**
     * Returns a const sequence containing the global indices of the
     * shape functions associated to this cell. Global indices start at
     * 0 and are consecutive.
     */
    GlobalIndexRange globalIndices(Cell const& cell) const
    {
    return globalIndices(indexSet.index(cell));
  }

  /**
   * Returns a const sequence containing the global indices of the
   * shape functions associated to this cell. Global indices start at
   * 0 and are consecutive.
   */
  GlobalIndexRange globalIndices(size_t n) const 
  {
      return GlobalIndexRange(globIdx[n].begin(),globIdx[n].end());
    }

    /**
     * \brief Returns an empty range just for initialization, since RangeView is not default constructible.
     */
    static SortedIndexRange initSortedIndexRange() 
    {
      static std::vector<IndexPair> dummy; // empty
      return SortedIndexRange(dummy.begin(),dummy.end());
    }
    
    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(Cell const& cell) const 
    {
      return sortedIndices(indexSet().index(cell));
    }
    
    /**
     * \brief Returns an immutable sequence of (global index, local index) pairs sorted in ascending global index order.
     */
    SortedIndexRange sortedIndices(size_t n) const 
    {
      return SortedIndexRange(sortedIdx[n].begin(),sortedIdx[n].end());
    }
    
    /**
     * Returns the number of global degrees of freedom managed.
     */
    size_t size() const { return indexSet.size(dim-1); }

    /**
     * \brief A class mapping local shape function values and derivatives to global shape function values and derivatives.
     */
    typedef HcurlConverter<Grid> Converter;
    

    /**
     * \brief A class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For edge elements, this is a diagonal matrix with entries 1 or -1
     * depending on the global vs. local orientation of the edge.
     */
    class Combiner {
    public:
      /// @param cell         the cell for which the combiner is to construct
      /// @param is           the index set used for global orientation of edges by global vertex indices
      /// @param vertexIndex  for each shape function, i.e. for each edge, a pair of local starting vertex index and local end vertex index
      Combiner(Cell const& cell, IndexSet const& is, std::vector<std::pair<int,int> > const& vertexIndex)
      : orient(vertexIndex.size())
      {
        // first compute the global indices of cell (i.e. simplex) vertices.
        int globalVertexIndex[dim+1];
        for (int i=0; i<=dim; ++i)
          globalVertexIndex[i] = is.subIndex(cell,i,dim);

        // induce a global orientation based on the numbering of vertices
        // edge points locally from locally lower to locally higher index
        // edge points globally from globally lower index to globally higher index
        for (int i=0; i<vertexIndex.size(); ++i) {
          int p0 = vertexIndex[i].first;
          int p1 = vertexIndex[i].second;
          orient[i] = globalVertexIndex[p0]<globalVertexIndex[p1]? 1: -1;
        }
      }

      /**
       * \brief In-place computation of \f$ A \leftarrow A K \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void rightTransform(Matrix& A) const {
        // multiply from right -> modify columns
        assert(A.M()==orient.size());
        for (int i=0; i<A.N(); ++i)
          for (int j=0; j<A.M(); ++j)
            A[i][j] *= orient[j];
      }

      /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
      template <int n, int m>
      void rightTransform(std::vector<VariationalArg<Scalar,n,m> >& v) const {
        assert(v.size()==orient.size());
        for (int i=0; i<v.size(); ++i) {
          v[i].value *= orient[i];
          v[i].gradient *= orient[i];
        }
      }

      /**
       * \brief In-place computation of \f$ A \leftarrow K^+ A \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void leftPseudoInverse(Matrix& A) const {
        // Since K^{-1} = K, this is simple...
        // multiply from left -> modify rows
        assert(A.N()==orient.size());
        for (int i=0; i<A.N(); ++i)
          for (int j=0; j<A.M(); ++j)
            A[i][j] *= orient[i];
      }

      /// Implicit conversion to a sparse matrix.
      /// This is just the diagonal.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >() const
        {
        int n = orient.size();
        Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > K(n,n,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::random);
        for (int i=0; i<n; ++i)
          K.incrementrowsize(i);
        K.endrowsizes();
        for (int i=0; i<n; ++i)
          K.addindex(i,i);
        K.endindices();
        for (int i=0; i<n; ++i)
          *K[i].begin() = orient[i];
        return K;
        }

    private:
      // contains the diagonal of K
      std::vector<Scalar> orient;
    };

    /**
     * Returns a combiner for the given cell.
     * \param cell the grid cell for which the combiner is requested
     * \param index the index of the cell
     */
    Combiner combiner(Cell const& cell, size_t index) const 
    { 
      assert(indexSet.index(cell)==index);
      return Combiner(cell,indexSet,vertexIndex_); 
    }

    /**
     * @brief (Re)computes the internal data.
     *
     * This has to be called after grid modifications.
     */
    void update() {
      // Precompute and cache all the global indices. First allocate the
      // memory needed to prevent frequent reallocation.
      globIdx.resize(indexSet.size(0));

      // Step through all cells.
      //    typedef typename IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;

      typedef typename GridView::template Codim<0>::Iterator CellIterator;

      CellIterator end = gridView.template end<0>() ;

      //    CellIterator end = indexSet.template end<0,Dune::All_Partition>();
      for (CellIterator ci=gridView.template begin<0>(); ci!=end; ++ci) {
      size_t const k = indexSet.index(*ci);
      std::vector<size_t>& gidx = globIdx[k];
        ShapeFunctionSet const& sf = shapefunctions(*ci);
        // Allocate needed memory.
        gidx.resize(sf.size());
      sortedIdx[k].resize(gidx.size());
        for (int i=0; i<gidx.size(); ++i)
      {
          gidx[i] = indexSet.subIndex(*ci,std::get<2>(sf[i].location()),dim-1);
        sortedIdx[k][i] = std::make_pair(gidx[i],i);
      }
      std::sort(sortedIdx[k].begin(),sortedIdx[k].end());
      }
    }

  private:
    GridView                   gridView;
    IndexSet const&            indexSet;

    // For each cell, this vector contains a collection with global
    // indices of each shape function on the cell. Note that up to their
    // sign, global shape functions and global ansatz functions
    // coincide.
    std::vector<std::vector<size_t> > globIdx;

    // In sortedIdx, for each cell there is a vector of (global index, local index) pairs, sorted ascendingly
    // by the global index. This could be computed outside on demand, but the sorting takes time and may be 
    // amortized over multiple assembly passes/matrix blocks if cached here.
    std::vector<std::vector<IndexPair> > sortedIdx;
  
    // For each edge in the reference simplex, this array contains the
    // local vertex indices of the start and end points of the edge.
    std::vector<std::pair<int,int> > vertexIndex_;
  };
} // end of namespace Kaskade



#endif
