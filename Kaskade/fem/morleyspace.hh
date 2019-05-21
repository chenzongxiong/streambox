/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MORLEYSPACE_HH
#define MORLEYSPACE_HH

#include <numeric>

#include "fem/converter.hh"
#include "fem/functionspace.hh"
#include "fem/lagrangeshapefunctions.hh"
#include "fem/views.hh"

//---------------------------------------------------------------------

namespace Kaskade
{
  /**
   * \ingroup fem
   * \brief Degrees of freedom manager for Morley nonconforming elements.
   * 
   * The Morley elements are piecewise quadratic, continuous at the mesh vertices, and have continuous 
   * edge-normal derivatives at the edge midpoints.
   *
   * \tparam ScalarType scalar type
   * \tparam GV         grid view
   */
  template <class ScalarType, class GV>
  class MorleyMapper
  {
  public:
    typedef ScalarType Scalar;


    typedef GV                          GridView;
    typedef typename GridView::Grid     Grid;
    typedef typename GridView::IndexSet IndexSet;

    static int const dim = Grid::dimension;
    static_assert(dim==2, "Morley elements implemented just for 2D problems.");
    
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
    typedef typename LagrangeShapeFunctionSetContainer<typename Grid::ctype,dim,Scalar>::value_type ShapeFunctionSet;

    typedef LagrangeSimplexShapeFunctionSet<typename Grid::ctype,dim,ScalarType> ShapeFunctionSetImplementation;

    /**
     * \brief Constructs a MorleyMapper for a given grid view. 
     */
    MorleyMapper(GridView const& gridView_):
      gridView(gridView_), indexSet(gridView_.indexSet())
    {
      update();
    }

    virtual ~MorleyMapper() {}

    /**
     * \brief Returns the maximal polynomial order of shape functions encountered in any cell.
     */
    int maxOrder() const
    {
      return 2;
    }

    ShapeFunctionSet const& shapefunctions(Cell const& cell) const
    {
      assert(cell.type().isSimplex());
      return lagrangeShapeFunctionSet<typename Grid::ctype,dim,Scalar>(cell.type(),2);
    }

    ShapeFunctionSet const& shapefunctions(size_t n) const
    {
      assert(cell.type().isSimplex());
      return lagrangeShapeFunctionSet<typename Grid::ctype,dim,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dim),2);
    }
    
    /**
     * \brief Returns an empty range just for initialization purposes, since RangeView is not default constructible.
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
     * \brief Returns a const sequence containing the global indices of the shape functions associated to this cell. 
     * 
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
     * 
     * This is just the number of vertices plus the number of edges.
     */
    size_t size() const { return indexSet.size(dim-1) + indexSet.size(dim); }

    typedef ScalarConverter<Cell> Converter;



    /**
     * \brief A combiner class implementing a matrix \f$ K \f$ mapping a subset
     * of global degrees of freedom (those given by globalIndices()) to
     * local degrees of freedom (shape functions).
     *
     * For Morley elements, this is a dense 6x6 matrix such that the global ansatz functions 
     * are continuous at vertices and their gradients are continuous at edge midpoints.
     * 
     * 
     */
    class Combiner {
    public:
      /// @param cell         the cell for which the combiner is to construct
      /// @param is           the index set used for global orientation of edges by global vertex indices
      /// @param vertexIndex  for each shape function, i.e. for each edge, a pair of local starting vertex index and local end vertex index
      Combiner(Cell const& cell, IndexSet const& is)//, std::vector<std::pair<int,int> > const& vertexIndex)
      {
        // Compute K by solving the interpolation problem. First evaluate the shape functions and their gradients at the 
        // vertices and the edge midpoints.
        Converter psi(cell);
        ShapeFunctionSet const& sfs = lagrangeShapeFunctionSet<typename Grid::ctype,dim,Scalar>(cell.type(),2);
        assert(sfs.size() == 6); // we expect 6 shape functions for 2D quadratic Morley elements
        //auto refTriangle = Dune::ReferenceElement::simplex<typename Grid::ctype,dim>();

        auto const& refTriangle = Dune::GenericReferenceElements<typename Grid::ctype,dim>::simplex();
        auto cellIndex = is.index(cell);
        
        // First create global ansatz functions on edges, i.e. they vanish on vertices and have
        // a unit normal derivative exactly on one edge. For global "continuity" of gradients, the 
        // sign of the outer normal derivative depends on the global orientation of the edge.
        for (auto edge=cell.ileafbegin(); edge!=cell.ileafend(); ++edge)
        {
          Dune::FieldVector<double,1> xiEdge(0.5);  // edge midpoint in edge local coordinates
          auto outerNormal = edge->unitOuterNormal(xiEdge);

          int const j = edge->indexInInside();
          auto xi = refTriangle.position(j,1); // edge midpoint in cell local coordinates
          psi.setLocalPosition(xi);
          
          Scalar sign = 1.0;
          if (edge->neighbor() && is.index(*edge->outside()) < cellIndex)
            sign = -1.0;
          
          for (int i=0; i<sfs.size(); ++i)
          {
            VariationalArg<Scalar,dim,1> dsf(0,sfs[i].evaluateDerivative(xi)); // gradients at edge midpoints, value are of no importance
            dsf = psi.global(dsf);                                             // convert them to the shape function on the actual element
            K[j][i] = sign * (dsf.derivative[0] * outerNormal);
          }
        }
        
        for (int j=0; j<3; ++j)
        {
          auto xi = refTriangle.position(j,2); // vertex position in local coordinates
          psi.setLocalPosition(xi);
          
          for (int i=0; i<sfs.size(); ++i)
          {
            VariationalArg<Scalar,dim,1> sf(sfs[i].evaluateFunction(xi)); // values at vertices, gradienst are of no importance
            sf = psi.global(sf,0);                                        // convert them to the shape function on the actual element
            K[j+3][i] = sf.value;
          }
        }
        
        K.invert();
      }

      /**
       * \brief In-place computation of \f$ A \leftarrow A K \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void rightTransform(Matrix& A) const {
        // multiply from right -> modify columns
        A.rightmultiply(K);
      }

      /// In-place computation of row vectors \f$ v \leftarrow v^T K \f$.
      template <int n, int m>
      void rightTransform(std::vector<VariationalArg<Scalar,n,m> >& vs) const {
        assert(vs.size()==6);
        std::vector<VariationalArg<Scalar,n,m> > phis(6);
        for (int i=0; i<6; ++i)
        {
          phis[i].value = 0;
          phis[i].derivative = 0;
          phis[i].hessian = 0;
          
          // form linear combinations of shape functions leading to global ansatz functions
          for (int j=0; j<6; ++j)
          {
            phis[i].value += K[j][i] * vs[j].value;
            phis[i].derivative += K[j][i] * vs[j].derivative;
            phis[i].hessian += K[j][i] * vs[j].hessian;
          }
        }
        
        // return result
        vs.swap(phis);
      }

      /**
       * \brief In-place computation of \f$ A \leftarrow K^+ A \f$.
       * \tparam Matrix A matrix class satisfying the Dune::DenseMatrix interface.
       */
      template <class Matrix>
      void leftPseudoInverse(Matrix& A) const {
        Dune::FieldMatrix<Scalar,6,6> Kinv = K;
        Kinv.invert();
        A.leftmultiply(Kinv);
      }

      /// Implicit conversion to a sparse matrix.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >() const
      {
        std::cerr << "not implemented yet\n";
        abort();
      }

    private:
      // contains the diagonal of K
      Dune::FieldMatrix<Scalar,6,6> K;
    };

    /**
     * Returns a combiner for the given cell.
     * \param cell the grid cell for which the combiner is requested
     * \param index the index of the cell
     */
    Combiner combiner(Cell const& cell, size_t index) const 
    { 
      assert(indexSet.index(cell)==index);
      return Combiner(cell,indexSet);
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
      sortedIdx.resize(globIdx.size());
      
      size_t const ne = indexSet.size(1); // number of edges

      // Step through all cells.
      auto end = gridView.template end<0>() ;
      for (auto ci=gridView.template begin<0>(); ci!=end; ++ci) {
        size_t const k = indexSet.index(*ci);        // cell number
        std::vector<size_t>& gidx = globIdx[k];      // global inices on the cell
        // Allocate needed memory.
        gidx.resize(6);                              // 2D Morley elements of second order have 6 local DOFs
        sortedIdx[k].resize(gidx.size());
        
        for (int i=0; i<3; ++i)
        {
          gidx[i] = indexSet.subIndex(*ci,i,1);         // edge ansatz functions
          sortedIdx[k][i] = std::make_pair(gidx[i],i);
          
          gidx[i+3] = indexSet.subIndex(*ci,i,2) + ne;  // vertex ansatz functions
          sortedIdx[k][i+3] = std::make_pair(gidx[i+3],i+3);
        }
        std::sort(sortedIdx[k].begin(),sortedIdx[k].end());
      }
    }

  private:
    GridView                   gridView;
    IndexSet const&            indexSet;

    // For each cell, this vector contains a collection with global
    // indices of each shape function on the cell. 
    std::vector<std::vector<size_t> > globIdx;

    // In sortedIdx, for each cell there is a vector of (global index, local index) pairs, sorted ascendingly
    // by the global index. This could be computed outside on demand, but the sorting takes time and may be 
    // amortized over multiple assembly passes/matrix blocks if cached here.
    std::vector<std::vector<IndexPair> > sortedIdx;
  };
} // end of namespace Kaskade



#endif
