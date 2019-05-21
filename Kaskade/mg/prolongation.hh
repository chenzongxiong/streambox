/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MG_PROLONGATION_HH
#define MG_PROLONGATION_HH

#include <boost/fusion/include/vector.hpp>
#include <boost/timer/timer.hpp>

#include "fem/barycentric.hh"
#include "fem/spaces.hh"
#include "linalg/dynamicMatrix.hh"
#include "linalg/localMatrices.hh"
#include "linalg/threadedMatrix.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  /**
   * \ingroup iterative
   * \brief Computes an interpolation-based prolongation matrix from a (supposedly) coarser space to a finer space.
   *
   * \tparam CoarseSpace a FEFunctionSpace type
   * \tparam FineSpace a FEFunctionSpace type
   *
   * \param coarseSpace the (supposedly) coarser space (domain)
   * \param fineSpace   the (supposedly) finer space (image)
   *
   * Both function spaces have to be defined on the very same grid view.
   */
  template <class CoarseSpace, class FineSpace>
  NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>> prolongation(CoarseSpace const& coarseSpace, FineSpace const& fineSpace)
  {
    namespace bf = boost::fusion;

    Timings& timer = Timings::instance();
    timer.start("space prolongation");

    std::vector<bf::vector<std::vector<size_t>,std::vector<size_t>,
                           DynamicMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>>>> interpolationData(fineSpace.gridView().size(0));

    // Intermediate shape function values, declared here to prevent frequent reallocations
    typename CoarseSpace::Mapper::ShapeFunctionSet::SfValueArray afValues;

    // Step through all the cells and perform local interpolation on each cell.
    // TODO: do this in parallel
    for (auto const& cell: elements(fineSpace.gridView()))
    {
      auto index = fineSpace.indexSet().index(cell);

      auto const& fineIndices = fineSpace.mapper().globalIndices(index);
      auto const& coarseIndices = coarseSpace.mapper().globalIndices(index);

      std::copy(fineIndices.begin(),  fineIndices.end(),  std::back_inserter(bf::at_c<0>(interpolationData[index])));
      std::copy(coarseIndices.begin(),coarseIndices.end(),std::back_inserter(bf::at_c<1>(interpolationData[index])));

      auto const& fineSfs = fineSpace.mapper().shapefunctions(cell);
      auto const& coarseSfs = coarseSpace.mapper().shapefunctions(cell);

      // Obtain interpolation nodes of target on this cell
      auto const& iNodes(fineSpace.mapper().shapefunctions(cell).interpolationNodes());

      // Evaluate coarse space global shape functions
      evaluateGlobalShapeFunctions(coarseSpace,cell,iNodes,afValues,coarseSfs);

      // Interpolate fine space global values
      approximateGlobalValues(fineSpace,cell,afValues,bf::at_c<2>(interpolationData[index]),fineSfs);
    }
    timer.stop("space prolongation");

    timer.start("matrix creation");
    // Create a sparse prolongation matrix sparsity pattern. Note that the local prolongation matrices
    // are often sparse (e.g. for Lagrangian elements). We filter out the zero entries in order not to
    // create too densely populated prolongations, which in turn lead to very expensive conjugation
    // computations.
    NumaCRSPatternCreator<size_t> creator(fineSpace.degreesOfFreedom(),coarseSpace.degreesOfFreedom());
    for (auto const& block: interpolationData)
      for (int i=0; i<bf::at_c<0>(block).size(); ++i)
        for (int j=0; j<bf::at_c<1>(block).size(); ++j)
        {
          auto entry = bf::at_c<2>(block)[i][j];
          if (std::abs(entry) > 1e-12)
            creator.addElement(bf::at_c<0>(block)[i],bf::at_c<1>(block)[j]);
        }

    // Now create and fill the matrix itself. Overwrite already existing entries
    // (this means the weighing is just "last one wins"
    NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>,size_t> p(creator);
    for (auto const& block: interpolationData)
      for (int i=0; i<bf::at_c<0>(block).size(); ++i)
        for (int j=0; j<bf::at_c<1>(block).size(); ++j)
        {
          auto entry = bf::at_c<2>(block)[i][j];
          if (std::abs(entry) > 1e-12)
            p[bf::at_c<0>(block)[i]][bf::at_c<1>(block)[j]] = entry;
        }

    timer.stop("matrix creation");
    return p;
  }

  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup iterative
   * \brief Computes a stack of prolongation matrices for higher order finite element spaces.
   *
   * This computes a stack of prolongation matrices for a scale of finite element spaces with
   * increasing polynomial ansatz order. From level to level, the ansatz order is doubled (maybe +1
   * for odd degrees).
   *
   * \tparam Mapper a finite element local to global mapper with scalar shape functions.
   */
  template <class Mapper>
  std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>> prolongationStack(FEFunctionSpace<Mapper> const& space)
  {
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>> stack;
    int p = space.mapper().maxOrder();

    if (p > 1) // there is some need for prolongation
    {
      H1Space<typename Mapper::Grid> coarseSpace(space.gridManager(),space.gridView(),p/2);  // a coarser space (lower order)
      stack = prolongationStack(coarseSpace);
      stack.push_back(prolongation(coarseSpace,space));
    }

    return stack;
  }

  // ---------------------------------------------------------------------------------------------------------

  /// \cond internals
  namespace ProlongationDetail
  {
    struct Node
    {
      size_t p, q; // global (fine grid) indices of parent vertices
      int level;   // level of the vertex
    };
  }
  /// \endcond

  /**
   * \ingroup iterative
   * \brief A prolongation operator for P1 finite elements from a coarser grid level to the next finer level.
   *
   * For P1 finite elements on simplicial grids, the prolongation from a coarser level to the next finer is a matrix \f$ P \f$
   * in which each row contains either one entry of value 1 (if the vertex is already contained in the coarser level grid)
   * or two entries that sum up to 1 (if the vertex is created on the finer level by bisecting the edge between two parent
   * nodes). In the latter case, the two entries default to 0.5 each.
   *
   * Storing this as a sparse matrix data structure is quite inefficient.
   *
   * This class is a specialized implementation of such prolongation matrices.
   */
  class MGProlongation
  {
  public:
    /**
     * \brief Constructor.
     * \param parents a sequence of node information. The parent indices for fine grid nodes shall refer to the fine grid numbering.
     * \param indexInCoarse a vector of length of fine grid vertices. For each coarse grid node, it contains the index of that node
     *                      in the coarse grid.
     * \param nc the number of coarse grid nodes (i.e. number of columns in P)
     * \param fineLevel the level of fine grid nodes
     */
    MGProlongation(std::vector<ProlongationDetail::Node> const& parents, std::vector<size_t> const& indexInCoarse, size_t nc, int fineLevel);

    /**
     * \brief Returns the parent nodes indices.
     *
     * For a fine grid node, this returns the parent node indices in the coarse grid. For a coarse
     * grid node, the "parent indices" coincide and refer to the node index in the coarse grid.
     */
    std::array<size_t,2> const& parents(size_t i) const
    {
      return entries[i];
    }

    /**
     * \brief Matrix vector multiplcation (update mode).
     *
     * This computes \f$ f \leftarrow f + Pc \f$.
     *
     * \param f the fine grid target vector. It has to have the correct size.
     */
    template <class Vector>
    void umv(Vector const& c, Vector& f) const
    {
      assert(f.size()==entries.size());
      assert(c.size()==nc);

      for (size_t row=0; row<entries.size(); ++row)
      {
        auto e = entries[row];
        f[row] += 0.5*(c[e[0]] + c[e[1]]);
      }
    }

    /**
     * \brief Matrix vector multiplcation.
     *
     * This computes \f$ f \leftarrow Pc \f$, and is faster than but equivalent to
     * \code
     * f = 0;
     * umv(c,f);
     * \endcode
     *
     * \param f the fine grid target vector. It has to have the correct size.
     */
    template <class Vector>
    void mv(Vector const& c, Vector& f) const
    {
      assert(f.size()==entries.size());
      assert(c.size()==nc);

      for (size_t row=0; row<entries.size(); ++row)
      {
        auto e = entries[row];
        f[row] = 0.5*(c[e[0]] + c[e[1]]);
      }
    }

    /**
     * \brief Transpose matrix vector multiplication
     *
     * This computes \f$ c \leftarrow P^T f \f$.
     *
     * \param c the coarse grid target vector. It is resized to contain the result.
     */
    template <class Vector>
    void mtv(Vector const& f, Vector& c) const
    {
      assert(f.size()==entries.size());
      if (c.size()!=nc)
        c.resize(nc);
      for (size_t i=0; i<c.N(); ++i)
        c[i] = 0;

      for (size_t row=0; row<entries.size(); ++row)
      {
        c[entries[row][0]] += 0.5*f[row];
        c[entries[row][1]] += 0.5*f[row];
      }
    }

    /**
     * \brief The number of rows.
     */
    size_t N() const
    {
      return entries.size();
    }

    /**
     * \brief The number of columns.
     */
    size_t M() const
    {
      return nc;
    }

    /**
     * \brief Galerkin projection of fine grid matrices.
     * This computes \f$ P^T A P \f$.
     *
     * \tparam Entry the entry type of Galerkin matrix to be projected, in general a quadratic Dune::FieldMatrix type
     * \tparam Index the Index type of the matrix to be projected
     *
     * \param A the quadratic fine level Galerkin matrix of size \f$ N\times N \f$ to be projected
     * \param onlyLowerTriangle if true, only the lower triangular part of symmetric A is touched, and only the lower triangular part of \f$ P^T A P\f$ is created
     */
    template <class Entry, class Index>
    NumaBCRSMatrix<Entry,Index> galerkinProjection(NumaBCRSMatrix<Entry,Index> const& A, bool onlyLowerTriangle = false) const
    {
      assert(A.N()==N() && A.M()==N());
      Timings& timer = Timings::instance();

      // First, create the sparsity pattern of the projected matrix.
      NumaCRSPatternCreator<Index> creator(M(),M(),onlyLowerTriangle);

      // If C = P^T A P, we have that C_{ij} = \sum_{k,l} P_{ki} P_{lj} A_{kl}. Hence, the entry A_{kl} contributes to
      // all C_{ij} for which there are nonzero entries P_{ki} and P_{lj} in the rows k and l of P. Thus we can simply
      // run through all nonzeros A_{kl} of A, look up the column indices i,j of rows k and l of P, and flag C_{ij}
      // as nonzero.

      // Step through all entries of A
      timer.start("h conjugation pattern");
      for (Index k=0; k<A.N(); ++k)
      {
        auto const& is = entries[k];                        // indices i for which Pki != 0
        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          Index const l = ca.index();
          auto const& js = entries[l];                      // indices j for which Plj != 0
          // add all combinations i,j
          if (is[0]==is[1] && js[0]==js[1])
            creator.addElement(is[0],js[0]);
          else
            creator.addElements(std::begin(is),std::end(is),std::begin(js),std::end(js));

          if (onlyLowerTriangle && k>l)                     // subdiagonal entry (k,l) of A -> entry (l,k) must be treated implicitly:
            creator.addElements(std::begin(js),std::end(js),std::begin(is),std::end(is)); // add all combinations (j,i)
        }
      }
      timer.stop("h conjugation pattern");

      // An alternative way of computing the sparsity pattern would be to use that the nonzero entries j in column i
      // of C are exactly those for which there is k with (nonzero P_{jk} and there is l with (nonzero P_{il} and A_{lk})).
      // Hence we can obtain the column index set J directly by the following steps:
      // (i) find all l with P_{li} nonzero -> L  [requires to access columns of P - compute the transpose patterns once]
      // (ii) find all k with A_{lk} nonzero for some l in L -> K  [probably sorting K and removing doubled entries would be a good idea here]
      // (iii) find all j with P_{kj} nonzero for some k in K -> J
      // Compared to the above implementation this would have the following (dis)advantages
      // + easy to do in parallel (since write operations are separated)
      // + fewer scattered write accesses to memory
      // - more complex implementation
      // - requires the transpose pattern of P

      // Create the sparse matrix.
      timer.start("matrix creation");
      NumaBCRSMatrix<Entry,Index> pap(creator);
      timer.stop("matrix creation");

      // Fill the sparse matrix PAP. This is done as before by stepping through all Akl entries and scatter
      // Pki*Plj*Akl into PAPij.
      timer.start("matrix conjugation");
      for (Index k=0; k<A.N(); ++k)
      {
        auto const& is = entries[k];
        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          Index const l = ca.index();
          auto const& js = entries[l];

          if (is[0]==is[1] && js[0]==js[1])               // a coarse grid node - this means four identical contributions
            pap[is[0]][js[0]] += *ca;                     // with factor 1/4. Substitute with one contribution with factor 1.
          else
          {
            auto pap0 = pap[is[0]];
            pap0[js[0]] += 0.25 * *ca;
            pap0[js[1]] += 0.25 * *ca;

            auto pap1 = pap[is[1]];
            pap1[js[0]] += 0.25 * *ca;
            pap1[js[1]] += 0.25 * *ca;
          }

          if (onlyLowerTriangle)
            abort();  // not yet implemented
        }
      }
      timer.stop("matrix conjugation");

      return pap;
    }

  private:
    // Internally we treat all rows equally: A row with one entry of value one is represented as
    // two entries which sum up to 1 (that happen to reference the same column index...). This allows a
    // very simple and uniform implementation. The vector entries contains the column indices of
    // the two entries in each row.
    std::vector<std::array<size_t,2>> entries;

    // Number of columns (i.e. coarse grid nodes).
    size_t nc;
  };

  std::ostream& operator<<(std::ostream& out, MGProlongation const& p);

  /**
   * \ingroup multigrid
   * \brief Creates a Galerkin projected Matrix \f$ P^T A P \f$ from a prolongation \f$ P \f$ and a
   *        symmetric matrix \f$ A \f$.
   * \param onlyLowerTriangle if true, A contains onl the lower triangular part of the symmetric matrix.
   */
  template <class Entry, class Index>
  auto conjugation(MGProlongation const& p, NumaBCRSMatrix<Entry,Index> const& a, bool onlyLowerTriangle)
  {
    return p.template galerkinProjection(a,onlyLowerTriangle);
  }

  // ---------------------------------------------------------------------------------------------------------

  /// \cond internals
  namespace ProlongationDetail
  {
    // For each corner of the given cell, if there are any of them which are not corners of the father cell
    // (i.e. created by bisection of an edge), enter the parent vertices into the parents vector.
    template <class LeafView, class Cell, class Parents>
    void computeParents(LeafView const& leafView, Cell const& cell, Parents& parents)
    {
      // If this cell is a coarse grid cell, none of its vertices have parents, and we're done.
      if (!cell.hasFather())
        return;

      int const dim = LeafView::dimension;

      // For each corner, check its position in the parent. If it has barycentric coordinates with one
      // entry approximately one (all others zero), it is a corner of the father cell and will be treated
      // later (or has already been treated). Otherwise there will be two entries 0.5 (all others zero),
      // and these entries denote father corners which are the parents.
      assert(cell.type().isSimplex());
      auto const& geo = cell.geometryInFather();
      int nCorners = geo.corners();
      for (int i=0; i<nCorners; ++i)
      {
        // TODO: We could first obtain the index of the corner and check whether its parents have
        //       already been determined, and skip the geometric considerations in that case.
        //       This should save some time.
        //       On the other hand, obtaining indices can also be expensive, and in the current
        //       implementation we only have to get them for corners that actually are no corners
        //       in the father cell. This saves some time, too.
        //       We should check which option is faster, and whether it makes a big difference in
        //       the first place.

        // barycentric coordinates of corner in the father cell
        auto b = barycentric(geo.corner(i));

        // find potential parent vertices (those with barycentric coordinates 0.5)
        int pcount = 0;
        int pc[2];
        for (int k=0; k<b.N(); ++k)                         // check all barycentric coordinates
          if (std::abs(b[k]-0.5) < 0.01)                    // if close to 0.5 accept
          {
            pc[pcount] = (k+1) % (dim+1);                   // map barycentric coordinate number to Dune corner number
            ++pcount;                                       // note down that we've found one more parent vertex
          }

        if (pcount == 2)                                    //  corner i is a father edge midpoint
        {
          auto const& is = leafView.indexSet();
          parents.push_back(std::make_pair(is.index(cell.template subEntity<dim>(i)),
                                           Node{is.index(cell.father().template subEntity<dim>(pc[0])),
                                                is.index(cell.father().template subEntity<dim>(pc[1])),
                                                cell.level()}));
        }
      }
    }

    /**
     * \ingroup multigrid
     * \brief Creates a stack of prolongations from parent-child relationships in grids
     * \param level
     */
    std::vector<MGProlongation> makeProlongationStack(std::vector<Node> parents, int level, size_t minNodes);
  };
  /// \endcond

  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Computes a sequence of prolongation matrices for P1 finite elements in hierarchical grids.
   *
   * We assume that each vertex not contained in the coarse grid has been created by bisecting an edge,
   * such that there are exactly two "parent vertices".
   * \tparam Grid the Dune grid type on which the P1 space is defined. The grid has to be a simplicial grid.
   * \param grid the grid itself
   * \param minNodes minimum number of nodes to keep as coarse grid
   */
  template <class GridMan>
  std::vector<MGProlongation> prolongationStack(GridMan const& gridman, size_t minNodes=0)
  {
    Timings& timer = Timings::instance();

    auto const& grid = gridman.grid();
    auto const& leafView = grid.leafGridView();
    auto const& cellRanges = gridman.cellRanges(grid.levelGridView(0));

    // In parallel compute the parent nodes for edge midpoints cell by cell
    timer.start("parent computation");
    std::vector<std::vector<std::pair<size_t,ProlongationDetail::Node>>> myParents(cellRanges.maxRanges());
    parallelFor([&](int k, int n)
    {
      for (auto const& coarseCell: cellRanges.range(n,k))
        for (auto const& cell: descendantElements(coarseCell,grid.maxLevel()))
          ProlongationDetail::computeParents(leafView,cell,myParents[k]);
    },myParents.size());

    // Now that all parent nodes have been identified, consolidate them in an easily
    // indexable array.
    std::vector<ProlongationDetail::Node> parents(leafView.size(leafView.dimension),ProlongationDetail::Node{0,0,-1});
    for (auto const& ps: myParents)
      for (auto const& p: ps)
        parents[p.first] = p.second;
    timer.stop("parent computation");


    // Now construct the prolongation matrices for all levels.
    timer.start("made stack parents");
    auto ps = ProlongationDetail::makeProlongationStack(std::move(parents),grid.maxLevel(),minNodes);
    timer.stop("made stack parents");
    return ps;
  }


  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Class for multigrid stacks.
   *
   * This provides the storage of and access to prolongations and projected Galerkin matrices,
   * but leaves the construction of these to derived classes.
   */
  template <class Prolongation, class Entry, class Index>
  class MultiGridStack
  {
  public:

    /**
     * \brief Constructor
     *
     * This takes both a stack of prolongations and a matching stack of projected
     */
    MultiGridStack(std::vector<Prolongation>&& ps, std::vector<NumaBCRSMatrix<Entry,Index>>&& as)
    : prolongations(std::move(ps)), galerkinMatrices(std::move(as))
    {
      assert(prolongations.size()+1==galerkinMatrices.size());
    }

    /**
     * \brief Constructor
     *
     * This takes a stack of prolongations and creates the stack of projected Galerkin matrices from the given
     * fine grid matrix. Most useful for geometric multigrid, where the prolongations are defined solely in terms of the grid.
     *
     * \see makeMultiGridStack
     */
    MultiGridStack(std::vector<Prolongation>&& ps, NumaBCRSMatrix<Entry,Index>&& A, bool onlyLowerTriangle)
    {
      Timings& timer = Timings::instance();

      prolongations = std::move(ps);
      galerkinMatrices.push_back(std::move(A));

      assert((prolongations.size() == 0) || (prolongations.back().N() == galerkinMatrices[0].N()));
      assert(galerkinMatrices[0].N() == galerkinMatrices[0].M());

      // Step through the prolongations from top level to bottom
      timer.start("matrix projection");
      for (auto pi=prolongations.crbegin(); pi!=prolongations.crend(); ++pi)
        galerkinMatrices.insert(galerkinMatrices.begin(),conjugation(*pi,galerkinMatrices.front(),onlyLowerTriangle));
      timer.stop("matrix projection");
    }

    MultiGridStack(MultiGridStack&& other) = default;

    /**
     * \brief The number of grid levels
     */
    int levels() const
    {
      return galerkinMatrices.size();
    }

    /**
     * \brief Returns the prolongation from given level to next higher one.
     * \param level precondition 0 <= level < levels()-1
     */
    Prolongation const& p(int level) const
    {
      return prolongations[level];
    }

    /**
     * \brief Returns the projected Galerkin matrix on the given level.
     * \param level precondition 1 <= level < levels()
     */
    NumaBCRSMatrix<Entry,Index> const& a(int level) const
    {
      return galerkinMatrices[level];
    }

    /**
     * \brief Returns the projected Galerkin matrix on the coarsest level.
     *
     * This is explicitly a mutable reference, such that the coarse grid matrix can be moved
     * from. This is useful, as in the multigrid, the coarsest level matrix is not referenced
     * (the coarse grid preconditioner has its own copy, maybe obtained by moving from here...).
     */
    NumaBCRSMatrix<Entry,Index>& coarseGridMatrix()
    {
      return galerkinMatrices[0];
    }

    void report(std::ostream& out) const
    {
      for (auto const& p: prolongations)
        out << "Prolongation:\n" << p;

      for (auto const& a: galerkinMatrices)
        out << "GalerkinMatrix: \n" << a;
    }

  private:
    std::vector<Prolongation>                 prolongations;    // prolongations grid i -> i+1, i=0,...,n-1
    std::vector<NumaBCRSMatrix<Entry,Index>>  galerkinMatrices; // Galerkin matrices on grid i, i=0,...,n
  };

  template <typename Prolongations, typename Entry, typename Index>
  std::ostream& operator<<(std::ostream& out, MultiGridStack<Prolongations,Entry,Index> const& mgStack) { mgStack.report(out); return out; }

  /**
   * \ingroup multigrid
   * \brief Convenience routine for creating multigrid stacks
   *
   * Given a stack of prolongations and the top level Galerkin matrix, this creates the complete stack
   * including all projected Galerkin matrices.
   *
   * \param ps a vector of prolongation matrices
   * \param A the top level (finest) Galerkin matrix to be projected down
   * \param onlyLowerTriangle if true, only the lower triangular part of A will be referenced
   */
  template <class Prolongation, class Entry, class Index>
  MultiGridStack<Prolongation,Entry,Index> makeMultiGridStack(std::vector<Prolongation>&& ps, NumaBCRSMatrix<Entry,Index>&& A,
                                                              bool onlyLowerTriangle)
  {
    return MultiGridStack<Prolongation,Entry,Index>(std::move(ps),std::move(A),onlyLowerTriangle);
  }

  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks based on geometric coarsening for P1 elements
   */
  template <class GridMan, class Entry, class Index>
  MultiGridStack<MGProlongation,Entry,Index> makeGeometricMultiGridStack(GridMan const& gridManager, NumaBCRSMatrix<Entry,Index>&& A,
                                                                         size_t minNodes=10000, bool onlyLowerTriangle=false)
  {
    return MultiGridStack<MGProlongation,Entry,Index>(prolongationStack(gridManager,minNodes),std::move(A),onlyLowerTriangle);
  }

  /**
   * \ingroup multigrid
   * \brief Creates stack of prolongations and projected Galerkin matrices.
   *
   * \param A the symmetric sparse matrix
   * \param n stop the coarsening if the number of rows/cols drops below this number
   * \param onlyLowerTriangle if true, only the lower triangular part of symmetric A is accessed
   */
  template <class Entry, class Index>
  MultiGridStack<MGProlongation,Entry,Index> makeAlgebraicMultigridStack(NumaBCRSMatrix<Entry,Index>&& A,
                                                                         Index n=0, bool onlyLowerTriangle=false);


  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks based on coarsening by reducing the ansatz order from P to P1
   */
  template <typename FineSpace, typename Matrix>
  auto makePMultiGridStack(FineSpace const& space, Matrix&& A, bool onlyLowerTriangle)
  {
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>> prolongations;

    if (space.mapper().maxOrder() > 1)
    {
      H1Space<typename FineSpace::Grid> coarseSpace(space.gridManager(),space.gridView(),1);  // a coarser space (lower order)
      prolongations.push_back(prolongation(coarseSpace,space));
    }

    return makeMultiGridStack(std::move(prolongations),std::move(A),onlyLowerTriangle);
  }

  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks between two spaces.
   *
   * An approximation of the projected Galerkin matrix is provided. Often, this can be assembled much more
   * efficiently than projected.
   */
  template <typename FineSpace, typename CoarseSpace, typename Matrix>
  auto makePMultiGridStack(FineSpace const& fineSpace, Matrix&& fA, CoarseSpace const& coarseSpace, Matrix&& cA)
  {
    using Prolongation = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
    std::vector<Prolongation> prolongations;
    prolongations.push_back(prolongation(coarseSpace,fineSpace));

    std::vector<Matrix> matrices{std::move(cA),std::move(fA)};

    return MultiGridStack<Prolongation,typename Matrix::block_type,
                                       typename Matrix::size_type>(std::move(prolongations),std::move(matrices));
  }

}

#endif
