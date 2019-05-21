/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef COARSENING_HH
#define COARSENING_HH


#include <boost/timer/timer.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

#include "fem/errorest.hh"
#include "fem/embedded_errorest.hh"
#include "fem/fetransfer.hh"
#include "fem/firstless.hh"
#include "fem/fixfusion.hh"
#include "fem/iterate_grid.hh"
#include "utilities/power.hh"

#include "linalg/dynamicMatrix.hh"

namespace Kaskade
{
  namespace CoarseningDetail
  {

    using namespace boost::fusion;
    
    // A class representing a local projector for a given FE space
    template <class Sp>
    struct Projector
    {
      typedef Sp Space;
      
      Space const*                                                  space;  // pointer to FE space
      std::vector<size_t>                                           gidx;   // global indices
      DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,1,1>>  q;      // stores a factor of the projection matrix Q Q^T.

      // these are declared here to be used in GetLocalTransferProjection - avoid frequent reallocations
      DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,1,1>> pLocal;
      std::vector<int>                                             idx; 
    };
    
    // A functor creating local projectors for a given FE space
    struct CreateProjection
    {
      template <class T> struct result {};
      
      template <class SpacePtr>
      struct result<CreateProjection(SpacePtr)> {
        typedef typename boost::remove_pointer<typename boost::remove_reference<SpacePtr>::type>::type Space;
        typedef Projector<Space> type;
      };
      
      template <class SpacePtr>
      typename result<CreateProjection(SpacePtr)>::type operator()(SpacePtr space) const
      {
        typename result<CreateProjection(SpacePtr)>::type p;
        p.space = space;
        return p;
      }
    };


    // A functor filling a local projector with actual data for a given cell
    template <class CellPointer>
    struct GetLocalTransferProjection
    {
      // father points to the father cell,
      // children contains pointer to its children
      GetLocalTransferProjection(CellPointer father_, std::vector<CellPointer> const& children_): father(father_), children(children_)
      { }

      template <class Projector>
      void operator()(Projector& p) const
      {
	typedef typename Projector::Space Space;
	// We could use LocalTransfer here, but this is conceptual and computational overkill. It supports different FE spaces (we need only one here)
	// and creates new shape function sets just in case the cell vanishes on coarsening. We do the projection here before coarsening...
	// Moreover we know here that all children of the father cell are leafs.

        // Following a direct implementation exploiting our a priori knowledge for simplicity and performance.
        // Obtain all global indices associated to leaf cells within this father cell
        p.gidx.clear();
        for (int i=0; i<children.size(); ++i) 
        {
          auto gi = p.space->mapper().globalIndices(*children[i]);
          p.gidx.insert(p.gidx.end(),gi.begin(),gi.end());
        }
        // global indices can be duplicate - make unique
        std::sort(p.gidx.begin(),p.gidx.end());
        p.gidx.erase(std::unique(p.gidx.begin(),p.gidx.end()),p.gidx.end());
        
        // Step through all children and obtain the prolongation matrix. Scatter this to the 
        // prolongation matrix P for all children's degrees of freedom.
        p.q.setSize(p.gidx.size(),p.space->mapper().shapefunctions(*father).size());
//        p.q = 0 // causes compilation error with dune-2.4.0 and clang++ on OS X (Darwin)
        p.q.fill(0);
        typedef typename Space::Grid Grid;
        typedef typename Space::Scalar   Scalar;
        int const sfComponents = Space::sfComponents;
        DynamicMatrix<Dune::FieldMatrix<Scalar,sfComponents,1>> afValues; // declare out of loop to prevent frequent reallocations
        std::vector<Dune::FieldVector<typename Grid::ctype,Grid::dimension> > iNodes; // declare out of loop to prevent frequent reallocations
        for (int i=0; i<children.size(); ++i) 
        {
          // Compute the mapping of global indices on the child to the set of global indices on all children.
          auto gi = p.space->mapper().globalIndices(*children[i]);
          p.idx.resize(gi.size());
          for (int j=0; j<gi.size(); ++j) 
            p.idx[j] = std::lower_bound(p.gidx.begin(),p.gidx.end(),gi[j]) - p.gidx.begin();
          
          // Compute local prolongation matrix p.pLocal. As father and child spaces are the same, we can use the same shape function set
          // (assuming that the type of cell doesn't change...).
          auto const& sfs = p.space->mapper().shapefunctions(*children[i]);
          localProlongationMatrix(*p.space,children[i],*p.space,father,p.pLocal,sfs,sfs);
          
          // Enter prolongation values on this child into the whole prolongation matrix for all children
          for(size_t k=0; k<p.idx.size(); ++k)
            for(int j=0; j<p.q.M(); ++j)
              p.q[p.idx[k]][j] = p.pLocal[k][j];
        }
        
        // Now the projection onto a coarser space is P P^+, where P^+, the pseudoinverse of P, is the restriction matrix.
        // Note that this formulation of the (not uniquely defined restriction) is different from the interpolation-based
        // restriction as obtained from LocalTransfer. Nevertheless, both make sense and both appear to work very similar.
        // In fact one might argue that the pseudoinverse formulation here makes even more sense, as it leads to a 
        // least-squares approximation of the fine grid solution by the coarse grid restriction. "Least squares" is, however,
        // a norm-dependent concept, and the euclidean norm of FE coefficient vectors as implicitly used here might not 
        // be the best suited.
        // With a QR decomposition of P = Q [R;0], the pseudoinverse is P^+ = [R^{-1} 0] Q^T and hence the projection
        // P P^+ = Q [I 0; 0 0] Q^T. If Q = [Q1 Q2], then P P^+ = Q1 Q1^T holds. The image space of Q1 is the same as that
        // of P. For computing Q1, we therefore do not need to compute a full QR decomposition, but only the first columns of
        // Q. The most simple method is Gram-Schmidt. Though it is numerically instable if the columns of P are almost linear
        // dependent, we use it here, as the columns of P are associated to the shape functions on the father cell and 
        // should be comfortably linearly independent, though not orthogonal.
        // Instead of forming Q1 Q1^T, we store just Q1, and apply the product if needed. As Q1 is high and tall, this is 
        // (slightly) more efficient.
        for (int j=0; j<p.q.M(); ++j) // step through all columns
        {
          Scalar tmp = 0;
          
          // Normalize column j.
          for (int i=0; i<p.q.N(); ++i)
            tmp += power<2>(p.q[i][j]);
          tmp = std::sqrt(tmp);
          for (int i=0; i<p.q.N(); ++i)
            p.q[i][j] /= tmp;
          
          // Orthogonalize remaining columns.
          for (int k=j+1; k<p.q.M(); ++k)
          {
            tmp = 0;
            for (int i=0; i<p.q.N(); ++i) 
              tmp += p.q[i][j] * p.q[i][k];
            for (int i=0; i<p.q.N(); ++i)
              p.q[i][k] -= tmp*p.q[i][j];
          }
        }
      }

    private:
      CellPointer father; 
      std::vector<CellPointer> const& children; 
    };


    // A functor that applies the appropriate projector provided on construction to the 
    // FE function that is given as argument.
    template <class Projectors>
    struct ProjectCoefficients
    {
      // Construct the functor, giving a list of (local) projectors indexed by the FE space index.
      ProjectCoefficients(Projectors const& projectors_): projectors(projectors_) {}

      // Apply the appropriate projector to the provided FE function. The pair consists of a variable description
      // (containing the space index) and the FE function itself.
      template <class Pair>
      void operator()(Pair const& pair) const
      {
        typedef typename boost::remove_reference<typename result_of::value_at_c<Pair,0>::type>::type VarDesc;
        typedef typename boost::remove_reference<typename result_of::value_at_c<Pair,1>::type>::type Function;

        int const sIdx = VarDesc::spaceIndex;
        auto const& proj = at_c<sIdx>(projectors);
        int const n = proj.gidx.size();

        // Multiplication of projector Q Q^T and coefficient vector. The
        // multiplication has to be done manually because of
        // type-incompatible entries. This results in the coefficients of
        // the projection residual.
        auto const& q = proj.q;
        int const m = q.M();
        typename Function::StorageType x(n);
        typename Function::StorageType y(m);

        // y <- Q^T v
        for (int i=0; i<m; ++i) {
          y[i] = 0;
          for (int j=0; j<n; ++j)
            y[i].axpy(q[j][i][0][0],(at_c<1>(pair).coefficients())[proj.gidx[j]]);
        }
        
        // x <- Q y
        for (int i=0; i<n; ++i) {
          typename Function::StorageValueType x(0);
          for (int j=0; j<m; ++j)
            x.axpy(q[i][j][0][0],y[j]);
          (at_c<1>(pair).coefficients())[proj.gidx[i]] = x;
        }
      }

    private:
      Projectors const& projectors;
    };
    
    // Convenience function for template type deduction.
    template <class Projectors>
    ProjectCoefficients<Projectors> getCoefficientProjectors(Projectors const& ps) { return ProjectCoefficients<Projectors>(ps); }

    
    
    class GroupByCell
    {
    public:
      GroupByCell(std::vector<size_t> const& groups_)
      : groups(&groups_), nGroups(*std::max_element(groups_.begin(),groups_.end())+1) {}
      
      size_t operator[](size_t idx) const { return (*groups)[idx]; }
      size_t nGroups;
      
    private:
      std::vector<size_t> const* groups;
    };


  } // End of namespace CoarseningDetail


  /** \ingroup adapt
   * In order to maintain compatibility with existing code,
   * this overload is needed, as non const reference parameters (coarsenedCells) can not
   * have default values. All it does is call the original method with a temporary
   * parameter.
   */
  template <class VariableSetDescription, class Scaling>
  void coarsening(VariableSetDescription const& varDesc,
                  typename VariableSetDescription::VariableSet const& sol,
                  Scaling const& scaling,
                  std::vector<std::pair<double,double> > const& tol,
                  GridManager<typename VariableSetDescription::Grid>& gridManager, int verbosity=1, int minRefLevel=0)
  {
    std::vector<bool> tmp;
    coarsening(varDesc, sol, scaling, tol, gridManager, tmp, verbosity, minRefLevel);
  }

  /** \ingroup adapt
   * \brief coarsening routine
   *
   * Perform a projection from a fine grid onto a coarse grid. If the error is small, locally
   * coarsen the grid there.
   *
   * \param varDesc
   * \param sol 
   * \param scaling
   * \param tol 
   * \param gridManager
   * \param coarsenedCells
   * \param verbosity
   * \param minRefLevel do not coarse cells that are on this level or below.
   */
  template <class VariableSetDescription, class Scaling>
  void coarsening(VariableSetDescription const& varDesc,
                  typename VariableSetDescription::VariableSet const& sol,
                  Scaling const& scaling,
                  std::vector<std::pair<double,double> > const& tol,
                  GridManager<typename VariableSetDescription::Grid>& gridManager,
                  std::vector<bool>& coarsenedCells, int verbosity=1, int minRefLevel=0)
  {
    using namespace boost::fusion;
    using namespace CoarseningDetail;

// boost::timer::cpu_timer timer;
// boost::timer::cpu_timer timerA, timerB;
// timerA.stop(); timerB.stop();


    typedef typename VariableSetDescription::Grid           Grid;
    typedef typename Grid::LocalIdSet    IdSet;
    typedef typename Grid::LeafIndexSet  IndexSet;
    typedef typename Grid::template Codim<0>::Entity        Cell;
    typedef typename Grid::template Codim<0>::EntityPointer CellPointer;
    typedef typename Cell::HierarchicIterator               HierarchicIterator;

    IndexSet const& indexSet = varDesc.indexSet;

    // A set of FE functions that will contain the projection.
    typename VariableSetDescription::VariableSet pSol(sol);

    // A container of local projectors - defined here to avoid multiple allocations inside the cell loop
    auto projectors = as_vector(transform(varDesc.spaces,CreateProjection()));

    // a container for caching computed cell pointers to children. Prevents doubled hierarchical grid traversal 
    // and frequent reallocation
    std::vector<CellPointer> children;
    
    // A mapping from leaf cells to their coarsening group number. A
    // coarsening group is defined as follows. If all the direct
    // children of a non-leaf cell are leafs (i.e. none of the children
    // has been refined), these children form a coarsening group. Either
    // all or none cells of a coarsening group are marked for
    // coarsening. The special group number 0 collects all cells which
    // are not element of a regular coarsening group (and hence cannot be
    // removed during refinement).
    std::vector<size_t> coarseningGroup(indexSet.size(0),0);
    size_t coarseningGroupCount = 0;

    // Step through all leaf cells and process their fathers.
    typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator CellIterator;
    CellIterator end = gridManager.grid().template leafend<0>();
    for (CellIterator ci=gridManager.grid().template leafbegin<0>(); ci!=end; ++ci) // TODO: parallelize this loop
      // Father has already been processed if cell belongs to a regular
      // cell group. No father exists if the level is 0.
      if (coarseningGroup[indexSet.index(*ci)]==0 && ci->level()>0) { // not yet processed, not on coarse grid
        // Only process fathers which are suitable for coarsening.
        CellPointer father = ci->father();
        bool canBeCoarsened = true;
        children.clear();
        HierarchicIterator first = father->hbegin(father->level()+1), last = father->hend(father->level()+1);
        for (HierarchicIterator hi=first; canBeCoarsened && hi!=last; ++hi) // early termination in case we find a non-leaf child
        {
          canBeCoarsened = canBeCoarsened && hi->isLeaf();
          children.push_back(CellPointer(hi));
        }
        
        if (canBeCoarsened) {
// timerA.resume();
          // Gather all the children into a new coarsening group.
          ++coarseningGroupCount;
          for (auto const& c: children)
            coarseningGroup[indexSet.index(*c)] = coarseningGroupCount; // not coarseningGroupCount-1: group number 0 is for not coarsenable cells
          // Compute the projections to the locally coarser subspace.
          for_each(projectors,GetLocalTransferProjection<CellPointer>(father,children));
// timerA.stop();
// timerB.resume();          

          // Apply the projectors to the FE functions.
          for_each2(typename VariableSetDescription::Variables(),pSol.data,getCoefficientProjectors(projectors));
// timerB.stop();
        }
      }

    // Create the projection residual (i.e. the error introduced by coarsening).
    pSol -= sol;

//     std::cerr << "Coarsening part Ia time: " << timer.format()
//               << "       A               : " << timerA.format()
//               << "       B               : " << timerB.format();
//     timer.start();

    // Compute scaled L2 norms of the function and difference,
    // attributing the cell contributions to its coarsening group.
    ErrorestDetail::GroupedSummationCollector<GroupByCell> sum{GroupByCell(coarseningGroup)};
    scaledTwoNormSquared(join(typename VariableSetDescription::Variables(),typename VariableSetDescription::Variables()),
                         join(pSol.data,sol.data),varDesc.spaces,scaling,sum);
//     std::cerr << "Coarsening part Ib time: " << timer.format();
//     timer.start();

    // For each variable, compute the (squared) relative error
    // contribution, that is e_i^2 = |p_i|^2/(atol^2+|f|^2*rtol^2),
    // where f is the function and p the difference to its hierarchic
    // projection.
    std::vector<double> norm2(varDesc.noOfVariables,0);
    for (int i=0; i<varDesc.noOfVariables; ++i)
      for (int j=0; j<=coarseningGroupCount; ++j)
        norm2[i] += sum.sums[j][i+varDesc.noOfVariables];
    if ( verbosity>0 )
    {
      std::cout << "coarsening solution norms2: "; 
      std::copy(norm2.begin(),norm2.end(),
      std::ostream_iterator<double>(std::cout," ")); std::cout  << '\n';
    };

    for (int i=0; i<varDesc.noOfVariables; ++i) {
      double toli = power(tol[i].first,2) + norm2[i]*power(tol[i].second,2);
      for (int j=0; j<=coarseningGroupCount; ++j)
        sum.sums[j][i] /= toli;
    }
    

    // Select a maximal subset C of the coarsening groups, such that for
    // each variable the sum of C's relative errors is below 1.
    //
    // A greedy algorithm may work as follows: Maintain a vector E_i (of
    // size the number of variables) that contains the remaining allowed
    // relative error. Of the remaining coarsening groups j with
    // relative errors e_ji select the one for which max_i e_ji / E_i is
    // minimal.
    //
    // However, this has quadratic complexity. Thus we resort to the
    // following simplification.  First we assign to each coarsening
    // group j the maximal relative error e_j = max_i e_ji encountered
    // in any variable. Subsequently, these are sorted in ascending
    // order. Then the first k groups for which sum_{j=0}^k e_j <= 1
    // hold are selected. Non-coarsenable cells (in the irregular
    // coarsening group 0) are ignored.
    //
    // WARNING: Despite the fact that in principle we can compute the
    // projection error exactly, the value we obtain is just an
    // estimate. This is not only due to the heuristic suboptimal
    // selection of C, but also due to the fact that the coarsening
    // errors are modeled strictly locally. In conforming meshes, mesh
    // adaptation can lead to nonlocal effects, e.g. the removal of no
    // longer needed green closures (the introduced error is not taken
    // into account), the introduction of green closures on newly
    // coarsened cells (in which case the error could be below our
    // estimate, depending on the actual transfer - currently the
    // transfer is suboptimal and our error estimate should be exact),
    // or a cell group can actually be retained for mesh topology
    // reasons even if marked for coarsening (in which case the actual
    // local error would be zero).
    std::vector<std::pair<double,size_t> > e(coarseningGroupCount);
    for (int j=0; j<coarseningGroupCount; ++j) {
      e[j].second = j+1;                            // omit the irregular coarsening group 0
      for (int i=0; i<varDesc.noOfVariables; ++i)
        e[j].first = std::max(e[j].first,sum.sums[j+1][i]); 
    }

    std::sort(e.begin(),e.end(),FirstLess());

    std::vector<size_t> selectedCoarseningGroups;
    double totalError = 0;
    for (int i=0; i<e.size(); ++i) {
      totalError += e[i].first;
      if (totalError<=1)
        selectedCoarseningGroups.push_back(e[i].second);
      else
        break;
    }

    // for compression: need to keep track of coarsening history
    coarsenedCells.clear() ;
    coarsenedCells.resize( indexSet.size(0), false ) ;
    // no coarsening: just stop here
    //   return ;

    // If no cell can be coarsened, we don't need to touch the grid.
    if (selectedCoarseningGroups.empty())
      return;

    // Now we sweep over the whole grid and mark all cells for
    // coarsening which belong to a selected coarsening group.
    std::sort(selectedCoarseningGroups.begin(),selectedCoarseningGroups.end());
    for (CellIterator ci=gridManager.grid().template leafbegin<0>(); ci!=gridManager.grid().template leafend<0>(); ++ci)
    {
      auto idx = indexSet.index(*ci);
      if (ci->level()>minRefLevel && std::binary_search(selectedCoarseningGroups.begin(),selectedCoarseningGroups.end(),
                                                        coarseningGroup[idx])) // for uniform refinement: just delete level
      {
        coarsenedCells[idx] = true ; // for compression: coarsening history
        gridManager.mark(-1,*ci);
      }
    }

//     std::cerr << "Coarsening part II time: " << timer.format();
//     timer.start();

    gridManager.adaptAtOnce();
    
//     std::cerr << "Coarsening part III time: " << timer.format();
  }
} /* end of namespace Kaskade */

#endif
