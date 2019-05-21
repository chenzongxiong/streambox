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

#ifndef FETRANSFER_HH
#define FETRANSFER_HH


  /**
   * @file
   * @brief  Tools for transfer of data between grids
   * @author Martin Weiser, Anton Schiela
   */
#include <limits>
#include <memory>       // std::unique_ptr
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <iterator>
#include <utility>      // std::move


#include "dune/grid/common/grid.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/common/fmatrix.hh"
#include "dune/grid/io/file/dgfparser/dgfparser.hh"
//#include "dune/grid/geometrygrid.hh" // why should this be needed?

#include "fem/fixdune.hh"
#include "fem/mllgeometry.hh"
#include "linalg/dynamicMatrix.hh"

namespace Kaskade
{
  template<class G, class T> class CellData;
  template <class Space, class CP> class TransferData;
  struct AdaptationCoarseningPolicy ;


  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \ingroup fetransfer
   * \brief Computes the values of global ansatz functions at the
   * given points (which are supposed to be local coordinates inside the
   * cell to which the evaluator is currently moved).
   *
   * This realizes the matrix product \f$\Psi\Phi \f$ (see
   * LocalToGlobalMapperConcept). On return, \arg sfValues[i][j]
   * contains the value of global shape function j at global(x[i]).
   *
   * See also \ref approximateGlobalValues.
   */
  template <class Space, class ShapeFunctionSet>
  void evaluateGlobalShapeFunctions(Space const& space,
                                    typename Space::Grid::template Codim<0>::Entity const& cell,
                                    std::vector<Dune::FieldVector<typename Space::Grid::ctype,Space::Grid::dimension> > const& x,
                                    DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,Space::sfComponents,1> >& sfValues,
                                    ShapeFunctionSet const& shapeFunctionSet)
  {

    sfValues.setSize(x.size(),shapeFunctionSet.size());

    // Evaluation and transformation of shape function values. This
    // realizes \Phi.
    shapeFunctionSet.evaluate(x,sfValues);

    // Convert the local to global values. This realizes \Psi
    typename Space::Mapper::Converter psi(cell);
    for (int i=0; i<sfValues.N(); ++i) 
    {
      psi.setLocalPosition(x[i]);
      for (int j=0; j<sfValues.M(); ++j)
        sfValues[i][j] = psi.global(sfValues[i][j]);
    }
  }

  template <class Space>
  void evaluateGlobalShapeFunctions(Space const& space,
                                    typename Space::Grid::template Codim<0>::Entity const& cell,
                                    std::vector<Dune::FieldVector<typename Space::Grid::ctype,Space::Grid::dimension> > const& x,
                                    DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,Space::sfComponents,1> >& sfValues)
  {
    evaluateGlobalShapeFunctions(space,cell,x,sfValues,space.mapper().shapefunctions(cell));
  }
  
  /**
   * \ingroup fetransfer
   * \brief Computes global shape function coefficients such
   * that the given set of global function values at the interpolation
   * nodes of the shape function set are approximated by the resulting
   * linear combinations of global shape functions.
   *
   * This realizes the matrix product \f$ \Phi^{+} \Psi^{-1} \f$
   * (see LocalToGlobalMapperConcept).
   *
   * On exit, the value of globalValues is undefined.
   *
   * The actual meaning of "approximate" (formally denoted by a
   * pseudoinverse symbol here) is defined by the specific FE
   * space. However, for any space it holds that after executing
   *
   * \code
   * evaluateGlobalShapeFunctions(space,cell,space.mapper().shapefunctions(cell).interpolationNodes(),sfValues);
   * approximateGlobalValues(space,cell,sfValues,coeff);
   * \endcode
   *
   * the matrix coeff is the identity matrix.
   * 
   * \param[in] globalValues global function values at interpolation nodes of the shape function set (will be overwritten and invalidated)
   * \param[out] coeff coefficients of shape function linear combination such that the global values at the interpolation points are approximated
   * \param sfs the space's shape function set on the given cell (note that the cell needs not be contained in the space's grid view at all, e.g. when
   *            creating multilevel prolongations)
   */
  template <class Space, class ShapeFunctionSet>
  void approximateGlobalValues(Space const& space,
                               typename Space::Grid::template Codim<0>::Entity const& cell,
                               DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,Space::sfComponents,1>>& globalValues,
                               DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,1,1> >& coeff,
                               ShapeFunctionSet const& sfs)
  {
    // Transform global values to local values. This is \Psi^{-1}.
    typename Space::Mapper::Converter psi(cell);
    
    auto const& iNodes = sfs.interpolationNodes();

    for (int i=0; i<globalValues.N(); ++i) 
    {
      psi.setLocalPosition(iNodes[i]);
      for (int j=0; j<globalValues.M(); ++j)
        globalValues[i][j] = psi.local(globalValues[i][j]);
    }
    // Interpolate local values, giving shape function coefficients. This is \Phi^+.
    sfs.interpolate(globalValues,coeff);
  }

  /**
   * \brief 
   * \param cell the cell on which to interpolate. This has to be contained in the space's grid view.
   */
  template <class Space>
  void approximateGlobalValues(Space const& space,
                               typename Space::Grid::template Codim<0>::Entity const& cell,
                               DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,Space::sfComponents,1>>& globalValues,
                               DynamicMatrix<Dune::FieldMatrix<typename Space::Scalar,1,1> >& coeff)
  {
    approximateGlobalValues(space,cell,globalValues,coeff,space.mapper().shapefunctions(cell));
  }

  /**
   * \brief Computes a local prolongation matrix from father cell to child cell. 
   * 
   * In the notation of LocalToGlobalMapperConcept, this
   * can be written as \f[ p \leftarrow \Phi_c^+ \Psi_c^{-1} \Psi_f \Phi_f. \f]
   *
   * Although suggestively called father and child, the only geometric
   * relation required is that all interpolation nodes of the child cell
   * are also contained inside the father cell.
   *
   * Note that the FE spaces on father and child may be different,
   * providing interpolation matrices for different FE
   * approximations. The requirement is that the scalar type and the
   * number of shape function components coincide.
   * 
   * \param childSfs a shape function set to use for the child space on the given child cell (which need not be contained in the child space's grid view)
   * \param fatherSfs a shape function set to use for the father space on the given father cell (which need not be contained in the fathers space's grid view)
   */
  template <class ChildSpace, class FatherSpace>
  void localProlongationMatrix(ChildSpace const& childSpace,
                               typename ChildSpace::Grid::template Codim<0>::EntityPointer const& child,
                               FatherSpace const& fatherSpace,
                               typename FatherSpace::Grid::template Codim<0>::EntityPointer const& father,
                               DynamicMatrix<Dune::FieldMatrix<typename ChildSpace::Scalar,1,1> >& prolongation,
                               typename ChildSpace::Mapper::ShapeFunctionSet const& childSfs,
                               typename FatherSpace::Mapper::ShapeFunctionSet const& fatherSfs)
  {
    typedef typename ChildSpace::Grid Grid;
    typedef typename ChildSpace::Scalar   Scalar;
    int const sfComponents = ChildSpace::sfComponents;

    // Obtain the child interpolation nodes and map them to the father.
    std::vector<Dune::FieldVector<typename Grid::ctype,Grid::dimension> > iNodes(childSfs.interpolationNodes());

    std::vector<int> dummy(iNodes.size());
    MultiLevelLocalGeometry<Grid> const mlGeometry(childSpace.grid(),child,father,MultiLevelLocalGeometry<Grid>::ChildIsGlobal);
    mlGeometry.local(iNodes,dummy);
    assert(dummy.size()==iNodes.size()); // \todo: we assume the child is contained in the father. Is this guaranteed?

    // Evaluate global shape functions on father cell.
    DynamicMatrix<Dune::FieldMatrix<Scalar,sfComponents,1>> afValues;
    evaluateGlobalShapeFunctions(fatherSpace,*father,iNodes,afValues,fatherSfs);

    // Interpolate global values on the child cell.
    approximateGlobalValues(childSpace,*child,afValues,prolongation,childSfs);
  }

  /**
   * \ingroup fetransfer
   * \brief Computes a local prolongation matrix from father cell to child cell. 
   * 
   * \param child a pointer to a cell contained in the child space's grid view
   * \param father a pointer to a cell contained in the father space's grid view
   */
  template <class ChildSpace, class FatherSpace>
  void localProlongationMatrix(ChildSpace const& childSpace,
                               typename ChildSpace::Grid::template Codim<0>::EntityPointer const& child,
                               FatherSpace const& fatherSpace,
                               typename FatherSpace::Grid::template Codim<0>::EntityPointer const& father,
                               DynamicMatrix<Dune::FieldMatrix<typename ChildSpace::Scalar,1,1> >& prolongation)
  {
    localProlongationMatrix(childSpace,child,fatherSpace,father,prolongation,
                            childSpace.mapper().shapefunctions(*child),fatherSpace.mapper().shapefunctions(*father));
  }
  
  /**
   * \ingroup fetransfer
   * These two structures are used for determining whether to perform a
   * "normal" grid transfer (default, AdaptionCoarseningPolicy, same
   * behaviour as before), or to calculate prolongation matrices for
   * a FE function space based on a HierarchicIndexSet (see mgtools.hh).
   */
  struct AdaptationCoarseningPolicy
  {
    void update( int ) { }

    template <class Cell>
    bool accept(Cell const& cell) const { return cell.isLeaf(); }

    template <class Cell>
    bool operator()(Cell const& cell) const { return cell.mightBeCoarsened(); }
  };

  struct MultilevelCoarseningPolicy
  {
    MultilevelCoarseningPolicy(int level_):
      level(level_)
    {}

    void update( int level_ ) { level = level_ ; }

    template <class Cell>
    bool accept(Cell const& cell) const { return cell.level()==level; }

    template <class Cell>
    bool operator()(Cell const& cell) const { /*return cell.level()==level;*/return false ; }

  private:
    int level;
  };

  /**
   * \ingroup fetransfer
   * \brief Class that stores for a cell ("this entity") a local restriction matrix.
   *
   * That is the matrix \f[ R = \Phi^+ \Psi^{-1} D^{-1}
   * \sum_{\textrm{leaf children } c} \Psi_c \Phi_c K_c. \f] 
   * Here, \f$ D \f$ is a diagonal matrix with \f$ D_{ii} \f$ equal to the number
   * of leaf child cells contributing to interpolation node \f$ i \f$.
   *
   * Note that this matrix maps global degrees of freedom (i.e. global
   * ansatz function coefficients) to local shape function
   * coefficients on the father cell.
   *
   * If this is a leaf cell, this reduces to \f$ R = K \f$.
   * 
   * \tparam ChildSpace  a FEFunctionSpace defined on the leaf grid view
   * \tparam FatherSpace a FEFunctionSpace 
   * \tparam CoarseningPolicy DOCME
   */
  template <class ChildSpace, class FatherSpace, class CoarseningPolicy=AdaptationCoarseningPolicy>
  class LocalTransfer
  {
  public:
    typedef typename FatherSpace::Grid            Grid;
    typedef typename FatherSpace::IndexSet        IndexSet;
    typedef typename FatherSpace::Scalar              Scalar;
    typedef typename Grid::GlobalIdSet::IdType Id;
    typedef Dune::FieldVector<Scalar,FatherSpace::sfComponents>   SfValue;
    typedef typename Grid::template Codim<0>::EntityPointer       CellPointer;

    typedef DynamicMatrix<Dune::FieldMatrix<Scalar,1,1> >          LocalTransferMatrix;

    enum { RESTRICTION=1, PROLONGATION=2 };

    LocalTransfer(LocalTransfer&& other)
    : globalIndices(std::forward<std::vector<size_t> >(other.globalIndices)), prolongationMatrix(std::forward<LocalTransferMatrix>(other.prolongationMatrix)),
      restrictionMatrix(std::forward<LocalTransferMatrix>(other.restrictionMatrix)), shapeFunctionSet(std::move(other.shapeFunctionSet))
    {}

    LocalTransfer(LocalTransfer const& other)
    : globalIndices(other.globalIndices), prolongationMatrix(other.prolongationMatrix), restrictionMatrix(other.restrictionMatrix),
      shapeFunctionSet(new typename FatherSpace::Mapper::ShapeFunctionSet(*other.shapeFunctionSet))
    {}

    /**
     * \param flags computes restriction only if flags|RESTRICTION==true, similarly for prolonagtion
     */
    LocalTransfer(ChildSpace const& childSpace, FatherSpace const& fatherSpace, CellPointer father, int flags,
                  CoarseningPolicy const& acceptancePolicy = CoarseningPolicy() )
    {
      // TODO: Higher performance solution (for special case, though) has been developed in coarsening.hh.
      //       See whether this can be used to improve the performance here.
      Grid const& grid = childSpace.grid();
      assert(&grid == &fatherSpace.grid());

      // Collect all leaf children descending from the father cell
      std::vector<CellPointer> children;
      if ( acceptancePolicy.accept(*father) )
        children.push_back(father);
      else {
        auto end = father->hend(grid.maxLevel());
        for (auto hi=father->hbegin(grid.maxLevel()); hi!=end; ++hi)
          if(acceptancePolicy.accept(*hi))
            children.push_back(CellPointer(hi));
      }

      // Compute the global indices associated to the leaf cells that
      // are descendants from the father cell.
      for (int i=0; i<children.size(); ++i) {
        auto gi = childSpace.mapper().globalIndices(*children[i]);
        globalIndices.insert(globalIndices.end(),gi.begin(),gi.end());
      }
      std::sort(globalIndices.begin(),globalIndices.end());
      globalIndices.erase(std::unique(globalIndices.begin(),globalIndices.end()),globalIndices.end());

      // Initialize local prolongation and restriction data.
      prolongationMatrix.setSize(globalIndices.size(),fatherSpace.mapper().shapefunctions(*father).size());
//     prolongationMatrix = 0 // causes compilation error with dune-2.4.0 and clang++ on OS X (Darwin)
      prolongationMatrix.fill(0);
      typedef typename FatherSpace::Mapper::ShapeFunctionSet::SfValueArray SfValueArray;
      SfValueArray restrictionData;
      restrictionData.setSize(fatherSpace.mapper().shapefunctions(*father ).interpolationNodes().size(),globalIndices.size());
//     restrictionData = 0 // causes compilation error with dune-2.4.0 and clang++ on OS X (Darwin)
      restrictionData.fill(0);


      // Compute local prolongation and restriction matrices.
      // First declare variables that allocate memory in front of the loop to prevent frequent reallocations.
      std::vector<int> fCount(restrictionData.N(),0);
      LocalTransferMatrix p; 
      SfValueArray afValues; 
      std::vector<size_t> globIdx;
      typedef Dune::FieldVector<typename Grid::ctype,Grid::dimension> Position;
      std::vector<Position> iNodes;
      std::vector<int> nodesInside;
      
      for (int i=0; i<children.size(); ++i) {

        // Compute the mapping of degrees of freedom on the child to the set of dofs on all children.
        typename ChildSpace::Mapper::GlobalIndexRange gi = childSpace.mapper().globalIndices(*children[i]);
        globIdx.resize(gi.size());
        for (int j=0; j<gi.size(); ++j) {
          globIdx[j] = std::lower_bound(globalIndices.begin(),globalIndices.end(),gi[j]) - globalIndices.begin();
          assert(globIdx[j]<globalIndices.size());
        }

        // Compute and scatter prolongation
        localProlongationMatrix(childSpace,children[i],fatherSpace,father,p,childSpace.mapper().shapefunctions(*children[i]),
                                                                            childSpace.mapper().shapefunctions(*children[i])); // this makes only sense if the father and child
        for(size_t k=0; k<globIdx.size(); ++k)                                                                                 // shape functions are the same...
          for(int j=0; j<p.M(); ++j)
            prolongationMatrix[globIdx[k]][j] = p[k][j];


        // Compute and scatter restriction
        // Create a geometry mapping from the child coordinate system (local) to myself (global).
        MultiLevelLocalGeometry<Grid> mlGeometry(grid,children[i],father,MultiLevelLocalGeometry<Grid>::FatherIsGlobal,false);
        shapeFunctionSet.reset(new typename FatherSpace::Mapper::ShapeFunctionSetImplementation(static_cast<typename FatherSpace::Mapper::ShapeFunctionSetImplementation const&>(fatherSpace.mapper().shapefunctions(*father))));
        iNodes = fatherSpace.mapper().shapefunctions(*father).interpolationNodes();
        nodesInside.resize(iNodes.size());
        for (int j=0; j<nodesInside.size(); ++j)
          nodesInside[j] = j;
        mlGeometry.local(iNodes,nodesInside);

        evaluateGlobalShapeFunctions(childSpace,*children[i],iNodes,afValues);

        childSpace.mapper().combiner(*children[i],childSpace.indexSet().index(*children[i])).rightTransform(afValues);

        for (int k=0; k<afValues.N(); ++k) {
          int rIdx = nodesInside[k];
          assert(0<=rIdx && rIdx<fCount.size());
          fCount[rIdx]++;
          for (int j=0; j<p.M(); ++j) {
            assert(!std::isnan(restrictionData[rIdx][globIdx[j]]));
            restrictionData[rIdx][globIdx[j]] += afValues[k][j];
            assert(!std::isnan(restrictionData[rIdx][globIdx[j]]));
          }
        }
      }

      // Do averaging of restriction of discontinuous FE functions.
      for (int i=0; i<restrictionData.N(); ++i)
        for (int j=0; j<restrictionData.M(); ++j) {
          assert(!std::isnan(restrictionData[i][j]));
          restrictionData[i][j] /= fCount[i];
          assert(!std::isnan(restrictionData[i][j]));
        }

      // Interpolate father shape functions.
      approximateGlobalValues(fatherSpace,*father,restrictionData,restrictionMatrix);
    }

    std::vector<size_t> globalIndices;
    LocalTransferMatrix prolongationMatrix;
    LocalTransferMatrix restrictionMatrix;
    std::unique_ptr<typename FatherSpace::Mapper::ShapeFunctionSetImplementation> shapeFunctionSet;
  };


  //---------------------------------------------------------------------


  /**
   * A class storing data collected before grid adaptation that is
   * necessary to transfer FE functions. This includes global indices of
   * shape functions on leaf cells and restriction matrices to coarser
   * cells in case some leaf cells might be removed during grid
   * adaptation.
   *
   * This class is used by GridManager to perform the data
   * transfer. Usually, this class is not needed by any other client
   *
   * \tparam Space FE function space for which the transfer data is to be gathered
   * \tparam CoarseningPolicy
   */
  template <class Space, class CoarseningPolicy=AdaptationCoarseningPolicy>
  class TransferData {

  public:
    typedef typename Space::Grid            Grid;
    typedef typename Space::GridView        GridView;
    typedef typename Space::IndexSet        IndexSet;
    typedef typename Space::Scalar              Scalar;
    typedef typename Grid::GlobalIdSet::IdType Id;
    typedef Dune::FieldVector<Scalar,Space::sfComponents>             SfValue;
    typedef typename Grid::template Codim<0>::EntityPointer       CellPointer;

  private:
    typedef DynamicMatrix<Dune::FieldMatrix<Scalar,1,1> >    RestrictionMatrix;
    typedef LocalTransfer<Space,Space,CoarseningPolicy> RestrictionData;
    typedef std::map<Id, RestrictionData>               HistoryMap;

  public:

    /**
     * Matrix that transforms a data vector v_1 corresponding to the old
     * grid to a data vector v_2 corresponding to the new grid, such
     * that the functions represented by these vectors coincide as much
     * as possible.
     */
    class TransferMatrix
    {
    public:
      /** Transforms oldCoeff, which lives on the old grid to an equivalent
       * vector that lives on the new grid*/
      template <class StorageValue>
      std::unique_ptr<Dune::BlockVector<StorageValue> > apply(Dune::BlockVector<StorageValue> const& oldCoeff) const
      {
        std::unique_ptr<Dune::BlockVector<StorageValue> > newCoeff(new Dune::BlockVector<StorageValue>(tm.size()));
        for (int i=0; i<tm.size(); ++i) {
          StorageValue tmp(0.0);
          for (int j=0; j<tm[i].first.size(); ++j)
            tmp += tm[i].first[j].second * oldCoeff[tm[i].first[j].first];
          (*newCoeff)[i] = tmp;
        }
        return newCoeff;
      }

    private:

      friend class TransferData;

      TransferMatrix(size_t rows, size_t cols):
        tm(rows), nrows(rows), ncols(cols)
      {
        for (int i=0; i<tm.size(); ++i){
          tm[i].second = std::numeric_limits<int>::max(); // sentinel: restriction level
          tm[i].first.clear();
        }
      }

      template <class ColIter, class EntryIter>
      void storeRow(int row, int restrictionLevel, ColIter firstCol, ColIter lastCol,
          EntryIter firstEntry)
      {
        if (tm[row].second > restrictionLevel) {
          tm[row].second = restrictionLevel;
          tm[row].first.clear();
          tm[row].first.reserve(std::distance(firstCol,lastCol));
          while (firstCol != lastCol) {
            if (std::abs(*firstEntry)>1e-14)
              tm[row].first.push_back(std::make_pair(*firstCol,*firstEntry));
            ++firstCol;
            ++firstEntry;
          }
        }
      }

      template <class ColIter, class EntryIter>
      void storeRowWeighed(int row, double weight, ColIter firstCol, ColIter lastCol,
          EntryIter firstEntry)
      {
        tm[row].first.reserve(std::distance(firstCol,lastCol)+tm[row].first.size());
        while (firstCol != lastCol) {
          //  Commented out: find out, if element already exists in row. If yes: add its value
          //                                                             If no: insert new element
          //  Current (probably more efficient) implementation: always insert new element
          //          this saves searching time, but leads to some memory overhead
          //
          //      int i=0;
          //      tm[row].first.reserve(lastCol-firstCol+tm[row].first.size());
          //       while(i<tm[row].first.size() && tm[row].first[i].first != *firstCol)
          //       {
          //         ++i;
          //       }
          //       if(i!=tm[row].first.size())
          //         tm[row].first[i].second+=(*firstEntry)*weight;
          //       else {
          if (std::abs(*firstEntry) > 1e-14)
            tm[row].first.push_back(std::make_pair(*firstCol,(*firstEntry)*weight));
          //      }

          ++firstCol;
          ++firstEntry;
        }
      }

      void scaleRow(int row, double weight)
      {
        for(int i=0; i<tm[row].first.size(); ++i)
          tm[row].first[i].second *= weight;
      }

    public:
      void out(std::ostream& o) const
      {
        o << "shape=[" << nrows << ' ' << ncols << "];\n";
        o << "data=[\n";
        for (int i=0; i<tm.size(); ++i)
          for (int j=0; j<tm[i].first.size(); ++j)
            o << i+1 << ' ' << tm[i].first[j].first+1 << ' ' << tm[i].first[j].second << '\n';
        o << "];\n";

      }

      int nRows() const { return nrows; }
      int nCols() const { return ncols; }

    private:
      // Data for TransferMatrix - a vector with one entry for each row
      // tm.first: sparse column-vector <index,value>  tm.second: globalIdx
      std::vector<std::pair<std::vector<std::pair<int,Scalar> >,int> > tm;
      int nrows,ncols;
    }; // end of class TransferMatrix



    /**
     * Gathers all information that might be lost after the refinement
     * process. To be called after preAdapt() and before adapt().
     */
    TransferData(Space const& space_, CoarseningPolicy const& mightBeCoarsened = CoarseningPolicy() ):
      space(space_), DOFcoarseSpace(space.degreesOfFreedom())
    {
      typename Grid::GlobalIdSet const& idSet = space.grid().globalIdSet();
      //     typedef typename IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;
      typedef typename GridView::template Codim<0>::Iterator CellIterator;

      typename Space::Evaluator evaluator(space);

      // Gather all information that is necessary to create local restriction matrices
      CellIterator cend = space.gridView().template end<0>();
      for (CellIterator ci=space.gridView().template begin<0>(); ci!=cend; ++ci) {

        // obtain and store only global indices of ansatz functions for all cells
        historyData.insert(std::make_pair(idSet.id(*ci),RestrictionData(space,space,CellPointer(ci),0,mightBeCoarsened)));


        // lower level cells might be important for grid transfer
        // todo@ !(entity->mightBeCoarsened) does not mean that this entity will exist after a grid change.
        // e.g. if an entity with the same father is marked for coarsening, then there are problems.
        CellPointer ancestor(ci);
        while (ancestor->mightVanish() && ancestor->level()>0) {
          ancestor = ancestor->father();

          typename HistoryMap::iterator rd = historyData.find(idSet.id(*ancestor));
          if (rd==historyData.end()) {
            // Restriction data for ancestor does not yet exist -> create
            rd = historyData.insert(historyData.begin(),
                std::make_pair(idSet.id(*ancestor),RestrictionData(space,space,ancestor,RestrictionData::RESTRICTION,mightBeCoarsened)));
          }
        }
      }
    }


    /**
     * \brief Create a TransferMatrix.
     *
     * Its application transfers a data vector from a coarse
     * grid to a fine grid. To be called after adapt() and
     * before postAdapt(). In this version we assume that we
     * have a transfer between two different grids, while
     * the type of space stays the same.
     *
     * For each new leaf cell \f$ c \f$, we find the closest ancestor
     * \f$ f \f$ with a restriction matrix \f$ R \f$ mapping global
     * ansatz function coefficients to ancestor shape function
     * coefficients. Building a prolongation matrix \f$ P \f$ mapping
     * ancestor shape function coefficients to leaf shape function
     * coefficients, the local transfer matrix can be written as \f[ T =
     * K_c^+ P R. \f]
     */
    std::unique_ptr<TransferMatrix> transferMatrix() const
      {
      typename Grid::GlobalIdSet const& idSet = space.grid().globalIdSet();
      std::unique_ptr<TransferMatrix> transfer(new TransferMatrix(space.degreesOfFreedom(),DOFcoarseSpace));

      typedef typename GridView::template Codim<0>::Iterator CellIterator;

      CellIterator cend = space.gridView().template end<0,Dune::All_Partition>();
      for (CellIterator ci=space.gridView().template begin<0,Dune::All_Partition>(); ci!=cend; ++ci) {
        // Among the ancestors of the current cell that have a
        // restriction matrix stored, we select the one with the highest
        // refinement level (i.e. the most exact representation of the
        // FE functions before mesh adaptation.
        CellPointer ancestor(ci);
        typename HistoryMap::const_iterator restriction;
        int colevel = 0;
        while((restriction=historyData.find(idSet.id(*ancestor))) == historyData.end() && ancestor->level()> 0){
          colevel++;
          ancestor = ancestor->father();
        }
        assert(restriction != historyData.end());
        RestrictionMatrix const& rm = restriction->second.restrictionMatrix;

        auto index = space.indexSet().index(*ci);

        RestrictionMatrix prolongation;
        // WARNING: having the same shape function set for child and father makes only sense if these are the same...
        //localProlongationMatrix(space,CellPointer(ci),space,ancestor,prolongation,*(restriction->second.shapeFunctionSet),*(restriction->second.shapeFunctionSet));
        // Then let's try the following.
        localProlongationMatrix(space,CellPointer(ci),space,ancestor,prolongation,space.mapper().shapefunctions(index),*(restriction->second.shapeFunctionSet));

        // compute matrix product prolongation * restriction
        RestrictionMatrix localTransfer;
        MatMult(localTransfer,prolongation,rm);

        // apply K^+
        space.mapper().combiner(*ci,index).leftPseudoInverse(localTransfer);

        // scatter the local transfer
        typename Space::Mapper::GlobalIndexRange rowIdx = space.mapper().globalIndices(index);
        std::vector<size_t> const& columnIdx = restriction->second.globalIndices;
        //assert(!columnIdx.empty());
        for (int i=0; i<localTransfer.N(); ++i)
        {
          transfer->storeRow(rowIdx[i], colevel,
              columnIdx.begin(), columnIdx.end(),
              localTransfer[i].begin());
        }
      }
      return transfer;
      }

    /// Create a TransferMatrix.
    /**
     * Its application transfers a Data Vector
     * from a coarse grid to a fine grid. To be called after adapt() and
     * before postAdapt(). In this version we assume that we have a transfer
     * between two different grids, and two different types of spaces, e.g. globally
     * continuous and discontinous spaces. A transfer from discontinuous to continuous
     * spaces is performed by local averaging.
     */
    template<class Space2>
    std::unique_ptr<TransferMatrix> transferMatrix(Space2 const& space2) const
    {
      typename Grid::GlobalIdSet const& idSet = space2.grid().globalIdSet();
      std::unique_ptr<TransferMatrix> transfer(new TransferMatrix(space2.degreesOfFreedom(),DOFcoarseSpace));

      std::vector<typename Grid::ctype> volume(space2.degreesOfFreedom(),0.0);
      typename Grid::ctype myvolume = 0.0;

      typedef typename GridView::template Codim<0>::Iterator CellIterator;

      CellIterator cend = space2.gridView().template end<0,Dune::All_Partition>();
      for (CellIterator ci=space2.gridView().template begin<0,Dune::All_Partition>();
          ci!=cend; ++ci) {
        auto index = space2.indexSet().index(*ci);
        //     for (CellIterator ci=space2.indexSet().template begin<0,Dune::All_Partition>();
        //   ci!=space2.indexSet().template end<0,Dune::All_Partition>(); ++ci) {
        CellPointer ancestor(ci);
        typename HistoryMap::const_iterator restriction;
        int colevel = 0;
        while((restriction=historyData.find(idSet.id(*ancestor))) == historyData.end() && ancestor->level()> 0){
          colevel++;
          ancestor = ancestor->father();
        }
        assert(restriction != historyData.end());
        RestrictionMatrix const& rm = restriction->second.restrictionMatrix;

        myvolume = 1.0/ci->geometry().volume();
        for (int i=0; i<space2.mapper().globalIndices(*ci).size(); ++i)
          volume[space2.mapper().globalIndices(*ci)[i]] += myvolume;

        RestrictionMatrix prolongation;
        localProlongationMatrix(space,CellPointer(ci),space,ancestor,prolongation,restriction.shapeFunctionSet,restriction.shapeFunctionSet);
        // WARNING: using the same shape function set for child and father spaces makes only sense if these are the same

        // compute matrix product prolongation * restriction
        RestrictionMatrix localTransfer;
        MatMult(localTransfer,prolongation,/*restriction*/rm);

        // apply K^+
        space.mapper().combiner(*ci,index).leftPseudoInverse(localTransfer);

        // scatter the local transfer
        typename Space2::Mapper::GlobalIndexRange rowIdx = space2.mapper().globalIndices(index);
        std::vector<size_t> const& columnIdx=restriction->second.globalIndices;
        for (int i=0; i<localTransfer.N(); ++i)
          transfer->storeRowWeighed(rowIdx[i],myvolume,
              columnIdx.begin(), columnIdx.end(),
              localTransfer[i].begin());
      }

      for(int i=0; i<volume.size(); ++i) transfer->scaleRow(i,1.0/volume[i]);

      return transfer;
    }

  private:
    Space const&  space;
    HistoryMap    historyData;
    int           DOFcoarseSpace;
  };


  /**
   * \ingroup fem
   * \brief Transfer a FunctionSpaceElement that lives on one space to a FunctionSpaceElement that lives on another space
   * The transfer from discontinuous to continuous spaces is performed
   * by volume weighed averaging
   */
  template <class Fu1, class Fu2>
  void spaceTransfer(Fu1& f1, Fu2 const& f2)
  {
    typedef TransferData<typename Fu2::Space> TrData;

    TrData trdata(f2.space());
    std::unique_ptr<typename TrData::TransferMatrix> trMat=trdata.transferMatrix(f1.space());


    (*f1) = *(trMat->apply(*f2));
  }

  
  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup fetransfer
   * \brief A Weighing that associates to each cell its inverse volume.
   */
  struct InverseVolume
  {
    template<class Entity>
    double operator()(Entity const& e) {return 1.0/e.geometry().volume();}
    double scalingFactor() {return 1.0;}
    double offset() {return 0.0;}
  };

  /**
   * \ingroup fetransfer
   * \brief A Weighing that associates to each cell its volume.
   */
  struct Volume
  {
    template<class Entity>
    double operator()(Entity const& e) {return e.geometry().volume();}
    double scalingFactor() {return 1.0;}
    double offset() {return 0.0;}
  };

  /**
   * \ingroup fetransfer
   * \brief A Weighing that associates to each cell a constant weight of 1.
   */
  struct PlainAverage
  {
    template<class Entity>
    double operator()(Entity const& e) {return 1.0;}
    double scalingFactor() {return 1.0;}
    double offset() {return 0.0;}
  };

  template<int scl=1>
  struct SumUp
  {
    template<class Entity>
    double operator()(Entity const& e) {return 1.0;}
    double scalingFactor() {return 0.0;}
    double offset() {return scl;}
  };


  /**
   * \ingroup fetransfer
   * \brief Interpolates FunctionSpaceElement to FunctionSpaceElement
   *
   * This is particulary useful together with FunctionViews and makes
   * it possible to transform a FunctionView to a FunctionSpaceElement
   *
   * \param fse function space element which is set to (approximately) the value of fu
   * \param fu function space element of a (probably different) FE space
   *
   * \tparam Weighing has to be default constructible
   * \tparam FSElement type of the target finite element function 
   * \tparam Function type of the source function 
   */
  template<typename Weighing=PlainAverage, typename FSElement, typename Function>
  void interpolateGlobally(FSElement& fse, Function const& fu)
  {
    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;

    DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar,ImageSpace::sfComponents,1> > globalValues;
    DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, 1, 1 > > localCoefficients;
    std::vector<typename Grid::ctype> sumOfWeights(fse.space().degreesOfFreedom(),0.0);

    // Initialize the target to zero.
    fse.coefficients() = typename FSElement::Scalar(0.0);

    typename ImageSpace::Evaluator isfs(fse.space());
    typename Function::Space::Evaluator dsfs(fu.space());
    Weighing w;
    
    using ValueType = decltype(fu.value(dsfs));
    std::vector<ValueType> fuvalue; // function values - declare here to prevent reallocations

    // Step through all the cells and perform local interpolation on each cell.
    for (auto const& cell: elements(fse.space().gridView())) 
    {
      auto index = fse.space().indexSet().index(cell);
      
      isfs.moveTo(cell);
      dsfs.moveTo(cell);

      typename Grid::ctype myWeight = w(cell);
      for (int i=0; i<isfs.size(); ++i) 
        sumOfWeights[isfs.globalIndices()[i]] += myWeight;

      // Obtain interpolation nodes of target on this cell
      auto const& iNodes(isfs.shapeFunctions().interpolationNodes());

      // Evaluate source
      globalValues.setSize(iNodes.size(),1);
      localCoefficients.setSize(iNodes.size(),1);
      fuvalue.resize(iNodes.size());
      for (int i=0; i<iNodes.size(); ++i) 
      {
        dsfs.evaluateAt(iNodes[i]);
        fuvalue[i] = fu.value(dsfs);
      }

      // Perform local interpolation
      for(int k=0; k<FSElement::components/ImageSpace::sfComponents; ++k) 
      {
        for (int i=0; i<iNodes.size(); ++i)
          for(int j=0; j<ImageSpace::sfComponents; ++j)
            globalValues[i][0][j][0] = myWeight*fuvalue[i][ImageSpace::sfComponents*k+j];

        approximateGlobalValues(fse.space(),cell,globalValues,localCoefficients);

        assert(localCoefficients.N() == isfs.globalIndices().size());

        // transform global values onto shape function coefficients
        fse.space().mapper().combiner(cell,index).leftPseudoInverse(localCoefficients);

        // add up the shape function coefficients to the ansatz function coefficients
        for (int i=0; i<localCoefficients.N(); ++i) 
          fse.coefficients()[isfs.globalIndices()[i]][k] += localCoefficients[i][0];
      }
    }

    // On vertices, edges, and faces, degrees of freedom may have received
    // updates from multiple cells. Average this out.
    for(int i=0; i<sumOfWeights.size(); ++i)
      fse.coefficients()[i] /= (sumOfWeights[i]*w.scalingFactor()+w.offset());
  }



  /**
   * \ingroup fetransfer
   * \brief Interpolates WeakFunctionViews to FunctionSpaceElement
   *
   * This is particulary useful together with WeakFunctionViews and makes
   * it possible to transform a WeakFunctionView to a FunctionSpaceElement.
   * 
   * Remember that WeakFunctionView s support the evaluation by cell and local coordinate.
   * If the source function is a finite element function, prefer interpolateGlobally.
   */
  template<typename Weighing=PlainAverage, typename FSElement, typename Function>
  void interpolateGloballyWeak(FSElement& fse, Function const& fu)
  {
    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;

    DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;
    DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, 1, 1>> localCoefficients;
    std::vector<typename Grid::ctype> sumOfWeights(fse.space().degreesOfFreedom(),0.0);
    fse.coefficients() = typename ImageSpace::Scalar(0.0);

    typename ImageSpace::Evaluator isfs(fse.space());
    Weighing w;


    auto const cend = fse.space().gridView().template end<0>();
    using ValueType = decltype(fu.value(*cend,Dune::FieldVector<typename Grid::ctype, ImageSpace::dim>()));
    std::vector<ValueType> fuvalue; // declare here to prevent reallocations
    for (auto ci=fse.space().gridView().template begin<0>(); ci!=cend; ++ci)
    {
      auto index = fse.space().indexSet().index(*ci);
      
      isfs.moveTo(*ci);
      typename Grid::ctype myWeight = w(*ci);
      for (int i=0; i<isfs.size(); ++i) 
        sumOfWeights[isfs.globalIndices()[i]] += myWeight;

      auto const& iNodes(isfs.shapeFunctions().interpolationNodes());

      globalValues.setSize(iNodes.size(),1);
      localCoefficients.setSize(iNodes.size(),1);
      fuvalue.resize(iNodes.size());
      for (int i=0; i<iNodes.size(); ++i)
        fuvalue[i] = fu.value(*ci,iNodes[i]);

      for(int k=0; k<FSElement::components/ImageSpace::sfComponents; ++k)
      {
        for (int i=0; i<iNodes.size(); ++i)
          for(int j=0; j<ImageSpace::sfComponents; ++j)
            globalValues[i][0][j][0] = myWeight*fuvalue[i][ImageSpace::sfComponents*k+j];

        approximateGlobalValues(fse.space(),*ci,globalValues,localCoefficients);

        assert(localCoefficients.N() == isfs.globalIndices().size());
        fse.space().mapper().combiner(*ci,index).leftPseudoInverse(localCoefficients);
        for (int i=0; i<localCoefficients.N(); ++i)
          fse.coefficients()[isfs.globalIndices()[i]][k]+= localCoefficients[i][0];
      }
    }
    
    for(int i=0; i<sumOfWeights.size(); ++i)
      fse.coefficients()[i] /= sumOfWeights[i]*w.scalingFactor()+w.offset();
  }


  
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup fetransfer
   * \brief An adaptor that allows to turn lambda functions into function views.
   * \tparam Space a FEFunctionSpace on which the underlying function lives
   * \tparam Functor a callable object (e.g. a lambda function) mapping space evaluators to values
   */
  template <class Space_, class Functor>
  class FunctionViewAdaptor
  {
  public:
    using Space = Space_;
    
    FunctionViewAdaptor(Space const& space, Functor const& g_): g(g_), sp(space)
    {}
    
    Space const& space() const
    {
      return sp;
    }
    
    template <class ...Args>
    auto value(Args... args) const 
    {
      return g(std::forward<Args>(args)...);
    }
    
    
  private:
    Functor g;
    Space const&  sp;
  };
  
  /**
   * \ingroup fetransfer
   * \brief A convenience functor that supports the easy creation of function view adaptors.
   * \tparam Space a FEFunctionSpace type 
   * \tparam Functor a callable object with operator()(Space::Evaluator eval)
   * \relates FunctionViewAdaptor
   */
  template <class Space, class Functor>
  auto makeFunctionView(Space const& space, Functor const& g)
  {
    return FunctionViewAdaptor<Space,Functor>(space,g);
  }
  
  /**
   * \ingroup fetransfer
   * \brief A convenience functor that supports the easy creation of function view adaptors.
   * \tparam Space a FEFunctionSpace type 
   * \tparam Functor a callable object with operator()(Cell, LocalPosition)
   * \relates FunctionViewAdaptor
   */
  template <class Functor>
  auto makeWeakFunctionView(Functor const& g)
  {
    return FunctionViewAdaptor<Functor,Functor>(g,g); // this is misusing FunctionViewAdaptor slightly...
  }
  
} /* end of namespace Kaskade */
#endif
