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

#ifndef MULTIVARIATIONAL_ASSEMBLER_HH
#define MULTIVARIATIONAL_ASSEMBLER_HH

#include <iostream>
#include <numeric>
#include <memory>
#include <tuple>
#include <new>
#include <type_traits>

#include <boost/version.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/timer/timer.hpp>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/find_if.hpp>

#include "dune/istl/bcrsmatrix.hh"
#include "dune/istl/bdmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/common/fmatrix.hh"
#include "dune/geometry/type.hh"
#include "dune/grid/common/capabilities.hh"

#include "fem/firstless.hh"
#include "fem/fixfusion.hh"
#include "fem/functional_aux.hh"
#include "fem/gridmanager.hh"
#include "fem/quadrature.hh"
#include "fem/variables.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/threadedMatrix.hh"

#include "utilities/threading.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/typeTraits.hh"


// this is omitted in boost/fusion/sequence/container/vector/detail/vector_n.hpp
#undef N

namespace Kaskade
{

  /**
   * @file
   * @brief  Compute local and global stiffness matrices for variational problems with multiple variables.
   * @author Martin Weiser
   */

  static int const timerStatistics = 0x1;       // gather global execution time statistics
  static int const scatterStatistics = 0x2;
  static int const lockStatistics = 0x4;
  static int const localTimerStatistics = 0x0;  // gather local (per element) time statistics
  /**
   * \brief flags for gathering performance and timing statistics
   */
  static int const statistics = 0;

  //---------------------------------------------------------------------

  /**
   * \cond internals
   *
   * This namespace contains template metaprogramming details used in
   * the VariationalFunctionalAssembler.
   */
  namespace AssemblyDetail {

    using namespace boost::fusion;


    /**
     * \brief A class that supports exclusive access to several row groups of a matrix by providing appropriate locking.
     */
    class RowGroupManager
    {
    public:
      /**
       * \brief Default constructor creates a (pretty useless) empty row group manager
       */
      RowGroupManager();

      /**
       * \brief Copy constructor.
       */
      RowGroupManager(RowGroupManager const& m);

      /**
       * \brief Destructor
       */
      virtual ~RowGroupManager();

      /**
       * \brief Assignment.
       */
      RowGroupManager& operator=(RowGroupManager const& m);

      /**
       * \brief initializes the row group manager to row groups of approximately the same number of rows
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] rows the total number of rows
       */
      void init(int nrg, size_t rows);

      /**
       * \brief initializes the row group manager to row groups of approximately the same number of matrix entries
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] colIndex a container of column indices for matrix entries
       */
      void init(int nrg, std::vector<std::vector<size_t> > const& colIndex);

      /**
       * \brief initializes the row group manager to row groups of approximately the same number of matrix entries
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] rowSize a container of row sizes
       */
      void init(int nrg, std::vector<size_t> const& rowSize);

      /**
       * \brief Requests exclusive access to row group \a n.
       * \return the half open range of row indices belonging to this group.
       */
      std::pair<size_t,size_t> lock(int n);

      /**
       * \brief Releases the exclusive access to row group \a n.
       */
      void unlock(int n);

      /**
       * \brief Returns the number of row groups in this matrix block.
       */
      int size() const;


    private:
      // start indices of row groups
      std::vector<size_t> rowGroupStart;

      // Support for simultaneous write access to different row groups.
#ifndef KASKADE_SEQUENTIAL
      boost::mutex* mutexes;
#endif
    };

    // ----------------------------------------------------------------------------------------------------

    /**
     * \brief A class that describes a matrix block and holds the associated global matrix.
     *
     * \tparam Policy \todo docme
     * \tparam RowVar \todo docme
     * \tparam ColVar \todo docme
     *
     * Despite holding the global stiffness matrix by a shared pointer,
     * the copy and assignment are deep (i.e. expensive).
     */
    template <class Policy, class RowVar, class ColVar>
    struct MatrixBlock
    {
      static int const rowId = RowVar::id;
      static int const colId = ColVar::id;
      static int const rowSpaceIndex = RowVar::spaceIndex;
      static int const colSpaceIndex = ColVar::spaceIndex;
      static bool const symmetric = Policy::template BlockInfo<rowId,colId>::symmetric;
      static bool const mirror    = Policy::template BlockInfo<rowId,colId>::mirror;
      static bool const lumped    = Policy::template BlockInfo<rowId,colId>::lumped;

      typedef MatrixBlock<Policy,RowVar,ColVar> Self;
      typedef typename Policy::Scalar Scalar;
      typedef typename Policy::Spaces Spaces;
      typedef Dune::FieldMatrix<Scalar,RowVar::m,ColVar::m> BlockType;
      typedef NumaBCRSMatrix<BlockType> Matrix;

      typedef typename SpaceType<Spaces,RowVar::spaceIndex>::type RowSpace;
      typedef typename SpaceType<Spaces,ColVar::spaceIndex>::type ColSpace;
      typedef typename RowSpace::GridView GridView;
      typedef typename GridView::template Codim<0>::Iterator CellIterator;
//       typedef typename GridView::IntersectionIterator FaceIterator;

      static int const dim = ColSpace::dim;

      /**
       * \brief Default constructor creates a (pretty useless) empty matrix block.
       */
      MatrixBlock() {}

      /**
       * \brief Deep copy constructor.
       */
      MatrixBlock(MatrixBlock const& mb): matrix(new Matrix(*mb.matrix)), isDense(mb.isDense), rowGroups(mb.rowGroups) {}

      /**
       * \brief Deep copy.
       */
      MatrixBlock& operator=(MatrixBlock const& mb)
      {
        *matrix = *mb.matrix;
        isDense = mb.isDense;
        rowGroups = mb.rowGroups;
      }

      MatrixBlock& operator=(Scalar x)
      {
        assert(matrix.get());
        *matrix = x;
        return *this;
      }

      MatrixBlock& operator+=(MatrixBlock const& mb)
      {
        *matrix += *mb.matrix;
        return *this;
      }


      /**
       * \brief Provides read access to the global matrix.
       */
      Matrix& globalMatrix() { assert(matrix.get()); return *matrix; }

      /**
       * \brief Provides write access to the global matrix.
       */
      Matrix const& globalMatrix() const { assert(matrix.get()); return *matrix; }

      /**
       * \brief Provides persistent read access to the matrix.
       *
       * The matrix remains
       * valid even on destruction of the assembler or calling flush(),
       * but may be overwritten by the next assembly call.
       */
      std::shared_ptr<Matrix const> globalMatrixPointer() const { return matrix; }


      /**
       * \brief Computes the sparsity pattern of the matrix block. For each cell
       * (codim 0 entity) of the mesh, all shape functions on the cell are
       * assumed to interact, hence their associated global indices induce
       * a scattering of nonzero entries into the matrix.
       *
       * For symmetric matrix blocks, only the lower triangular part is stored.
       *
       * \param nrg number of row groups (should be in the order of magnitude of the number of threads)
       */
      void init(Spaces const& spaces, CellIterator first, CellIterator last, int nrg)
      {
        RowSpace const& rowSpace = *at_c<RowVar::spaceIndex>(spaces);
        ColSpace const& colSpace = *at_c<ColVar::spaceIndex>(spaces);

        // @todo: respect *dynamic* presence flag. static one is respected
        // on construction of sequence of matrix blocks. Respect symmetric flag.
        if (!Policy::template BlockInfo<rowId,colId>::present) {
          assert("Aieee! Nonpresent matrix block detected!\n"==0);
          abort();
        }

        boost::timer::cpu_timer timer;

        size_t nnz = 0;
        size_t const rows = rowSpace.degreesOfFreedom();
        size_t const cols = colSpace.degreesOfFreedom();

        // If mass lumping is desired, we just create a diagonal matrix.
        // Otherwise, we need to establish the connectivity pattern of
        // the ansatz functions.
        if (Policy::template BlockInfo<rowId,colId>::lumped)
        {
          assert(RowVar::spaceIndex==ColVar::spaceIndex && rows==cols);
          NumaCRSPatternCreator<> creator(rows,cols,false,1);
          for (size_t i=0; i<rows; ++i)
            creator.addElements(&i,&i+1,&i,&i+1);
          matrix.reset(new NumaBCRSMatrix<BlockType>(creator)); // zero init
          nnz = rows;

          // Initialize row groups with roughly the same number of rows in each group
          rowGroups.init(nrg,rows);
        }
        else
        {
          // In case of global support of shape functions of one of the involved spaces (usually ConstantSpace),
          // the resulting matrix is dense. Use a shortcut to define the sparsity structure.
          bool const dense = RowSpace::Mapper::globalSupport || ColSpace::Mapper::globalSupport;


          // Allocate sparse matrix and define the sparsity structure. Estimate the required number of entries per row
          // for efficient preallocation.
          // Ws 2016-01-22: Zero nnz init is actually faster on P1 pattern, probably because otherwise the
          // chunks put a high burden on the memory system by requesting memory at the same time. This appears
          // to parallelize quite badly and hence is turned off in parallel mode.
          int preallocateEntriesPerRow = 0;
          if (NumaThreadPool::instance().nodes()==1)                       // sequential allocation in NumaCRSPatternCreator
          {
            int localSize = colSpace.mapper().globalIndices(0).size();     // # entries per row of local matrices (rough estimate)
            preallocateEntriesPerRow = dense? 0:                           // for dense matrices we'll use a specialized version
                                       dim<=1? 2*localSize:                // make room for (roughly) all elemental matrices appended...
                                       dim==2? 8*localSize:                // ... this is overprovisioning, but faster due to less ...
                                               20*localSize;               // ... allocations and data movement.
          }
          NumaCRSPatternCreator<> creator(rows,cols,symmetric,preallocateEntriesPerRow);

          if (statistics & timerStatistics)
            std::cout << "initial creator construction (" << rowId << "," << colId << "):  " << timer.format();
          timer.start();


          if (dense)
          {
            creator.addAllElements();
            if (statistics & timerStatistics)
              std::cout << "entering indices for (" << rowId << "," << colId << "):        " << timer.format();
            timer.start();
          }
          else
          {
            // Traditional local FE spaces - that's the hard case. We have to traverse the grid and pick up all
            // interactions of (local) basis functions.

            // Iteration over all cells. Note that the grid views of row and
            // column spaces need not cover the whole grid (e.g., in case of
            // local support). It would be best to iterate over the grid
            // view with smaller support, but for simplicity we choose the
            // row space.

            typename RowSpace::Evaluator rsfs(rowSpace);
            typename ColSpace::Evaluator csfs(colSpace);

            std::vector<std::vector<size_t>> rowIndices, colIndices;

            size_t const nCells = rowSpace.indexSet().size(0);

            std::vector<typename RowSpace::Mapper::GlobalIndexRange> rIndices(nCells,rowSpace.mapper().initGlobalIndexRange());
            std::vector<typename ColSpace::Mapper::GlobalIndexRange> cIndices(nCells,colSpace.mapper().initGlobalIndexRange());
            for (size_t i=0; i<nCells; ++i)
            {
              rIndices[i] = rowSpace.mapper().globalIndices(i);
              cIndices[i] = colSpace.mapper().globalIndices(i);
            }
            if (statistics & timerStatistics)
              std::cout << "gathering indices for (" << rowId << "," << colId << "):         " << timer.format();
            timer.start();

            creator.addElements(rIndices,cIndices);

            if(Policy::considerInnerFaces)
            {
              std::cout << "considering inner faces" << std::endl;
              typename RowSpace::Evaluator rsfs(rowSpace);
              typename ColSpace::Evaluator csfs(colSpace);

              for (CellIterator ci=first; ci!=last; ++ci)
              {
                // obtain the shape functions
                rsfs.moveTo(*ci);
                csfs.moveTo(*ci);

                typename RowSpace::Evaluator rsfs2(rowSpace);
                typename ColSpace::Evaluator csfs2(colSpace);

                GridView const& gridView = rowSpace.gridView();
                auto fend = gridView.iend(*ci);
                for(auto fit = gridView.ibegin(*ci); fit!=fend; ++fit)
                {
                  if(fit->boundary() || !fit->neighbor()) continue;
                  if(gridView.indexSet().index(*ci) > gridView.indexSet().index(*fit->outside())) continue;

                  csfs2.moveTo(*fit->outside());
                  // Append column indices to affected rows.
                  creator.addElements(rsfs.globalIndices().begin(), rsfs.globalIndices().end(),
                                      csfs2.globalIndices().begin(),csfs2.globalIndices().end());

                  rsfs2.moveTo(*fit->outside());
                  // Append column indices to affected rows.
                  creator.addElements(rsfs2.globalIndices().begin(),rsfs2.globalIndices().end(),
                                      csfs.globalIndices().begin(), csfs.globalIndices().end());
                } // end iteration over faces
              } // end iteration over cells
            } // end Policy::considerInnerFaces

            if (statistics & timerStatistics)
              std::cout << "entering indices for (" << rowId << "," << colId << "):          " << timer.format();
            timer.start();

            creator.balance();
            if (statistics & timerStatistics)
              std::cout << "balancing for (" << rowId << "," << colId << "):                 " << timer.format();
            timer.start();
          }

          // Now that the pattern creator is filled, create pattern and matrix.

          // Count number of nonzeros.
          nnz = creator.nonzeroes();

          // Create the matrix (zero initialization of entries).
          matrix.reset(new Matrix(creator));
          if (statistics & timerStatistics)
            std::cout << "creating pattern & matrix for (" << rowId << "," << colId << "): " << timer.format();
          timer.start();

          // Initialize row groups with roughly the same number of entries in each group
          std::vector<size_t> colIndex(matrix->N());
          for (auto ri=matrix->begin(); ri!=matrix->end(); ++ri)
            colIndex[ri.index()] = ri->size();
          rowGroups.init(nrg,colIndex);
          if (statistics & timerStatistics)
            std::cout << "creating row groups for (" << rowId << "," << colId << "):       " << timer.format();
          timer.start();
        }

        // Check density (nonzeros>(rows*cols/2). Be careful not to get
        // overflow for very large and sparse matrices!
        if (std::min(rows,cols)>0)
          isDense = nnz/std::min(rows,cols) > std::max(rows,cols)/2;
        else
          isDense = false; // well, empty matrices are as sparse as one can hope for ;)

        if (statistics & timerStatistics)
            std::cout << "init cleanup for (" << rowId << "," << colId << "):              " << timer.format() << "\n";
      }

      // Returns true if the matrix is essentially dense (more than 50% nonzeros).
      bool dense() const { return isDense; }

      /**
       * \brief Provides access to the row group management
       */
      RowGroupManager& rowGroup() { return rowGroups; }

      // forward declaration
      struct LocalMatrices;

      // A structure for combining local matrices with its local and global row and column indices
      struct LocalMatrix
      {
        static bool const lumped = MatrixBlock<Policy,RowVar,ColVar>::lumped;

        // just because the ranges are not default constructible
        LocalMatrix(): ridx_(RowSpace::Mapper::initSortedIndexRange()), cidx_(ColSpace::Mapper::initSortedIndexRange()), data(nullptr) {}

        typedef typename RowSpace::Mapper::SortedIndexRange SortedRowIdx;
        typedef typename ColSpace::Mapper::SortedIndexRange SortedColIdx;

        /**
         * \brief A sequence of (global row, local row) indices.
         */
        SortedRowIdx ridx() const { return ridx_; }

        /**
         * \brief A sequence of (global col, local col) indices.
         */
        SortedColIdx cidx() const { return cidx_; }

        BlockType& operator()(int row, int col)
        {
          // For lumped matrix blocks, only the diagonal of local matrices are computed
          // and accessed. Store these not as a diagonal of a full matrix, but in a
          // contiguous vector.
          if (lumped)
          {
            assert(row==col);
            return data[row];
          }
          else
            return data[row*cidx().size()+col];
        }

        BlockType const& operator()(int row, int col) const
        {
          // For lumped matrix blocks, only the diagonal of local matrices are computed
          // and accessed. Store these not as a diagonal of a full matrix, but in a
          // contiguous vector.
          if (lumped)
          {
            assert(row==col);
            return data[row];
          }
          else
            return data[row*cidx().size()+col];
        }

      private:
        friend class LocalMatrices;

        template <class RowEvaluator, class ColEvaluator>
        void init(RowEvaluator const& reval, ColEvaluator const& ceval, BlockType* data_)
        {
          ridx_ = reval.sortedIndices();
          cidx_ = ceval.sortedIndices();
          data = data_;
        }

        SortedRowIdx ridx_;
	SortedColIdx cidx_;

        // a pointer to the raw storage for the matrix entries
        BlockType* data;
      };

      // A structure for holding a sequence of several local matrices to be filled sequentially
      // and to be scattered together.
      struct LocalMatrices
      {
        typedef Self MatrixBlock;
        typedef DynamicMatrix<BlockType> LocalMatrixType;

        /**
         * \param n the maximal number of local matrices to store in this chunk
         * \param maxStorage_ the desired maximum size of the occupied memory. This is no hard limit, merely a hint. Defaults to 256kB.
         */
        LocalMatrices(int n, MatrixBlock& mb_, size_t maxStorage_=256*1024)
        : localMatrices(n), currentLocalMatrix(0), lastLocalMatrix(-1), mb(&mb_), maxStorage(maxStorage_)
        {
          // prevent frequent reallocation at startup. A (conservative) guess is (dim+1) dofs per cell (P1 on
          // triangles). Except for constant or at least piecewise constant functions it's hard to
          // imagine smaller local matrices.
          localData.reserve(std::min( static_cast<size_t>(lumped?(dim+1):(dim+1)*(dim+1))*n, 11*maxStorage/10  ));
        }

        // Storage for local stiffness matrices.
        std::vector<LocalMatrix> localMatrices;

        // Which is the current local matrix to be filled?
        int currentLocalMatrix, lastLocalMatrix;

        // provides access to the global matrix and its properties.
        MatrixBlock* mb;


        /**
         * \brief initializes the current local matrix
         */
        template <class RowEvaluator, class ColEvaluator>
        void init(RowEvaluator const& reval, ColEvaluator const& ceval)
        {
          ++lastLocalMatrix;
          assert(localMatrices[lastLocalMatrix].data == nullptr); // make sure it's not yet initialized

          // for lumped matrices, only the diagonal is stored in a contiguous vector
          size_t s = lumped ? reval.size() : reval.size()*ceval.size();
          assert(reval.size()==ceval.size() || !lumped);

          //          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (localData.capacity()<localData.size()+s)
          {
            // not enough space available: reserve more. As the pointers into the vector
            // are invalidated by the reallocation, we have to adjust the data pointers
            // of the previous blocks.
            localData.reserve(localData.size()+s);

            size_t offset = 0;
            for (int i=0; i<lastLocalMatrix; ++i)
            {
              localMatrices[i].data = &localData[offset];
              offset += lumped ? localMatrices[i].ridx().size() : localMatrices[i].ridx().size()*localMatrices[i].cidx().size();
            }
            assert(offset==localData.size());
          }

          // Now we're sure there is enough space in the data array. Extend it by zeros.
          size_t offset = localData.size();

          localData.insert(localData.end(),s,BlockType(0));
          localMatrices[lastLocalMatrix].init(reval,ceval,&localData[offset]);
        }

        /**
         * \brief resets the local matrices to an empty state
         */
        void clear()
        {
          localData.clear();
          assert(lastLocalMatrix<(int)localMatrices.size());
          for (int i=0; i<=lastLocalMatrix; ++i)
            localMatrices[i].data = nullptr;
          currentLocalMatrix = 0;
          lastLocalMatrix = -1;
        }

        /**
         * \brief reports the size of the local matrices storage in bytes
         *
         * This can be used to limit the number of local matrices such that their memory fits into the CPU cache.
         */
        size_t storageSize() const
        {
          return localData.size() * sizeof(BlockType);
        }

        /**
         * \brief reports the maximal desired local matrices storage size
         */
        size_t storageSizeLimit() const
        {
          return maxStorage;
        }

      private:
        // This allows to hold a single array into which the local matrices are scattered, and
        // thus improves memory locality (for better cache hit rates and less false sharing).
        // Contains the local matrices one after the other, for lumped matrices, only storage
        // for the diagonal is provided.
        std::vector<BlockType>                  localData;      // entries of elemental matrices
        size_t maxStorage;
      };

    private:
      std::shared_ptr<Matrix> matrix;
      bool                    isDense;
      RowGroupManager         rowGroups;
    };

    // ---------------------------------------------------------------------------------------

    template <int row, int col>
    struct IsBlock
    {
      template <class MatrixBlock>
      struct apply
      {
        using type = std::conditional_t<MatrixBlock::rowId==row && MatrixBlock::colId==col,boost::mpl::true_,boost::mpl::false_>;
      };
    };

    // ---------------------------------------------------------------------------------------

    /**
     * \brief defines a flat sequence of matrix blocks
     */
    template <class Policy, class AnsatzVariables, class TestVariables>
    class BlockArray
    {
      // Retain only those blocks which are present in the functional
      struct PresenceFilter
      {
        template <class Block>
        struct apply
        {
          typedef boost::mpl::bool_<Policy::template BlockInfo<Block::rowId,Block::colId>::present> type;
        };
      };

      static auto create(AnsatzVariables a, TestVariables t)
      {
        auto joiner = [=](auto blocks, auto rowVar)
        {
          using RowVar = decltype(rowVar);
          auto thisRow = transform(a,[](auto avar){ return MatrixBlock<Policy,RowVar,decltype(avar)>(); });
          return join(blocks,thisRow);
        };
        auto allBlocks = accumulate(t,vector<>(),joiner);
        return as_vector(filter_if<PresenceFilter>(allBlocks));
      }

    public:
      typedef decltype(create(AnsatzVariables(),TestVariables())) type;
    };


    //---------------------------------------------------------------------

    /**
     * \brief A class that defines one block of the right hand side
     */
    template <class Policy, class RV>
    struct RhsBlock
    {
      typedef typename std::remove_reference<RV>::type RowVar;

      static int const rowId = RowVar::id;

      typedef typename Policy::Scalar Scalar;
      typedef typename Policy::Spaces Spaces;
      typedef Dune::FieldVector<Scalar,RowVar::m> BlockType;
      typedef typename Dune::BlockVector<BlockType> Rhs;
      typedef RhsBlock<Policy,RowVar> Self;

      static int const rowSpaceIndex = RowVar::spaceIndex;
      typedef typename SpaceType<Spaces,RowVar::spaceIndex>::type RowSpace;
      typedef typename RowSpace::GridView::template Codim<0>::Iterator CellIterator;


      /**
       * \brief Initializes the rhs block info.
       */
      void init(Spaces const& spaces, CellIterator /*first*/, CellIterator /*last*/, int nrg)
      {
        RowSpace const& rowSpace = *at_c<RowVar::spaceIndex>(spaces);
        size_t rows = rowSpace.degreesOfFreedom();

        rowGroupManager.init(nrg,rows);
      }

      RowGroupManager& rowGroup() { return rowGroupManager; }

      // A class storing a local element rhs vector and the corresponding indices
      struct LocalVector
      {
        // just because the ranges are not default constructible
        LocalVector(): ridx(RowSpace::Mapper::initSortedIndexRange()) {}

        std::vector<BlockType> vector;

        typedef typename RowSpace::Mapper::SortedIndexRange SortedRowIdx;

        /**
         * \brief A container of (global row, local row) indices.
         */
        SortedRowIdx ridx;
      };

      // A structure for holding a sequence of several local vectors to be filled sequentially
      // and to be scattered together.
      struct LocalVectors
      {
        typedef Self RhsBlock;
        typedef std::vector<BlockType> LocalVectorType;

        LocalVectors(int n, RhsBlock& rb_): localVectors(n), currentLocalVector(0), rb(&rb_) {}

        // Storage for local rhs vectors.
        std::vector<LocalVector> localVectors;

        // Which is the current local vector to be filled?
        int currentLocalVector;

        // Pointing back to its rhs block
        RhsBlock* rb;
      };

    private:
      RowGroupManager rowGroupManager;
    };

    /**
     * \brief Defines a flat sequence of present rhs blocks
     */
    template <class Policy, class TestVariables>
    class RhsBlockArray
    {
      // Define a mapping from Variable to RhsBlock for Variable as row variable.
      // Retain only those blocks which are present in the functional
      struct PresenceFilter
      {
        template <class Block>
        struct apply
        {
          typedef boost::mpl::bool_<Policy::template RhsBlockInfo<Block::rowId>::present> type;
        };
      };

      static auto create(TestVariables t)
      {
        auto allBlocks = transform(t,[](auto tvar) { return RhsBlock<Policy,decltype(tvar)>(); });
        return as_vector(filter_if<PresenceFilter>(allBlocks));
      }

    public:
      typedef decltype(create(TestVariables())) type;
    };

    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that for a given test variable defines the type of local rhs blocks.
     * \TODO rewrite this to use multiple local rhs vectors to be scattered simultaneously (as already done for matrix blocks)
     */
    template <class Policy>
    struct RhsLocalData
    {
      /**
       * \param[in] n_ desired number of local vectors to complete before scattering
       */
      RhsLocalData(int n_): n(n_) {}

      template <class RhsBlock>
      typename RhsBlock::LocalVectors operator()(RhsBlock const& rb) const
      {
        // due to boost::fusion::transform providing only const references to the source sequence elements, we
        // need a const cast here :(
        return typename RhsBlock::LocalVectors(n,removeConst(rb));
      }

    private:
      int n;
    };

    /**
     * \brief A boost::fusion functor that defines for each matrix block the local matrices storage for assembly.
     */
    struct MatrixLocalData
    {
      /**
       * \brief Constructor.
       * \param n_ The number of local element stiffness matrices in each block.
       * \param s_ The memory size in bytes the local matrices should not significantly exceed.
       */
      MatrixLocalData(int n_, size_t s_): n(n_), s(s_) {}

      template <class MatrixBlock>
      typename MatrixBlock::LocalMatrices operator()(MatrixBlock const& mb) const
      {
        // When switching from boost 1.51 to 1.57, it seems that the filter_view only presents
        // const references to the elements (why?), but we need non-const ones. Taking non-const
        // reference as arguments leads to failure of overload resolution, omitting this
        // functor, and ultimately to error. As a workaround, we take const reference and cast
        // it away.
        return typename MatrixBlock::LocalMatrices(n,removeConst(mb),s);
      }

    private:
      int n;
      size_t s;
    };

    //---------------------------------------------------------------------

    template <class Spaces>
    class GetMaxDerivativeBase
    {
    public:
      GetMaxDerivativeBase(Spaces const& spaces_, std::map<void const*,int>& deriv_): spaces(spaces_), deriv(deriv_) {}

    private:
      Spaces const&        spaces;
      std::map<void const*,int>& deriv;

    protected:
      template <int spaceIndex>
      void enterSpace(int d) const
      {
        void const* space = at_c<spaceIndex>(spaces);
        auto i = deriv.find(space);
        if (i==deriv.end())
          deriv[space] = d;
        else
          i->second = std::max(i->second,d);
      }
    };

    template <class Functional, class Spaces>
    class GetMaxDerivativeMatrix: public GetMaxDerivativeBase<Spaces>
    {
    public:
      GetMaxDerivativeMatrix(Spaces const& spaces, std::map<void const*,int>& deriv): GetMaxDerivativeBase<Spaces>(spaces,deriv)  {}

      template <class MBlock>
      void operator()(MBlock const&) const
      {
        int d = Functional::template D2<MBlock::MatrixBlock::rowId,MBlock::MatrixBlock::colId>::derivatives;
        this->template enterSpace<MBlock::MatrixBlock::rowSpaceIndex>(d);
        this->template enterSpace<MBlock::MatrixBlock::colSpaceIndex>(d);
      }
    };


    template <class Functional, class Spaces>
    class GetMaxDerivativeRhs: public GetMaxDerivativeBase<Spaces>
    {
    public:
      GetMaxDerivativeRhs(Spaces const& spaces, std::map<void const*,int>& deriv): GetMaxDerivativeBase<Spaces>(spaces,deriv)  {}

      template <class RBlock>
      void operator()(RBlock const&) const
      {
        int d = Functional::template D1<RBlock::RhsBlock::rowId>::derivatives;
        this->template enterSpace<RBlock::RhsBlock::rowSpaceIndex>(d);
      }
    };


    /**
     * \brief Computes the maximum occuring derivative for each space.
     *
     * \returns a map from space address to the highest occuring derivative order
     */
    template <class Functional, class Spaces, class MBlocks, class RBlocks>
    std::map<void const*,int> derivatives(Functional const& f, Spaces const& spaces, MBlocks const& mblocks, RBlocks const& rblocks)
    {
      std::map<void const*,int> deriv; // map from space address to maximum occuring derivative

      for_each(mblocks,GetMaxDerivativeMatrix<Functional,Spaces>(spaces,deriv));
      for_each(rblocks,GetMaxDerivativeRhs<Functional,Spaces>(spaces,deriv));
      return deriv;
    }

    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that, given a test variable, initializes the current local rhs vector tor correct size and 0 values.
     */
    template <class Evaluators, class TestVariables>
    struct ClearLocalRhs
    {
      ClearLocalRhs(Evaluators const& evaluators_):
        eval(evaluators_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        // Deletes all entries from vector and inserts the needed number
        // of zero elements. Thus, finite elements with different number
        // of shape functions are supported. Memory reallocations do
        // rarely occur.
        assert(lv.currentLocalVector <  lv.localVectors.size());

        int const rSpaceIdx = result_of::value_at_c<TestVariables,LocalVectors::RhsBlock::rowId>::type::spaceIndex;

        auto& reval = at_c<rSpaceIdx>(eval);                    // evaluator

        auto& rc = lv.localVectors[lv.currentLocalVector];        // current local vector structure

        // set current local vector to correct size and value 0
        rc.vector.resize(reval.size());
        std::fill(rc.vector.begin(),rc.vector.end(),0);

        // retrieve sorted indices
        rc.ridx = reval.sortedIndices();
      }

      Evaluators const&  eval;
    };

    /**
     * \brief A boost::fusion functor that initializes the current local stiffness matrix to zero.
     */
    template <class Evaluators, class AnsatzVariables, class TestVariables>
    struct ClearLocalMatrix
    {
      ClearLocalMatrix(Evaluators const& reval_, Evaluators const& ceval_)
        : reval(reval_), ceval(ceval_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        assert(lm.lastLocalMatrix+1<lm.localMatrices.size()+1);

        typedef typename LocalMatrices::MatrixBlock MatrixBlock;

        int const cSpaceIdx = result_of::value_at_c<AnsatzVariables,MatrixBlock::colId>::type::spaceIndex;
        int const rSpaceIdx = result_of::value_at_c<TestVariables,MatrixBlock::rowId>::type::spaceIndex;

        lm.init(at_c<rSpaceIdx>(reval),at_c<cSpaceIdx>(ceval));
      }

      Evaluators const& reval;
      Evaluators const& ceval;
    };


    /**
     * \brief A boost::fusion functor that updates the variable currentLocalMatrix in BlockMatrix::LocalMatrices.
     *
     * When all local matrices corresponding to a cell have been assembled and scattered we have to set
     * the variable currentLocalMatrix to the next entry in the buffer that will be filled, i.e.
     * currentLocalMatrix = lastLocalMatrix+1
     */
    struct SetCurrentLocalMatrix
    {
      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        lm.currentLocalMatrix = lm.lastLocalMatrix+1;
      }
    };

    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that evaluates local right hand side contributions
     */
    template <class TestVariables, class Evaluators, class Real, class Cache>
    struct UpdateLocalRhs
    {
      UpdateLocalRhs(Evaluators const& evaluators_, Real integrationFactor_, Cache const& cache_):
        evaluators(evaluators_), integrationFactor(integrationFactor_), cache(cache_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        typedef typename LocalVectors::RhsBlock RhsBlock;

        int const rSpaceIndex = result_of::value_at_c<TestVariables,  RhsBlock::rowId>::type::spaceIndex;
        auto& rc = lv.localVectors[lv.currentLocalVector].vector; // current local vector
        size_t const rows = rc.size();
        for (size_t i=0; i<rows; ++i)
        {
          rc[i] += integrationFactor * cache.template d1<RhsBlock::rowId>(at_c<rSpaceIndex>(evaluators).evalData[i]);
          assert(!std::isnan(rc[i][0])); // check that first entry in rhs vector is valid
        }
      }

      Evaluators const&    evaluators;
      Real                 integrationFactor;
      Cache const&         cache;
    };


    /**
     * \brief A boost::fusion functor which adds up contributions to the local stiffness matrices during the integration
     */
    template <class AnsatzVariables, class TestVariables, class Evaluators, class Real, class Cache>
    struct UpdateLocalMatrix
    {
      UpdateLocalMatrix(Evaluators const& evaluators_, Real integrationFactor_, Cache const& cache_):
        evaluators(evaluators_), integrationFactor(integrationFactor_), cache(cache_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        typedef typename LocalMatrices::MatrixBlock MatrixBlock;

        int const rSpaceIndex = result_of::value_at_c<TestVariables,  MatrixBlock::rowId>::type::spaceIndex;
        int const cSpaceIndex = result_of::value_at_c<AnsatzVariables,MatrixBlock::colId>::type::spaceIndex;


        assert(lm.currentLocalMatrix<=lm.lastLocalMatrix);
        auto& A = lm.localMatrices[lm.currentLocalMatrix];  // the current local Galerkin matrix to contribute to

        // retrieve number of localized ansatz functions
        int const rows = A.ridx().size();
        int const cols = A.cidx().size();

        // loop over all combinations of localized ansatz and test functions
        for (int i=0; i<rows; ++i)
        {
          int begin = MatrixBlock::lumped? i: 0;
          int end = (MatrixBlock::symmetric||MatrixBlock::lumped)? i+1: cols;
          for (int j=begin; j<end; ++j)
          {
            A(i,j) += integrationFactor * cache.template d2<MatrixBlock::rowId,MatrixBlock::colId>(
                at_c<rSpaceIndex>(evaluators).evalData[i], at_c<cSpaceIndex>(evaluators).evalData[j]);
            assert(!std::isnan(A(i,j)[0][0])); // check that top left entry is valid
          }
        }
      }

      Evaluators const& evaluators;
      Real              integrationFactor;
      Cache const&      cache;
    };


    template <class AnsatzVariables, class TestVariables, class Evaluators, class Real, class Cache>
    struct UpdateLocalMatrixFromInnerBoundaryCache
    {
      UpdateLocalMatrixFromInnerBoundaryCache(Evaluators const& insideEval, Evaluators const& outsideEval, Real integrationFactor_, Cache const& cache_, size_t offset_=0):
        insideEvaluator(insideEval), outsideEvaluator(outsideEval), integrationFactor(integrationFactor_), cache(cache_), offset(offset_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        typedef typename LocalMatrices::MatrixBlock MatrixBlock;

        int const rSpaceIndex = result_of::value_at_c<TestVariables,  MatrixBlock::rowId>::type::spaceIndex;
        int const cSpaceIndex = result_of::value_at_c<AnsatzVariables,MatrixBlock::colId>::type::spaceIndex;


        assert(lm.currentLocalMatrix+offset<=lm.lastLocalMatrix);
        auto& A = lm.localMatrices[lm.currentLocalMatrix+offset];

        // retrieve number of localized ansatz functions
        int const rows = A.ridx().size();
        int const cols = A.cidx().size();

        // loop over all combinations of localized ansatz and test functions
        for (int i=0; i<rows; ++i)
        {
          int begin = MatrixBlock::lumped? i: 0;
          int end = (MatrixBlock::symmetric||MatrixBlock::lumped)? i+1: cols;
          for (int j=begin; j<end; ++j)
          {
            A(i,j) += integrationFactor * cache.template d2<MatrixBlock::rowId,MatrixBlock::colId>(
                at_c<rSpaceIndex>(insideEvaluator).evalData[i], at_c<cSpaceIndex>(outsideEvaluator).evalData[j], (offset==0));
            assert(!std::isnan(A(i,j)[0][0])); // check that top left entry is valid
          }
        }
      }

      Evaluators const& insideEvaluator, outsideEvaluator;
      Real              integrationFactor;
      Cache const&      cache;
      size_t offset;
    };
    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion that scatters a set of local rhs vectors into the global ones
     */
    template <class GlobalRhs>
    struct ScatterLocalRhs
    {
      ScatterLocalRhs(GlobalRhs& globalRhs_, bool immediate_ = false):
        globalRhs(globalRhs_), immediate(immediate_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        typedef typename LocalVectors::RhsBlock RhsBlock;

        if (immediate || lv.currentLocalVector+1==lv.localVectors.size())
        {
          int const last = immediate ? lv.currentLocalVector : lv.localVectors.size();

          // step through all row groups and scatter all the included local vectors' components
          for (int rg=0; rg<lv.rb->rowGroup().size(); ++rg)
          {
            // obtain the row group
            auto rowRange = lv.rb->rowGroup().lock(rg);

            // step through all local vectors
            for (int i=0; i<last; ++i)
            {
              auto& ri = lv.localVectors[i].vector;
              auto& ridx = lv.localVectors[i].ridx;

              assert(ri.size()==ridx.size());

              for (int r=0; r<ridx.size(); ++r)
              {
                size_t const rGlobalIndex = ridx[r].first;

                // check which rows are inside the row group
                if (rowRange.first<=rGlobalIndex && rGlobalIndex<rowRange.second)
                {
                  size_t const rLocalIndex = ridx[r].second;
                  at_c<RhsBlock::rowId>(globalRhs.data)[rGlobalIndex] += ri[rLocalIndex];
                }
              }
            }

            // release row group
            lv.rb->rowGroup().unlock(rg);
          }

          // everything's scattered - clean up
          lv.currentLocalVector = 0;
        }
        else
          ++lv.currentLocalVector; // prepare next local vector for updating
      }

      GlobalRhs&        globalRhs;
      bool              immediate;
    };

    //---------------------------------------------------------------------

    // This has to be done via partial specialization, as transpose returns
    // an invalid type if the matrix entries themselves are not symmetric.
    template <bool sym, class MatrixBlock>
    struct Symmetrize { static void symmetrize(typename MatrixBlock::LocalMatrix&) {} };

    /**
     * \brief Copies the lower triangular part into the upper triangular part, transposing each element.
     */
    template <class MatrixBlock>
    struct Symmetrize<true,MatrixBlock>
    {
      static void symmetrize(typename MatrixBlock::LocalMatrix& lm)
      {
        int const N = lm.ridx().size();
        assert(N==lm.cidx().size());

        for (int i=0; i<N-1; ++i)
          for (int j=i+1; j<N; ++j)
            lm(i,j) = transpose(lm(j,i));
      }
    };


    /**
     * \brief A boost::fusion functor for symmetrizing local Galerkin matrices.
     *
     * This does symmetrize the local Galerkin matrix in case it should be
     * symmetric (in which case only the lower triangular part is
     * assembled).
     *
     * This functor does only work on local Galerkin
     * matrices and is therefore thread-safe.
     */
    struct SymmetrizeLocalMatrix {

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        typedef typename LocalMatrices::MatrixBlock MatrixBlock;

        // If a block is symmetric, only the lower part of the local
        // matrix is filled. Mirror to the upper part for easy access in
        // the following loop. For lumped (purely diagonal) blocks, of course,
        // nothing is to do.
        if (!MatrixBlock::lumped)
          Symmetrize<MatrixBlock::symmetric,MatrixBlock>::symmetrize(lm.localMatrices[lm.currentLocalMatrix]);
      }
    };


    /**
     * \brief A boost::fusion functor for scattering local stiffness matrices into the global stiffness matrix.
     */
    class ScatterLocalMatrix
    {
    public:
      ScatterLocalMatrix(bool immediate_ = false): immediate(immediate_), nextMatrices(1) {}

      void setNextMatrices(size_t n)
      {
        nextMatrices = n;
      }

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        // Only if all local matrices are used up, or the allowed memory is exhausted (or if it is explicitly desired)
        // scatter them simultaneously into the global matrix.
        if (immediate || lm.lastLocalMatrix+nextMatrices >= lm.localMatrices.size() || lm.storageSize()>=lm.storageSizeLimit())
        {
          lm.mb->globalMatrix().scatter(lm.localMatrices.begin(),lm.localMatrices.begin()+lm.lastLocalMatrix+1);

          // everything's scattered: clean up
          lm.clear();
        }
      }

    private:

      bool immediate;
      size_t nextMatrices;
    };


    struct PrintResult
    {
      template <class Rhs>
      void operator()(Rhs const& r) const { std::cout << "rhs:\n" << r << "\n"; }
    };

    //---------------------------------------------------------------------
    //---------------------------------------------------------------------

    /**
     * \brief provides some meta information for assembly problems (here for variational functionals)
     *
     * In contrast to weak formulations, variational functionals provide a method d0() to evaluate the
     * scalar value of the functional. The matrix is always symmetric (it's a second derivative).
     */
    template <class VariationalFunctional>
    struct VariationalFunctionalPolicy
    {
      typedef typename VariationalFunctional::Scalar                Scalar;
      typedef typename VariationalFunctional::AnsatzVars::Spaces    Spaces;
      static constexpr bool considerInnerFaces = hasInnerBoundaryCache<VariationalFunctional>();

      template <int rowId>
      struct RhsBlockInfo
      {
        static bool const present = VariationalFunctional::template D1<rowId>::present;
      };


      template <int rowId, int colId>
      struct BlockInfo
      {
        static bool const present = VariationalFunctional::template D2<rowId,colId>::present && rowId>=colId;
        static bool const symmetric = (rowId==colId || VariationalFunctional::template D2<rowId,colId>::symmetric)
                      && (result_of::value_at_c<typename VariationalFunctional::TestVars::Variables,rowId>::type::spaceIndex
                          == result_of::value_at_c<typename VariationalFunctional::AnsatzVars::Variables,colId>::type::spaceIndex);
        static bool const mirror = rowId>colId;
        static bool const lumped = symmetric && VariationalFunctional::template D2<rowId,colId>::lumped;
      };

      /**
       * Returns the functional's value that has been assembled in a
       * previous call to assemble() with flags|VALUE==true. If no such
       * call to assemble() occured before, the value is undefined.
       */
      Scalar functional() const { return fValue; }

      //    protected:
    public:
      static constexpr bool hasValue = true;

      Scalar fValue;

      template <class Cache>
      static Scalar d0(Cache const& cache) { return cache.d0(); }
    };

    /**
     * \brief provides some meta information for assembly problems (here for weak formulations)
     *
     * In contrast to variational functional problems, weak formulations are not induced by real-valued
     * functionals. As a consequence, there is no functional value d0() to be evaluated, and the matrix
     * need not be symmetric.
     */
    template <class Functional>
    struct WeakFormulationPolicy
    {
      typedef typename Functional::Scalar                Scalar;
      typedef typename Functional::AnsatzVars::Spaces    Spaces;
      static constexpr bool considerInnerFaces = hasInnerBoundaryCache<Functional>();

      template <int rowId>
      struct RhsBlockInfo
      {
        static bool const present = Functional::template D1<rowId>::present;
      };


      template <int rowId, int colId>
      struct BlockInfo
      {
        static bool const present = Functional::template D2<rowId,colId>::present;
        static bool const symmetric = Functional::template D2<rowId,colId>::symmetric
            && (result_of::value_at_c<typename Functional::TestVars::Variables,rowId>::type::spaceIndex
                == result_of::value_at_c<typename Functional::AnsatzVars::Variables,colId>::type::spaceIndex);
        static bool const mirror = false;
        static bool const lumped = symmetric && Functional::template D2<rowId,colId>::lumped;
      };

      //    protected:
    public:
      static constexpr bool hasValue = false;

      Scalar fValue;

      template <class Cache>
      static Scalar d0(Cache const& cache) { return 0; }
    };


    template <class Problem>
    struct PolicySelector
    {
      typedef typename boost::mpl::if_c<Problem::type==VariationalFunctional,
          VariationalFunctionalPolicy<Problem>,
          WeakFormulationPolicy<Problem> >::type type;
    };


    /**
     * \ingroup assembly
     * \brief A trivial block filter: take all blocks
     */
    struct TakeAllBlocks
    {
      template <int row>          struct D1 { static bool const assemble = true; };
      template <int row, int col> struct D2 { static bool const assemble = true; };
    };

    // Mapping a block filter to an MPL predicate
    template <class BlockFilter>
    struct MatrixBlockFilter
    {
      template <class X>
      struct apply
      {
        typedef boost::mpl::bool_<BlockFilter::template D2<X::rowId,X::colId>::assemble> type;
      };
    };

    template <class> struct Fill;
    template <class> struct FillPointer;
    template <class,int,int> struct FillPointerWithBlock;

    template <class Block, bool, bool>
    struct CopyBlock
    {
      template <class OtherBlock>
      static void apply(OtherBlock const&, std::unique_ptr<Block>&, bool, size_t){}
    };

    template <class Scalar, int n, class Allocator, bool symmetric>
    struct CopyBlock<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Allocator>,symmetric,true>
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Allocator> Block;

      static void apply(Block const& from, std::unique_ptr<Block>& to, bool extractOnlyLowerTriangle, size_t nnz)
      {
        // copy
        if(extractOnlyLowerTriangle == symmetric)
        {
          to.reset(new Block(from));
          return;
        }

        // only lower triangle is stored in the assembler -> restore full matrix
        if( (!extractOnlyLowerTriangle && symmetric) || (extractOnlyLowerTriangle && !symmetric) )
        {
          std::vector<std::vector<size_t> > bcrsPattern(from.N());
          for(size_t i=0; i<from.N(); ++i)
            for(size_t j=0; j<=i; ++j)
            {
              if(from.exists(i,j))
              {
                bcrsPattern[i].push_back(j);
                if(i!=j && !extractOnlyLowerTriangle)
                  bcrsPattern[j].push_back(i);
              }
            }
          for(std::vector<size_t>& row : bcrsPattern) std::sort(row.begin(),row.end());

          to.reset(new Block(from.N(),from.M(),nnz,Block::row_wise));

          // init sparsity pattern
          for (typename Block::CreateIterator row=to->createbegin(); row!=to->createend(); ++row)
            for(size_t col : bcrsPattern[row.index()]) row.insert(col);

          // read data
          for(size_t row=0; row<bcrsPattern.size(); ++row)
            for(size_t col : bcrsPattern[row])
            {
              if(row >= col) (*to)[row][col] = from[row][col];
              else if(!extractOnlyLowerTriangle)
              {
                for(size_t i=0; i<n; ++i)
                  for(size_t j=0; j<n; ++j)
                    (*to)[row][col][i][j] = from[col][row][j][i];
              }
            }
        }
      }
    };


    template <class BCRSMat, int rowIndex, int columnIndex>
    struct BlockToBCRS
    {
      explicit BlockToBCRS(std::unique_ptr<BCRSMat>& m, bool extractOnlyLowerTriangle_, size_t nnz_) : mat(m), extractOnlyLowerTriangle(extractOnlyLowerTriangle_), nnz(nnz_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        CopyBlock<BCRSMat,MatrixBlock::symmetric,MatrixBlock::rowId==rowIndex && MatrixBlock::colId==columnIndex>::apply(mb.globalMatrix(),mat,extractOnlyLowerTriangle,nnz);
      }

    private:
      std::unique_ptr<BCRSMat>& mat;
      bool extractOnlyLowerTriangle;
      size_t nnz;
    };


    template <class Scalar, int n, int m, class Allocator, int row, int column>
    struct FillPointerWithBlock<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator>,row,column>
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> Matrix;

      template <class MatrixBlockArray>
      static std::unique_ptr<Matrix> apply(MatrixBlockArray const& block, bool extractOnlyLowerTriangle, size_t nnz)
        {
          std::unique_ptr<Matrix> result;

          for_each(block,BlockToBCRS<Matrix,row,column>(result,extractOnlyLowerTriangle,nnz));
          return result;
        }
    };
  } // End of namespace AssemblyDetail
  /// \endcond
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector checks the existence
   * of boundary intersections on construction and simply looks it up
   * when questioned.
   *
   * If multiple assembly passes are performed on the same mesh, this
   * tends to be faster than the ForwardingBoundaryDetector or even the
   * TrivialBoundaryDetector.
   */
  template <class GridView>
  class CachingBoundaryDetector
  {
    typedef typename GridView::Grid Grid;
    typedef typename Grid::template Codim<0>::Entity Cell;
    typedef CachingBoundaryDetector<GridView> Self;

  public:
    CachingBoundaryDetector(GridSignals& signals, GridView const& gridView_) :
      gridView(gridView_), indexSet(gridView.indexSet())
    {
      refConnection = signals.informAboutRefinement.connect(int(GridSignals::allocResources),
                                                            [=](GridSignals::Status when){this->update(when);});
      update(GridSignals::AfterRefinement);
    }

    CachingBoundaryDetector(GridView const& gridView_) :
      gridView(gridView_), indexSet(gridView.indexSet())
    {
      update(GridSignals::AfterRefinement);
    }

    bool hasBoundaryIntersections(Cell const& cell) const
    {
      assert(indexSet.size(0)==boundaryFlag.size()); // check for grid adaptation...
      return boundaryFlag[indexSet.index(cell)];
    }

  private:
    GridView const&                    gridView;
    typename GridView::IndexSet const& indexSet;
    std::vector<bool>                  boundaryFlag;
    boost::signals2::scoped_connection  refConnection;

    void update(int const status)
    {
      using namespace Dune;

      if (status == GridSignals::AfterRefinement)
      {
        boundaryFlag.resize(indexSet.size(0));
        // TODO: run this loop in parallel (no so easy because of concurrent access to the bool vector...)
        for (auto const& cell: elements(gridView))
          boundaryFlag[indexSet.index(cell)] = cell.hasBoundaryIntersections();
      }
    }
  };

  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector simply forwards the
   * questioning to the entity.
   */
  class ForwardingBoundaryDetector
  {
  public:
    template <class GridView>
    ForwardingBoundaryDetector(GridSignals const&, GridView const&)
    {}

    template <class GridView>
    ForwardingBoundaryDetector(GridView const&)
    {}

    template <class Entity>
    bool hasBoundaryIntersections(Entity const& cell) const
    {
      return cell.hasBoundaryIntersections();
    }
  };

  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector simply says every
   * cell might have boundary intersections.
   */
  class TrivialBoundaryDetector
  {
  public:
    template <class GridView>
    TrivialBoundaryDetector(GridSignals const&, GridView const&)
    {}

    template <class GridView>
    TrivialBoundaryDetector(GridView const&)
    {}

    template <class Entity>
    bool hasBoundaryIntersections(Entity const& cell) const
    {
      return true;
    }
  };


  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  extern bool doDebug;

  namespace AssemblyDetail
  {
    /**
     * \brief Policy that controls the optional assembly over inner intersections.
     */
    template <class Functional, bool considerInnerFaces> class InnerBoundaryAssemblerImpl;

    /**
     * \brief Specialization for the case that inner intersections will be ignored.
     */
    template <class Functional>
    class InnerBoundaryAssemblerImpl<Functional,false>
    {
      typedef typename Functional::Scalar Scalar;
      typedef typename Functional::AnsatzVars AnsatzVariables;
      typedef typename Functional::TestVars TestVariables;
      typedef typename AnsatzVariables::Spaces Spaces;
      typedef typename AnsatzVariables::Grid Grid;
      typedef typename Grid::LeafGridView GridView;
      typedef typename GridView::template Codim<0>::Iterator CellIterator;
      typedef typename CellIterator::Entity Cell;
      typedef typename GridView::IntersectionIterator FaceIterator;
      typedef ShapeFunctionCache<Grid,Scalar> SfCache;

    public:
      InnerBoundaryAssemblerImpl(Spaces const&, GridView const&){}

      void createInnerBoundaryCache(Functional const&, unsigned int){}

      template <class LocalRhs, class LocalMatrices>
      void assembleInnerBoundaryCache(FaceIterator const&, int const, Evaluators<Spaces> const&, bool,
                                      Scalar&, LocalRhs&, LocalMatrices&, unsigned int, size_t&){}

      template <class LocalRhs, class LocalMatrices, class Evaluators>
      void initLocalStorageForInnerBoundaryCache(LocalRhs&, LocalMatrices&, Cell const&, Evaluators const&, unsigned int){}

      size_t nNewMatrices(Cell const&) const { return 1; }
    };


    /**
     * \brief Specialization for the case that inner intersections must be considered.
     * i.e. Discontinuous Galerkin Methods, ...
     */
    template <class Functional>
    class InnerBoundaryAssemblerImpl<Functional,true>
    {
      typedef typename Functional::Scalar Scalar;
      typedef typename Functional::AnsatzVars::Variables AnsatzVariables;
      typedef typename Functional::TestVars::Variables TestVariables;
      typedef typename Functional::AnsatzVars::Spaces Spaces;
      typedef typename Functional::AnsatzVars::Grid Grid;
      typedef typename Grid::LeafGridView GridView;
      typedef typename GridView::template Codim<0>::Iterator CellIterator;
      typedef typename CellIterator::Entity Cell;
      typedef typename Grid::ctype CoordType;
      typedef typename GridView::IntersectionIterator FaceIterator;
      typedef typename Functional::InnerBoundaryCache InnerBoundaryCache;
      typedef Dune::QuadratureRule<CoordType,Grid::dimension-1> QuadRule;

      enum { VALUE=0x1, RHS=0x2, MATRIX=0x4 };
      static constexpr int dim = Grid::dimension;

    public:
      InnerBoundaryAssemblerImpl(Spaces const& spaces, GridView const& gridView_)
      : sfCache(), neighbourEvaluators(getEvaluators(spaces,&sfCache)),
        gridView(gridView_), innerBoundaryCache(nullptr)
      {}

      void createInnerBoundaryCache(Functional const& f, unsigned int flags)
      {
        innerBoundaryCache.reset(new InnerBoundaryCache(f.createInnerBoundaryCache(flags)));
      }

      template <class LocalRhs, class LocalMatrices>
      void assembleInnerBoundaryCache(FaceIterator const& face, int const p, Evaluators<Spaces>& evaluators, bool hasValue,
                                      Scalar& floc, LocalRhs& localRhs, LocalMatrices& localMatrices, unsigned int flags, size_t& offset)
      {
        // only consider faces that do not lie on the boundary and have a neighbour
        // the latter may not be of interest
        if( face->neighbor() && !face->boundary() )
        {
          // update offset to get access to the correct local matrix
          ++offset;
          innerBoundaryCache->moveTo(face);

          // move evaluator for the neighbouring cell
          moveEvaluatorsToCell(neighbourEvaluators,*(face->outside()));

          // init evaluator on inside cell
          Dune::GeometryType gt = face->geometryInInside().type();
          assert(gt == face->geometryInOutside().type());

          QuadratureTraits<QuadRule>  faceQuadratureCache;
          QuadRule const qr = faceQuadratureCache.rule(gt,p);
          useQuadratureRuleInEvaluators(evaluators,qr,face->indexInInside());
          useQuadratureRuleInEvaluators(neighbourEvaluators,qr,face->indexInOutside());
          size_t nQuadPos = qr.size();

          for(size_t g=0; g<nQuadPos; ++g)
          {
            // position of integration point
            const Dune::FieldVector<CoordType,dim-1>& quadPos = qr[g].position();
            double weightTimesDetJac = qr[g].weight() * face->geometry().integrationElement(quadPos);

            // Move evaluators to integration point
            // TODO: enable caching, currently caching is not in use as we have to match the local ids of the quadrature points in the
            //       local coordinates of the neighboring codim-0-entities
            moveEvaluatorsToIntegrationPoint(evaluators,face->geometryInInside().global(quadPos));
            moveEvaluatorsToIntegrationPoint(neighbourEvaluators,face->geometryInOutside().global(quadPos));

            innerBoundaryCache->evaluateAt(quadPos,evaluators,neighbourEvaluators);

            // update value
            if (hasValue && (flags & VALUE))
	      floc += weightTimesDetJac*PolicySelector<Functional>::type::d0(*innerBoundaryCache);
            // step through all rhs blocks and update their local rhs
	    if (flags & RHS)
	      boost::fusion::for_each(localRhs,UpdateLocalRhs<TestVariables,Evaluators<Spaces>,Scalar,
				      InnerBoundaryCache>(evaluators,weightTimesDetJac,*innerBoundaryCache));
            // update local matrices
            if (flags & MATRIX)
            {
              // update diagonal block, contribution of coupling terms on cell
              boost::fusion::for_each(localMatrices,UpdateLocalMatrixFromInnerBoundaryCache<AnsatzVariables,TestVariables,Evaluators<Spaces>,Scalar,
				      InnerBoundaryCache>(evaluators,evaluators,weightTimesDetJac,*innerBoundaryCache,0));
              // update off-diagonal block, contribution of coupling terms with neighbouring cell
              boost::fusion::for_each(localMatrices,UpdateLocalMatrixFromInnerBoundaryCache<AnsatzVariables,TestVariables,Evaluators<Spaces>,Scalar,
				      InnerBoundaryCache>(evaluators,neighbourEvaluators,weightTimesDetJac,*innerBoundaryCache,offset));
            }
          }
        }
      }

      template <class LocalRhs, class LocalMatrices, class Evaluators>
      void initLocalStorageForInnerBoundaryCache(LocalRhs& localRhs, LocalMatrices& localMatrices, Cell const& cell, Evaluators const& evaluators, unsigned int flags)
      {
        if (flags & MATRIX)
        {
          FaceIterator fend = gridView.iend(cell);
          for(FaceIterator face = gridView.ibegin(cell); face!=fend; ++face)
            if(face->neighbor() && !face->boundary())
            {
              moveEvaluatorsToCell(neighbourEvaluators,*(face->outside()) );
              for_each(localMatrices,ClearLocalMatrix<Evaluators,AnsatzVariables,TestVariables>(evaluators,neighbourEvaluators));
            }
        }
      }

      /**
       * \brief Computes the number of local matrices affected by assembly on the given cell.
       *
       * If inner boundary terms are to be assembled, degrees of freedom living on this cell
       * are connected to those living on the neighbouring cells. Thus in one step we assemble
       * multiple local matrices.
       */
      size_t nNewMatrices(Cell const& cell)
      {
        size_t nextMatrices = 1;
        auto fend = gridView.iend(cell);
        for (auto face = gridView.ibegin(cell); face!=fend; ++face)
          if(face->neighbor() && !face->boundary())
	    ++nextMatrices;
        return nextMatrices;
      }

    private:
      ShapeFunctionCache<Grid,Scalar> sfCache;
      Evaluators<Spaces> neighbourEvaluators;
      GridView const& gridView;
      std::unique_ptr<InnerBoundaryCache> innerBoundaryCache;
    };

    /**
     * \brief Policy that checks for the existence of 'InnerBoundaryCache'.
     *
     * If InnerBoundaryCache is detected, assembly also takes into account
     * inner intersections.
     */
    template <class Functional> using InnerBoundaryAssembler = InnerBoundaryAssemblerImpl<Functional,hasInnerBoundaryCache<Functional>()>;

  } // end of namespace AssemblyDetail


  namespace Assembler
  {
    /**
     * \ingroup assembly
     * \brief What to assemble.
     *
     * Binary flags to be or'ed for supplying assembly requests.
     */
    enum { VALUE=0x1, RHS=0x2, MATRIX=0x4, EVERYTHING=0x7 };
  }


  /**
   * \ingroup assembly
   * \brief A class for assembling Galerkin representation matrices and right hand sides for variational functionals depending on multiple variables.
   *
   * Template parameters:
   * \tparam F: The variational functional. This has to be a model of the VariationalFunctional concept.
   * \tparam BoundaryDetector: A policy class type for detecting cells
   *   incident to the boundary. Currently CachingBoundaryDetector,
   *   ForwardingBoundaryDetector, and TrivialBoundaryDetector may be used.
   * \tparam QuadRule: allows choice of quadrature formula. Default: provided by DUNE, some special formulas: special_quadrature.hh
   *
   * After a grid refinement an assembler is not valid anymore and has to be
   * constructed again. This is, because the assembler copies the index set
   * at construction, which changes, if the grid changes
   *
   */
  template <class F,
            class BoundaryDetector = CachingBoundaryDetector<typename F::AnsatzVars::GridView>,
            class QuadRule =  Dune::QuadratureRule<typename F::AnsatzVars::Grid::ctype, F::AnsatzVars::Grid::dimension> >
  class VariationalFunctionalAssembler : public AssemblyDetail::PolicySelector<F>::type, public AssemblyDetail::InnerBoundaryAssembler<F>
  {
    typedef typename AssemblyDetail::PolicySelector<F>::type            Policy;
    typedef AssemblyDetail::InnerBoundaryAssembler<F>                   InnerBoundary;

  public:
    typedef VariationalFunctionalAssembler<F,BoundaryDetector,QuadRule> Self;

    /// functional
    typedef          F                            Functional;
    /// ansatz variables DEPRECATED - use AnsatzVariableSetDescription instead. Will be removed after 2014-05-01.
    typedef typename Functional::AnsatzVars       AnsatzVariableSet  __attribute((deprecated));

    typedef typename Functional::AnsatzVars       AnsatzVariableSetDescription;

    /// test variables DEPRECATED - use TestVariableSetDescription instead. Will be removed after 2014-05-01.
    typedef typename Functional::TestVars         TestVariableSet __attribute((deprecated));

    typedef typename Functional::TestVars         TestVariableSetDescription;

    static_assert(std::is_same<typename AnsatzVariableSetDescription::Spaces, typename TestVariableSetDescription::Spaces>::value,
		  "VariableSetDescriptions for ansatz spaces and test spaces must provide the same space list.");
    /// grid
    typedef typename AnsatzVariableSetDescription::Grid      Grid;

    /**
     * \brief The grid view on which the variables live.
     */
    typedef typename AnsatzVariableSetDescription::GridView  GridView;

    /// spaces
    typedef typename AnsatzVariableSetDescription::Spaces    Spaces;
    /// index set
    typedef typename AnsatzVariableSetDescription::IndexSet  IndexSet;

    /// deprecated, will be removed after 2013-12-31
    typedef typename Functional::Scalar           Scalar;

    /**
     * \brief Underlying field type.
     */
    typedef typename Functional::Scalar           field_type;

  private:
    typedef typename AnsatzVariableSetDescription::Variables AnsatzVariables;
    typedef typename TestVariableSetDescription::Variables   TestVariables;

    // TODO: document me
    static bool const consistentSymmetry = Functional::type==WeakFormulation
                                           || std::is_same<AnsatzVariableSetDescription,TestVariableSetDescription>::value;

  public:

    /**
     * \brief A boost::fusion sequence of AssemblyDetail::MatrixBlock elements for present matrix blocks.
     *
     * The elements of this sequence type include both block info and the global matrices.
     */
    typedef typename AssemblyDetail::BlockArray<Policy,AnsatzVariables,TestVariables>::type MatrixBlockArray;

    /**
     * \brief A boost::fusion sequence of AssemblyDetail::RhsBock elements for present rhs blocks.
     *
     * The elements of this sequence type include the rhs block infos.
     */
    typedef typename AssemblyDetail::RhsBlockArray<Policy,TestVariables>::type RhsBlockArray;

    /**
     * \brief A LinearProductSpace type of right hand side coefficient vectors.
     */
    typedef typename TestVariableSetDescription::template CoefficientVectorRepresentation<>::type RhsArray;

    /// Construct an empty assembler
    VariationalFunctionalAssembler(GridManagerBase<Grid>& gridManager_, Spaces const& spaces)
    : InnerBoundary(spaces, boost::fusion::at_c<0>(spaces)->gridView()),
      spaces_(spaces),
      gridManager(gridManager_),
      gridView(boost::fusion::at_c<0>(spaces)->gridView()),
      indexSet(gridView.indexSet()),
      boundaryDetector(gridManager_.signals,gridView),
      nSimultaneousBlocks(80),
      rowBlockFactor(2.0),
      localStorageSize(128000) // roughly 128 kB
    {
      assert(consistentSymmetry);
      refConnection = gridManager_.signals.informAboutRefinement.connect(int(GridSignals::freeResources),
                                                                         [this](GridSignals::Status when){this->reactToRefinement(when);});
    }

    /// Construct an empty assembler, gridManager is implicitly passed via spaces.
    explicit VariationalFunctionalAssembler(Spaces const& spaces)
    : VariationalFunctionalAssembler(boost::fusion::deref(boost::fusion::begin(spaces))->gridManager(),spaces)
    {
    }

    // injected here for backward compatibility and convenience. TODO: Remove in the long run?
    static unsigned int const VALUE  = Assembler::VALUE;
    static unsigned int const RHS    = Assembler::RHS;
    static unsigned int const MATRIX = Assembler::MATRIX;
    static unsigned int const EVERYTHING = Assembler::EVERYTHING;


    /**
     * \brief Defines how many cells are assembled locally before scattering them together into the global data structures.
     *
     * Higher numbers reduce the overhead of locking for mutually exclusive write access to global matrix/rhs and improve
     * memory access locality, but increase the amount of thread-local computation and storage. The default value is 42, but the
     * performance appears to be not very sensitive to this parameter in the range [20,60].
     */
    Self& setNSimultaneousBlocks(int n)
    {
      nSimultaneousBlocks = n;
      return *this;
    }

    /**
     * \brief Defines how many more row blocks in each matrix are used compared to the number of threads.
     *
     * When scattering the local matrices/rhs into the global data structures, blocks of sequential rows are
     * processed in turn. Each block of rows is equipped with a lock in order to guarantee exclusive write access.
     *
     * A high number of row blocks allows fine granularity and improves cache locality and simultaneous scattering into different
     * parts of the global data structures, but incurs some inefficiency due to more frequent locking and more
     * thread-local processing. The default value is to use two times the number of threads.
     */
    Self& setRowBlockFactor(double a)
    {
      rowBlockFactor = a;
      return *this;
    }

    /**
     * \brief Defines how many memory the locally assembled matrices may occupy before they are scattered
     *
     * The amount of memory (in bytes) should scale with CPU cache, roughly a quarter of the second level cache size should be
     * a reasonable value. Note that this is merely a hint, not a hard limit. The actual local matrix storage
     * can be slightly larger.
     */
    Self& setLocalStorageSize(size_t s)
    {
      localStorageSize = s;
      return *this;
    }

    /// Assembly without block filter
    void assemble(F const& f, unsigned int flags=Assembler::EVERYTHING, int nThreads=0, bool verbose=false)
    {
      if(flags)
        assemble<AssemblyDetail::TakeAllBlocks>(f,flags,nThreads,verbose);
    }

  private:
    // Define iterator and entity type for stepping through all cells of the grid.
    typedef typename GridView::template Codim<0>::Iterator CellIterator;
    typedef typename CellIterator::Entity Entity;


  public:

    /**
     * \brief Create data in assembler
     *
     * Assembles the value and/or derivative and/or Hessian (according
     * to the given flags, default is to assemble all three) of the
     * FE-discretized variational functional given by f. The block
     * filter can be used for partial assembly of specified blocks in
     * the matrix and rhs.
     *
     * The data that is not concerned by the given flags and block
     * filter is left unmodified. Be careful: This allows to drive the
     * Galerkin operator representation into a semantically inconsistent
     * state. This need not be a problem (instead it can be useful and more efficient in
     * certain situations), but one has to be aware of that fact.
     *
     * The assembled data can afterwards be accessed via the functional,
     * rhs, and matrix methods.
     *
     * \param f the variational functional or weak formulation to be assembled
     * \param flags a bitset determining what to assemble
     * \param nTasks number of tasks to submit to the thread pool for parallel execution. 0 (default) means to use the
     *                 default hardware concurrency reported by the system. Sequential computation is performed if the grid is reported
     *                 not to be thread-safe. You can call \ref GridManagerBase::enforceConcurrentReads
     *                 to lie about thread-safety in case you know (or believe) the implementation to be
     *                 thread safe.
     * \param verbose whether to output status messages to standard output
     *
     * \see NumaThreadPool
     */
    template <class BlockFilter>
    void assemble(F const& f, unsigned int flags=Assembler::EVERYTHING, int nTasks=0, bool verbose=false)
    {
      using namespace Dune;
      using namespace boost::fusion;
      using namespace AssemblyDetail;

      boost::timer::cpu_timer globalTimer;

      auto const& cellRanges = gridManager.cellRanges(gridView);
      if (statistics & timerStatistics)
        std::cout << "computed cell ranges: " << globalTimer.format();
      globalTimer.start();

#ifndef KASKADE_SEQUENTIAL
      // Choose the number of tasks
      if (nTasks < 1)
        nTasks = NumaThreadPool::instance().cpus();

      if (!gridManager.gridIsThreadSafe())
      {
        if (verbose) std::cout << "Grid is not thread safe. Reducing number of used threads to 1." << std::endl;
        nTasks = 1;
      }
      else if(verbose) std::cout << "#tasks:" << nTasks << " " << std::endl;
      assert(nTasks>=1);

      nTasks = std::min(nTasks,cellRanges.maxRanges());
#else
      nTasks = 1;
#endif
      // set entries to zero
      auto clearGlobalData = [](auto& x) { x = 0; };

      if (flags & RHS)
        for_each(getRhs().first.data,clearGlobalData);

      if (flags & MATRIX)
      {
        typedef filter_view<MatrixBlockArray,MatrixBlockFilter<BlockFilter> > MatrixBlockFilter;
        for_each(MatrixBlockFilter(getMatrix()),clearGlobalData);
      }

      if (statistics & timerStatistics)
        std::cerr << "Data structure creation: " << globalTimer.format();
      globalTimer.start();

      // If there is only one task, things are somewhat simpler and more efficient...
      if (nTasks == 1) {
        Scalar retval =
            assembleCellRange<BlockFilter>(f,flags,gridView.template begin<0>(),gridView.template end<0>(),0);
        if (flags & VALUE)
          this->fValue = retval;
      }
      else
      {
        // Partition the cells in ranges of almost equal size and do the assembly in parallel.
        std::vector<Scalar> values(nTasks,0.0);

        parallelFor([&cellRanges,this,&f,flags,&values](int i, int n) {
                      auto range = cellRanges.range(n,i);
                      values[i] = this->assembleCellRange<BlockFilter>(f,flags,range.begin(),range.end(),i);
                    },
                    nTasks);

        // At the end, add all the functional values up
        if (flags & VALUE)
          this->fValue = std::accumulate(values.begin(),values.end(),0.0);
      }

      // announce validity of assembled data
      validparts |= flags;

      if (statistics & timerStatistics)
        std::cerr << "assembly time: " << globalTimer.format();
    }

  private:


    /**
     * \brief Assembles stiffness matrix and/or right hand side on a given part of the cells.
     *
     * \param[in] f the variational functional (or weak formulation) to be used
     * \param[in] flags what to assemble
     * \param[in] first start of cell range
     * \param[in] last one behind end of cell range
     * \param[in] threadNo thread identifier (for debugging/timing output)
     *
     * This method is intended to be thread-safe.
     */
    template <class BlockFilter>
    Scalar assembleCellRange(F const& f, unsigned int flags, CellIterator cell, CellIterator last, int threadNo)
    {
      using namespace Dune;
      using namespace boost::fusion;
      using namespace AssemblyDetail;

      // timer for gathering performance statistics
      boost::timer::cpu_timer scatterTimer, evalTimer, localTimer, setupTimer, totalTimer;
      if (statistics & localTimerStatistics)
      {
        evalTimer.stop();
        scatterTimer.stop();
        localTimer.stop();
      }
      // Obtain global data structures. Note that here is no race condition even if the
      // data is constructed on demand, as the required data structures have been created
      // and initialized before.
      MatrixBlockArray* globalMatrix    = flags & MATRIX? & getMatrix()    : nullptr;
      RhsArray*         globalRhs       = flags & RHS   ? & getRhs().first : nullptr;
      RhsBlockArray*    globalRhsBlocks = flags & RHS   ? & getRhs().second: nullptr;

      typedef filter_view<MatrixBlockArray,AssemblyDetail::MatrixBlockFilter<BlockFilter> > MatrixBlockFilter;

      // Local contribution to global value.
      double fValue = 0;

      int const dim = Grid::dimension;
      typedef typename Grid::ctype CoordType;

      // Allocate local rhs vectors that can be reused on each
      // element. Since all variables may have different element types
      // (i.e. different number of components), this has to be one local
      // rhs per variable.
      auto localRhs = boost::fusion::as_vector(boost::fusion::transform(*globalRhsBlocks,RhsLocalData<Policy>(nSimultaneousBlocks))); // assemble several cells en bloc before scattering

      auto localMatrices = boost::fusion::as_vector(boost::fusion::transform(MatrixBlockFilter(*globalMatrix),MatrixLocalData(nSimultaneousBlocks,localStorageSize))); // assemble several cells before scattering

      // element-local functional value
      Scalar floc = 0;

      // Shape function cache. Remember that every thread has to use its own cache, thus we create our own here.
      typedef ShapeFunctionCache<Grid,Scalar> SfCache;
      SfCache sfCache;

      // Quadrature rule caches. Remember that every thread has to use its own quadrature caches, thus we create our own here.
      QuadratureTraits<QuadRule>                          quadratureCache;
      QuadratureTraits<QuadratureRule<CoordType,dim-1> >  faceQuadratureCache;

      // Evaluators for shape functions of all FE spaces. Remember that
      // every thread has to use its own collection of evaluators!
      auto evaluators = getEvaluators(spaces(),&sfCache,AssemblyDetail::derivatives(f,spaces(),localMatrices,localRhs));
      using Evaluators = decltype(evaluators);

      auto domainCache = f.createDomainCache(flags);
      auto boundaryCache = f.createBoundaryCache(flags);
      this->createInnerBoundaryCache(f,flags);

      // The symmetrizer functor.
      SymmetrizeLocalMatrix symmetrizeLocalMatrix;

      // The scatter functor. This takes care of locking the global data structures properly.
      ScatterLocalMatrix scatterLocalMatrix;

      // Dummy variable sets to be fed to boost::fusion algorithms.
      typename AnsatzVariableSetDescription::Variables ansatzVars;
      typename TestVariableSetDescription::Variables testVars;

      if (statistics & localTimerStatistics) setupTimer.stop();

      // Assemble representation: step through all elements, assemble
      // the local matrix and rhs, and scatter their entries into the
      // global data structures.
      if (statistics & localTimerStatistics) evalTimer.resume();

      while(cell!=last)
      {
        // tell the problem on which cell we are
        domainCache.moveTo(*cell);

        // for all spaces involved, compute the shape functions and
        // their global indices, which are needed for evaluating the
        // functional's derivative.
        moveEvaluatorsToCell(evaluators,*cell,indexSet.index(*cell));
        int const shapeFunctionMaxOrder = maxOrder(evaluators);

        if (statistics & localTimerStatistics) localTimer.resume();
        if (statistics & localTimerStatistics) evalTimer.stop();

        // clear the local matrices and right hand sides.
        floc = 0;
        if (flags & RHS)
          for_each(localRhs,ClearLocalRhs<Evaluators,TestVariables>(evaluators));

        if (flags & MATRIX)
        {
          for_each(localMatrices,SetCurrentLocalMatrix());
          for_each(localMatrices,ClearLocalMatrix<Evaluators,AnsatzVariables,TestVariables>(evaluators,evaluators));
        }
        this->initLocalStorageForInnerBoundaryCache(localRhs,localMatrices,*cell,evaluators,flags);

        // loop over all quadrature points and integrate local stiffness matrix and rhs
        int const integrationOrderOnCell = f.integrationOrder(*cell,shapeFunctionMaxOrder,false);
        GeometryType gt = cell->type();

        QuadRule const qr = quadratureCache.rule(gt,integrationOrderOnCell);
        useQuadratureRuleInEvaluators(evaluators,qr,0);

        size_t nQuadPos = qr.size();
        for (size_t g=0; g<nQuadPos; ++g)
        {
          // pos of integration point
          FieldVector<CoordType,dim> const& quadPos = qr[g].position();

          // for all spaces involved, update the evaluators associated
          // to this quadrature point
          moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qr,g,0);

          // prepare evaluation of functional
          domainCache.evaluateAt(quadPos,evaluators);

          CoordType weightTimesDetJac(cell->geometry().integrationElement(quadPos)); // determinant of jacobian
          assert(weightTimesDetJac > 0);
          // weight of quadrature point
          weightTimesDetJac *= qr[g].weight();

          // step through all local matrix and rhs blocks and update their data
          if (this->hasValue && (flags & VALUE))
            floc += weightTimesDetJac*Policy::d0(domainCache);
          if (flags & RHS)
            for_each(localRhs,UpdateLocalRhs<TestVariables,Evaluators,Scalar,typename F::DomainCache>(evaluators,weightTimesDetJac,domainCache));
          if (flags & MATRIX)
            for_each(localMatrices,UpdateLocalMatrix<AnsatzVariables,TestVariables,Evaluators,Scalar,typename F::DomainCache>(evaluators,weightTimesDetJac,domainCache));
        }


        // Handle Robin (Cauchy) type boundary conditions and interior face terms. Loop over all
        // faces (codim 1) that are exterior boundaries or are relevant due to interior face terms and integrate
        // the contribution. We assume that intersections on the domain boundary always consist of the whole face (codim 1
        // subentity).
        if(boundaryDetector.hasBoundaryIntersections(*cell) || hasInnerBoundaryCache<Functional>())
        {
          auto faceEnd = gridView.iend(*cell);
          size_t localEdgeNum = 0;
          for (auto face=gridView.ibegin(*cell); face!=faceEnd; ++face)
          {
            // assemble on inner faces if InnerBoundaryCache exists
            int const integrationOrderOnFace = f.integrationOrder(*cell,shapeFunctionMaxOrder,true);
            this->assembleInnerBoundaryCache(face, integrationOrderOnFace, evaluators, this->hasValue, floc, localRhs, localMatrices, flags, localEdgeNum);

            // assemble on domain boundary
            if (face->boundary())
            {
              boundaryCache.moveTo(face);
              GeometryType gt = face->geometryInInside().type();
              QuadratureRule<CoordType,dim-1> const qr = faceQuadratureCache.rule(gt,integrationOrderOnFace);
              useQuadratureRuleInEvaluators(evaluators,qr,face->indexInInside());

              size_t nQuadPos = qr.size();
              for (size_t g=0; g<nQuadPos; ++g)
              {
                // pos of integration point
                const FieldVector<CoordType,dim-1>& quadPos = qr[g].position();

                // weight of quadrature point
                double weight = qr[g].weight();
                CoordType const detjac = face->geometry().integrationElement(quadPos); // determinant of jacobian

                // evaluate values at integration points for all spaces
                // involved, update the evaluators associated to this
                // quadrature point
                //
                // TODO: check whether this call is
                // thread-safe. (currently as of 2010-02-11 this is not
                // guaranteed by the Dune interface for UG). Copying the
                // quadrature position instead of obtaining a const
                // reference minimizes the time window for concurrent
                // access.
                FieldVector<CoordType,dim> const quadPosInCell = face->geometryInInside().global(quadPos);

                // @warning: make SURE the following does always hold: (i)
                // On boundary faces, there is exactly one intersection, and
                // this covers the whole face. (ii) When an intersection
                // covers the whole codim 1 subentity, its geometry in the
                // cell is the same as specified in the reference
                // element. Currently this is not guaranteed by the Dune
                // interface. Maybe there is a relation to the ominous
                // TwistUtility from Freiburg?
#if !defined(NDEBUG) && defined(KASKADE_SEQUENTIAL)
                // It seems as if this access to reference elements is
                // not thread-safe. Probably this is due to singleton
                // object construction in reference elements.
                //               Dune::GenericReferenceElement<CoordType,dim> const& refElem =
                //                 Dune::GenericReferenceElements<CoordType,dim>::general(ci->type());
                //               FieldVector<CoordType,dim> diff = refElem.template global<1>(quadPos,face->indexInInside(),1) - quadPosInCell;

                // @warning refElem.template global<1>(quadPos,face->indexInInside(),1) gives the global coordinate _in the reference element_, while quadPosInCell is the coordinate in the actual cell! With Dune 2.0 they differ for some edges, so this test needs to be skipped. Warning (ii) above is violated, but there seems to be no problem.

                //               if(diff.two_norm() > 1e-12)
                //               {
                //                 std::cout << "Geometry Error: " << diff.two_norm()
                //                           << " Quadpos: " << quadPos << std::endl
                //                           << " SubEntityID: " << face->indexInInside()
                //                           << " at QPoint: " << quadPosInCell
                //                           << " which should be identical to: "
                //        << refElem.template global<1>(quadPos,face->indexInInside(),1)
                //        << std::endl;
                //                 assert(0);
                //               }
#endif

                // quadrature points on faces are not consistently indexed (if even symmetric...) - hence evaluating
                // shape functions based on quadrature points is impossible
                // moveEvaluatorsToIntegrationPoint(evaluators,quadPosInCell,qr,g,face->indexInInside());
                moveEvaluatorsToIntegrationPoint(evaluators,quadPosInCell);

                boundaryCache.evaluateAt(quadPos,evaluators);

                // step through all rhs blocks and update their local rhs
                if (this->hasValue && (flags & VALUE))
                  floc += weight*detjac*Policy::d0(boundaryCache);
                if (flags & RHS)
                  for_each(localRhs,UpdateLocalRhs<TestVariables,Evaluators,Scalar,
                                                   typename F::BoundaryCache>(evaluators,weight*detjac,boundaryCache));
                if (flags & MATRIX)
                  for_each(localMatrices,UpdateLocalMatrix<AnsatzVariables,TestVariables,
                                                           Evaluators,Scalar,typename F::BoundaryCache>(evaluators,weight*detjac,boundaryCache));
              } // done integration point loop
            } // done boundary faces
          } // done face loop
        } // done boundary conditions

        // symmetrize local matrix
        if (flags & MATRIX)
          for_each(localMatrices,symmetrizeLocalMatrix);

        if (statistics & localTimerStatistics) scatterTimer.resume();
        if (statistics & localTimerStatistics) localTimer.stop();

        // Move to next cell. This happens before scattering as we need to know the number of local matrices that
        // will be created for the next cell in order to guarantee that the allocated storage of the buffer can
        // hold all local matrices of the next cell.
        ++cell;
        // Scatter local data into global data (note that this occurs blockwise, i.e.
        // the actual scatter may be delayed)
        fValue += floc;
        if (flags & RHS)
          for_each(localRhs,ScatterLocalRhs<RhsArray>(*globalRhs));
        if (flags & MATRIX)
        {
          // Compute the number of local matrices that will be created for the next cell.
          scatterLocalMatrix.setNextMatrices( (cell!=last) ? this->nNewMatrices(*cell) : 0 );
          for_each(localMatrices,scatterLocalMatrix);
        }

        if (statistics & localTimerStatistics) evalTimer.resume(); // count cell iterator increment as well
        if (statistics & localTimerStatistics) scatterTimer.stop();
      } // end iteration over cells

      if (statistics & localTimerStatistics) scatterTimer.resume();
      if (statistics & localTimerStatistics) evalTimer.stop();

      // scattering of local matrices and vectors may not be complete (there may be local data left over in the buffers).
      // flush them.
      if (flags & RHS)
        for_each(localRhs,ScatterLocalRhs<RhsArray>(*globalRhs,true));

      if (flags & MATRIX)
        for_each(localMatrices,ScatterLocalMatrix(true));
      if (statistics & localTimerStatistics) scatterTimer.stop();

      // report statistics
      if (statistics & localTimerStatistics)
      {
        std::cerr << "--------------\nThread " << threadNo << ":\n";
        double totalTime = setupTimer.elapsed().wall + evalTimer.elapsed().wall + scatterTimer.elapsed().wall + localTimer.elapsed().wall;
        std::cerr.precision(1);
        std::cerr.setf(std::ios_base::fixed,std::ios_base::floatfield);
        std::cerr << "Total time:    " << 100* totalTimer.elapsed().wall / totalTime << "% -- " << totalTimer.format()
        << "Setup time:    " << 100* setupTimer.elapsed().wall / totalTime << "% -- " << setupTimer.format()
        << "Eval times:    " << 100* evalTimer.elapsed().wall / totalTime << "% -- " << evalTimer.format()
        << "Scatter times: " << 100* scatterTimer.elapsed().wall / totalTime << "% -- " << scatterTimer.format()
        << "Local times:   " << 100 * localTimer.elapsed().wall / totalTime << "% -- " << localTimer.format();
        std::cerr << "--------------\n";
      }

      return fValue;
    }

  public:
    /**
     * \brief Returns a bitfield specifyign which of the parts have been assembled since
     * construction or flush() according to the format (VALUE|RHS|MATRX).
     */
    int valid() { return validparts; }

    /**
     * \brief Destructs parts of the assembled quatities according to the format (VALUE|RHS|MATRX)
     *
     * The default flag leads to destruction of all three values.
     */
    void flush(int flags=(VALUE|RHS|MATRIX))
    {
      validparts &= !flags;
      if (flags & VALUE)  this->fValue = 0;
      if (flags & RHS)
      {
        rhss.reset();
        rhsBlocks.reset();
      }
      if (flags & MATRIX) matrixBlocks.reset();
    }

    /**
     * \brief the size of a matrix block (in terms of scalar rows/columns)
     * \return a pair (nrows,ncols)
     */
    std::pair<size_t,size_t> size(int row, int col) const
    {
      return std::make_pair(TestVariableSetDescription::degreesOfFreedom(spaces(),row,row+1),
                            AnsatzVariableSetDescription::degreesOfFreedom(spaces(),col,col+1));
    }


    /**
     * \brief Writes the submatrix defined by the half-open block ranges
     * [rbegin,rend), [cbegin,cend) to a newly created std::unique_ptr<MatrixType>.
     *
     * Use this method only if your matrix type does not support move-semantics.
     *
     * If extractOnlyLowerTriangle is true, the Galerkin operator is treated as
     * if its upper triangle was structurally zero.
     *
     * In order to be able to read the assembled matrix into the smart pointer
     * a specialization of AssemblyDetail::FillPointer, i.e. AssemblyDetail::FillPointer<MatrixType>
     * must be available!!!
     *
     * Kaskade provides spezializations of AssemblyDetail::Fill for
     *  - MatrixType = MatrixAsTriplet<Scalar>
     *  - MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >
     *
     */
    template <class MatrixType, int rbegin=0, int rend=TestVariableSetDescription::noOfVariables,
              int cbegin=0, int cend=AnsatzVariableSetDescription::noOfVariables>
    std::unique_ptr<MatrixType> getPointer(bool extractOnlyLowerTriangle) const __attribute__((deprecated)) /*definition at end of file*/;

    /**
     * \brief Extracts the submatrix defined by the half-open block ranges
     * [rbegin,rend), [cbegin,cend).
     *
     * If extractOnlyLowerTriangle is true, the Galerkin operator is treated as
     * if its upper triangle was structurally zero.
     *
     * In order to be able to read the assembled matrix into the smart pointer
     * a specialization of AssemblyDetail::Fill, i.e. AssemblyDetail::Fill<MatrixType>
     * must be available!!!
     *
     * Kaskade provides spezializations of AssemblyDetail::Fill for
     *  - MatrixType = MatrixAsTriplet<Scalar>
     *  - MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >
     *
     * \param extractOnlyLowerTriangle if true, and the matrix is symmetric, only the lower triangular part is extracted
     * \param rbegin number first block row to store
     * \param rend   number of one past the last block row to store
     * \param cbegin number first block column to store
     * \param cend   number of one past the last block column to store
     */
    template <class Matrix>
    Matrix get(bool extractOnlyLowerTriangle=false, int rbegin=0, int rend=TestVariableSetDescription::noOfVariables,
                                                    int cbegin=0, int cend=AnsatzVariableSetDescription::noOfVariables) const
    {
      // block sizes
      auto cSizes = AnsatzVariableSetDescription::variableDimensions(spaces());
      auto rSizes = TestVariableSetDescription::variableDimensions(spaces());

      // row and column offsets
      std::vector<size_t> rowOff(rend-rbegin,0), colOff(cend-cbegin,0);
      std::partial_sum(rSizes.begin()+rbegin,rSizes.begin()+rend-1,rowOff.begin()+1);
      std::partial_sum(cSizes.begin()+cbegin,cSizes.begin()+cend-1,colOff.begin()+1);

      return AssemblyDetail::Fill<Matrix>::apply(getMatrix(), rbegin, cbegin, rowOff, colOff, extractOnlyLowerTriangle,
                                                 nnz(rbegin, rend, cbegin, cend, extractOnlyLowerTriangle), nrows(rbegin,rend), ncols(cbegin,cend));
    }

    template <int row, int col>
    auto const& get() const
    {
      return (*boost::fusion::find_if<AssemblyDetail::IsBlock<row,col>>(getMatrix())).globalMatrix();
    }

    /**
     * \brief Extracts the submatrix defined by the half-open block ranges given as template parameter.
     *
     * Provide the block information via IstlInterface_Detail::BlockInfo.
     */
    template <class MatrixType, class BlockInformation>
    MatrixType get(bool extractOnlyLowerTriangle) const
    {
      return get<MatrixType,BlockInformation::firstRow,BlockInformation::lastRow,BlockInformation::firstCol,BlockInformation::lastCol>(extractOnlyLowerTriangle);
    }

    template <class Matrix, int row, int column>
    std::unique_ptr<Matrix> getBlockPointer(bool extractOnlyLowerTriangle=false) const
    {
      return AssemblyDetail::FillPointerWithBlock<Matrix,row,column>::apply(getMatrix(),extractOnlyLowerTriangle,nnz(row,row+1,column,column+1));
    }

    /**
     * \brief Writes the subvector defined by the half-open block range [rbegin,rend) to the output iterator xi.
     */
    template <class DataOutIter>
    void toSequence(int rbegin, int rend, DataOutIter xi) const
    {
      // WARNING: this assumes that for_each processes the rhsArray in
      // correct order from front to back!
      for_each(getRhs().first.data,BlockToSequence<DataOutIter>(rbegin,rend,xi));
    }


    /**
     * \brief Returns the number of structurally nonzero entries in the
     * submatrix defined by the half-open block ranges [rbegin,rend),
     * [cbegin,cend).
     *
     * This is the number of elements written by the
     * corresponding call of toTriplet().
     */
    size_t nnz(size_t rbegin=0, size_t rend=TestVariableSetDescription::noOfVariables,
               size_t cbegin=0, size_t cend=AnsatzVariableSetDescription::noOfVariables, bool extractOnlyLowerTriangle=false) const
    {
      size_t n = 0;
      boost::fusion::for_each(getMatrix(),CountNonzeros(n,rbegin,rend,cbegin,cend,extractOnlyLowerTriangle));
      return n;
    }

    /**
     * \brief Returns the number of scalar rows in the row block range [firstBlock,lastBlock[.
     */
    size_t nrows(int firstBlock=0, int lastBlock=TestVariableSetDescription::noOfVariables) const
    {
      return TestVariableSetDescription::degreesOfFreedom(spaces(),firstBlock,lastBlock);
    }

    /**
     * \brief Returns the number of scalar cols in the column block range [firstBlock,lastBlock[.
     */
    size_t ncols(int firstBlock=0, int lastBlock=AnsatzVariableSetDescription::noOfVariables) const
    {
      return AnsatzVariableSetDescription::degreesOfFreedom(spaces(),firstBlock,lastBlock);
    }

    /**
     * \brief Returns the list of spaces used.
     */
    Spaces const& spaces() const { return spaces_; }

    /**
     * \brief Returns a mutable reference to the sequence of matrix blocks.
     *
     * \todo: check const correctness!
     */
    MatrixBlockArray& getMatrix() const {
      if(!matrixBlocks) {
        matrixBlocks.reset(new MatrixBlockArray());
        // Choose the number of row groups in the matrix blocks. The obvious choice would be 1 in sequential mode, but it
        // appears that a somewhat larger number is beneficial, probably due to better memory access locality. For multithreaded
        // mode, we choose the number of hardware threads as the default, since then all threads can (in principle) be active
        // simultaneously.
        int nrg = 4;
#ifndef KASKADE_SEQUENTIAL
        nrg = std::max((int)(rowBlockFactor*boost::thread::hardware_concurrency()),4);
#endif
//         boost::fusion::for_each(*matrixBlocks,InitializeBlock(spaces(),gridView.template begin<0>(),gridView.template end<0>(),nrg));
        boost::fusion::for_each(*matrixBlocks,[&](auto& block)
        {
          block.init(this->spaces(),gridView.template begin<0>(),gridView.template end<0>(),nrg);
        });
      }

      return *matrixBlocks;
    }


    /**
     * \brief Returns a contiguous subrange of the rhs coefficient vectors.
     *
     * This returns a boost::fusion view on a sequence of right hand side coefficient vectors.
     * Use it to initialize subvectors of rhs coefficient vectors:
     * \code
     * TestVariableSetDescription::CoefficientVectorRepresentation<0,3>::type r(assembler.rhs<0,3>());
     * \endcode
     */
    template <int first=0, int last=TestVariableSetDescription::noOfVariables>
    auto rhs() const
    {
      return make_range<first,last>(getRhs().first.data);
    }

  private:
    std::pair<RhsArray&,RhsBlockArray&> getRhs() const
    {
      if(rhss.get()==0)
      {
        // construct global vectors
// old version, should be equivalent to the following line.
//        rhss.reset(new RhsArray(Variables_Detail::VariableRangeCreator<TestVariables,ConstructCoefficientVector<Spaces>>
//            ::apply(ConstructCoefficientVector<Spaces>(spaces()))));
        rhss.reset(new RhsArray(TestVariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces())));
        // construct block infos
        rhsBlocks.reset(new RhsBlockArray());
        // Choose the number of row groups in the matrix blocks. The obvious choice would be 1 in sequential mode, but it
        // appears that a somewhat larger number is beneficial, probably due to better memory access locality. For multithreaded
        // mode, we choose the number of hardware threads as the default, since then all threads can (in principle) be active
        // simultaneously.
        int nrg = 4;
#ifndef KASKADE_SEQUENTIAL
        nrg = std::max((int)(rowBlockFactor*boost::thread::hardware_concurrency()),4);
#endif
        boost::fusion::for_each(*rhsBlocks,[&](auto& block)
        {
          block.init(this->spaces(),gridView.template begin<0>(),gridView.template end<0>(),nrg);
        });
      }

      return std::pair<RhsArray&,RhsBlockArray&>(*rhss,*rhsBlocks);
    }


    // Resource management: delete matrix and rhs data structures as well as cell ranges on
    // refinement since they are no longer meaningful on the new grid.
    void reactToRefinement(GridSignals::Status const status)
    {
      if (status == GridSignals::BeforeRefinement)
        flush();
    }


    // resource management
    boost::signals2::scoped_connection refConnection;

    Spaces const&                      spaces_;
    GridManagerBase<Grid> const&       gridManager;
    GridView const&                    gridView;
    IndexSet const&                    indexSet;

    mutable std::shared_ptr<MatrixBlockArray> matrixBlocks; // matrix block info including global matrices
    mutable std::shared_ptr<RhsArray>         rhss;         // global rhs vectors
    mutable std::shared_ptr<RhsBlockArray>    rhsBlocks;    // rhs block info
    BoundaryDetector                          boundaryDetector;
    int                                       validparts;


    int nSimultaneousBlocks;
    double rowBlockFactor;
    size_t localStorageSize;


    template <class IdxOutIter, class DataOutIter>
    struct BlockToTriplet
    {
      BlockToTriplet(size_t rbegin_, size_t cbegin_, std::vector<size_t> const& rowOff_, std::vector<size_t> const& colOff_,
          IdxOutIter& ri_, IdxOutIter& ci_, DataOutIter& xi_, bool extractOnlyLowerTriangle_):
            rbegin(rbegin_), cbegin(cbegin_), rowOff(rowOff_), colOff(colOff_),
            ri(ri_), ci(ci_), xi(xi_), extractOnlyLowerTriangle(extractOnlyLowerTriangle_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        // Check if block is in requested range
        if (inRange(mb.rowId,mb.colId))
          Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ri,ci,xi,
              rowOff[MatrixBlock::rowId-rbegin],colOff[MatrixBlock::colId-cbegin],
              mb.rowId==mb.colId && extractOnlyLowerTriangle,
              mb.symmetric);
        // For mirror blocks, the transposed block needs to be written
        // if the transposed block is in the requested
        // range. Transposition is implicitly achieved by swapping
        // column and row index output iterators.
        if (MatrixBlock::mirror && inRange(mb.colId,mb.rowId) && extractOnlyLowerTriangle==false)
          Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ci,ri,xi,
              colOff[MatrixBlock::rowId-cbegin],rowOff[MatrixBlock::colId-rbegin],
              false,
              mb.symmetric);
      }

      bool inRange(size_t r, size_t c) const
      {
        return r>=rbegin && r<rbegin+rowOff.size() && c>=cbegin && c<cbegin+colOff.size();
      }


      size_t rbegin, cbegin;
      std::vector<size_t> const& rowOff;
      std::vector<size_t> const& colOff;
      IdxOutIter& ri;
      IdxOutIter& ci;
      DataOutIter& xi;
      bool extractOnlyLowerTriangle;
    };

    struct CountNonzeros
    {
      CountNonzeros(size_t& n_, size_t rbegin_, size_t rend_, size_t cbegin_, size_t cend_, bool onlyLowerTriangle_):
        n(n_), rbegin(rbegin_), rend(rend_), cbegin(cbegin_), cend(cend_), onlyLowerTriangle(onlyLowerTriangle_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        size_t myN = Matrix_to_Triplet<typename MatrixBlock::Matrix>::nnz(mb.globalMatrix(),
                                                                          (mb.rowId==mb.colId) && onlyLowerTriangle,
                                                                          mb.symmetric);
        // Check if block is in requested range
        if (inRange(MatrixBlock::rowId,MatrixBlock::colId))
          n += myN;

        // For subdiagonal blocks, the transposed block needs to be
        // counted if onlyLowerTriangle is false and the transposed
        // block is in the requested range and we need to mirror this
        // block.
        if (MatrixBlock::mirror && inRange(mb.colId,mb.rowId) && onlyLowerTriangle==false)
          n += myN;
      }

      bool inRange(size_t r, size_t c) const
      {
        return r>=rbegin && r<rend && c>=cbegin && c<cend;
      }



      size_t& n;
      size_t rbegin, rend, cbegin, cend;
      bool onlyLowerTriangle;
    };

    template <class DataOutIter>
    struct BlockToSequence
    {
      BlockToSequence(int& rbegin_, int& rend_, DataOutIter& out_): rbegin(rbegin_), rend(rend_), out(out_) {}
      template <class VectorBlock> void operator()(VectorBlock const& v) const
      {
        if (rbegin<=0 && rend>0)
          out = vectorToSequence(v,out);
        --rbegin;
        --rend;
      }
    private:
      int& rbegin;
      int& rend;
      DataOutIter& out;
    };

  };


  template <class F, class BoundaryDetector,class QuadRule>
  template <class MatrixType, int rbegin, int rend, int cbegin, int cend>
  std::unique_ptr<MatrixType> VariationalFunctionalAssembler<F,BoundaryDetector,QuadRule>::getPointer(bool extractOnlyLowerTriangle) const
  {
    using namespace boost::fusion;

    // block sizes
    auto cSizes = AnsatzVariableSetDescription::variableDimensions(spaces());
    auto rSizes = TestVariableSetDescription::variableDimensions(spaces());

    // row and column offsets
    std::vector<size_t> rowOff(rend-rbegin,0), colOff(cend-cbegin,0);
    std::partial_sum(rSizes.begin()+rbegin,rSizes.begin()+rend-1,rowOff.begin()+1);
    std::partial_sum(cSizes.begin()+cbegin,cSizes.begin()+cend-1,colOff.begin()+1);

    return AssemblyDetail::FillPointer<MatrixType>::apply(getMatrix(), rbegin, cbegin, rowOff, colOff, extractOnlyLowerTriangle,
        nnz(rbegin, rend, cbegin, cend, extractOnlyLowerTriangle), nrows(rbegin,rend), ncols(cbegin,cend));
  }
  // End of VariationalFunctionalAssembler
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
} /* end of namespace Kaskade */

#endif
