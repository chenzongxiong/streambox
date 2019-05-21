/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef JACOBI_PRECONDITIONER_HH
#define JACOBI_PRECONDITIONER_HH

#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvercategory.hh"

#include "fem/fixdune.hh"
#include "fem/istlinterface.hh"

#include "linalg/matrixTraits.hh"
#include "linalg/symmetricOperators.hh"

#include "utilities/threading.hh"

namespace Kaskade
{
  /// \internal
  namespace JacobiPreconditionerDetail
  {
    template <class Entry, int row=0, int col=0>
    struct DiagonalBlock  // TODO: implement move constructor for efficient return by ExtractDiagonal functor
    {
      using field_type = typename Entry::field_type;
      
      DiagonalBlock() = default;
      DiagonalBlock(DiagonalBlock&& other) = default;
      DiagonalBlock& operator=(DiagonalBlock&& other) = default;
      
      template <class Matrix, class Weight>
      DiagonalBlock(Matrix const& mat, Weight w)
      : size(mat.N()), dps(NumaThreadPool::instance().nodes())
      {
        assert(mat.N() == mat.M()); // make sure it's quadratic
        size_t const nodes = NumaThreadPool::instance().nodes();
        for (size_t i=0; i<nodes; ++i)
          diag.push_back(std::vector<Entry,NumaAllocator<Entry>>(uniformWeightRangeStart(i+1,nodes,size)-uniformWeightRangeStart(i,nodes,size),
                                                                 Entry(),NumaAllocator<Entry>(i)));

        // copy and invert diagonal
        for (auto ri=mat.begin(); ri!=mat.end(); ++ri)
        {
          auto const i = ri.index();
          size_t const node = uniformWeightRange(i,nodes,size);     // which NUMA node
          size_t const k = i - uniformWeightRangeStart(node,nodes,size); // which offset in node block
          
          if (k<0 || k>=diag[node].size()) 
          {
            std::cerr << "i=" << i << " nodes=" << nodes << " node=" << node << " k=" << k << " => i=" << k + (node*size)/nodes << " diag[node].size=" << diag[node].size() << " n=" << size << "\n";
            abort();
          }
          
          auto ci = ri->begin();
          auto cend = ri->end();
          for (; ci!=cend && ci.index()!=i; ++ci)      // find diagonal entry TODO: more efficient search
            ;
          Entry& a = diag[node][k];
          if (ci.index() == i) 
          {
            a = *ci;
            a.invert(); // TODO: this calls DynamicMatrix::invert and performs memory allocation even for fixed size 3x3 matrices via std::vector...
            a *= w;
          }
          else
            a = 0; // ignore zero diagonal elements
        }
      }
      
      /**
       * \brief Computes \f$ x \leftarrow D^{-1} y \f$ and returns \f$ y^T x \f$.
       * \tparam X a vector class with element access
       * \tparam Y a vector class with element access
       */
      template <class X, class Y>
      field_type apply(X& x, Y const& y) const
      {
        // Perform one Jacobi step for Ax=y, i.e. for A = L+D+R we compute x = D^{-1}y 
        // Note that the Dune::SeqJac preconditioner as of Dune 2.2.1 allocates and frees a new vector in every
        // iteration, and iterates each time over the whole matrix, which is quite costly. 
        // This implementation is roughly five times faster.
       std::fill(dps.begin(),dps.end(),0.0);
       parallelForNodes([&x,&y,this](size_t node, size_t nodes)
       {
         assert(nodes==NumaThreadPool::instance().nodes());
         size_t const start = uniformWeightRangeStart(node,nodes,this->size);
         auto const& d = this->diag[node];
         size_t const n = d.size();
         field_type dp = 0;
         for (size_t k=0; k<n; ++k)
         {
           size_t const i = k + start;
           x[i] = d[k] * y[i];
           dp += x[i] * y[i];
         }
         this->dps[node] = dp;
       });
       return std::accumulate(dps.begin(),dps.end(),0.0);
      }
       
      
      static int const rowId = row;
      static int const colId = col;

      size_t size;
      
      // diag[i] contains entries [i*n/nodes, (i+1)*n/nodes[
      std::vector<std::vector<Entry,NumaAllocator<Entry>>> diag;
      
    private:
      mutable std::vector<field_type> dps;   // for storing duality pairings of parallel blocks (could be local to apply, but is here to prevent reallocation
    };
    
    // --------------------------------------------------------------------------------------------
    
    // boost::fusion functor extracting, inverting, and scaling the diagonal of a matrix
    struct ExtractDiagonal
    {
      ExtractDiagonal(double w_): w(w_) {}
      
      template <class T> struct result {};
      
      template <class Block>
      struct result<ExtractDiagonal(Block)> {
        typedef typename std::remove_reference<Block>::type B;
        typedef DiagonalBlock<typename B::BlockType,B::rowId,B::colId> type;
      };
      
      template <class MB>
      typename result<ExtractDiagonal(MB)>::type operator()(MB const& mb) const 
      {
        return typename result<ExtractDiagonal(MB)>::type(*mb.matrix,w);
      }
      
    private:
      double w;
    };
    
    // --------------------------------------------------------------------------------------------

    template <class Matrix, class Domain=typename MatrixTraits<Matrix>::NaturalDomain, class Range=typename MatrixTraits<Matrix>::NaturalRange>
    class JacobiPreconditionerBase: public SymmetricPreconditioner<Domain,Range>
    {
    public:
      static int const category = Dune::SolverCategory::sequential;
      
      using matrix_type = Matrix;
      using domain_type = Domain;
      using range_type = Range;
      using field_type = typename Matrix::field_type;
      
      /**
       * \brief Default constructor.
       * 
       * Specifying nothing, this leaves the preconditioner in a pretty useless state. Assign a reasonable value before using.
       */
      JacobiPreconditionerBase() = default;
      
      JacobiPreconditionerBase(JacobiPreconditionerBase&&) = default;
      JacobiPreconditionerBase& operator=(JacobiPreconditionerBase&&) = default;

      /**
       * \brief Constructor.
       * \param mat the matrix
       */
      JacobiPreconditionerBase(Matrix const& mat)
      : diag(mat,1.0)
      {  }
      
      virtual void apply(domain_type& x, range_type const& b) { applyDp(x,b); }
      
      virtual field_type applyDp (domain_type& x, range_type const& y)
      {
        return diag.apply(x,y);
      }
      
      virtual bool requiresInitializedInput() const { return false; }
      
    private:
      DiagonalBlock<typename matrix_type::block_type> diag;
    };
  }
  /// \endinternal
  
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup iterative
   * \brief A simple single-level Jacobi preconditioner working on AssembledGalerkinOperator.
   * \tparam Operator satisfying the Dune::AssembledLinearOperator concept.
   */
  template <class Operator>
  class JacobiPreconditioner: public JacobiPreconditionerDetail::JacobiPreconditionerBase<typename Operator::matrix_type,
                                                                                          typename Operator::domain_type,
                                                                                          typename Operator::range_type>
  {
    using Base = JacobiPreconditionerDetail::JacobiPreconditionerBase<typename Operator::matrix_type,typename Operator::domain_type,typename Operator::range_type>;
    using field_type = typename Base::field_type;
    
  public:
    /**
     * \brief Constructor.
     * \param op the linear operator to be preconditioned
     */
    JacobiPreconditioner(Operator const& op)
    : Base(op.getmat())
    { }
  };
  
  /**
   * \ingroup iterative
   * \brief A simple single-level Jacobi preconditioner working on sparse matrices.
   */
  template <class Entry, class Index>
  class JacobiPreconditioner<NumaBCRSMatrix<Entry,Index>>: public JacobiPreconditionerDetail::JacobiPreconditionerBase<NumaBCRSMatrix<Entry,Index>>
  {
    using Base = JacobiPreconditionerDetail::JacobiPreconditionerBase<NumaBCRSMatrix<Entry,Index>>;
  public:
    using field_type = typename Base::field_type;
    
    /**
     * \brief Default constructor.
     * 
     * Specifying nothing, this leaves the preconditioner in a pretty useless state. Assign a reasonable value before using.
     */
    JacobiPreconditioner() = default;
    
    JacobiPreconditioner(JacobiPreconditioner&&) = default;
    
    JacobiPreconditioner& operator=(JacobiPreconditioner&&) = default;

    /**
     * \brief Constructor.
     * \param mat the matrix to be preconditioned
     */
    JacobiPreconditioner(NumaBCRSMatrix<Entry,Index> const& mat)
    : Base(mat)
    {}
  };
  
  /**
   * \ingroup iterative
   * \brief Convenience function for creating Jacobi preconditioners.
   */
  template <class Operator>
  auto makeJacobiPreconditioner(Operator const& op)
  {
    return JacobiPreconditioner<Operator>(op);
  }
  
  /**
   * \ingroup iterative
   * \brief A recursive Jacobi preconditioner for heterogeneous square block matrices from a variational assembler.
   *
   * This is more memory efficient than the Dune::SeqJac implementation working on the
   * AssembledGalerkinOperator matrix, since the extraction of the matrices from
   * the VariationalFunctionalAssembler is avoided here. Only the diagonal elements
   * are extracted (and inverted once on construction).
   *
   * As a (small) drawback, the current implementation is limited to a single Jacobi iteration, as 
   * it cannot compute residual vectors.
   */
  template <class Assembler,
            int firstRow, int lastRow,
            int firstCol, int lastCol,
            class BlockFilter,
            bool symmetric>
  class JacobiPreconditioner<AssembledGalerkinOperator<Assembler,firstRow,lastRow,firstCol,lastCol,BlockFilter,symmetric>> 
  : public SymmetricPreconditioner<typename AssembledGalerkinOperator<Assembler,firstRow,lastRow,firstCol,lastCol,BlockFilter,symmetric>::Domain, 
                                   typename AssembledGalerkinOperator<Assembler,firstRow,lastRow,firstCol,lastCol,BlockFilter,symmetric>::Range>
  {
    typedef AssembledGalerkinOperator<Assembler,firstRow,lastRow,firstCol,lastCol,BlockFilter,symmetric> Operator;
    typedef typename Operator::Domain Domain;
    typedef typename Operator::Range Range;
    typedef typename Operator::Scalar Scalar;

    // A block filter that selects only the diagonal blocks of a heterogeneous block matrix.
    struct DiagonalFilter {
      template <class Block> struct apply { typedef boost::mpl::bool_<Block::rowId==Block::colId> type; };
    };


  public:
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \brief Constructs a Jacobi preconditioner.
     *
     * \arg op the assembled galerkin operator
     * \arg w_ the relaxation factor
     */
    explicit JacobiPreconditioner(Operator const& op, Scalar w=1.0):
      blocks(boost::fusion::transform(
        boost::fusion::filter_if<DiagonalFilter>(
          boost::fusion::transform(boost::fusion::filter_if<BlockFilter>(IstlInterfaceDetail::AllBlocks<typename Operator::Assembler>::apply(op.getAssembler())),
                                   IstlInterfaceDetail::Translate<firstRow,firstCol>())),
        JacobiPreconditionerDetail::ExtractDiagonal(w)))
    {}

    /**
     * \brief Computes Jacobi step \f$ x = D^{-1} y \f$ for \f$ A=L+D+R \f$ and linear equation \f$ Ax=y \f$.
     * 
     * Note that one Jacobi step for each block on the diagonal of the heterogeneous block operator is performed.
     */
    virtual void apply (Domain& x, Range const& y)
    {
      applyDp(x,y);
    }

    /**
     * \brief Computes Jacobi step \f$ x = D^{-1} y \f$ and \f$ x^Ty \f$ for \f$ A=L+D+R \f$ and linear equation \f$ Ax=y \f$.
     * 
     * Note that one Jacobi step for each block on the diagonal of the heterogeneous block operator is performed.
     */
    virtual Scalar applyDp (Domain& x, Range const& y)
    {
      Scalar dp = 0;
      boost::fusion::for_each(blocks,ApplyJacobi(x,y,dp));
      return dp;
    }

    /** 
     * \brief Returns true if the target vector x has to be initialized to zero before calling apply or applyDp
     * 
     * Does not need init to zero.
     */
    virtual bool requiresInitializedInput() const { return false; }

  private:
    // a boost::fusion list of all the matrix blocks to consider. Select all blocks in the subranges covered by the
    // operator, shift the indices to start at 0, then select the ones on the diagonal and extract, invert, and store
    // the diagonal.
    typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::transform<
        typename boost::fusion::result_of::filter_if<
          typename boost::fusion::result_of::transform<
            typename boost::fusion::result_of::filter_if<typename IstlInterfaceDetail::AllBlocks<typename Operator::Assembler>::result,BlockFilter>::type,
            IstlInterfaceDetail::Translate<firstRow,firstCol> >::type,
          DiagonalFilter>::type,
        JacobiPreconditionerDetail::ExtractDiagonal>::type>::type blocks;

    struct ApplyJacobi
    {
      ApplyJacobi(Domain& x_, Range const& y_, Scalar& dp): x(x_), y(y_), dualPairing(dp) {}

      template <class Block>
      void operator()(Block const& mb) const 
      {
        using namespace boost::fusion;
        dualPairing += mb.apply(at_c<Block::colId>(x.data),at_c<Block::rowId>(y.data));
      }

    private:
      Domain&      x;
      Range const& y;
      Scalar&      dualPairing; 
    };
    
  };
} // namespace Kaskade
#endif
