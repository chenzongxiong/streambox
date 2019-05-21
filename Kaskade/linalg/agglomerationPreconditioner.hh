/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AGGLOMERATIONPRECONDITIONER_HH
#define AGGLOMERATIONPRECONDITIONER_HH

#include <fem/views.hh>
#include <linalg/jacobiPreconditioner.hh>
#include <linalg/symmetricOperators.hh>
#include <linalg/threadedMatrix.hh>


namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace AgglomerationPreconditionerDetails
  {
    template <class Entry>
    class AgglomerationPreconditionerEngine
    {
      using field_type = typename Entry::field_type;
      static int const blocksize = Entry::rows;        
      
    public:
      // Matrix has to satisfy the Dune::NumaBCRSMatrix interface
      template <class Matrix, class Space>
      AgglomerationPreconditionerEngine(Matrix const& A, Space const& space)
      {
        using std::begin; using std::end;
        
        // The coarse grid projection is B = P^T A P, and its diagonal that we use here is
        // B_{ii} = \sum_{k,l} (P^T)_{ik}^T A_{kl} P_{li}. Entries in the prolongation P are
        // just ones (or identities in case of multicomponent systems), such that we can 
        // drop the transposition of the entries. Accessing P row-wise is more efficient in 
        // the CRS data structures, hence we write B_{ii} = \sum_{k,l} (P^T)_{ik} A_{kl} (P^T)_{il}.
        // Let K_j denote the set of column indices of nonzero entries in the i-th row of P^T.
        // Then we can write B_{ii} = \sum_{k,l\in K_i} A_{kl}.
        
        // create prolongation
        size_t const n = space.degreesOfFreedom();  // rows: total degrees of freedom 
        size_t const m = space.indexSet().size(0);  // columns: number of cells
        NumaCRSPatternCreator<> creator(n,m);
        for (auto ci=space.gridView().template begin<0>(); ci!=space.gridView().template end<0>(); ++ci)
        {
          auto cidx = space.indexSet().index(*ci);                       // cell index is coarse (column) index
          auto const& gidx = space.mapper().globalIndices(cidx);         // FE indices are fine (row) indices
          creator.addElements(begin(gidx),end(gidx),&cidx,&cidx+1,true); // scatter into the one column (col index is trivially sorted)
        }
        P = NumaBCRSMatrix<Entry>(std::make_shared<NumaCRSPattern<>>(creator),          // create prolongation matrix with all
                                  unitMatrix<typename Entry::field_type,blocksize>());  // nonzero entries equal to 1 (identity
                                                                                        // for multicomponent systems).
        Pt = NumaBCRSMatrix<Entry>(P,false,true,false);  // create transposed matrix
        
        // create and invert diagonal of projected Galerkin matrix
        inverseDiagonalPAP.resize(m,0);
        for (size_t i=0; i<m; ++i)
        {
          auto const& row = P[i];
          auto ki = rangeView(row.getindexptr(),row.getindexptr()+row.size()); // column indices of i-th row of P (sorted)
          
          // sum up A_{kl} for k,l \in K_i
          for (auto k: ki)
          {
            auto const& Ak = A[k];
            auto al = begin(Ak);  // iterator of row A_k
            auto l  = begin(ki);  // iterator of K_i
            while (al!=end(Ak) && l!=end(ki))
            {
              if (al.index()==*l) // match found A_{kl} with l \in K_i
              {
                inverseDiagonalPAP[i] += *al;
                ++al;
                ++l;
              }
              else if (*l < al.index()) // no match: advance smaller index
                ++l;
              else
                ++al;
            }
            // invert (P^TAP)_ii
            inverseDiagonalPAP[i].invert();
          }
        }
      }
      
      // adds P (P^TAP)^{-1} P^T y to x
      template <class Domain, class Range>
      void applyAgglomerated(Domain& x, Range const& y) 
      {
        // Compute P (P^TAP)^{-1} P^T y
        size_t const m = inverseDiagonalPAP.size();
        Dune::BlockVector<Dune::FieldVector<field_type,blocksize>> z(m);
        
        Pt.mv(y,z);                             // z = P^T y
        for (size_t i=0; i<m; ++i)              // z = (P^TAP)^{-1} z
          z[i] = inverseDiagonalPAP[i] * z[i];  //   = (P^TAP)^{-1}P^T y
        P.umv(z,x);                             // x <- P z = P(P^TAP)^{-1}P^T y
      }
      
       // adds P (P^TAP)^{-1} P^T y to x
      template <class Domain, class Range>
      field_type applyAgglomeratedDp(Domain& x, Range const& y) 
      {
        applyAgglomerated(x,y);
        
        field_type dp = 0;                 // compute dp = x^T y
        for (size_t i=0; i<P.N(); ++i)
          dp += x[i] * y[i];
        return dp;
      }
      
    private:
      std::vector<Entry> inverseDiagonalPAP;  // diag(P^TAP)^{-1}
      NumaBCRSMatrix<Entry> P, Pt;            // prolongation P and restriction Pt, stored separately here because application of transpose is inefficient      
    };
  }
  
  /**
   * \endcond
   */
  
  
  /**
   * \ingroup linalg
   * \brief An agglomeration-based preconditioner for elliptic problems discretized with
   *        higher-order Lagrange finite elements.
   * 
   * This implements an additive two-level subspace correction preconditioner. The first (fine) level
   * consists of the individual ansatz functions, i.e. a Jacobi preconditioner. The second (coarse) 
   * level consists of functions constant one on a single cell, zero on cells not incident to vertices
   * of that cell, and something else in between. A drawback of the approach is that 
   * the projected Galerkin matrix on the coarse space is denser than a usual linear FE stiffness
   * matrix. Thus, the coarse space functions span individual subspaces (which realizes a Jacobi smoother 
   * on the coarse space), such that we only need to store the diagonal of the projected Galerkin matrix.
   * 
   * For higher order ansatz spaces, the dimension of the coarse space is quite small. 
   * 
   * \tparam Operator the Galerkin operator
   * 
   * The Galerkin operator type can satisfy either the Dune::AssembledLinearOperator interface or, by specialization,
   * be an AssembledGalerkinOperator of block size 1 x 1. For the Dune::AssembledLinearOperator, the matrix type
   * has to implement the Dune::BCRSMatrix interface, the domain and range types have to implement the Dune::BlockVector
   * interface (e.g. FunctionSpaceElement).
   */
  template <class Operator>
  class AgglomerationPreconditioner: public SymmetricPreconditioner<typename Operator::domain_type,typename Operator::range_type>
  {
    using Base = SymmetricPreconditioner<typename Operator::domain_type,typename Operator::range_type>;
    using typename Base::field_type;
    using typename Base::domain_type;
    using typename Base::range_type;
    
    static int const blocksize = range_type::block_type::dimension;
    static_assert(blocksize==domain_type::block_type::dimension,"Operator matrix entries must be square for Jacobi preconditioning.");
    using Entry = Dune::FieldMatrix<field_type,blocksize,blocksize>;
    
    
  public:
    static int const category = Dune::SolverCategory::sequential;
    
    /**
     * \brief Constructor.
     * \tparam Space Lagrangian finite element space
     * \param opA the assembled symmetric Galerkin operator
     * \param space the Lagrangian FE space with which the operator A is discretized
     */
    template <class Space>
    AgglomerationPreconditioner(Operator const& opA, Space const& space)
    : jac1(opA), // fine grid Jacobi preconditioner
      engine(opA.getmat(),space)
    {}
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ 
     */
    virtual void apply(domain_type& x, range_type const& y) 
    {
      jac1.applyDp(x,y);          // apply fine grid Jacobi
      engine.apply(x,y);          // apply agglomerated Jacobi
    }
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ and returns \f$ \langle By, y \rangle \f$.
     */
    virtual field_type applyDp(domain_type& x, range_type const& y) 
    {
      jac1.applyDp(x,y);          // apply fine grid Jacobi
      return engine.applyDp(x,y); // apply agglomerated Jacobi
    }
    
    /** 
     * \brief Returns true if the target vector x has to be initialized to zero before calling apply or applyDp
     */
    virtual bool requiresInitializedInput() const 
    {
      // we apply the fine grid Jacobi first, hence rely on that.
      return jac1.requiresInitializedInput();
    }
    
  private:
    JacobiPreconditioner<Operator> jac1;    // diag(A)^{-1}
    AgglomerationPreconditionerDetails::AgglomerationPreconditionerEngine<Entry> engine;
  };
  
  
  
  template <class Assembler,
            int firstRow, 
            int firstCol>
  class AgglomerationPreconditioner<AssembledGalerkinOperator<Assembler,firstRow,firstRow+1,firstCol,firstCol+1,
                                                              IstlInterfaceDetail::RangeBlockSelector<firstRow,firstRow+1,firstCol,firstCol+1>,
                                                              true>>
  : public SymmetricPreconditioner<typename AssembledGalerkinOperator<Assembler,firstRow,firstRow+1,firstCol,firstCol+1,
                                                              IstlInterfaceDetail::RangeBlockSelector<firstRow,firstRow+1,firstCol,firstCol+1>,
                                                              true>::domain_type,typename AssembledGalerkinOperator<Assembler,firstRow,firstRow+1,firstCol,firstCol+1,
                                                              IstlInterfaceDetail::RangeBlockSelector<firstRow,firstRow+1,firstCol,firstCol+1>,
                                                              true>::range_type>
  {
    using Operator = AssembledGalerkinOperator<Assembler,firstRow,firstRow+1,firstCol,firstCol+1,
                                                              IstlInterfaceDetail::RangeBlockSelector<firstRow,firstRow+1,firstCol,firstCol+1>,
                                                              true>;
    using Base = SymmetricPreconditioner<typename Operator::domain_type,typename Operator::range_type>;
    using typename Base::field_type;
    using typename Base::domain_type;
    using typename Base::range_type;
    
    using RangeEntry = typename boost::fusion::result_of::value_at_c<typename range_type::Sequence,0>::type::block_type;
    using DomainEntry = typename boost::fusion::result_of::value_at_c<typename domain_type::Sequence,0>::type::block_type;
    static int const blocksize = RangeEntry::dimension;
    static_assert(blocksize==DomainEntry::dimension,"Matrixe entries have to be square for Jacobi preconditioning");
    using Entry = Dune::FieldMatrix<field_type,blocksize,blocksize>;
    
    
  public:
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \brief Constructor.
     * \tparam Space Lagrangian finite element space
     * \param opA the assembled symmetric Galerkin operator
     * \param space the Lagrangian FE space with which the operator A is discretized
     */
    template <class Space>
    AgglomerationPreconditioner(Operator const& opA, Space const& space)
    : jac1(opA), // fine grid Jacobi preconditioner
      engine(opA.template get<NumaBCRSMatrix<Entry>>(),space)
    {}
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ 
     */
    virtual void apply(domain_type& x, range_type const& y) 
    {
      jac1.applyDp(x,y);   // apply fine grid Jacobi
      auto const& y0 = boost::fusion::at_c<0>(y.data);
      auto& x0 = boost::fusion::at_c<0>(x.data);
      engine.applyDp(x0,y0); // apply agglomerated Jacobi
    }
    
    /**
     * \brief Computes \f$ x \leftarrow By \f$ and returns \f$ \langle By, y \rangle \f$.
     */
    virtual field_type applyDp(domain_type& x, range_type const& y) 
    {
      jac1.applyDp(x,y);   // apply fine grid Jacobi
      auto const& y0 = boost::fusion::at_c<0>(y.data);
      auto& x0 = boost::fusion::at_c<0>(x.data);
      return engine.applyDp(x0,y0); // apply agglomerated Jacobi
    }
    
    /** 
     * \brief Returns true if the target vector x has to be initialized to zero before calling apply or applyDp
     */
    virtual bool requiresInitializedInput() const 
    {
      // we apply the fine grid Jacobi first, hence rely on that.
      return jac1.requiresInitializedInput();
    }
    
  private:
    JacobiPreconditioner<Operator> jac1;    // diag(A)^{-1}
    AgglomerationPreconditionerDetails::AgglomerationPreconditionerEngine<Entry> engine;
  };
  
  
}

#endif