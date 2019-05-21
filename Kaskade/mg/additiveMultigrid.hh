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

#ifndef MG_AMULTIGRID_HH
#define MG_AMULTIGRID_HH

#include <future>
#include <memory>
#include <type_traits>

#include <dune/istl/preconditioner.hh>

#include "fem/spaces.hh"
#include "linalg/conjugation.hh"
#include "linalg/direct.hh"
#include "linalg/domainDecompositionPreconditioner.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "mg/prolongation.hh"

namespace Kaskade
{
  /**  
   * \ingroup multigrid
   * \brief An additive multigrid preconditioner for P1 elements.
   * 
   * This realizes an additive multigrid preconditioner. With a Jacobi smoother, the BPX preconditioner is obtained.
   * 
   * \tparam Smoother a smoother type satisfying the SymmetricPreconditioner interface. The smoother must \em not rely 
   * on the matrix provided on construction to be available subsequently. Hence it should copy all required parts.
   * 
   * Convenience functions for creating additive multigrid preconditioners are available as makeAdditiveMultiGrid and makeBPX.
   * 
   * \todo: consider inheriting from Kaskade::SymmetricPreconditioner
   */
  template <class Smoother, class Prolongation, class CoarsePreconditioner=Smoother>
  class AdditiveMultiGrid: public SymmetricPreconditioner<typename Smoother::domain_type,typename Smoother::range_type>
  {
  public:
    using domain_type = typename Smoother::domain_type;
    using range_type = typename Smoother::range_type;
    using Base = SymmetricPreconditioner<domain_type,range_type>;
    
    /**
     * \brief Default constructor.
     */
    AdditiveMultiGrid() = default;
    
    /**
     * \brief Constructor.
     * \tparam Entry the type of Galerkin matrix entries
     * \tparam Index the index type of Galerkin matrices (usually size_t)
     */
    template <class Entry, class Index, class MakeSmoother, class MakeCoarsePreconditioner>
    AdditiveMultiGrid(NumaBCRSMatrix<Entry,Index> A, std::vector<Prolongation> Ps, MakeSmoother const& makeSmoother, 
                      MakeCoarsePreconditioner const& makeCoarsePreconditioner, bool onlyLowerTriangle=false)
    : prolongations(std::move(Ps))
    {
      for (int l=prolongations.size(); l>0; --l)
      {
        smoothers.insert(smoothers.begin(),makeSmoother(A,l));
        A = conjugation(prolongations[l-1],A,onlyLowerTriangle);
      }
      coarsePreconditioner = std::make_unique<CoarsePreconditioner>(makeCoarsePreconditioner(A,0));
    }
    
    AdditiveMultiGrid(AdditiveMultiGrid&& other) = default;
    
    AdditiveMultiGrid& operator=(AdditiveMultiGrid&& other) = default;
    
    virtual void apply(domain_type& x, range_type const& r)
    {
      range_type b = r;
      runMG(x,b,prolongations.size());
    }
    
    virtual typename Base::field_type applyDp(domain_type& x, range_type const& r)
    {
      apply(x,r);
      return x*r;
    }

    virtual bool requiresInitializedInput() const
    {
      for (auto const& s: smoothers)
        if (s.requiresInitializedInput())
          return true;
      return false;
    }
    
  private:
    void runMG(domain_type& x, range_type& r, int const level) 
    {
      if (level == 0)
      {
        x = 0;
        coarsePreconditioner->pre(x,r);
        coarsePreconditioner->apply(x,r);
        coarsePreconditioner->post(x);
      }
      else
      {
        // first the restriction -- the smoother may modify the residual r, but the coarse grid correction won't...
        auto const& p = prolongations[level-1];
        range_type cr(p.M());
        p.mtv(r,cr);                            // cr = P^T r
        domain_type cx;                         // for the coarse grid correction
        
        // Run the coarse grid correction in parallel. Depending on the prolongation structure and
        // the machine this code runs on, this appears to have a minor to moderate benefit. We use 
        // std::async instead of the thread pool in order to prevent deadlocks due to recursive
        // submission of tasks.
        auto coarseGridCorrectionTicket = std::async(std::launch::async,[&] ()
        { 
          cx.resize(p.M()); cx = 0;
          runMG(cx,cr,level-1);                 // cx ~ (P^T A P)^{-1} cr
        });

        // Now the smoothing while the coarse grid correction is on the way
        auto& s = smoothers[level-1];             // Smoother S
        s.pre(x,r);
        s.apply(x,r);
        s.post(x);
        
        // now wait for the coarse grid correction to be available and add it up
        coarseGridCorrectionTicket.wait();
        p.umv(cx,x);                            // dx = dx + P*cx
      }
    }
    
    std::vector<Prolongation>             prolongations;
    std::vector<Smoother>                 smoothers;
    std::unique_ptr<CoarsePreconditioner> coarsePreconditioner;
  };
  
  // ---------------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup multigrid
   * \brief Convenience function creating additive multigrid preconditioners for P1 elements.
   * 
   * \param A the sparse Galerkin matrix for P1 finite elements on the leaf view of the given simplicial grid
   * \param grid a simplicial grid
   * \param makeSmoother a callable object taking a (projected) Galerkin matrix of type NumaBCRSMatrix<Entry,Index> and a level, returning a preconditioner
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   * 
   * \see AdditiveMultiGrid
   */
  template <class Entry, class Index, class Prolongation, class MakeSmoother, class MakeCoarsePreconditioner>
  auto makeAdditiveMultiGrid(NumaBCRSMatrix<Entry,Index> const& A, std::vector<Prolongation> const& prolongations, 
                             MakeSmoother const& makeSmoother, MakeCoarsePreconditioner const& makeCoarsePreconditioner,
                             bool onlyLowerTriangle=false)
  {
    using Smoother = std::result_of_t<MakeSmoother(NumaBCRSMatrix<Entry,Index>,int)>;
    using CoarsePreconditioner = std::result_of_t<MakeCoarsePreconditioner(NumaBCRSMatrix<Entry,Index>,int)>;
    return AdditiveMultiGrid<Smoother,Prolongation,CoarsePreconditioner>(A,prolongations,makeSmoother,makeCoarsePreconditioner,onlyLowerTriangle);
  }
 
  
  /**
   * \ingroup multigrid
   * \brief Convenience function creating a BPX preconditioner for P1 elements.
   * 
   * \param A the sparse Galerkin matrix for P1 finite elements on the leaf view of the given simplicial grid
   * \param grid a simplicial grid
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   * 
   * \see AdditiveMultiGrid
   */
  template <class Entry, class Index, class GridMan>
  auto makeBPX(NumaBCRSMatrix<Entry,Index> const& A, GridMan const& gridman, bool onlyLowerTriangle=false)
  {
    using Traits = MatrixTraits<NumaBCRSMatrix<Entry,Index>>;
    auto smoother = [](auto const& a, int ){ return makeJacobiPreconditioner(a); };
    auto makeCoarsePreconditioner = [](auto const& a, int ) { return DirectSolver<typename Traits::NaturalDomain,typename Traits::NaturalRange>(a); };    
    return makeAdditiveMultiGrid(A,prolongationStack(gridman),smoother,makeCoarsePreconditioner,onlyLowerTriangle);
  }
 
  // ---------------------------------------------------------------------------------------------------------
  
  
  /**
   * \ingroup multigrid
   * \brief Convenience function creating additive multigrid preconditioner of V-cycle type for higher order elements.
   * 
   * This is realized as an additive h-multigrid nested within a two-grid scheme. 
   * The outer two-level algorithm uses the smoother obtained from \arg makeSmoother on the p-elements defined by the \arg space, 
   * and uses the h-multigrid as a coarse level preconditioner/solver.
   * 
   * Note that a pure Jacobi preconditioner becomes quickly inefficient when the polynomial ansatz order increases. From degree
   * 3 on, using an overlapping domain decomposition smoother such as PatchDomainDecompositionPreconditioner is probably more efficient.
   * 
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param p1Space a linear finite element space
   * \param makeSmoother a callable object taking a (projected) Galerkin matrix of type NumaBCRSMatrix<Entry,Index> and a level, returning a preconditioner
   * \param makeCoarseSmoother a callable object taking a (projected) Galerkin matrix of type NumaBCRSMatrix<Entry,Index> and a level, returning a preconditioner
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   * 
   * \see AdditivePMultiGrid
   * \see PatchDomainDecompositionPreconditioner
   */
  template <class Entry, class Index, class FineSpace, class MakeSmoother, class MakeCoarseSmoother>
  auto makeAdditivePMultiGrid(NumaBCRSMatrix<Entry,Index> A, FineSpace const& space, H1Space<typename FineSpace::Grid> const& p1Space, 
                               MakeSmoother const& makeSmoother, MakeCoarseSmoother const& makeCoarseSmoother, bool onlyLowerTriangle=false)  
  {
    using Traits = MatrixTraits<NumaBCRSMatrix<Entry,Index>>;
    auto makeCoarsePreconditioner = [&](NumaBCRSMatrix<Entry,Index> const& a, int)
    {
      auto makeCoarseSolver = [](auto const& a, int ) { return DirectSolver<typename Traits::NaturalDomain,typename Traits::NaturalRange>(a); }; 
      return makeAdditiveMultiGrid(a,prolongationStack(space.gridManager()),makeCoarseSmoother,makeCoarseSolver,onlyLowerTriangle);
    };
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>> prolongations(1,prolongation(p1Space,space));
    return makeAdditiveMultiGrid(A,prolongations,makeSmoother,makeCoarsePreconditioner,onlyLowerTriangle);
  }
  
  /**
   * \ingroup multigrid
   * \brief Convenience function creating additive multigrid preconditioner of V-cycle type with domain decomposition and 
   * Jacobi smoother for higher order elements.
   * 
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param p1Space a linear finite element space
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   * 
   * \see AdditivePMultiGrid
   */
  template <class Entry, class Index, class FineSpace>
  auto makePBPX(NumaBCRSMatrix<Entry,Index> A, FineSpace const& space, H1Space<typename FineSpace::Grid> const& p1Space, 
                bool onlyLowerTriangle=false)  
  {
    auto jacobi = [](auto const& a, int level){ return JacobiPreconditioner<NumaBCRSMatrix<Entry,Index>>(a); };
    auto patchJacobi = [&space](auto const& a, int level){ return PatchDomainDecompositionPreconditioner<FineSpace,Entry::rows>(space,a); };
    return makeAdditivePMultiGrid(A,space,p1Space,patchJacobi,jacobi,onlyLowerTriangle);
  }
}

#endif