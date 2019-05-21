
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

#ifndef MG_MULTIPLICATIVEMULTIGRID_HH
#define MG_MULTIPLICATIVEMULTIGRID_HH

#include <type_traits>

#include <boost/timer/timer.hpp>

#include <dune/istl/preconditioner.hh>

#include "fem/spaces.hh"
#include "linalg/conjugation.hh"
#include "linalg/direct.hh"
#include "linalg/domainDecompositionPreconditioner.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/symmetricOperators.hh"
#include "mg/prolongation.hh"
#include "utilities/memory.hh"
#include "utilities/timing.hh"

namespace Kaskade
{


  // ----------------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief A general multiplicative multigrid preconditioner.
   *
   * This realizes a classical V-cycle with a specified number of pre- and post-smoothings. Two (possibly different)
   * preconditioners are used for (i) smoothing on the levels and (ii) approximately solving the coarse grid system.
   *
   * Note that this is only a preconditioner, no solver. In fact, the computed step is not guaranteed to decrease the energy
   * of the error (depending on the smoother - Gauss-Seidel does, but Jacobi need not).
   *
   * \tparam Entry the type of projected Galerkin matrix entries
   * \tparam Index the index type of Galerkin matrices (usually size_t)
   * \tparam Smoother a smoother type
   * \tparam Prolongation the type of prolongations (conceptually a sparse matrix type)
   *
   * The smoother may rely on the matrix provided on construction remaining available during any
   * method calls besides the smoother destruction. Hence it is perfectly suitable for the smoother
   * just to reference parts of the matrix.
   *
   * A convenient construction of multiplicative multigrid preconditioners is provided by the functions
   * makeMultiplicativeMultiGrid and makeJacobiMultiGrid for h-multigrid in grid hierarchies, and
   * makeMultiplicativePMultiGrid makeJacobiPMultiGrid for p-h-multigrid exploiting both polynomial order
   * and grid hierarchies.
   *
   * Multigrid preconditioners can be put together easily from the three main ingredients multigrid stack,
   * coarse solver, and smoother as follows:
   * \code
   * auto mgStack = makeGeometricMultiGridStack(gridman,duplicate(A),onlyLowerTriangle);
   * auto coarsePreconditioner = makeDirectPreconditioner(std::move(mgStack.coarseGridMatrix()));
   * auto mg = makeMultiplicativeMultiGrid(std::move(mgStack),MakeJacobiSmoother(),moveUnique(std::move(coarsePreconditioner)),nPre,nPost);
   * mg.setSmootherStepSize(0.5);
   * \endcode
   */
  template <class Entry, class Index, class Smoother, class Prolongation>
  class MultiplicativeMultiGrid: public SymmetricPreconditioner<typename Smoother::domain_type,typename Smoother::range_type>
  {
  public:
    using field_type = typename Smoother::field_type;
    using domain_type = typename Smoother::domain_type;
    using range_type = typename Smoother::range_type;
    using CoarsePreconditioner = Dune::Preconditioner<domain_type,range_type>; // symmetric would be better? but elaborate implementation for direct solvers

    /**
     * \brief Default constructor.
     */
    MultiplicativeMultiGrid() = default;

    /**
     * \brief Constructor.
     * \param A the Galerkin matrix to be preconditioned
     * \param Ps the stack of multigrid prolongations
     * \param makeSmoother a callable object with arguments (Matrix const&, int level) that creates a smoother of type Smoother
     *
     * \see makeMultiplicativeMultiGrid
     * \see makeJacobiMultiGrid
     */
    template <class MakeSmoother>
    MultiplicativeMultiGrid(MultiGridStack<Prolongation,Entry,Index>&& mgStack_, MakeSmoother const& makeSmoother,
                            std::unique_ptr<CoarsePreconditioner>&& coarsePreconditioner_, int nPre_=3, int nPost_=3)
    : mgStack(std::move(mgStack_)), nPre(nPre_), nPost(nPost_), coarsePreconditioner(std::move(coarsePreconditioner_)),
      linesearch(false), smootherStepSize(1.0)
    {
      // Create the smoothers for all levels.
      Timings& timer = Timings::instance();
      timer.start("smoother construction");
      for (int l=1; l<mgStack.levels(); ++l)
        smoothers.push_back(makeSmoother(mgStack.a(l)));
      timer.stop("smoother construction");
    }

    MultiplicativeMultiGrid(MultiplicativeMultiGrid&& other) = default;

    MultiplicativeMultiGrid& operator=(MultiplicativeMultiGrid&& other) = default;


    /**
     * \brief Application of preconditioner.
     *
     * Precondition: x = 0
     */
    virtual void apply(domain_type& x, range_type const& r)
    {
      assert(x.two_norm()==0);
      range_type b = r;
      runMG(x,b,mgStack.levels()-1);
    }

    virtual field_type applyDp(domain_type& x, range_type const& r)
    {
      apply(x,r);
      return x*r;
    }

    virtual bool requiresInitializedInput() const
    {
      return true;
    }

    /**
     * \brief Sets the number of pre- and post-smoothing iterations to perform.
     *
     * Note that the sum of pre- and postsmoothings must be positive.
     *
     * \param nPre number of pre-smoothings (nonnegative)
     * \param nPost number of post-smoothings (nonnegative)
     */
    void setSmoothings(int nPre_, int nPost_)
    {
      nPre = nPre_;
      nPost = nPost_;
      assert(nPre>=0 && nPost>=0 && nPre+nPost>0);
    }

    /**
     * \brief Define the smoother step size.
     *
     * The default step length is 1.
     */
    void setSmootherStepSize(double w)
    {
      assert(w>0);
      smootherStepSize = w;
    }

    double getSmootherStepSize() const { return smootherStepSize; }

    /**
     * \brief Enable or disable the line search option.
     *
     * With line search option, at the end of each level in the V-cycle a line search is performed for
     * an optimal scaling of the correction. This may improve the contraction (but need not, for simple
     * Laplace type problems it does not) and guarantee convergence of the multigrid as a fixed point
     * iteration (removing the need for an outer stepsize loop), but incurs one more matrix-vector product
     * (an overhead of up to 30%) and renders the V-cycle a nonlinear scheme, not suited as preconditioner
     * in CG.
     *
     * The default value on construction is off.
     */
    void setLinesearch(bool ls)
    {
      linesearch = ls;
    }

    /**
     * \brief Provides access to the coarse grid preconditioner.
     */
    CoarsePreconditioner& getCoarsePreconditioner()
    {
      return *coarsePreconditioner;
    }

    /**
     * \brief Returns a pair of number of pre- and post-smoothing iterations to perform.
     */
    std::pair<int,int> getSmoothings()
    {
      return std::make_pair(nPre,nPost);
    }

  private:
    // precondition: x==0
    void runMG(domain_type& x, range_type& r, int level)
    {
      if (level==0)
      {
        coarsePreconditioner->pre(x,r);
        coarsePreconditioner->apply(x,r);
        coarsePreconditioner->post(x);
      }
      else
      {
        auto& s = smoothers[level-1];             // Smoother S
        auto const& a = mgStack.a(level);         // Galerkin matrix A
        double w = smootherStepSize;
        domain_type dx(x);                        // temporary vector
        range_type adx(r);

        // pre-smoothing
        bool const sRequiresZeroInput = s.requiresInitializedInput();
        s.pre(x,r);
        for (int i=0; i<nPre; ++i)
        {
          if (sRequiresZeroInput) dx = 0;
          s.apply(dx,r);                          // dx = B^{-1} r
          if (linesearch)
          {
            double dxadx = a.mv(dx,adx);          // adx = A*dx           residual update direction
            w = dx*r / (dxadx);                   // w = dx*r / dx*A*dx   "optimal" step length
            r.axpy(-w,adx);                       // r = r - w*adx        update residual
          }
          else
            a.usmv(-w,dx,r);                      // r = r - w*A*dx       update residual (without intermediate adx access)
          x.axpy(w,dx);                           // x = x + w*dx         update iterate
        }
        s.post(x);

        // coarse grid correction
        auto const& p = mgStack.p(level-1);
        range_type cr; cr.resize(p.M());
        p.mtv(r,cr);                             // cr = P^T r
        domain_type cx; cx.resize(p.M()); cx = 0;
        runMG(cx,cr,level-1);                    // cx ~ (P^T A P)^{-1} cr
        p.mv(cx,dx);                             // dx = P*cx
        x += dx;                                 // x = x+dx
        a.usmv(-1,dx,r);                         // r = r-A*dx


        // post-smoothing
        s.pre(x,r);
        for (int i=0; i<nPost; ++i)
        {
          if (sRequiresZeroInput) dx = 0;
          s.apply(dx,r);                          // dx = B^{-1} r
          if (i+1<nPost || linesearch)            // update residual only if needed lateron
          {
            if (linesearch)
            {
              double dxadx = a.mv(dx,adx);        // adx = A*dx
              w = dx*r / (dxadx);                 // w = dx*r / dx*A*dx
              r.axpy(-w,adx);                     // r = r - w*adx
            }
            else
              a.usmv(-w,dx,r);                    // r = r - w*A*dx
          }
          x.axpy(w,dx);                           // x = x + w*dx
        }
        s.post(x);

        // optional line search
        if (linesearch)
        {
          a.mv(x,dx);                             // dx = Ax
          double omega = (x*r) / (dx*x);          // omega = r^T x / x^T A x
          x *= 1+omega;                           // x = x + omega * x
        }
      }
    }

    MultiGridStack<Prolongation,Entry,Index>  mgStack;
    std::vector<Smoother>                     smoothers;
    int                                       nPre, nPost;
    std::unique_ptr<CoarsePreconditioner>     coarsePreconditioner;
    bool                                      linesearch;
    double                                    smootherStepSize;
  };

  // ----------------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative multigrid preconditioner of V-cycle type for P1 elements.
   *
   * \param mgStack a stack of prolongations and matching projected Galerkin matrices
   * \param makeSmoother a callable object taking a (projected) Galerkin matrix of type NumaBCRSMatrix<Entry,Index>
   *                     and a level, returning a preconditioner
   * \param coarsePreconditioner a symmetric preconditioner to be used for "solving" on the coarsest level
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class Prolongation, class MakeSmoother, class CoarsePreconditioner>
  auto makeMultiplicativeMultiGrid(MultiGridStack<Prolongation,Entry,Index>&& mgStack,
                                   MakeSmoother const& makeSmoother,
                                   std::unique_ptr<CoarsePreconditioner>&& coarsePreconditioner,
                                   int nPre=3, int nPost=3)
  {
    using Smoother = std::result_of_t<MakeSmoother(NumaBCRSMatrix<Entry,Index>)>;
    return MultiplicativeMultiGrid<Entry,Index,Smoother,Prolongation>(
                std::move(mgStack),makeSmoother,std::move(coarsePreconditioner),nPre,nPost);
  }

  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup smoothers
   * \brief Functor for creating Jacobi smoothers.
   */
  class MakeJacobiSmoother
  {
  public:
    template <typename Matrix>
    auto operator()(Matrix const& a) const
    {
      return makeJacobiPreconditioner(a);
    }
  };

  /**
   * \ingroup smoothers
   * \brief Functor for creating overlapping Schwarz smoothers.
   */
  template <typename Space>
  class MakeAdditiveSchwarzSmoother
  {
  public:
    MakeAdditiveSchwarzSmoother(Space const& space_): space(space_) {}

    template <typename Entry, typename Index>
    auto operator()(NumaBCRSMatrix<Entry,Index> const& a) const
    {
      return PatchDomainDecompositionPreconditioner<Space,Entry::rows>(space,a);
    }

  private:
    Space const& space;
  };

  /**
   * \ingroup multigrid
   * \brief Creates a direct solver for the given matrix.
   */
  template <typename Matrix>
  auto makeDirectPreconditioner(Matrix&& A, DirectType directType=DirectType::MUMPS)
  {
    return DirectSolver<typename MatrixTraits<Matrix>::NaturalDomain,typename MatrixTraits<Matrix>::NaturalRange>(A,directType);
  }

  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative Jacobi multigrid for P1 elements.
   *
   * This creates a V-cycle multigrid preconditioner with Jacobi smoother. A direct solver is used for grid level 0.
   *
   * \param A the sparse Galerkin matrix for P1 finite elements on the leaf view of the given simplicial grid
   * \param gridman a grid manager base of a simplicial grid
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class GridMan>
  auto makeJacobiMultiGrid(NumaBCRSMatrix<Entry,Index> const& A, GridMan const& gridman,
                           int nPre=3, int nPost=3, bool onlyLowerTriangle=false)
  {
    Timings& timer = Timings::instance();

    timer.start("MG stack creation");
    auto mgStack = makeGeometricMultiGridStack(gridman,duplicate(A),onlyLowerTriangle);
    timer.stop("MG stack creation");

    timer.start("direct solver creation");
    auto coarsePreconditioner = makeDirectPreconditioner(std::move(mgStack.coarseGridMatrix()));
    timer.stop("direct solver creation");

    timer.start("MG creation");
    auto mg = makeMultiplicativeMultiGrid(std::move(mgStack),MakeJacobiSmoother(),moveUnique(std::move(coarsePreconditioner)),nPre,nPost);
    timer.stop("MG creation");
    mg.setSmootherStepSize(0.5);

    return mg;
  }

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative multigrid preconditioner of V-cycle type for higher order elements.
   *
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param p1Space a linear finite element space
   * \param makeSmoother a callable object taking a (projected) Galerkin matrix of type NumaBCRSMatrix<Entry,Index> and a level, returning a preconditioner
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class Space, class MakeSmoother>
  auto makeMultiplicativePMultiGrid(NumaBCRSMatrix<Entry,Index>&& A, Space const& space,
                                    MakeSmoother const& makeSmoother, int nPre=3, int nPost=3, bool onlyLowerTriangle=false)
  {
    Timings& timer = Timings::instance();
    timer.start("create P multigrid stack");
    auto mgStack = makePMultiGridStack(space,std::move(A),onlyLowerTriangle);
    timer.stop("create P multigrid stack");

    timer.start("create MG as coarse solver");
    auto coarsePreconditioner = makeJacobiMultiGrid(std::move(mgStack.coarseGridMatrix()),space.gridManager(),nPre,nPost,onlyLowerTriangle);
    coarsePreconditioner.setSmootherStepSize(0.6);
    timer.stop("create MG as coarse solver");

    timer.start("create multigrid");
    auto mg = makeMultiplicativeMultiGrid(std::move(mgStack),makeSmoother,moveUnique(std::move(coarsePreconditioner)),nPre,nPost);
    timer.stop("create multigrid");

    return mg;
  }

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative multigrid preconditioner of V-cycle type for higher order elements.
   *
   * This is realized as a multiplicative h-multigrid nested within a two-grid scheme.
   * The outer two-level algorithm uses an overlapping domain decomposition smoother on the patches around grid vertices.
   * Classical h-multigrid with point Jacobi smoother is used as a coarse level preconditioner/solver.
   *
   * The default smoother step sizes are 0.5 for the domain decomposition smoother and 0.8 for the Jacobi smoother. In a
   * preliminary test (Poisson equation on the unit square, 2016-01) these values turned out to be quite good.
   *
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class Space, class CoarseSpace>
  auto makeBlockJacobiPMultiGrid(NumaBCRSMatrix<Entry,Index> A, Space const& space,
                                 NumaBCRSMatrix<Entry,Index> coarseA, CoarseSpace const& coarseSpace,
                                 int nPre=3, int nPost=3)
  {
    Timings& timer = Timings::instance();
    timer.start("create P multigrid stack");
    auto mgStack = makePMultiGridStack(space,std::move(A),coarseSpace,std::move(coarseA));
    timer.stop("create P multigrid stack");

    timer.start("create MG as coarse solver");
    auto coarsePreconditioner = makeJacobiMultiGrid(std::move(mgStack.coarseGridMatrix()),space.gridManager(),nPre,nPost);
    coarsePreconditioner.setSmootherStepSize(0.6);
    timer.stop("create MG as coarse solver");

    timer.start("create multigrid");
    auto mg = makeMultiplicativeMultiGrid(std::move(mgStack),MakeAdditiveSchwarzSmoother<Space>(space),
                                          moveUnique(std::move(coarsePreconditioner)),nPre,nPost);
    mg.setSmootherStepSize(0.5);
    timer.stop("create multigrid");

    return mg;
  }

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative multigrid preconditioner of V-cycle type for higher order elements.
   *
   * This is realized as a multiplicative h-multigrid nested within a two-grid scheme.
   * The outer two-level algorithm uses an overlapping domain decomposition smoother on the patches around grid vertices.
   * Classical h-multigrid with point Jacobi smoother is used as a coarse level preconditioner/solver.
   *
   * The default smoother step sizes are 0.5 for the domain decomposition smoother and 0.8 for the Jacobi smoother. In a
   * preliminary test (Poisson equation on the unit square, 2016-01) these values turned out to be quite good.
   *
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class Space>
  auto makeBlockJacobiPMultiGrid(NumaBCRSMatrix<Entry,Index> A, Space const& space,
                                 int nPre=3, int nPost=3, bool onlyLowerTriangle=false)
  {
    auto mg = makeMultiplicativePMultiGrid(std::move(A),space,MakeAdditiveSchwarzSmoother<Space>(space),
                                           nPre,nPost,onlyLowerTriangle);
    mg.setSmootherStepSize(0.5);
    return mg;
  }

  /**
   * \ingroup multigrid
   * \brief Convenience function creating multiplicative multigrid preconditioner of V-cycle type with Jacobi smoother for higher order elements.
   *
   * \param A the sparse Galerkin matrix for higher order finite elements on the leaf view of the given simplicial grid
   * \param space a higher order finite element space
   * \param nPre the number of pre-smoothings
   * \param nPost the number of post-smoothings
   * \param onlyLowerTriangle if true, A is assumed to be symmetric and only its lower triangular part is accessed
   *
   * \relates MultiplicativeMultiGrid
   */
  template <class Entry, class Index, class FineSpace>
  auto makeJacobiPMultiGrid(NumaBCRSMatrix<Entry,Index> A, FineSpace const& space,
                            int nPre=3, int nPost=3, bool onlyLowerTriangle=false)
  {
    auto mg = makeMultiplicativePMultiGrid(std::move(A),space,MakeJacobiSmoother(),nPre,nPost,onlyLowerTriangle);
    mg.setSmootherStepSize(0.5);
    return mg;
  }




  // ---------------------------------------------------------------------------------------------------------


}

#endif
