

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SDC_HH
#define SDC_HH

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/istl/operators.hh"
#include "dune/istl/solvers.hh"

#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

#include "linalg/dynamicMatrix.hh"


namespace Kaskade {
  
  /**
   * \ingroup timestepping
   * \brief Abstract base class of time grids for (block) spectral defect correction methods.
   * 
   * This class represents a time grid on \f$ [t_0,t_n]\f$ with \f$ n \f$ subintervals and \f$ n+1 \f$
   * time grid points (including the end points).
   */
  class SDCTimeGrid
  {
  public:
    /**
     * \brief The type used for real vectors.
     */
    typedef Dune::DynamicVector<double> RealVector;
    
    /**
     * \brief The type used for real (dense) matrices.
     */
    typedef Dune::DynamicMatrix<double> RealMatrix;
    
    /**
     * \brief Time points in the time step
     * 
     * The time step \f$ [t_0, t_n] \f$ contains \f$ n+1 \f$ time points \f$ t_i \f$,
     * including the end points. Those are provided here. The time points are stored in increasing
     * order.
     */
    virtual RealVector const& points() const = 0;
    
    /**
     * \brief Integration matrix \f$ S \f$
     * 
     * On the time interval, the time grid defines an interpolation scheme, such that given \f$ u(t_i) \f$
     * for all time points, a continuous function \f$ u(t) \f$ can be evaluated in the whole interval.
     * 
     * The Lagrangian interpolation functions \f$ L_k \f$ are defined by \f$ L_(t_i) = \delta_{ik} \f$.
     * The matrix \f$ S \in \mathbb{R}^{n\times n+1}\f$ contains the values 
     * \f[ S_{ik} = \int_{\tau=t_i}^{t_{i+1}} L_k(\tau) \, d\tau. \f]
     * This way, if \f$ u \f$ is defined in terms of its function values \f$ v_i = u(t_i) \f$, the 
     * integrals can be evaluated by a matrix-vector multiplication:
     * \f[ \int_{\tau=t_i}^{t_{i+1}} u(\tau) \, d\tau = (Sv)_i \f]
     * 
     * Note that not necessarily all time points \f$ t_i \f$ are used in formulating the Lagrangian interpolation
     * functions. E.g., on Radau points, the first point \f$ t_0 \f$ is omitted (leading to a zero column in \f$ S \f$).
     */
    virtual RealMatrix const& integrationMatrix() const = 0;
    
    /**
     * \brief Differentiation matrix \f$ D \f$
     *
     * On the time interval, the time grid defines an interpolation scheme, such that given \f$ u(t_i) \f$
     * for all time points, a continuous function \f$ u(t) \f$ can be evaluated in the whole interval.
     * 
     * The Lagrangian interpolation functions \f$ L_k \f$ are defined by \f$ L_k(t_i) = \delta_{ik} \f$.
     * The matrix \f$ D \in \mathbb{R}^{n+1\times n+1}\f$ contains the values 
     * \f[ D_{ik} = \dot L_k(\tau_i)  \f]
     * This way, if \f$ u \f$ is defined in terms of its function values \f$ v_i = u(t_i) \f$, its 
     * derivatives can be evaluated by a matrix-vector multiplication:
     * \f[ \dot u(\tau_i) = (Dv)_i \f]
     * 
     * Note that not necessarily all time points \f$ t_i \f$ are used in formulating the Lagrangian interpolation
     * functions. E.g., on Radau points, the first point \f$ t_0 \f$ is omitted (leading to a zero column in \f$ D \f$).
     */
    virtual RealMatrix const& differentiationMatrix() const = 0;
    
    /**
     * \brief Perform refinement of the grid, filling the prolongation matrix.
     * 
     * If the function representation is not sufficiently accurate, a finer grid of time points can be tried.
     * This method refines the grid to \f$ m+1 \f$ time points \f$ s_i \f$, \f$ m>n \f$, 
     * and fills a prolongation matrix \f$ P\in \mathbb{R}^{m+1\times n+1} \f$ such that with \f$ v_i = u(t_i) \f$
     * and \f$ w=Pv \f$ it holds that \f$ w_i = u(s_i) \f$.
     * 
     * Different derived classes may implement this in different ways, or provide their own, more flexible
     * ways of refining the grid.
     * 
     * \param[out] p the prolongation matrix
     */
    virtual void refine(RealMatrix& p) = 0;
    
    /**
     * \brief Compute interpolation coefficients.
     * 
     * Returns a matrix \f$ w \in \mathbb{R}^{m+1\times n+1} \f$, such that the interpolation polynomial \f$ p \f$ to the values \f$ y_i \f$ at 
     * grid points \f$ t_i \f$ can be evaluated as 
     * \f[ p(x_i) = \sum_{j=0}^n w_{ij} y_j, \quad i=0,\dots,m. \f]
     */
    virtual RealMatrix interpolate(RealVector const& x) const = 0;
  };
  
  // ---------------------------------------------------------------------------------------------------
  
  /** 
   * \ingroup timestepping
   * \brief Triangular approximate integration matrix \f$ \hat S \f$ corresponding to Euler integrator
   * 
   * This computes a lower "triangular" matrix \f$ \hat{S} \in \mathbb{R}^{n\times n+1} \f$ with \f$ \hat{S}_{ij} = 0 \f$ 
   * for \f$ i \le j \f$ that can be used to formulate SDC sweeps.
   * 
   * This is \f$ \hat{S}_{i,i+1} = t_{i+1}-t_i \f$, all other entries zero. In an SDC sweep,
   * this corresponds to an implicit Euler scheme for integrating the defect equation.
   */
  void eulerIntegrationMatrix(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix& Shat);
  
  /** 
   * \ingroup timestepping
   * \brief Triangular approximate integration matrix \f$ \hat S \f$ corresponding to specialized multi-step integrators
   * 
   * This computes a lower "triangular" matrix \f$ \hat{S} \in \mathbb{R}^{n\times n+1} \f$ with \f$ \hat{S}_{ij} = 0 \f$ 
   * for \f$ i \le j \f$ that can be used to formulate SDC sweeps.
   * 
   * If \f$ S^T = LU \f$, this is \f$ \hat S = U^T \f$, which, in an SDC sweep
   * corresponds to a cascade of tailor-made multi-step integrators. The consequence is that on the Dahlquist
   * test equation \f$ \dot u = \lambda u \f$, in the limit \f$ \lambda \to \infty \f$ the SDC iteration matrix 
   * \f$ G = I - \hat S^{-1} S \f$ is nilpotent and hence has a spectral radius of 0, which yields fast convergence.
   */
  void luIntegrationMatrix(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix& Shat);
  
  /**
   * \ingroup timestepping
   * \brief Types of weight functions for optimized integration matrices.
   */
  enum IntegrationMatrixOptimizationWeight { OW_FLAT };
  
  /**
   * \ingroup timestepping
   * \brief Triangular approximate differentiation and integration matrix \f$ \hat D, \hat S \f$ corresponding to specialized multi-step integrators
   * 
   * \param grid the time grid
   * \param k    the sweep number
   * \param w    the type of weight function for the optimization
   * 
   * This returns one of a set of precomputed, optimized integration matrices, if available for the given grid. If no precomputed
   * matrix is available, a lookup exception is thrown.
   */
  std::pair<SDCTimeGrid::RealMatrix,SDCTimeGrid::RealMatrix> 
  optimizedMatrices(SDCTimeGrid const& grid, int k, IntegrationMatrixOptimizationWeight w=OW_FLAT);

  /**
   * \ingroup timestepping
   * \brief Triangular approximate differentiation and integration matrix \f$ \hat D, \hat S \f$ corresponding to specialized multi-step integrators
   * 
   * \param grid     the time grid
   * \param k        the sweep number
   * \param Df       the differentiation matrix returned in case no precomputed optimized matrix is available
   * \param Sf       the integration matrix returned in case no precomputed optimized matrix is available
   * \param w        the type of weight function for the optimization
   * 
   * This returns one of a set of precomputed, optimized integration matrices, if available for the given grid. If no precomputed
   * matrix is available, the provided fallback matrix is returned. This is a convenience overload.
   */
  std::pair<SDCTimeGrid::RealMatrix,SDCTimeGrid::RealMatrix> 
  optimizedIntegrationMatrix(SDCTimeGrid const& grid, int k, SDCTimeGrid::RealMatrix const& Df, SDCTimeGrid::RealMatrix const& Sf, 
                             IntegrationMatrixOptimizationWeight w=OW_FLAT);

  
  
  /**
   * \ingroup timestepping
   * \brief spectral time grid for defect correction methods with Lobatto points
   * 
   * Note that collocation with Lobatto points is A-stable but not L-stable, and therefore SDC methods
   * on Lobatto grids may be a suboptimal choice for highly stiff problems (e.g., parabolic equations or Dirichlet 
   * b.c. realized as quadratic penalty). Consider using RadauTimeGrid in these cases.
   */
  class LobattoTimeGrid: public SDCTimeGrid
  {
  public:
    /**
     * \brief constructs a Lobatto grid with \f$ n+1 \f$ points on \f$ [a,b] \f$
     * 
     * This may throw LinearAlgebraException.
     */
    LobattoTimeGrid(int n, double a, double b);
    
    virtual RealVector const& points() const 
    {
      return pts;
    }
    
    virtual RealMatrix const& integrationMatrix() const
    {
      return integ;
    }
    
    virtual RealMatrix const& differentiationMatrix() const
    {
      return diff;
    }
    
    /**
     * \brief perform refinement of the grid, filling the prolongation matrix
     * 
     * This may throw LinearAlgebraException.
     */
    virtual void refine(RealMatrix& p);
    
    virtual RealMatrix interpolate(RealVector const& x) const;

  private:
    RealVector pts;
    RealMatrix integ;
    RealMatrix diff;
  };
  
  /**
   * \ingroup timestepping
   * \brief spectral time grid for defect correction methods with Radau points
   * 
   * Note that collocation with Radau points is L-stable, and therefore SDC methods
   * on Radau grids are best for highly stiff problems (e.g., parabolic equations or Dirichlet 
   * b.c. realized as quadratic penalty). Consider using LobattoTimeGrid in other cases.
   */
  class RadauTimeGrid: public SDCTimeGrid
  {
  public:
    /**
     * \brief constructs a Radau grid with \f$ n+1 \f$ points on \f$ [a,b] \f$
     * 
     * This may throw LinearAlgebraException.
     */
    RadauTimeGrid(int n, double a, double b);
    
    virtual RealVector const& points() const 
    {
      return pts;
    }
    
    virtual RealMatrix const& integrationMatrix() const
    {
      return integ;
    }
    
    virtual RealMatrix const& differentiationMatrix() const
    {
      return diff;
    }
    
    /**
     * \brief perform refinement of the grid, filling the prolongation matrix
     * 
     * This may throw LinearAlgebraException.
     */
    virtual void refine(RealMatrix& p);
    
    virtual RealMatrix interpolate(RealVector const& x) const;

  private:
    RealVector pts;
    RealMatrix integ;
    RealMatrix diff;
    
    void computeMatrices();
  };
  
  
  
  /**
   * \ingroup timestepping
   * \brief A single spectral defect correction iteration sweep
   * 
   * This function performs one spectral defect correction (SDC) iteration for the abstract reaction
   * diffusion equation
   * \f[ M \dot u = A u + M f(u) \f]
   * using the linearly implicit Euler method on the given time grid. The solution interpolation
   * vector \arg u contains the approximate values of \f$ u \f$ at the time nodes \f$ t_i \f$ in a
   * Lagrangian FE basis.
   * Here, \arg A and \arg M are fixed matrices, such that locally we perform a method of lines
   * in this time step.
   * 
   * SDC iterations are inexact Newton iterations for the Fredholm collocation time discretization
   * of an ODE on a time grid \f$ t_0, \dots, t_n \f$:
   * \f[ M (u_{i+1} - u_i) = \int_{\tau=t_i}^{t_{i+1}} p(\tau) \,d\tau \quad i=0,\dots,n-1, \f]
   * where \f$ p(t_i) = r(u_i) = Au_i + M f(u_i) \f$ represents the (probably polynomial) interpolant
   * of the right hand side. By linearity, the integral can be expressed as a linear combination of
   * the right hand side values:
   * \f[ M (u_{i+1} - u_i) = \sum_{j=0}^n S_{ij} r(u_j) \f]
   * Applying Newton's method to \f$ F(u) = M (u_{i+1} - u_i) - \sum_{j=0}^n S_{ij} r(u_j) \f$ yields
   * \f$ F'(u^k)\delta u^k = -F(u^k), \quad \delta u_i^k = u_i^{k+1} - u_i^k \f$, or, more elaborate,
   * \f[ M(\delta u_{i+1}^k - \delta u_i^k) - \sum_{j=0}^n S_{ij} r'(u_j^k) \delta u_j^k
   *     = -M (u_{i+1}^k - u_i^k) + \sum_{j=0}^n S_{ij} r(u_j^k) =: R_i, \quad i=0,\dots,n-1. \f]
   * The time-global coupling makes this system difficult to solve. An approximation of the
   * quadrature coefficients \f$ S_{ij} \f$ by a triangular \f$ \hat S_{ij} \f$  yields the simpler system
   * \f[ M(\delta u_{i+1}^k - \delta u_i^k) - \hat{S}_{i,i+1} r'(u_{i+1}^k) \delta u_{i+1}^k
   *     =  R_i + \sum_{j=0}^i \hat{S}_{i,j} r'(u_{j}^k) \delta u_j^k, \quad i=0,\dots,n-1, \f]
   * starting with \f$ \delta u_0^k = 0 \f$, or, equivalently,
   * \f[ (M- \hat{S}_{i,i+1} r'(u_{i+1}^k))(\delta u_{i+1}^k - \delta u_i^k)  
   *     =  \sum_{j=0}^{i}\hat{S}_{i,j}r'(u_{j}^k) \delta u_j^k +  \hat{S}_{i,i+1}r'(u_{i+1}^k) \delta u_i^k + R_i, \quad i=0,\dots,n-1, \f]
   * 
   * The norm of the correction \f$ \delta u^k \f$ is returned, i.e.
   * \f[ \left( (t_n-t_0)^{-1}\sum_{i=0}^{n-1} (t_{i+1}-t_i)(\delta u_i^k)^T (M- \hat{S}_{i,i+1} r'(u_{i+1}^k)) \delta u_i^k \right)^{1/2}. \f]
   * 
   * \tparam Matrix a sparse matrix type, usually Dune::BCRSMatrix or a NumaBCRSMatrix
   * \tparam Vectors a container type with elements from the domain type of Matrix
   * \tparam Solver
   * \tparam ReactionDerivatives
   * 
   * \param[in] grid the collocation time grid 
   * \param[in] Shat the triangular quadrature matrix approximation \f$ \hat S \f$
   * \param[in] solve a callable that supports solve(A,x,b) giving an approximative solution of \f$ Ax = b \f$, where A is of type Matrix
   * \param[in] M the mass matrix 
   * \param[in] Stiff the stiffness matrix \f$ A \f$ (usually negative semidefinite)
   * \param[in] rUi the right hand sides at the collocation times
   * \param[in] rDu the reaction term derivatives (a collection of vectors representing the diagonals of the derivatives)
   * \param[in] u the current iterate (values of \f$ u \f$ at the collocation times)
   * \param[out] du the approximate Newton correction \f$ \delta u \f$
   */
  template <class Matrix, class Vectors, class ReactionDerivatives, class Solver>
  typename Matrix::field_type sdcIterationStep(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix const& Shat, Solver const& solve, 
                                               Matrix const& M, Matrix const& Stiff,
                                               Vectors const& rUi, ReactionDerivatives const& rDu, Vectors const& u, Vectors& du)
  {
    auto const& pts = grid.points();
    int const n = pts.size()-1;      // number of subintervals

    assert(u.size()==n+1); // including start point
    typedef typename Vectors::value_type Vector;
    
    Vector tmp = u[0];
    std::vector<Vector> Mdu(n,tmp);
    for (int i=0; i<n; ++i)
    {
      // compute M (u_{i}-u_{i+1})  TODO: loop fusion
      tmp = u[i];
      tmp.axpy(-1.0,u[i+1]);
      M.mv(tmp,Mdu[i]);
    }
    return sdcIterationStep2(grid,Shat,solve,M,Stiff,rUi,rDu,Mdu,du);
  } 
  
 /**
   * \ingroup timestepping
   * \brief A single spectral defect correction iteration sweep
   * 
   * This is an overload with slightly different interface (required for some more or less arcane algorithmic variants). The only interface difference 
   * is that instead of values \f$ u_i \f$, the actually required differences \f$ M (u_{i+1}-u_i) \f$ are provided. Those are computed explicitly in
   * \ref sdcIterationStep.
   * 
   * \tparam Matrix a sparse matrix type, usually Dune::BCRSMatrix or a NumaBCRSMatrix
   * \tparam Vectors a container type with elements from the domain type of Matrix
   * \tparam ReactionDerivatives a container type with sparse matrix elements (usually BCRSMatrix or NumaBCRSMatrix)
   * \tparam Solver
   * 
   * \param[in] grid the collocation time grid with n+1 points
   * \param[in] Shat the triangular quadrature matrix approximation \f$ \hat S \f$
   * \param[in] solve a callable that supports solve(A,x,b) giving an approximative solution of \f$ Ax = b \f$, where A is of type Matrix
   * \param[in] M the mass matrix 
   * \param[in] Stiff the stiffness matrix (usually negative semidefinite)
   * \param[in] rUi the right hand sides at the collocation times (including interval start point), size n+1
   * \param[in] rDu the reaction term derivatives (a collection of sparse matrices, the sparsity pattern of which is a subset of that of A and M), size n+1
   * \param[in] Mdu the current iterate differences (values of \f$ M(u_{i}-u_{i+1}) \f$ for \f$ i=0,\dots,n-1 \f$)
   * \param[out] du the approximate Newton correction \f$ \delta u \f$
   * 
   * \return the energy norm of the correction du
   */
  template <class Matrix, class Vectors, class ReactionDerivatives, class Solver>
  typename Matrix::field_type sdcIterationStep2(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix const& Shat, Solver const& solve, 
                                                Matrix const& M, Matrix const& Stiff,
                                                Vectors const& rUi, ReactionDerivatives const& rDu, Vectors const& Mdu, Vectors& du)
  {
    auto const& pts = grid.points();
    int const n = pts.size()-1;      // number of subintervals

    assert(Mdu.size()>=n); 
    assert(du.size()>=n+1);  // including start point
    
    size_t const dofs = Mdu[0].size();

    // compute exact integration matrix
    auto const& S = grid.integrationMatrix();


    // initialize correction at starting point to zero
    du[0] = 0.0; 

    typedef typename Vectors::value_type Vector;
    Vector rhs(dofs), tmp(dofs);  // declare here to prevent frequent reallocation
    
    // perform n Euler steps
    typename Matrix::field_type norm = 0;
    Matrix J = M;
    for (int i=1; i<=n; i++)
    {
      // matrix J = M - Shat_i-1,i*(A+f_u)
      for (size_t row=0; row<J.N(); ++row)
      {
        auto colJ = J[row].begin(); 
        auto end = J[row].end();
        auto colM = M[row].begin();
        auto colA = Stiff[row].begin();
        auto colR = rDu[i][row].begin();
        auto endR = rDu[i][row].end();
        
        while (colJ != end)
        {
          *colJ = *colM - Shat[i-1][i] * *colA;
          
          if (colR != endR && colJ.index() == colR.index()) // f_u can have subset of sparsity pattern
          {
            *colJ -= std::min(0.5* *colM, Shat[i-1][i] * *colR); // guarantee M - Shat*fu is nonnegative -- reduce by at most 50%
            ++colR;
          }
          ++colJ; ++colM; ++colA;           
        }
      }
      

      //  right-hand side for linear system

      // M * ( u_i^{k} - u_{i+1}^k + du_i)
      rhs = Mdu[i-1];
      M.umv(du[i-1],rhs);

      // add sum_j S_ij r_j to right hand side
      for (int j=0; j<=n; ++j)
        rhs.axpy(S[i-1][j],rUi[j]);
         
      // add sum_j Shat_ij r'(u_j) du_j with r' = A + f_u
      tmp = 0;
      for (int j=0; j<i; ++j) // TODO: start at 1 instead of 0? du[0] is zero anyway...
      {
        rDu[j].usmv(Shat[i-1][j],du[j],rhs);
        tmp.axpy(Shat[i-1][j],du[j]);
      }
      Stiff.umv(tmp,rhs);
        
      // solve linear system
      du[i] = du[i-1]; // previous increment is probably a good starting value
      solve(J,du[i],rhs);
      
      // evaluate norm of correction
      norm += (pts[i]-pts[i-1]) * (du[i]*rhs);
    }   // end i - loop
    
    return std::sqrt(norm/(pts[n]-pts[0]));
  } 
  
  template <class Matrix, class Vectors, class Reaction, class Solver>
  typename Matrix::field_type sdcIterationStep3(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix const& Shat, Solver const& solve, 
                                                Matrix const& M, Matrix const& Stiff,
                                                Vectors const& rUi, Reaction const& fAt, Vectors const& Mdu, Vectors& du,
                                                int nReactionSweeps)
  {
    auto const& pts = grid.points();
    int const n = pts.size()-1;      // number of subintervals

    assert(Mdu.size()>=n); 
    assert(du.size()>=n+1);  // including start point
    
    size_t const dofs = Mdu[0].size();

    // compute exact integration matrix
    auto const& S = grid.integrationMatrix();
// auto const& S = Shat;

    // initialize correction at starting point to zero
    du[0] = 0.0; 

    typedef typename Vectors::value_type Vector;
    Vector rhs(dofs), tmp(dofs);  // declare here to prevent frequent reallocation
    
    
    // perform n basic steps
    typename Matrix::field_type norm = 0;
    Matrix J = M;
    for (int i=1; i<=n; i++)
    {
      // Each basic step consists of a splitting method, separating a linearly implicit Euler step
      // from the pointwise nonlinearity of the right hand side. First we perform the linearly implicit
      // Euler step.
      
      //  right-hand side for linear system

      // M * ( u_{i-1}^{k} - u_{i}^k + du_{i-1})
      rhs = Mdu[i-1];
      M.umv(du[i-1],rhs);

      // add sum_j S_ij r_j to right hand side
      for (int j=0; j<=n; ++j)
        rhs.axpy(S[i-1][j],rUi[j]);
         
      // add sum_j Shat_ij r'(u_j) du_j with r' = A + f_u, i.e. A (sum_j Shat_ij du_j) + sum_j Shat_ij f_uj du_j
      // addition of f_u is postponed to matrix loop
      tmp = 0;
      for (int j=0; j<i; ++j) // TODO: start at 1 instead of 0? du[0] is zero anyway...
        tmp.axpy(Shat[i-1][j],du[j]);
      Stiff.umv(tmp,rhs);
        
      // matrix J = M - Shat_i-1,i*(A+f_u)
      auto const f = fAt( pts[i] );
      for (size_t row=0; row<J.N(); ++row)
      {
        auto colJ = J[row].begin(); 
        auto end = J[row].end();
        auto colM = M[row].begin();
        auto colA = Stiff[row].begin();
        
        while (colJ != end)
        {
          auto fu = *colM * (f(row,0.0,1) + f(colJ.index(),0.0,1)) / 2;
//           *colJ = *colM - Shat[i-1][i] * (*colA + fu); // exact Newton, but for large step sizes the Jacobi matrix becomes singular...
//           *colJ = *colM - Shat[i-1][i] * (*colA + std::min(0.0,fu));  // this modification guarantees an invertible Jacobian
          *colJ = *colM - Shat[i-1][i] * (*colA);         // explicit reaction
          if (colJ.index()==row)
            *colJ -= Shat[i-1][i] * fu; // only diagonal of reaction mass matrix

//           // add f_u contribution to right hand side (no-op for Euler SDC as Shat[i-1][j]=0 for j<i)
//           for (int j=0; j<i; ++j) // TODO: start at 1 instead of 0? du[0] is zero anyway... 
//           {
//             auto f = fAt(pts[j]);
//             rhs[row] += *colM*Shat[i-1][j]*(f(row,0.0,1)+f(colJ.index(),0.0,1))/2 * du[j][colJ.index()];
//           }
          
          ++colJ; ++colM; ++colA;
        }
      }
      
      // solve linear system
      du[i] = du[i-1]; // previous increment is probably a good starting value
      solve(J,du[i],rhs);

std::cerr << "du=" << du[i] << "\n";      
// std::cerr << "|du| = " << du[i].two_norm2() << "  du[i]=" << du[i][0] << "  du[i-1]=" << du[i-1][0] << "\n";      
      // Second part is the remaining nonlinearity, which we address by a subsequent run of a couple of Euler steps
      double const tau = pts[i]-pts[i-1];
      double snorm = 0;
      for (size_t row=0; row<J.N(); ++row)
      {
        auto fend = fAt(pts[i]);
        double dfdu = std::min(0.0,fend(row,0.0,1)) * du[i][row];
        
        double s = 0;
        int l = nReactionSweeps;
        for (int j=0; j<l; ++j)
        {
          double theta = 0.0; // where in the subinterval to evaluate (theta=0 explicit Euler, theta=0.5 lin. impl. midpoint, theta=1 lin. impl Euler)
          auto f = fAt(pts[i-1]+(j+theta)*tau/l);
          double dut = ((l-j-theta)*du[i-1][row]+(j+theta)*du[i][row])/l;
          
// if (row==0)          
// std::cerr << "j=" << j << " du+s=" << dut+s << "   f(u+du+s)=" <<  f(row,dut+s,0) << "  f(u)=" <<  f(row,0.0,0) << "  dfdu=" <<  dfdu  <<  "  sum rhs=" << (f(row,dut+s,0) - f(row,0.0,0) - dfdu);
//           s += tau/l * (f(row,dut+s,0) - f(row,0.0,0) - dfdu) / (1-1.0*tau/l*std::min(0.0,f(row,dut+s,1)));
if (row==0)          
std::cerr << "j=" << j << " du+s=" << dut+s << "   f(u+du+s)=" <<  f(row,dut+s,0) << "  f(u)=" <<  f(row,0.0,0) << "  dfdu=" <<  f(row,dut+s,1)  <<  "  sum rhs=" << (f(row,dut+s,0) - f(row,0.0,0));
          s += tau/l * (f(row,dut+s,0) - f(row,0.0,0) - f(row,0.0,1)*(dut+s)) / (1-theta*tau/l*std::min(0.0,f(row,dut+s,1)));
if (row==0) std::cerr << "  -> s=" << s << "\n";          
        }
        du[i][row] += s;
        snorm += s*s;
      }
// std::cerr << "|s| = " << snorm << "\n";      
      
      // evaluate norm of correction
      norm += (pts[i]-pts[i-1]) * (du[i]*rhs);
    }   // end i - loop
    
    return std::sqrt(norm/(pts[n]-pts[0]));
  } 
  
  /**
   * \ingroup timestepping
   * \brief A single waveform relaxation iteration
   * 
   * This function performs one linearized waveform relaxation iteration for 
   * \f[ M \dot u = A u + M f(u) \f]
   * using the Jacobi method. It approximately solves the linear defect system
   * \f[ M \dot{\delta u} = A \delta u + M f'(u)\delta u + r, \quad r=Au+f(u)-\dot u \f]
   * for the correction \f$ \delta u \f$ to be added to the provided approximation \f$ u \f$. The time evolution
   * is discretized by a collocation method on the given time grid. Thus, for each entry \f$ u_i \f$, the collocation
   * solution of the scalar ODE 
   * \f[ M_{ii}\dot{\delta u_i} = A_{ii} \delta u_i + M_{ii}f'(u_i) \delta u_i + r_i. \f]
   *
   * The solution interpolation
   * vector \arg u contains the approximate values of \f$ u \f$ at the time nodes \f$ t_i \f$ in a
   * Lagrangian FE basis.
   * Here, \arg A and \arg M are fixed matrices, such that locally we perform a method of lines
   * in this time step.
   * 
   * Linearized waveform relaxation works here as
   * 
   * 
   * 
   * The norm of the correction \f$ \delta u^k \f$ is returned, i.e.
   * \f[  \f]
   * 
   * \tparam Matrix a sparse matrix type, usually Dune::BCRSMatrix or a NumaBCRSMatrix
   * \tparam Vectors a container type with elements from the domain type of Matrix
   * \tparam ReactionDerivatives
   * 
   * \param[in] grid the collocation time grid 
   * \param[in] Shat the triangular quadrature matrix approximation \f$ \hat S \f$
   * \param[in] solve a callable that supports solve(A,x,b) giving an approximative solution of \f$ Ax = b \f$, where A is of type Matrix
   * \param[in] M the mass matrix 
   * \param[in] A the stiffness matrix (usually negative semidefinite)
   * \param[in] rUi the right hand sides at the collocation times
   * \param[in] rDu the reaction term derivatives (a container of sparse matrices, the sparsity pattern of which is a subset of that of A and M), size n+1
   * \param[in] Mdu the current iterate differences (values of \f$ M(u_{i}-u_{i+1}) \f$ for \f$ i=0,\dots,n-1 \f$)
   * \param[out] du the approximate Newton correction \f$ \delta u \f$
   * 
   * \return the energy norm of the correction du
   */
  template <class Matrix, class Vectors, class ReactionDerivatives>
  double waveformRelaxationStep2(SDCTimeGrid const& grid, 
                                 Matrix const& M, Matrix const& A,
                                 Vectors const& rUi, ReactionDerivatives const& rDu, Vectors& du)
  {
    auto const& pts = grid.points();
    int const n = pts.size()-1;      // number of subintervals -- left interval boundary always included in points
    

    assert(rUi.size()>=n); 
    assert(du.size()>=n+1);  // including start point - usually 0
    
    size_t const dofs = rUi[0].size();
    
    // We need to solve the scalar ODE
    // Mii y_t = Aii y + Rii y + ri, y(0) = 0
    // We approximate y_t by the spectral differentiation matrix D as Dy and obtain, as system in the 
    // collocation points (excluding the initial value 0)
    // (MiiD - diag(Aii+Rii)) y = r. 
    // The system matrix is called S here.
    
    Dune::DynamicMatrix<double> S(n,n);  // n x n system matrix TODO: use Kaskade::DynamicMatrix?
    Dune::DynamicVector<double> r(n), y(n);
    auto const& D = grid.differentiationMatrix();
    
    std::cerr << "using D = \n" << D << "\n";
    
    // a function to extract diagonal entries 
    auto diagonal = [](Matrix const& M, size_t i) -> double
    {
      auto const& arow = M[i];
      auto first = arow.begin();
      auto last = arow.end();
      while (first!=last && first.index()<i)
        ++first;
      if (first==last)
        std::cerr << "first==last\n"; std::cerr.flush();
      return *first;
    };
    
    double retval = 0;
    
    for (size_t i=0; i<dofs; ++i)
    {
      std::cerr.flush();
      auto mii = diagonal(M,i);
      auto aii = diagonal(A,i);
      
      for (int j=0; j<n; ++j)
      {
        for (int k=0; k<n; ++k)
          S[j][k] = mii*D[j+1][k+1]; // remember D starts with left interval point
        S[j][j] -= aii+diagonal(rDu[j+1],i);
        r[j] = rUi[j+1][i];
      }
      S.solve(y,r);
      for (int j=0; j<n; ++j)
        du[j+1][i] = 0.5*y[j]; // Jacobi damping (otherwise we get spatially high-frequent errors)
      retval += 0.25 * y.two_norm2();
    }
    
    return std::sqrt(retval);
  }

}


#endif
