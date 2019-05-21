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

/**
 * \file
 */

#ifndef MULTI_GRID_SOLVER_HH_
#define MULTI_GRID_SOLVER_HH_

#include <memory> // std::unique_ptr
#include <vector>
#include <cmath>
#include <fstream>

#include <iostream> // std::cout, std::endl

#include "dune/grid/common/grid.hh"
#include "dune/istl/bcrsmatrix.hh"
#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"

#include "fem/fixdune.hh"
#include "linalg/conjugation.hh"
#include "linalg/triplet.hh"
#include "linalg/jacobiSolver.hh"

namespace Kaskade
{
  namespace MultiGridSolver_Detail
  {
    // convenient template alias
    template <class Scalar, int n, class Allocator> using BlockVector = Dune::BlockVector<Dune::FieldVector<Scalar,n>,Allocator>;

    /// Compute \f$ y = \alpha Px+y \f$. If resetSolution=true computes \f$ y = \alpha Px \f$.
    template <class Scalar, int n, bool resetSolution=false>
    void axpy(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> > const& x,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> >& y,
        Scalar alpha = 1.0)
    {
      if(resetSolution) y *= 0.0;
      assert(P.M()==x.N());
      assert(P.N()==y.N());

      auto riter = P.begin(), rend = P.end();
      for(;riter!=rend; ++riter)
      {
        auto citer = riter->begin(), cend = riter->end();
        for(;citer!=cend; ++citer) y[riter.index()] += alpha * (*citer)[0][0] * x[citer.index()];
      }
    }

    /// Compute \f$ y = \alpha P^T x+y \f$. If resetSolution=true computes \f$ y = \alpha P^T x \f$.
    template <class Scalar, int n, bool resetSolution=false>
    void atxpy(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> > const& x,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> >& y,
        Scalar alpha = 1.0)
    {
      if(resetSolution) y *= 0.0;
      assert(P.M()==y.N());
      assert(P.N()==x.N());

      auto riter = P.begin(), rend = P.end();
      for(;riter!=rend; ++riter)
      {
        auto citer = riter->begin(), cend = riter->end();
        for(;citer!=cend; ++citer) y[citer.index()] += alpha * (*citer)[0][0] * x[riter.index()];
      }
    }

    template <class Scalar, int n>
    void applyProlongation(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> > const& x,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> >& y)
    {
      axpy<Scalar,n,true>(P,x,y);
    }

    template <class Scalar, int n>
    void applyTransposedProlongation(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> > const& x,
        Dune::BlockVector<Dune::FieldVector<Scalar,n> >& y)
    {
      atxpy<Scalar,n,true>(P,x,y);
    }


    /**
     * Creates a vector of the prolongation matrices (base transformation..) and a
     * vector of the restricted stiffness matrices for each gridlevel
     *
     *
     * Access to BCRS prolongation matrix from level k to k+1:
     * ProlongationStack<Dune::UGGrid<dim> > ps( gridManager.grid() );
     * BCRSMatrix<FieldMatrix<double,1,1> >  prolongation = ps.stack[k];
     *
     * Access to k times restricted BCRS stiffness matrix with assembled Matrix A:
     * ProlongationStack<Dune::UGGrid<dim> > ps( gridManager.grid() );
     * MultiLevelStack<Dune::UGGrid<dim> > myStack( ps , A );
     * BCRSMatrix<FieldMatrix<double,1,1> >  levelStiff = myStack.levelMatrixStack[k];
     *
     *
     * We have three classes:  ParentalNodes
     *			   ProlongationStack
     *			   MultiLevelStack
     */

    /**
     * \brief Finds the parent nodes and their interpolation weight for each node in the grid.
     * usage eg.:  ParentalNodes<Grid> parentalNodes( gridManager.grid() );
     *
     * creates a "level-nested" vector "levelParentsVector" whose components are vectors, containing
     * the indices of vertices and barycentric coordinates corresponding to the parent cells

     * creates also the "level-nested" Vector "levelDisappear" whose entries will be set to 1 if
     * corresponding node disappears in current level and to 0 else
     */
    template < class Grid >
    class ParentalNodes
    {
      typedef typename Grid::LevelGridView::IndexSet::IndexType 	IndexType;
      typedef typename std::pair<IndexType,float> 			FatherIndexValue;
      typedef typename std::vector<FatherIndexValue> 		FatherSimplexPairs;

    public:
      std::vector< std::vector<FatherSimplexPairs> >  	levelParentsVector;
      std::vector< std::vector<int> > 			levelDisappear;

      ParentalNodes( Grid const& grd ): levelParentsVector(grd.maxLevel()), levelDisappear(grd.maxLevel())
      {
        for(int levelNr=1; levelNr<=grd.maxLevel(); levelNr++)
          findParents(grd,levelNr);
      }

    private:
      void findParents( const Grid& grd, int levelNr )
      {
        // First we iterate through all elements in
        // the current level. In each element we iterate through each corner node. First we check if we
        // visited it before (vector "processed"). If not we collect the barycentric coordinates in
        // the father element and the corresponding level indices of the father corner vertex in the
        // level above ( using vector<pair<int index,double barycentricCoordinate> > >, each vector
        // consisting of maximal four pairs ). We also collect the nodes, which disappear from one to
        // another level, we need them (!) to build up consistent prolongations.

        typedef typename Grid::template Codim<0>::LevelIterator       ElementLevelIterator;
        typedef typename Grid::LevelGridView::IndexSet::IndexType     IndexType;

        std::vector<FatherSimplexPairs> levelParents;
        levelParents.resize( grd.size(levelNr,Grid::dimension) );

        // helpvector to avoid multiple processing of identical vertices
        std::vector<bool> processed(grd.size(Grid::dimension),false);

        // disappearance monitor: necessary to keep track of disappearing nodes
        std::vector<bool> disappears( grd.size(levelNr-1,Grid::dimension) , true );

        typename Grid::LevelGridView::IndexSet const& indexSet = grd.levelView(levelNr).indexSet() ;
        typename Grid::LevelGridView::IndexSet const& fatherIndexSet = grd.levelView(levelNr-1).indexSet();

        // iterate through all cells on current level
        ElementLevelIterator itEnd = grd.template lend<0>(levelNr) ;
        for (ElementLevelIterator it = grd.template lbegin<0>(levelNr) ; it != itEnd; ++it)
        {
          // for storing the father level indices of the father's corner points
          std::vector<IndexType> fatherIndices;

          // VertexPointer to father Corners
          typename Grid::template Codim<0>::EntityPointer itFather = (*it).father();

          for(int i=0 ; i<=Grid::dimension ; i++ )
          {
            fatherIndices.push_back( fatherIndexSet.subIndex(*itFather, i, Grid::dimension) );
            // corners of the father cell won't disappear on coarser level
            disappears[ indexSet.subIndex( *itFather,i,Grid::dimension ) ] = false;
          }

          // now regard corners of subentity
          for( int i=0 ; i<=Grid::dimension ; i++ )
          {
            IndexType idx_child = indexSet.subIndex(*it, i, Grid::dimension);
            if (!processed[idx_child])
            {
              processed[idx_child] = true;

              std::vector<float> relativeCoord;
              relativeCoord.push_back(1);
              for( int component = 0 ; component < Grid::dimension ; component++ )
              {
                relativeCoord.push_back( it->geometryInFather().corner(i)[component] );
                relativeCoord.front() -=  relativeCoord.back();
              }

              for( int component = 0 ; component <= Grid::dimension ; component++ )
              {
                FatherIndexValue couple;
                couple.first =  fatherIndices[component];
                couple.second = relativeCoord[component];
                levelParents[idx_child].push_back(couple);
              }

              // TODO: instead use barycentric() from fem/barycentric.hh, something like this:
              // Dune::FieldVector<typename Grid::ctype,Grid::dimension+1> relativeCoord
              //                    = barycentric(it->geometryInFather().corner(i));
              // ...for() {...
              //       couple.second = relativeCoord[(component+1)%(Grid::dimension+1)];
              // ... }
              //
              // somehow barycentric coordinates are in different order than mine

            }
          }
        }

        // store which nodes will not turn up on the coarser level
        for (int k = 0 ; k < disappears.size(); k++ )
          if(disappears[k])
            levelDisappear[levelNr-1].push_back(k);

        // move data to the global store.
        levelParentsVector[levelNr-1].swap(levelParents);
      }
    };



    /// Compute prolongation matrices between consecutive grid levels.
    template <class Grid, class Scalar=double>
    std::vector<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > > computeProlongations( const Grid& grd )
    {
      std::vector<Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>> stack;

      int maxiLevel = grd.maxLevel();
      //  creates levelParentsVector
      ParentalNodes<Grid> parentalNodes( grd );

      std::vector<std::vector<int> > disappear(parentalNodes.levelDisappear);

      int disappearSumBefore = 0, disappearSumAfter  = 0;

      for(int levelNr = 1 ; levelNr <= maxiLevel ; levelNr++)
      {
        disappearSumBefore = disappearSumAfter;
        disappearSumAfter += disappear[levelNr-1].size();

        int prolCols = grd.size( levelNr-1 , Grid::dimension ) + disappearSumBefore ;
        int prolRows = grd.size( levelNr , Grid::dimension ) + disappearSumAfter ;

        // In order to create a BCRS Prolongation matrix we first compute an upper bound for the
        // number of nonzeros. The nodes on lower levels just retain their values and are not
        // interpolated, hence there is exactly one nonzero entry in the row.
        // Each child node's value is determined by the dim+1 values of its father simplex's corners, such that
        // each row has at most dim+1 nonzeros. Most refinement schemes lead to even sparser prolongation matrices,
        // as not all father's corners contribute. We start with the worst case here.
        int nz = disappearSumAfter + grd.size(levelNr,Grid::dimension)*(Grid::dimension+1);
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> prolMatrix(prolRows,prolCols,nz,Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>::row_wise);

        // Now we create the sparsity pattern. Here we can be more precise and allocate only those entries
        // which are actually nonzero. Remember that "very close to zero" barycentric coordinates have been
        // set to exactly zero.

        auto row = prolMatrix.createbegin(), rend = prolMatrix.createend();
        for (;row!=rend; ++row)
        {
          if (row.index() < disappearSumBefore)
            row.insert( row.index() );
          else if (row.index() >= disappearSumBefore && row.index() < disappearSumAfter  )
            row.insert(  disappearSumBefore + disappear[levelNr-1][ row.index() - disappearSumBefore ] );
          else
            for( int k = 0 ; k <= Grid::dimension ; k++ )
              if (parentalNodes.levelParentsVector[levelNr-1][row.index()-disappearSumAfter][k].second > 0)
                row.insert( parentalNodes.levelParentsVector[levelNr-1][row.index()-disappearSumAfter][k].first + disappearSumBefore );
        }

        // now enter the values
        for(int row = 0 ; row < disappearSumBefore ; row++)
          prolMatrix[row][row] = 1;
        for(int indx = 0 ; indx < disappear[levelNr-1].size(); indx++)
          prolMatrix[ disappearSumBefore + indx ][ disappear[levelNr-1][indx] + disappearSumBefore ] = 1;
        for(int row = disappearSumAfter ; row < prolRows ; row++)
        {
          for( int k = 0 ; k <= Grid::dimension ; k++ )
            if (parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter /*row.index()*/][k].second > 0)
            {
              prolMatrix[row][ parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter][k].first + disappearSumBefore ]
                               = parentalNodes.levelParentsVector[levelNr-1][row-disappearSumAfter][k].second;
            }
        }

        stack.push_back( prolMatrix );
      }

      return stack;
    }

    // ------------------------------------------------------------------------------------------------------

    /**
     * (III) MultiLevelStack:
     *
     * The obtained prolongationmatrices are used to build up a vector of n times restricted
     * stiffness matrices A[k], k=0,...,n-1
     * I used conjugation.hh to multiply P^T A P recursively.
     *
     * In conjugation.hh parameter onlyLowerTriangle can be set true, so that only the lower triangle of A will be visited
     * and only lower triangle of the level stiffness matrices will be created.
     *
     */
    template <class Grid, class Scalar, int nComponents=1, class Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,nComponents,nComponents> > >
    class MultiLevelStack
    {
      typedef Dune::FieldMatrix<Scalar,nComponents,nComponents> MatrixBlock;

    public:
      // sets the member levelMatrixStack, a vector of level stiffness matrices
      MultiLevelStack( Grid const& grid , Matrix const& A , bool onlyLowerTriangle = false)
      : prolongations(computeProlongations(grid)), levelMatrixStack(prolongations.size()+1)
      {
        // number of prolongations ( = # gridlevels - 1 )
        size_t prolongationStackSize = prolongations.size();
        // sets the first entry of levelMatrixStack to A (stiffness matrix on leaflevel)
//        levelMatrixStack.push_back(bcrsA);
        levelMatrixStack.back() = A;

        for(int level = prolongationStackSize; level>0; --level)
        {
          // format needed for toolbox conjugation.hh --> C=C+P^T*A*P
          std::unique_ptr<Dune::MatrixIndexSet> C = conjugationPattern( prolongations[level-1] , levelMatrixStack[level] , onlyLowerTriangle);
          Matrix Cnew;
          // exportIdx from Dune::MatrixIndexSet
          C->exportIdx(Cnew);
          // initialize Cnew
          Cnew = 0.0;
          // now Cnew = Cnew + P^T * A * P , where P and A are level corresponding
          // prolongation- and stiffness matrices
          conjugation( Cnew , prolongations[level-1] , levelMatrixStack[level] , onlyLowerTriangle );
          // append new stiffness matrix
          levelMatrixStack[level-1] = Cnew;
        }
      }

      /// Read-only access to level-th entry of matrix stack, i.e. to the stiffness matrix corresponding to level
      Dune::BCRSMatrix<MatrixBlock> const& stiffnessMatrix(size_t level) const { return levelMatrixStack[level]; }

      /// Number of grid levels.
      size_t size() const { return levelMatrixStack.size(); }

      Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& prolongation(size_t level) const { return prolongations[level]; }

    //private:
    private:
      std::vector<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > > prolongations;
      std::vector<Matrix> levelMatrixStack;
     };
  }
  // ------------------------------------------------------------------------------------------------------

  /**
   * \ingroup linalgsolution
   * MultigridSolver realizes a classical v-cycle multigrid solver
   *
   * truncation criterion is relative defect tolerance ||Ax-b||_2/||b||_2 < maxTol
   */

  template < class Grid, int nComponents=1>
  class MultigridSolver
  {
    typedef double Scalar;
    typedef Dune::BlockVector<Dune::FieldVector<double,nComponents> > CoeffVector ;
    typedef Dune::FieldMatrix<double,nComponents,nComponents> MatrixBlock;

  public:
    //forward declaration
    struct Parameter;

    /**
     * \param A const reference to the stiffness matrix on leaf level.
     */
    MultigridSolver(Dune::BCRSMatrix<MatrixBlock> const& A, Grid const& grid, Parameter params=Parameter())
    : parameter(params), levelStack(grid,A),
      solver(MatrixRepresentedOperator<Dune::BCRSMatrix<MatrixBlock>, CoeffVector, CoeffVector>(levelStack.stiffnessMatrix(0)),DirectType::MUMPS,MatrixProperties::GENERAL),
      mgPre("MULTIGRIDSOLVER: ")
    {}

    template <class Assembler, int row, int col, class BlockFilter>
    MultigridSolver(AssembledGalerkinOperator<Assembler,row,row+1,col,col+1,BlockFilter> const& A, Grid const& grid, Parameter params=Parameter(), bool transposed = false)
      : MultigridSolver( (transposed)
                         ? *(A.getAssembler().template get<MatrixAsTriplet<Scalar> >(false,row,row+1,col,col+1).transpose().template toBCRS<nComponents>())
                         : *(A.getAssembler().template get<MatrixAsTriplet<Scalar> >(false,row,row+1,col,col+1).template toBCRS<nComponents>()),
                        grid,
                        params)
//      : MultigridSolver(*(A.getAssembler().template getBlockPointer<Dune::BCRSMatrix<MatrixBlock>,row,col>()), grid, params)
    {}

    void apply(CoeffVector& solution, CoeffVector& rightHand)
    {
      int level = levelStack.size() - 1;
      double relDefect = 1;

      for(int i=0;i<parameter.maxSteps;i++)
      {
        mgSolveRecursive( solution , rightHand , level );
        // calculate relative Defect
        CoeffVector defect(rightHand);
        levelStack.stiffnessMatrix(level).mmv(solution,defect);
        relDefect = defect.two_norm()/rightHand.two_norm();
        if( relDefect < parameter.maxTol)
        {
          if(parameter.verbose) std::cout << mgPre << "converged after " << i << " iterations (relative defect = " << relDefect << ")." << std::endl;
          return;
        }
      }
      if(parameter.verbose) std::cout << mgPre << "NOT converged after " << parameter.maxSteps << " iterations (relative defect = " << relDefect << ")." << std::endl;
    }

    void setParameter(Parameter p) { parameter = p; }

  private:
    void mgSolveRecursive(CoeffVector& solution, CoeffVector& rightHand, int level)
    {
      if(level==0)
      {
        solver.apply(solution,rightHand);
      }
      else
      {
        // pre-smoothing
        applySmoothing( solution , rightHand , level );
        CoeffVector levelRightHand( levelStack.prolongation(level-1).M() ), tmpRhs(rightHand);
        levelRightHand = 0.0;
        levelStack.stiffnessMatrix(level).mmv( solution , tmpRhs );// mmv(const X &x, Y &y) ... y -= A x
        // restriction residual:
        MultiGridSolver_Detail::applyTransposedProlongation(levelStack.prolongation(level-1), tmpRhs, levelRightHand);
        // initialize level-error
        CoeffVector levelError(levelRightHand);
        levelError = 0;

        mgSolveRecursive( levelError , levelRightHand , level - 1 );
        // prolongate correction
        MultiGridSolver_Detail::axpy(levelStack.prolongation(level-1), levelError, solution);
        // post-smoothing
        applySmoothing( solution , rightHand , level );
      }
    }

    void applySmoothing( CoeffVector& solution , CoeffVector const& rhs , int level )
    {
      JacobiSolver<Scalar,nComponents> jacobi(levelStack.stiffnessMatrix(level), parameter.maxSmoothingSteps, 0.5);
      jacobi.apply(solution,rhs);
    }

    Parameter parameter;
    const  MultiGridSolver_Detail::MultiLevelStack<Grid,Scalar,nComponents> levelStack;
    DirectSolver<CoeffVector,CoeffVector> solver;
    std::string mgPre;

  public:
    struct Parameter
    {
      explicit Parameter(size_t maxSteps_, size_t maxSmoothingSteps_, double maxTol_, double relaxation_=0.5, bool verbose_=false)
      : maxSteps(maxSteps_),maxSmoothingSteps(maxSmoothingSteps_), maxTol(maxTol_), relaxation(relaxation_), verbose(verbose_)
      {}

      explicit Parameter(size_t maxSteps_=100, size_t maxSmoothingSteps_=10, size_t maxCGSteps_=1000, size_t lookahead_=25,
                         double maxTol_=1e-9, double relaxation_=0.5, bool verbose_=false, bool useDirectLevel0Solver_ = true)
      : maxSteps(maxSteps_),maxSmoothingSteps(maxSmoothingSteps_), maxCGSteps(maxCGSteps_), lookahead(lookahead_),
        maxTol(maxTol_), relaxation(relaxation_), verbose(verbose_), useDirectLevel0Solver(useDirectLevel0Solver_)
      {}

      Parameter(Parameter const&) = default;
      Parameter& operator=(Parameter const&) = default;

      size_t maxSteps, maxSmoothingSteps, maxCGSteps = 1000, lookahead = 25;
      double maxTol, relaxation, cgTol = 1e-9, eps = 1e-15;
      bool verbose, useDirectLevel0Solver = true;
    };
  };

  template <class Grid, int nComponents=1> using MultiGridSolver = MultigridSolver<Grid,nComponents>;
  
  template < class Grid, int nComponents=1>
  class MultiGridPreconditioner : public MultiGridSolver<Grid,nComponents>
  {
    typedef double Scalar;
    typedef Dune::FieldMatrix<double,nComponents,nComponents> MatrixBlock;
  public:
    typedef typename MultiGridSolver<Grid,nComponents>::Parameter Parameter;

    MultiGridPreconditioner(Dune::BCRSMatrix<MatrixBlock> const& A, Grid const& grid, Parameter p=Parameter())
     : MultiGridSolver<Grid,nComponents>(A,grid,Parameter(p.maxSteps,p.maxSmoothingSteps,0,p.relaxation,p.verbose))
    {}
    
    template <class Assembler, int row, int col, class BlockFilter>
    MultiGridPreconditioner(AssembledGalerkinOperator<Assembler,row,row+1,col,col+1,BlockFilter> const& A, Grid const& grid, Parameter p=Parameter())
     : MultiGridSolver<Grid,nComponents>(A,grid,Parameter(p.maxSteps,p.maxSmoothingSteps,0,p.relaxation,p.verbose))
    {}
  };
}

#endif
