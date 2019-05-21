/*
 * jacobiSolver.hh
 *
 *  Created on: 21.08.2013
 *      Author: bzflubko
 */

#ifndef JACOBISOLVER_HH_
#define JACOBISOLVER_HH_

namespace Kaskade
{
  /**
   * Jacobi iteration working on block matrices and vectors.
   *
   * \tparam nComponents number of components in each block entry (i.e. number of components of ansatz functions)
   */
  template <class Scalar, int nComponents=1>
  class JacobiSolver : public Dune::Preconditioner<Dune::BlockVector<Dune::FieldVector<Scalar,nComponents> >,Dune::BlockVector<Dune::FieldVector<Scalar,nComponents> > >
  {
  public:
    typedef Dune::FieldMatrix<Scalar,nComponents,nComponents> MatrixBlock;
    typedef Dune::FieldVector<Scalar,nComponents> VectorBlock;

    static int const category = Dune::SolverCategory::sequential;

    /// Constructor taking a reference to a bcrs matrix.
    /**
     * @TODO: constructor taking triplet matrix or operators???
     * This is difficult if we want to avoid unnecessary matrix copies. On the other hand this would allow a unique interface.
     */
    explicit JacobiSolver(Dune::BCRSMatrix<MatrixBlock> const& A_, size_t maxSteps_ = 100, double relaxation_=1.0)
    : A(A_), diag(A.N()), maxSteps(maxSteps_), relaxation(relaxation_)
    {
      extractDiagonal();
    }

    virtual ~JacobiSolver(){}

    void pre(Dune::BlockVector<VectorBlock>&,Dune::BlockVector<VectorBlock>&){}
    void post(Dune::BlockVector<VectorBlock>&){}

    void apply(Dune::BlockVector<VectorBlock>& solution, Dune::BlockVector<VectorBlock> const& rhs)
    {
      assert(solution.N()==rhs.N());
      Dune::BlockVector<VectorBlock> tmp(rhs);

      for(size_t i=0; i<maxSteps; ++i)
      {
        tmp = rhs;
        // apply non diagonal
        auto riter = A.begin(), rend = A.end();
        for(;riter!=rend; ++riter)      // iterate over rows
        {
          size_t row = riter.index();
          auto citer = riter->begin(), cend = riter->end();
          for(;citer!=cend; ++citer)    // iterate over cols
          {
            size_t col = citer.index();
            if(row!=col) tmp[row] -= *citer * solution[col];
          }
        }

        solution *= (1.-relaxation);
        for(size_t i=0; i<solution.N(); ++i)
        {
          solution[i] += relaxation * diag[i] * tmp[i];
        }
      }
    }

  private:
    /// Get diagonal blocks and invert them.
    void extractDiagonal()
    {
      auto riter = A.begin(), rend = A.end();
      for(;riter!=rend; ++riter)
      {
        size_t row = riter.index();
        auto citer = riter->begin(), cend = riter->end();
        for(;citer!=cend && row!=citer.index(); ++citer)
          ;
        if(row==citer.index())
        {
          diag[row] = *citer;
          diag[row].invert();
        }
        else
        {
          std::cout << "WARNING! Singular diagonal block found in jacobi solver" << std::endl;
          diag[row] = 0;
        }
      }
    }

    Dune::BCRSMatrix<MatrixBlock> const& A;
    std::vector<MatrixBlock> diag;
    size_t maxSteps;
    double relaxation;
  };


  /**
   * Jacobi iteration working on block matrices and vectors.
   *
   * \tparam nComponents number of components in each block entry (i.e. number of components of ansatz functions)
   */
  template <class Scalar, class Domain, class Range, class SparseIndexInt = int>
  class JacobiSolverForTriplets : public Dune::Preconditioner<Domain,Range>
  {
  public:
    typedef Dune::FieldMatrix<Scalar,1,1> MatrixBlock;
    typedef Dune::FieldVector<Scalar,1> VectorBlock;

    static int const category = Dune::SolverCategory::sequential;
    explicit JacobiSolverForTriplets(MatrixAsTriplet<Scalar,SparseIndexInt> const& A_, size_t maxSteps_ = 100, double relaxation_=1.0)
    : A(A_), invDiag(A.N()), maxSteps(maxSteps_), relaxation(relaxation_)
    {
      extractDiagonal();
    }

    template <class Assembler, int row, int rend, int col, int cend>
    explicit JacobiSolverForTriplets(AssembledGalerkinOperator<Assembler,row,rend,col,cend> const& A_, size_t maxSteps_ = 100, double relaxation_=1.0)
      : A(A_.getTriplet()), invDiag(A.N()), maxSteps(maxSteps_), relaxation(relaxation_)
    {}

    template <class Assembler, int row, int rend, int col, int cend>
    explicit JacobiSolverForTriplets(TransposedOperator<AssembledGalerkinOperator<Assembler,row,rend,col,cend> > const& A_, size_t maxSteps_ = 100, double relaxation_=1.0)
      : A(A_.getTriplet()), invDiag(A.N()), maxSteps(maxSteps_), relaxation(relaxation_)
    {}

    virtual ~JacobiSolverForTriplets(){}

    void pre(Domain&,Range&){}
    void post(Domain&){}

    void apply(Domain& solution, Range const& rhs)
    {
      assert(solution.dim()==rhs.dim());
      std::vector<Scalar> rhsVec, solVec;
      IstlInterfaceDetail::toVector(rhs,rhsVec);
      IstlInterfaceDetail::toVector(solution,solVec);
      //      for( size_t i=0; i<solVec.size(); ++i)
      //std::cout << i << ": " << solVec[i] << std::endl;

      for(size_t step=0; step<maxSteps; ++step)
      {
        auto tmp = rhsVec;
        // apply non-diagonal
        for(size_t i=0; i<A.data.size(); ++i)
        {
          if(A.ridx[i] != A.cidx[i])
          {
	    //              std::cout << "index: " << A.ridx[i] << ", " << A.cidx[i] << std::endl;
	    //std::cout << "subtracting: " << solVec[A.cidx[i]] << " * " << A.data[i] << std::endl;
              tmp[A.ridx[i]] -= solVec[A.cidx[i]] * A.data[i];
          }
        }
        for( double& d : solVec ) d *= (1. - relaxation);
//        solVec *= (1.-relaxation);
        for(size_t i=0; i<invDiag.size(); ++i)
        {
            solVec[i] += relaxation * invDiag[i] * tmp[i];
        }
      }
      IstlInterfaceDetail::fromVector(solVec,solution);
    }

  private:
    /// Get diagonal blocks and invert them.
    void extractDiagonal()
    {
      for(size_t i=0; i<A.data.size(); ++i)
        if(A.ridx[i] == A.cidx[i])
        {
          if(std::fabs(A.data[i]) > 1e-12)
            invDiag[A.ridx[i]] = 1./A.data[i];
          else
            invDiag[A.ridx[i]] = 0;
        }
    }

    MatrixAsTriplet<Scalar,SparseIndexInt> const& A;
    std::vector<Scalar> invDiag;
    size_t maxSteps;
    double relaxation;
  };

  template <class Scalar, class Domain, class Range, class SparseIndexInt = int>
  class JacobiPreconditionerForTriplets : public Dune::Preconditioner<Domain,Range>
  {
  public:
    typedef Dune::FieldMatrix<Scalar,1,1> MatrixBlock;
    typedef Dune::FieldVector<Scalar,1> VectorBlock;

    static int const category = Dune::SolverCategory::sequential;
    explicit JacobiPreconditionerForTriplets(MatrixAsTriplet<Scalar,SparseIndexInt> const& A)
    : invDiag(A.N())
    {
      extractDiagonal(A);
    }

    template <class Assembler, int row, int rend, int col, int cend>
    explicit JacobiPreconditionerForTriplets(AssembledGalerkinOperator<Assembler,row,rend,col,cend> const& A)
      : JacobiPreconditionerForTriplets(A.template getTriplet())
    {}

    template <class Assembler, int row, int rend, int col, int cend>
    explicit JacobiPreconditionerForTriplets(TransposedOperator<AssembledGalerkinOperator<Assembler,row,rend,col,cend> > const& A)
      : JacobiPreconditionerForTriplets(A.template getTriplet())
    {}

    virtual ~JacobiPreconditionerForTriplets(){}

    void pre(Domain&,Range&){}
    void post(Domain&){}

    void apply(Domain& solution, Range const& rhs)
    {
      assert(solution.dim()==rhs.dim());
      std::vector<Scalar> rhsVec, solVec;
      IstlInterfaceDetail::toVector(rhs,rhsVec);
      IstlInterfaceDetail::toVector(solution,solVec);

      for(size_t i=0; i<invDiag.size(); ++i)
        solVec[i] +=  invDiag[i] * rhsVec[i];

      IstlInterfaceDetail::fromVector(solVec,solution);
    }

  private:
    /// Get diagonal blocks and invert them.
    void extractDiagonal(MatrixAsTriplet<Scalar,SparseIndexInt> const& A)
    {
      for(size_t i=0; i<A.data.size(); ++i)
        if(A.ridx[i] == A.cidx[i])
        {
          if(std::fabs(A.data[i]) > 1e-12)
            invDiag[A.ridx[i]] = 1./A.data[i];
          else
            invDiag[A.ridx[i]] = 0;
        }
    }

    std::vector<Scalar> invDiag;
  };
}


#endif /* JACOBISOLVER_HH_ */
