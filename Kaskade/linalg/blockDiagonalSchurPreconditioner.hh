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

#ifndef BLOCK_DIAGONAL_SCHUR_PRECONDITIONER_HH
#define BLOCK_DIAGONAL_SCHUR_PRECONDITIONER_HH

#include <cmath>
#include <vector>

#include <boost/fusion/include/at_c.hpp>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvercategory.hh>

#include "fem/istlinterface.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/iccprecond.hh"
#include "linalg/iluprecond.hh"
#include "utilities/linalg/scalarproducts.hh"


namespace Kaskade
{
  namespace SchurPreconditionerDetail
  {
    template <class Operator, class Block>
    struct ExtractDomainAndRange
    {
      typedef typename Operator::Assembler::AnsatzVariableSet::template CoefficientVectorRepresentation<Block::firstCol,Block::lastCol>::type Domain;
      typedef typename Operator::Assembler::TestVariableSet::template   CoefficientVectorRepresentation<Block::firstRow,Block::lastRow>::type Range;
    };


    template <class Scalar_>
    class JacobiIteration
    {
      typedef Scalar_ Scalar;
      typedef MatrixAsTriplet<Scalar> Matrix;

    public:
      JacobiIteration() = default;

      template <class Operator>
      explicit JacobiIteration(Operator const& A_)
       : A(A_.template get<Matrix>(false)), diag(A.nrows(),0)
      {
        for(size_t i=0; i<A.ridx.size(); ++i) if(A.ridx[i] == A.cidx[i]) diag[A.ridx[i]] = A.data[i];
      }

      explicit JacobiIteration(Matrix const& A_)
       : A(A_), diag(A.nrows(),0)
      {
        for(size_t i=0; i<A.ridx.size(); ++i) if(A.ridx[i] == A.cidx[i]) diag[A.ridx[i]] = A.data[i];
      }

      template <class Domain, class Range>
      void apply(Domain& x, Range const& y, size_t nIter = 1) const
      {
        std::vector<Scalar> tmpx(x.dim()), tmpy;
        IstlInterfaceDetail::toVector(y,tmpy);
        for(size_t i=0; i<diag.size(); ++i) tmpx[i] = tmpy[i]/diag[i];

        for(size_t iter=0; iter < nIter; ++iter)
        {
          for(size_t i=0; i<A.ridx.size(); ++i) if(A.ridx[i] != A.cidx[i]) tmpy[A.ridx[i]] -= A.data[i]*tmpx[A.cidx[i]];
          for(size_t i=0; i<diag.size(); ++i) tmpx[i] += tmpy[i]/diag[i];
        }
        IstlInterfaceDetail::fromVector(tmpx,x);
      }

    private:
      Matrix const& A;
      std::vector<Scalar> diag;
    };

    template <class Scalar_>
    class InvertLumpedMatrix
    {
      typedef Scalar_ Scalar;
      typedef MatrixAsTriplet<Scalar> Matrix;

    public:
      InvertLumpedMatrix() = default;

      template <class Operator, typename... Args>
      explicit InvertLumpedMatrix(Operator const& A_, Args...)
      {
        Matrix A = A_.template get<Matrix>();
        diag.resize(A.nrows());
        for(size_t i=0; i<A.ridx.size(); ++i) if(A.ridx[i] == A.cidx[i]) diag[A.ridx[i]] = A.data[i];
      }

      template <typename... Args>
      explicit InvertLumpedMatrix(Matrix const& A, Args...)
       : diag(A.nrows(),0)
      {
        for(size_t i=0; i<A.ridx.size(); ++i) if(A.ridx[i] == A.cidx[i]) diag[A.ridx[i]] = A.data[i];
      }

      template <class Domain, class Range>
      void apply(Domain& x, Range const& y) const
      {
        tmpx.resize(x.dim());
        IstlInterfaceDetail::toVector(y,tmpy);
        for(size_t i=0; i<diag.size(); ++i) tmpx[i] = tmpy[i]/diag[i];
        IstlInterfaceDetail::fromVector(tmpx,x);
      }

    private:
      std::vector<Scalar> diag;
      mutable std::vector<Scalar> tmpx, tmpy;
    };

    
    template <class Domain, class Range>
    class ApplyDirectSolver
    {

    public:
      ApplyDirectSolver() = default;

      template <class Operator>
      explicit ApplyDirectSolver(Operator const& A, DirectType directType = DirectType::MUMPS, MatrixProperties property = MatrixProperties::GENERAL)
       : solver(DirectSolver<Domain,Range>(A,directType,property))
      {}

      explicit ApplyDirectSolver(MatrixAsTriplet<typename Domain::Scalar> const& A, DirectType directType = DirectType::MUMPS, MatrixProperties property = MatrixProperties::GENERAL)
       : solver(DirectSolver<Domain,Range>(A,directType,property))
      {}

      void apply(Domain& x, Range const& y) const 
      {
        solver.apply(y,x);
      }
      
    private:
      InverseLinearOperator<DirectSolver<Domain,Range> > solver;
    };  
  } /* end of namespace PreconditionerDetail */

  
  /// Approximation of the schur complement according to Pearson/Wathen '10
  /**
   *  \param BlockK type: IstlInterfaceDetail::BlockInfo, information on the block containing the stiffness matrix
   *  \param BlockM type: IstlInterfaceDetail::BlockInfo, information on the block containing the mass matrix
   */
  template <class Scalar_, class Domain_, class Range_,
             class BlockK = IstlInterfaceDetail::BlockInfo<2,3,0,1>,
             class BlockM = IstlInterfaceDetail::BlockInfo<0,1,0,1>,
             class LinearSolver = SchurPreconditionerDetail::ApplyDirectSolver<Domain_,Range_>
             >
  class ApproximateSchurComplement
  {
  public:
    typedef Scalar_ Scalar;
    typedef Domain_ Domain;
    typedef Range_ Range;
    typedef MatrixAsTriplet<Scalar> Matrix;

    /**
     * \param alpha Tikhonov regularization parameter
     * \param symmetricK_ true if stiffness matrix is symmetric, else false
     * \param symmetricM_ true if mass matrix is symmetric, else false
     */
    template <class Assembler>
    ApproximateSchurComplement(Assembler const& assembler, Scalar alpha, bool symmetricK_ = false, bool symmetricM_ = false)
    : symmetricK(symmetricK_ && symmetricM_), symmetricM(symmetricM_), propertyK(symmetricK ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL), propertyM(symmetricM ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL)
    {
      updateMatrices(assembler,alpha);
    }
    

    template <class Matrix>
    ApproximateSchurComplement(Matrix const& M_, Matrix const& K_, Scalar alpha, bool symmetricK_ = false, bool symmetricM_ = false)
    : sqrtOfAlpha(sqrt(alpha)), symmetricK(symmetricK_ && symmetricM_), symmetricM(symmetricM_),
      propertyK(symmetricK ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL), propertyM(symmetricM ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL),
      M(M_), K(K_), N(M_)
    {
      N /= sqrtOfAlpha;
      N += K;
      initSolver();
    }

    /**
     * \param alpha Tikhonov regularization parameter
     */
    template <class Assembler>
    void updateMatrices(Assembler const& assembler, Scalar alpha)
    {
      sqrtOfAlpha = sqrt(alpha);
      M = assembler.template get<Matrix,BlockM::firstRow,BlockM::lastRow,BlockM::firstCol,BlockM::lastCol>(symmetricM);
      N = M.get();
      N /= sqrtOfAlpha;
      K = assembler.template get<Matrix,BlockK::firstRow,BlockK::lastRow,BlockK::firstCol,BlockK::lastCol>(symmetricK);
      N += K;
      
      initSolver();
    }
    
    template <class Matrix>
    void updateMatrices(Matrix const& M_, Matrix const& K_, Scalar alpha)
    {
      M = M_;
      K = K_;
      N = M.get();
      N /= sqrtOfAlpha;
      N += K;
    }

    /// compute \f$x=S^{-1}y\f$, overrides y
    void solve(Domain& x, Range& y) const
    {         
      solver.apply(x,y);
      M.apply(x,y);
      solver.apply(x,y);
    }
    
    Matrix const& getMassMatrix() const
    {
      return M.get();
    }

    Matrix const& getStiffnessMatrix() const
    {
      return K;
    }

  private:
    void initSolver()
    {
      DirectType directType = DirectType::MUMPS;
      MatrixProperties property = (symmetricK && symmetricM) ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL;
      solver = LinearSolver(N,directType,property);
    }

    /// Tikhonov regularization parameter
    Scalar sqrtOfAlpha;
    bool symmetricK, symmetricM;
    MatrixProperties propertyK, propertyM;
    /// mass matrix
    MatrixRepresentedOperator<Matrix,Domain,Range> M;
    /// stiffness matrix
    MatrixAsTriplet<Scalar> K;
    /// sum of stiffness matrix and scaled mass matrix
    /**
     * Let \f$K\f$ be the stiffness matrix, \f$M\f$ be the mass matrix
     * and let \f$\alpha\f$ be the Tikhonov regularization parameter. Then
     * \f$N=K+\frac{1}{\sqrt{\alpha}}M\f$.
     */
    MatrixRepresentedOperator<Matrix,Domain,Range> N;

    LinearSolver solver;
  };
  

  /// Approximation of the schur complement
  /**
   *  \param BlockK type: IstlInterfaceDetail::BlockInfo, information on the block containing the stiffness matrix
   *  \param BlockM type: IstlInterfaceDetail::BlockInfo, information on the block containing the mass matrix
   */
  template <class Scalar_,
             class BlockK = IstlInterfaceDetail::BlockInfo<2,3,0,1>,
             class BlockM = IstlInterfaceDetail::BlockInfo<0,1,0,1>
             >
  class ApproximateSchurComplement2
  {
  public:
    typedef Scalar_ Scalar;
    typedef MatrixAsTriplet<Scalar> Matrix;

    /**
     * \param alpha Tikhonov regularization parameter
     * \param symmetricK_ true if stiffness matrix is symmetric, else false
     * \param symmetricM_ true if mass matrix is symmetric, else false
     */
    template <class Assembler>
    ApproximateSchurComplement2(Assembler const& assembler, Scalar /*alpha*/, bool symmetricK_ = false, bool symmetricM_ = false)
    : symmetricK(symmetricK_ && symmetricM_), symmetricM(symmetricM_), propertyK(symmetricK ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL), propertyM(symmetricM ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL)
    {
      updateMatrices(assembler);
    }

    /**
     * \param alpha Tikhonov regularization parameter
     */
    template <class Assembler>
    void updateMatrices(Assembler const& assembler)
    {
      M = assembler.template get<Matrix,BlockM::firstRow,BlockM::lastRow,BlockM::firstCol,BlockM::lastCol>(symmetricM);
      L = M;
      K = assembler.template get<Matrix,BlockK::firstRow,BlockK::lastRow,BlockK::firstCol,BlockK::lastCol>(symmetricK);
      L += K;
    }


    template <class Domain, class Range>
    void solve(Domain& x, Range const& y, size_t nIter) const
    {
      InverseLinearOperator<DirectSolver<Domain,Range> > Kinv(DirectSolver<Domain,Range>(K,DirectType::MUMPS,propertyK));

      Kinv.apply(y,x);

      std::vector<Scalar> tmpx, tmpy;
      IstlInterfaceDetail::toVector(x,tmpx);
      L.mv(tmpx,tmpy);
      Range tmp(y);
      IstlInterfaceDetail::fromVector(tmpy,tmp);

      Kinv.apply(tmp,x);
    }

    MatrixAsTriplet<Scalar> const& massMatrix() const
    {
      return M;
    }

    MatrixAsTriplet<Scalar> const& stiffnessMatrix() const
    {
      return K;
    }

  private:
    /// sum of stiffness matrix and mass matrix
    MatrixAsTriplet<Scalar> L;
    /// mass matrix
    MatrixAsTriplet<Scalar> M;
    /// stiffness matrix
    MatrixAsTriplet<Scalar> K;
    bool symmetricK, symmetricM;
    MatrixProperties propertyK, propertyM;
  };


  /**
   * Preconditioner à la Pearson/Wathen '10
   * 
   * based on an approximation of the schur complement
   * mesh independent and independent of the Tikhonov regularization parameter
   */
  template <class Operator, template <class,class,class,class,class,class> class SchurComplement = ApproximateSchurComplement,
             class BlockK = IstlInterfaceDetail::BlockInfo<2,3,0,1>,
             class BlockM = IstlInterfaceDetail::BlockInfo<0,1,0,1>,
             class LinearSolver = SchurPreconditionerDetail::ApplyDirectSolver<typename SchurPreconditionerDetail::ExtractDomainAndRange<Operator,BlockM>::Domain, typename SchurPreconditionerDetail::ExtractDomainAndRange<Operator,BlockM>::Range>,
             class LinearSolver_SchurComplement = LinearSolver
            >
  class BlockDiagonalSchurPreconditioner: public Dune::Preconditioner<typename Operator::Domain,typename Operator::Range>
  {
    typedef typename SchurPreconditionerDetail::ExtractDomainAndRange<Operator,BlockM>::Domain Domain0;
    typedef typename SchurPreconditionerDetail::ExtractDomainAndRange<Operator,BlockM>::Range Range0;
    typedef typename Operator::Assembler Assembler;
    typedef MatrixAsTriplet<typename Operator::Scalar> Matrix;

  public:
    typedef typename Operator::Scalar Scalar;
    typedef typename Operator::Domain Domain;
    typedef typename Operator::Range Range;
      
    static int const category = Dune::SolverCategory::sequential;
    
    explicit BlockDiagonalSchurPreconditioner(Operator const& A_,Scalar alpha_, bool symmetricK = false, bool symmetricM = false)
     : A(A_), alpha(alpha_), schurComplement(A.getAssembler(), alpha, symmetricK, symmetricM),
       propertyM(symmetricM ? MatrixProperties::SYMMETRIC : MatrixProperties::GENERAL), solver(schurComplement.getMassMatrix(),DirectType::MUMPS,propertyM)
    {}
    
    virtual ~BlockDiagonalSchurPreconditioner(){}

    virtual void apply (Domain& x, Range const& y)
    {
      using namespace boost::fusion;
      
      Range0 y0(at_c<0>(y.data));
      Domain0 x0(at_c<0>(x.data));

      std::cout << at_c<0>(x.data).dim() << std::endl;
      solver.apply(x0,y0);
      x0 *= 1./alpha;
      at_c<0>(x.data) = at_c<0>(x0.data);
      
      at_c<0>(y0.data) = at_c<1>(y.data);
      solver.apply(x0,y0);
//      x0 *= 1./alpha;
      at_c<1>(x.data) = at_c<0>(x0.data);

      at_c<0>(y0.data) = at_c<2>(y.data);
      schurComplement.solve(x0,y0);
      at_c<2>(x.data) = at_c<0>(x0.data);
    }

    virtual void pre (Domain&, Range&) {}
    virtual void post (Domain&) {}

  private:
    Operator const& A;
    Scalar alpha;
    SchurComplement<Scalar,Domain0,Range0,BlockK,BlockM,LinearSolver_SchurComplement> schurComplement;
    MatrixProperties propertyM;
    LinearSolver solver;
  };
} // namespace Kaskade
#endif
