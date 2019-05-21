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

#ifndef ISTLINTERFACE_HH
#define ISTLINTERFACE_HH

#include <memory>
#include <type_traits>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

#include <dune/istl/operators.hh>
#include "dune/istl/preconditioners.hh"

#include "fem/fixdune.hh"
#include "fem/functionspace.hh"

#include "linalg/triplet.hh"
#include "linalg/crsutil.hh"
#include "linalg/factorization.hh"
#include "linalg/symmetricOperators.hh"
#include "linalg/threadedMatrix.hh"

#include "utilities/duneInterface.hh"

//---------------------------------------------------------------------

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace IstlInterfaceDetail {

    template <class T>
    typename T::iterator begin(T& t)
    {
      return t.begin();
    }

    template <class T>
    typename T::const_iterator begin(T const& t)
    {
      return t.begin();
    }

    template <class Scalar>
    typename std::vector<Scalar>::iterator begin(std::vector<Scalar>& vec) { return vec.begin(); }

    template <class Scalar>
    typename std::vector<Scalar>::const_iterator begin(std::vector<Scalar> const& vec) { return vec.begin(); }

    using namespace boost::fusion;


    /** 
     * \brief This is a lightweight matrix block class, holding a reference to the actual matrix that is owned by the matrix block in the
     *        assembler as well as a NUMA-aware copy.
     * \tparam ElementType          the type of block matrix elements
     * \tparam rid the row index of the block
     * \tparam cid the column index of the block
     * \tparam sym whether the block is symmetric or not
     * \tparam trans whether the block is mirrored
     * 
     * If the block is mirrored, the assembler matrix may have a different element type than this block (m x n instead of n x m).
     */
    template <class ElementType, int rid, int cid, bool sym, bool trans>
    struct Block
    {
      static int const rowId = rid;          // row index of the block
      static int const colId = cid;          // col index of the block
      static bool const symmetric = sym;     // symmetry
      static bool const transposed = trans;  // whether only its mirror image is available in the assembler

      typedef ElementType BlockType;
      typedef typename BlockType::field_type Scalar;
      typedef NumaBCRSMatrix<typename Transpose<BlockType,transposed>::type> Matrix; // Matrix type in assembler
      typedef Dune::BlockVector<Dune::FieldVector<Scalar,BlockType::cols>> Domain; // take care of transposition
      typedef Dune::BlockVector<Dune::FieldVector<Scalar,BlockType::rows>> Range;
      typedef NumaBCRSMatrix<BlockType> TMatrix;

      explicit Block(std::shared_ptr<Matrix const> const& m): matrix(m) {}

      std::shared_ptr<Matrix const> matrix;
      std::shared_ptr<TMatrix> threadedMatrix;
    };

    // Provides a heterogeneous sequence of all matrix blocks appearing
    // conceptually in the Galerkin operator, even if they are not stored
    // explicitly in the assembler because of symmetry.
    template <class GOP>
    class AllBlocks
    {
      struct Copy
      {
        template <class MBlock>
        auto operator()(MBlock const& mb) const {
          typedef std::remove_reference_t<MBlock> MB;
          typedef Block<typename MB::BlockType,MB::rowId,MB::colId,MB::symmetric,false> type;
          return type(mb.globalMatrixPointer());
        }
      };

      // Filters out flagged blocks with row index -1
      struct Cleanup {
        template <class MB> struct apply { typedef boost::mpl::bool_<MB::rowId>=0> type; };
      };

      // Blocks not to be mirrored are flagged with row index -1 to be filtered out later by Cleanup.
      struct Mirror
      {        
        template <class T> struct result {};

        template <class MBlock>
        struct result<Mirror(MBlock)> {
          typedef typename std::remove_reference<MBlock>::type MB;
          typedef typename Transpose<typename MB::BlockType,true>::type ElementType;
          typedef Block<ElementType,
              MB::mirror?MB::colId:-1, MB::rowId, MB::symmetric, true> type;
        };

        template <class MB>
        typename result<Mirror(MB)>::type operator()(MB const& mb) const {
          return typename result<Mirror(MB)>::type(mb.globalMatrixPointer());
        }
      };

      // It would be more natural to filter the mirror blocks before
      // actually mirroring them, but this does frequently result in empty
      // filter lists. For some reasons (don't know why), the filter view
      // produces wrong type results when empty lists occur.

      typedef typename result_of::transform<typename GOP::MatrixBlockArray const,Copy>::type CopiedBlocks;
      typedef typename result_of::transform<typename GOP::MatrixBlockArray const,Mirror>::type MirroredBlocks;
      typedef typename result_of::join<CopiedBlocks const,MirroredBlocks const>::type JointBlocks;


    public:
      typedef typename result_of::filter_if<JointBlocks const,Cleanup>::type const A;
      typedef typename result_of::as_vector<A>::type result;

      static result apply(GOP const& op) {
        return filter_if<Cleanup>(join(transform(op.getMatrix(),Copy()),transform(op.getMatrix(),Mirror())));
      }
    };

    // A boost fusion predicate that is true for some matrix block if it
    // lies within the half-open subrange given by the template
    // parameters.
    template <int firstRow, int lastRow, int firstCol, int lastCol>
    struct RangeBlockSelector
    {
      template <class Block> struct apply
      {
        typedef boost::mpl::bool_<firstRow<=Block::rowId && Block::rowId<lastRow && firstCol<=Block::colId && Block::colId<lastCol> type;
      };
    };

    // Shifts row and column indices of matrix blocks downwards and initializes the threaded matrix held by this block.
    template <int firstRow, int firstCol>
    struct Translate {
      // Constructor takes the number of threads onto which the threaded matrix shall be distributed.
      // As a conservative default, this is 1 (sequential).
      Translate(int nThreads_=1): nThreads(nThreads_)
      {
      }

      template <class T> struct result {};

      template <class Blk>
      struct result<Translate(Blk)> {
        typedef typename std::remove_reference<Blk>::type Bl;
        typedef Block<typename Bl::BlockType,Bl::rowId-firstRow,Bl::colId-firstCol,Bl::symmetric,Bl::transposed> type;
      };

      template <class Bl>
      auto operator()(Bl const& bl) const 
      { 
        Block<typename Bl::BlockType,Bl::rowId-firstRow,Bl::colId-firstCol,Bl::symmetric,Bl::transposed>  blk(bl.matrix); 
        blk.threadedMatrix.reset(new typename Bl::TMatrix(*blk.matrix,bl.symmetric,bl.transposed,false));
        return blk;
      }

    private:
      int nThreads;
    };

    /// Coordinate mapping.
    /**
     * Enabled/disabled via std::enable_if on return type.  
     */
    template <class FSElement, class Vector>
    typename std::enable_if<!std::is_same<FSElement,Vector>::value,void>::type
    toVector(FSElement const& x, Vector& coefficients)
    {
      coefficients.resize(x.dim());
      vectorToSequence(x,begin(coefficients));
    }

    /// Coordinate mapping is the identity if range_type = vector_type
    template <class FSElement, class Vector>
    typename std::enable_if<std::is_same<FSElement,Vector>::value,void>::type
    toVector(FSElement const& x, Vector& coefficients)
    {
      coefficients = x;
    }

    /// Coordinate mapping
    template <class FSElement, class Vector>
    typename std::enable_if<!std::is_same<FSElement,Vector>::value,void>::type
    fromVector(Vector const& coefficients, FSElement& x)
    {
      assert(coefficients.size()>=x.dim());
      vectorFromSequence(x,begin(coefficients));
    }

    /// Coordinate mapping is the identity if range_type = vector_type
    template <class FSElement, class Vector>
    typename std::enable_if<std::is_same<FSElement,Vector>::value,void>::type
    fromVector(Vector const& coefficients, FSElement& x)
    {
      x = coefficients;
    }


    template <class T>
    void toFile(std::vector<T> const& v, std::string const& filename, size_t precision = 10)
    {
      std::ofstream os;
      os.precision(precision);
      os.open(filename);

      for(T t : v) os << t << "\t";
      os << std::endl;

      os.close();
    }


    // Container storing information on blocks
    template <int firstRow_, int lastRow_, int firstCol_,int lastCol_>
    struct BlockInfo
    {
      static int const firstRow = firstRow_;
      static int const lastRow = lastRow_;
      static int const firstCol = firstCol_;
      static int const lastCol = lastCol_;
    };
  } // End of namespace IstlInterfaceDetail

  /**
   * \endcond
   */

  //---------------------------------------------------------------------

/*  template <class Operator>
  class TransposedOperator
  {
  public:
    typedef typename Operator::Domain Range;
    typedef typename Operator::Range Domain;
    typedef typename Operator::Scalar Scalar;

    explicit TransposedOperator(Operator const& A_) : A(A_) {}

    void apply(const Domain &x, Range &b) const { return A.applyTransposed(x,b); }

    void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const { return A.applyscaleaddTransposed(alpha,x,b); }

    Operator const& getOperator() const { return A; }

  private:
    Operator const& A;
  };
*/


  /**
   * \ingroup linalg
   *
   * \brief Operator with matrix-representation and coordinate mappings.
   *
   * Note that the operator creates a copy of the passed matrix (if not passed as rvalue).
   * This allows the implementation of a rvalue-constructor and assures that further manipulation
   * of the matrix data does not change the operator's behaviour.
   *
   * Implements the MatrixRepresentedOperator concept.
   * 
   * \todo docme: How does this relate to Dune::MatrixAdapter and Dune::AssembledLinearOperator?
   */
  template <class Matrix_, class Domain_, class Range_> class MatrixRepresentedOperator;

  template <class Scalar_, class SparseInt, class Domain_, class Range_>
  class MatrixRepresentedOperator<MatrixAsTriplet<Scalar_,SparseInt>,Domain_,Range_>
      : public Dune::LinearOperator<Domain_,Range_>
  {
  public:
    typedef MatrixAsTriplet<Scalar_,SparseInt> Matrix;
    typedef Matrix matrix_type; // for compatibility with Dune::MatrixAdapter
    typedef Domain_ Domain;
    typedef Range_ Range;
    typedef Scalar_ Scalar;

    static int const category = Dune::SolverCategory::sequential;

    MatrixRepresentedOperator() = default;

    explicit MatrixRepresentedOperator(Matrix const& A_) : A(A_)
    {}

    explicit MatrixRepresentedOperator(Matrix&& A_) : A(A_)
    {}

    MatrixAsTriplet<Scalar,SparseInt> getTriplet() const { return get<MatrixAsTriplet<Scalar,SparseInt> >(); }

    /// Access matrix
    template <class OtherMatrix=Matrix>
    typename std::enable_if<std::is_same<Matrix,OtherMatrix>::value, Matrix const&>::type get() const { return A; }

    Matrix& get_non_const() { return A; }

    /// Access matrix in another matrix format specified by OtherMatrix.
    /**
     * Note that OtherMatrix must be constructible from Matrix.
     */
    template <class OtherMatrix>
    typename std::enable_if<!std::is_same<Matrix,OtherMatrix>::value, OtherMatrix>::type get() const { return OtherMatrix(A); }

    /// compute \f$b=Ax\f$
    void apply(Domain const& x, Range& b) const
    {
//      A.mv(x,b);
            std::vector<Scalar> tmpx, tmpb;
            domainToVector(x,tmpx);

            A.mv(tmpx,tmpb);

            vectorToRange(tmpb,b);
    }

    void applyTransposed(Range const& x, Domain& b) const
    {
//      A.mv(x,b);
            std::vector<Scalar> tmpx, tmpb;
            domainToVector(x,tmpx);

            A.mtv(tmpx,tmpb);

            vectorToRange(tmpb,b);
    }

    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    template <class OtherMatrix=Matrix> 
    typename std::enable_if<std::is_same<OtherMatrix,Matrix>::value,std::unique_ptr<Matrix> >::type getPointer() const
    {
      return std::unique_ptr<Matrix>(new Matrix(A));
    }

    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    template <class OtherMatrix> 
    typename std::enable_if<std::is_same<OtherMatrix,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > >::value,std::unique_ptr<OtherMatrix> >::type getPointer() const
    {
      return A.toBCRS();
    }

    /// compute \f$b=b+\alpha Ax\f$
    void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const
    {
//      A.usmv(alpha,x,b);
            std::vector<Scalar> tmpx, tmpb;
            domainToVector(x,tmpx);
            rangeToVector(b,tmpb);

            A.usmv(alpha,tmpx, tmpb);

            vectorToRange(tmpb,b);
    }

    void applyscaleaddTransposed(Scalar alpha, Range const& x, Domain& b) const
    {
//      A.usmv(alpha,x,b);
            std::vector<Scalar> tmpx, tmpb;
            domainToVector(x,tmpx);
            rangeToVector(b,tmpb);

            A.usmtv(alpha,tmpx, tmpb);

            vectorToRange(tmpb,b);
    }
    /// Get coefficient vector \f$ \mathrm{coefficients}\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
    /**
     * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$\mathrm{coefficients}=S_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> 
    void rangeToVector(Range const& y, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(y,coefficients);
    }

    /// Get coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$ from \f$x\in X\f$.
    /**
     * Apply \f$S_X\f$  to \f$x\in X\f$: \f$\mathrm{coefficients}=S_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> 
    void domainToVector(Domain const& x, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(x, coefficients);
    }

    /// Get \f$x\in X\f$ from coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$
    /**
     * Apply \f$S^{-1}_X\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToDomain(Vector const& coefficients, Domain& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /// Get \f$y\in Y\f$ from coefficient vector \f$\mathrm{coefficients}\in\mathbb{K}^m\f$
    /**
     * Apply \f$S^{-1}_Y\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToRange(Vector const& coefficients, Range& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    MatrixRepresentedOperator& operator=(MatrixRepresentedOperator const& other) = default;

    MatrixRepresentedOperator& operator=(Matrix const& B)
    {
      A = B;   return *this;
    }

    MatrixRepresentedOperator& operator*= (Scalar alpha)
    {
      A *= alpha;   return *this;
    }
    
    MatrixRepresentedOperator& operator/= (Scalar alpha)
    {
      A /= alpha;   return *this;
    }
    
    MatrixRepresentedOperator& operator+= (MatrixRepresentedOperator const& B)
    {
      A += B;   return *this;
    }
    
    MatrixRepresentedOperator& operator+= (Matrix const& B)
    {
      A += B;   return *this;
    }
    
    MatrixRepresentedOperator& operator-= (MatrixRepresentedOperator const& B)
    {
      A -= B;   return *this;
    }
    
    MatrixRepresentedOperator& operator-= (Matrix const& B)
    {
      A -= B;   return *this;
    }

  private:
    Matrix A;
  };

  template <class MatrixBlock, class Allocator, class Domain_, class Range_>
  class MatrixRepresentedOperator<Dune::BCRSMatrix<MatrixBlock,Allocator>,Domain_,Range_>
      : public Dune::LinearOperator<Domain_,Range_>
  {
  public:
    typedef Dune::BCRSMatrix<MatrixBlock,Allocator> Matrix;
    typedef Domain_ Domain;
    typedef Range_ Range;
    typedef typename GetScalar<Domain>::type Scalar;

    static int const category = Dune::SolverCategory::sequential;

    MatrixRepresentedOperator() = default;

    explicit MatrixRepresentedOperator(Matrix const& A_) : A(A_)
    {}

    explicit MatrixRepresentedOperator(Matrix&& A_) : A(A_)
    {}

    MatrixAsTriplet<Scalar> getTriplet() const { return get<MatrixAsTriplet<Scalar> >(); }

    /// Access matrix
    template <class OtherMatrix=Matrix>
    typename std::enable_if<std::is_same<Matrix,OtherMatrix>::value, Matrix const&>::type get() const { return A; }

    Matrix& get_non_const() { return A; }

    /// Access matrix in another matrix format specified by OtherMatrix.
    /**
     * Note that OtherMatrix must be constructible from Matrix.
     */
    template <class OtherMatrix>
    typename std::enable_if<!std::is_same<Matrix,OtherMatrix>::value, OtherMatrix>::type get() const { return OtherMatrix(A); }

    /// compute \f$b=Ax\f$
    void apply(Domain const& x, Range& b) const
    {
      A.mv(x,b);
    }

    void applyTransposed(Range const& x, Domain& b) const
    {
      A.mtv(x,b);
    }

    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    template <class OtherMatrix=Matrix>
    typename std::enable_if<std::is_same<OtherMatrix,Matrix>::value,std::unique_ptr<Matrix> >::type getPointer() const
    {
      return std::unique_ptr<Matrix>(new Matrix(A));
    }

    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    template <class OtherMatrix>
    typename std::enable_if<std::is_same<OtherMatrix,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > >::value,std::unique_ptr<OtherMatrix> >::type getPointer() const
    {
      return A.toBCRS();
    }

    /// compute \f$b=b+\alpha Ax\f$
    void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const
    {
      A.usmv(alpha,x,b);
    }

    void applyscaleaddTransposed(Scalar alpha, Range const& x, Domain& b) const
    {
      A.usmtv(alpha,x,b);
    }

    /// Get coefficient vector \f$ \mathrm{coefficients}\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
    /**
     * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$\mathrm{coefficients}=S_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> void rangeToVector(Range const& y, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(y,coefficients);
    }

    /// Get coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$ from \f$x\in X\f$.
    /**
     * Apply \f$S_X\f$  to \f$x\in X\f$: \f$\mathrm{coefficients}=S_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> void domainToVector(Domain const& x, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(x, coefficients);
    }

    /// Get \f$x\in X\f$ from coefficients vector \f$\mathrm{coefficients}\in\mathbb{K}^n\f$
    /**
     * Apply \f$S^{-1}_X\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToDomain(Vector const& coefficients, Domain& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /// Get \f$y\in Y\f$ from coefficient vector \f$\mathrm{coefficients}\in\mathbb{K}^m\f$
    /**
     * Apply \f$S^{-1}_Y\f$ to \f$\mathrm{coefficients}\f$: \f$x=S^{-1}_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToRange(Vector const& coefficients, Range& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    MatrixRepresentedOperator& operator=(MatrixRepresentedOperator const& other) = default;

    MatrixRepresentedOperator& operator=(Matrix const& B)
    {
      A = B;   return *this;
    }

    MatrixRepresentedOperator& operator*= (Scalar alpha)
    {
      A *= alpha;   return *this;
    }
    
    MatrixRepresentedOperator& operator/= (Scalar alpha)
    {
      A /= alpha;   return *this;
    }
    
    MatrixRepresentedOperator& operator+= (MatrixRepresentedOperator const& B)
    {
      A += B;   return *this;
    }
    
    MatrixRepresentedOperator& operator+= (Matrix const& B)
    {
      A += B;   return *this;
    }
    
    MatrixRepresentedOperator& operator-= (MatrixRepresentedOperator const& B)
    {
      A -= B;   return *this;
    }
    
    MatrixRepresentedOperator& operator-= (Matrix const& B)
    {
      A -= B;   return *this;
    }
    
  private:
    Matrix A;
  };

  //--------------------------------------------------------------------------------------------

  /// \cond internals
  namespace IstlInterfaceDetail {


    template <class Assembler,
              int firstRow, int lastRow,
              int firstCol, int lastCol,
              bool symmetric=false>
    struct Base
    {
      typedef typename Assembler::AnsatzVariableSet::template CoefficientVectorRepresentation<firstCol,lastCol>::type Domain;
      typedef typename Assembler::TestVariableSet::template CoefficientVectorRepresentation<firstRow,lastRow>::type   Range;

      static bool const inferredSymmetry = (lastRow-firstRow == lastCol-firstCol)              // make sure the number of domain/range variables is the same
                                                &&
                                                ((Assembler::Functional::type==VariationalFunctional       // a diagonally-centered block of a symmetric block matrix
                                                    && firstRow==firstCol && lastRow==lastCol)
                                                    || (firstRow+1==lastRow && firstCol+1==lastCol && false)); /* true if this one block is symmetric */  // TODO: implement the symmetry check

      typedef typename std::conditional<symmetric,
          SymmetricLinearOperator<Domain,Range>,
          Dune::LinearOperator<Domain,Range>>::type type;
    };
  } // end of namespace IstlInterfaceDetail
  /// \endcond



  /**
   * \ingroup linalg
   * \brief An adaptor presenting a Dune::LinearOperator <domain_type,range_type> interface to a
   * contiguous sub-blockmatrix of an assembled heterogeneous Galerkin operator.
   * 
   * \tparam Assembler the \ref VariationalFunctionalAssembler containing the Galerkin matrix
   * \tparam firstRow the first block row to be included
   * \tparam lastRow one behind the final row to be included (half-open range)
   * \tparam firstCol the first block column to be included
   * \tparam lastCol one behind the final column to be included (half-open range)
   * \tparam BlockFilter a boost::fusion predicate working on the assembler's matrix blocks, defining which blocks to include
   * \tparam symmetric whether the resulting operator is symmetric. The default is a conservative check whether the operator is symmetric or not.
   *
   * Due to data sharing, the GalerkinOperator is still valid
   * even if the referenced assembler is deleted or otherwise
   * invalidated. However, reassembling without invalidation of the
   * assembler's data structures also modifies the values in the
   * GalerkinOperator.
   */
  template <class Assembler_,
            int firstRow=0, int lastRow=Assembler_::TestVariableSet::noOfVariables,
            int firstCol=0, int lastCol=Assembler_::AnsatzVariableSet::noOfVariables,
            class BlockFilter=IstlInterfaceDetail::RangeBlockSelector<firstRow,lastRow,firstCol,lastCol>,
            bool symmetric=IstlInterfaceDetail::Base<Assembler_,firstRow,lastRow,firstCol,lastCol>::inferredSymmetry>
  class AssembledGalerkinOperator: public IstlInterfaceDetail::Base<Assembler_,firstRow,lastRow,firstCol,lastCol,symmetric>::type
  {
  public:
    typedef typename IstlInterfaceDetail::Base<Assembler_,firstRow,lastRow,firstCol,lastCol,symmetric>::type Base;
    typedef Assembler_ Assembler;
    typedef typename Assembler::AnsatzVariableSet::template CoefficientVectorRepresentation<firstCol,lastCol>::type Domain;
    typedef typename Assembler::TestVariableSet::template CoefficientVectorRepresentation<firstRow,lastRow>::type   Range;
    typedef typename Assembler::Scalar                                                                              Scalar;

    typedef IstlInterfaceDetail::BlockInfo<firstRow,lastRow,firstCol,lastCol> BlockInfo;

    // do we need this?
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \brief Constructor.
     * \param op the variational assembler.
     * \param onlyLowerTriangle if true, only the lower triangle is returned if a matrix is requested
     * \param nThreads use given number of threads (or a machine-dependent default if nThreads<=0).
     */
    explicit AssembledGalerkinOperator(Assembler const& assembler_, bool onlyLowerTriangle_=false, int nThreads_=0):
        onlyLowerTriangle(onlyLowerTriangle_),
        assembler(assembler_),
        blocks(boost::fusion::transform(boost::fusion::filter_if<BlockFilter>(IstlInterfaceDetail::AllBlocks<Assembler>::apply(assembler)),
                                        IstlInterfaceDetail::Translate<firstRow,firstCol>(nThreads_))),
        nThreads(nThreads_)
    { }

    virtual ~AssembledGalerkinOperator(){}

    /// update operator if grid has changed or assemble(...) has been called.
    void update()
    {
      blocks = boost::fusion::transform(boost::fusion::filter_if<BlockFilter>(IstlInterfaceDetail::AllBlocks<Assembler>::apply(assembler)),
                                        IstlInterfaceDetail::Translate<firstRow,firstCol>(nThreads));
    }

    virtual MatrixAsTriplet<Scalar> getTriplet() const { return get<MatrixAsTriplet<Scalar> >(); }

    /// Access matrix. Use if Matrix does support move-semantics.
    template <class Matrix> 
    Matrix get() const
    {
      return assembler.template get<Matrix>(onlyLowerTriangle, firstRow, lastRow, firstCol, lastCol);
    }

    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    template <class Matrix> 
    std::unique_ptr<Matrix> getPointer() const
    {
      return assembler.template getPointer<Matrix,firstRow,lastRow,firstCol,lastCol>(onlyLowerTriangle);
    }

    /**
     * \brief returns a reference to the matrix
     */
    MatrixAsTriplet<Scalar> getmat() const __attribute__((deprecated));

    /// Get coefficient vector \f$coefficients\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
    /**
     * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$coefficients=S_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> 
    void rangeToVector(Range const& y, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(y, coefficients);
    }

    /// Get coefficients vector \f$coefficients\in\mathbb{K}^n\f$ from \f$x\in X\f$.
    /**
     * Apply \f$S_X\f$  to \f$x\in X\f$: \f$coefficients=S_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> 
    void domainToVector(Domain const& x, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(x, coefficients);
    }

    /// Get \f$x\in X\f$ from coefficients vector \f$coefficients\in\mathbb{K}^n\f$
    /**
     * Apply \f$S^{-1}_X\f$ to \f$coefficients\f$: \f$x=S^{-1}_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToDomain(Vector const& coefficients, Domain& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /**
     * \brief Get \f$ y\in Y \f$ from coefficient vector \f$coefficients\in\mathbb{K}^m\f$
     * Apply \f$S^{-1}_Y\f$ to \f$coefficients\f$: \f$x=S^{-1}_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToRange(Vector const& coefficients, Range& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /// compute \f$ b \leftarrow Ax \f$
    virtual void apply(Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(blocks,ApplyScaleAdd(assembler,1.0,x,b,true));
    }

    void applyTransposed(Range const& x, Domain& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(blocks,TransposedApplyScaleAdd(assembler,1.0,x,b,true));
    }

    /**
     * \brief Computes \f$ b \leftarrow Ax \f$ and, in case \f$ A \f$ is symmetric, also \f$ \langle x, b \rangle \f$
     * 
     * If \f$ A \f$ is not symmetric, zero is returned.
     */
    virtual Scalar applyDp(Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      Scalar dp = 0;
      boost::fusion::for_each(blocks,ApplyScaleAdd(assembler,1.0,x,b,true,&dp));
      return symmetric? dp: 0;
    }

    virtual Scalar dp(Domain const& x, Range const& y) const
    {
      std::cerr << "AssembledGalerkinOperator::dp not yet implemented!\n"; abort();
      return 0; // TODO: fixme
    }

    /**
     * \brief  Compute \f$b=b+\alpha Ax\f$
     * Note that x and b must not refer to the same memory locations (in case Domain==Range).
     */
    virtual void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(blocks,ApplyScaleAdd(assembler,alpha,x,b,false));
    }

    virtual void applyscaleaddTransposed(Scalar alpha, Range const& x, Domain& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(blocks,TransposedApplyScaleAdd(assembler,alpha,x,b,false));
    }

    /**
     * \brief Provides access to the underlying assembler.
     */
    Assembler const& getAssembler() const 
    { 
      return assembler; 
    }

    bool getOnlyLowerTriangle() const { return onlyLowerTriangle; }

  protected:
    typedef typename boost::fusion::result_of::as_vector<
        typename IstlInterfaceDetail::AllBlocks<Assembler>::result
        >::type allblocks;

    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::filter_if<allblocks,BlockFilter>::type
        >::type filtered;

    bool onlyLowerTriangle;
    Assembler const& assembler;

    typename boost::fusion::result_of::as_vector<
    typename boost::fusion::result_of::transform<filtered,
    IstlInterfaceDetail::Translate<firstRow,firstCol>>::type
    >::type blocks;
    bool isSymmetric;


    struct ApplyScaleAdd 
    {
      Assembler const& assembler;
      // Constructs the functor for matrix-vector multiplication and addition. If init is true, the first block
      // encountered in each row is called with apply (using alpha=1 and initializing the rhs to zero), otherwise
      // applyscaleadd is used. This prevents a double access of the right hand side, first for setting to zero,
      // then for adding the MV product.
      ApplyScaleAdd(Assembler const& as, Scalar a_, Domain const& x_, Range& y_, bool init, Scalar* dp_=0): assembler(as), a(a_), x(x_), y(y_), dp(dp_? dp_: &dummy)
      {
        for (int i=0; i<lastRow; ++i)
          doInit[i] = init;
        *dp = 0;
      }

      template <class Block>
      void operator()(Block const& b) const {
        using namespace boost::fusion;

        typename result_of::at_c<typename Domain::Functions const,Block::colId>::type xi = at_c<Block::colId>(x.data);
        typename result_of::at_c<typename Range::Functions,Block::rowId>::type  yi = at_c<Block::rowId>(y.data);
        if (!b.threadedMatrix.get())
          abort();

        if (doInit[Block::rowId])
        {
          *dp += b.threadedMatrix->mv(xi,yi);
          doInit[Block::rowId] = false;
        }
        else
          *dp += b.threadedMatrix->usmv(a,xi,yi);      

        return;
      }

    private:
      Scalar a;
      Domain const& x;
      Range& y;
      Scalar* dp;
      Scalar dummy;                 // something to point to in case no target pointer has been provided
      mutable bool doInit[lastRow]; // ugly, but boost::fusion::for_each only works with const functors
    };

    struct TransposedApplyScaleAdd
    {
      Assembler const& assembler;
      // Constructs the functor for matrix-vector multiplication and addition. If init is true, the first block
      // encountered in each row is called with apply (using alpha=1 and initializing the rhs to zero), otherwise
      // applyscaleadd is used. This prevents a double access of the right hand side, first for setting to zero,
      // then for adding the MV product.
      TransposedApplyScaleAdd(Assembler const& as, Scalar a_, Range const& x_, Domain& y_, bool init, Scalar* dp_=0): assembler(as), a(a_), x(x_), y(y_), dp(dp_? dp_: &dummy)
      {
        for (int i=0; i<lastCol; ++i)
          doInit[i] = init;
        *dp = 0;
      }

      template <class Block>
      void operator()(Block const& b) const {
        using namespace boost::fusion;

        typename result_of::at_c<typename Range::Functions const,Block::rowId>::type xi = at_c<Block::rowId>(x.data);
        typename result_of::at_c<typename Domain::Functions,Block::colId>::type  yi = at_c<Block::colId>(y.data);
        if (!b.threadedMatrix.get())
          abort();

        if (doInit[Block::colId])
        {
          *dp += b.threadedMatrix->mv(xi,yi);
          doInit[Block::colId] = false;
        }
        else
          *dp += b.threadedMatrix->usmv(a,xi,yi);

        return;
      }

    private:
      Scalar a;
      Range const& x;
      Domain& y;
      Scalar* dp;
      Scalar dummy;                 // something to point to in case no target pointer has been provided
      mutable bool doInit[lastCol]; // ugly, but boost::fusion::for_each only works with const functors
    };

    int nThreads;
  };


  /**
   * \todo docme
   */
  template <class NormalStepAssembler, class TangentialStepAssembler, int stateId=1, int controlId=0, int adjointId=2>
  class LagrangeOperator: public Dune::LinearOperator<typename NormalStepAssembler::AnsatzVariableSet::template CoefficientVectorRepresentation<>::type,typename NormalStepAssembler::TestVariableSet::template CoefficientVectorRepresentation<>::type>
  {
    typedef typename NormalStepAssembler::AnsatzVariableSet::template CoefficientVectorRepresentation<>::type Domain;
    typedef typename NormalStepAssembler::TestVariableSet::template CoefficientVectorRepresentation<>::type Range;
    static constexpr int firstRow = 0;
    static constexpr int firstCol = 0;
    static constexpr int lastRow = NormalStepAssembler::TestVariableSet::noOfVariables;
    static constexpr int lastCol = NormalStepAssembler::AnsatzVariableSet::noOfVariables;
    typedef IstlInterfaceDetail::RangeBlockSelector<firstRow,lastRow,firstCol,lastCol> NormalStepBlockFilter;
    typedef IstlInterfaceDetail::RangeBlockSelector<firstRow,adjointId,firstCol,adjointId> TangentialStepBlockFilter;

  public:
    typedef typename NormalStepAssembler::Scalar Scalar;
    typedef MatrixAsTriplet<Scalar> Triplet;

    typedef IstlInterfaceDetail::BlockInfo<firstRow,lastRow,firstCol,lastCol> BlockInfo;

    static int const category = Dune::SolverCategory::sequential;

    /**
     * \param op the variational assembler.
     * \param onlyLowerTriangle
     * \param nThreads use given number of threads (or a machine-dependent default if nThreads<=0).
     */
    LagrangeOperator(NormalStepAssembler const& normalStepAssembler_, TangentialStepAssembler const& tangentialStepAssembler_,  bool onlyLowerTriangle_=false, int nThreads_=0):
      onlyLowerTriangle(onlyLowerTriangle_),
      normalStepAssembler(normalStepAssembler_), tangentialStepAssembler(tangentialStepAssembler_),
      normalStepBlocks(boost::fusion::transform(boost::fusion::filter_if<NormalStepBlockFilter>(IstlInterfaceDetail::AllBlocks<NormalStepAssembler>::apply(normalStepAssembler)),IstlInterfaceDetail::Translate<firstRow,firstCol>(nThreads_))),
      tangentialStepBlocks(boost::fusion::transform(boost::fusion::filter_if<TangentialStepBlockFilter>(IstlInterfaceDetail::AllBlocks<TangentialStepAssembler>::apply(tangentialStepAssembler)),IstlInterfaceDetail::Translate<firstRow,firstCol>(nThreads_))),
      nThreads(nThreads_)
    {}

    virtual ~LagrangeOperator(){}

    /// update operator if grid has changed or assemble(...) has been called.
    void update()
    {
      assert(false);
      //      blocks = boost::fusion::transform(boost::fusion::filter_if<BlockFilter>(IstlInterfaceDetail::AllBlocks<Assembler>::apply(assembler)),
      //          IstlInterfaceDetail::Translate<firstRow,firstCol>(nThreads));
    }

    virtual MatrixAsTriplet<Scalar> getTriplet() const
    {
      std::vector<size_t> offset(3,0);
      Triplet result, tmp;
      for(size_t i=0; i<3; ++i)
      {
        for(size_t j=0; j<3; ++j)
        {
          if(i==adjointId || j== adjointId || (i==controlId && j==controlId)) tmp = normalStepAssembler.template get<Triplet>(onlyLowerTriangle,i,i+1,j,j+1);
          else tmp = tangentialStepAssembler.template get<Triplet>(onlyLowerTriangle,i,i+1,j,j+1);
          tmp.shiftIndices(offset[i],offset[j]);
          result += tmp;
          if(i==0)
          {
            if(j==0) offset[1] = tmp.nrows();
            if(j==1) offset[2] = tmp.nrows();
          }
        }
      }

      return result;
    }

    //    /// Access matrix. Use if Matrix does support move-semantics.
    //    template <class Matrix> Matrix get() const
    //    {
    //      return assembler.template get<Matrix>(onlyLowerTriangle, firstRow, lastRow, firstCol, lastCol);
    //    }

    //    /// Access matrix via unique_ptr. Use if Matrix does not support move-semantics.
    //    template <class Matrix> std::unique_ptr<Matrix> getPointer() const
    //        {
    //      return assembler.template getPointer<Matrix,firstRow,lastRow,firstCol,lastCol>(onlyLowerTriangle);
    //  }

    /**
     * \brief returns a reference to the matrix
     */
    //    MatrixAsTriplet<Scalar> getmat() const __attribute__((deprecated));

    /// Get coefficient vector \f$coefficients\in\mathbb{K}^m\f$ from \f$y\in Y\f$.
    /**
     * Apply \f$S_Y\f$  to \f$y\in Y\f$: \f$coefficients=S_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> void rangeToVector(Range const& y, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(y, coefficients);
    }

    /// Get coefficients vector \f$coefficients\in\mathbb{K}^n\f$ from \f$x\in X\f$.
    /**
     * Apply \f$S_X\f$  to \f$x\in X\f$: \f$coefficients=S_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::iterator std::begin(Vector&);
     */
    template <class Vector> void domainToVector(Domain const& x, Vector& coefficients) const
    {
      IstlInterfaceDetail::toVector(x, coefficients);
    }

    /// Get \f$x\in X\f$ from coefficients vector \f$coefficients\in\mathbb{K}^n\f$
    /**
     * Apply \f$S^{-1}_X\f$ to \f$coefficients\f$: \f$x=S^{-1}_X(x)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToDomain(Vector const& coefficients, Domain& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /**
     * \brief Get \f$ y\in Y \f$ from coefficient vector \f$coefficients\in\mathbb{K}^m\f$
     * Apply \f$S^{-1}_Y\f$ to \f$coefficients\f$: \f$x=S^{-1}_Y(y)\f$.
     *
     * The used vector type Vector must provide:
     * - its iterator type via typename Vector::const_iterator
     * - (possibly overloads of) the free functions:
     *    - typename Vector::const_iterator std::begin(Vector const&);
     */
    template <class Vector>
    void vectorToRange(Vector const& coefficients, Range& x) const
    {
      IstlInterfaceDetail::fromVector(coefficients, x);
    }

    /// compute \f$ b \leftarrow Ax \f$
    virtual void apply(Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(normalStepBlocks,NormalStepApplyScaleAdd(normalStepAssembler,1.0,x,b,true));
      boost::fusion::for_each(tangentialStepBlocks,TangentialStepApplyScaleAdd(tangentialStepAssembler,1.0,x,b,false));
    }

    /**
     * \brief Computes \f$ b \leftarrow Ax \f$ and, in case \f$ A \f$ is symmetric, also \f$ \langle x, b \rangle \f$
     *
     * If \f$ A \f$ is not symmetric, zero is returned.
     */
    virtual Scalar applyDp(Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      std::cerr << "AssembledGalerkinOperator::dp not yet implemented!\n"; abort();
      return 0;
    }

    virtual Scalar dp(Domain const& x, Range const& y) const
    {
      std::cerr << "AssembledGalerkinOperator::dp not yet implemented!\n"; abort();
      return 0; // TODO: fixme
    }

    /**
     * \brief  Compute \f$b=b+\alpha Ax\f$
     * Note that x and b must not refer to the same memory locations (in case Domain==Range).
     */
    virtual void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const
    {
      assert(static_cast<void const*>(&x) != static_cast<void const*>(&b));
      boost::fusion::for_each(normalStepBlocks,NormalStepApplyScaleAdd(normalStepAssembler,1.0,x,b,false));
      boost::fusion::for_each(tangentialStepBlocks,TangentialStepApplyScaleAdd(tangentialStepAssembler,1.0,x,b,false));
    }

    /**
     * \brief Provides access to the underlying assembler.
     */
    //    Assembler const& getAssembler() const
    //    {
    //      return assembler;
    //    }
    //
    //    bool getOnlyLowerTriangle() const { return onlyLowerTriangle; }

  protected:
    typedef typename IstlInterfaceDetail::AllBlocks<NormalStepAssembler>::result AllNormalStepBlocksSeq;
    typedef typename IstlInterfaceDetail::AllBlocks<TangentialStepAssembler>::result AllTangentialStepBlocksSeq;
    typedef typename boost::fusion::result_of::filter_if<AllNormalStepBlocksSeq,NormalStepBlockFilter>::type FilteredNormalStepBlocksSeq;
    typedef typename boost::fusion::result_of::filter_if<AllTangentialStepBlocksSeq,TangentialStepBlockFilter>::type FilteredTangentialStepBlocksSeq;

    typedef typename boost::fusion::result_of::as_vector<AllNormalStepBlocksSeq>::type AllNormalStepBlocks;
    typedef typename boost::fusion::result_of::as_vector<AllTangentialStepBlocksSeq>::type AllTangentialStepBlocks;
    typedef typename boost::fusion::result_of::as_vector<FilteredNormalStepBlocksSeq>::type FilteredNormalStepBlocks;
    typedef typename boost::fusion::result_of::as_vector<FilteredTangentialStepBlocksSeq>::type FilteredTangentialStepBlocks;

    typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<FilteredNormalStepBlocks,IstlInterfaceDetail::Translate<firstRow,firstCol>>::type>::type NormalBlocks;
    typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::transform<FilteredTangentialStepBlocks,IstlInterfaceDetail::Translate<0,0>>::type>::type TangentialBlocks;

    bool onlyLowerTriangle;
    NormalStepAssembler const& normalStepAssembler;
    TangentialStepAssembler const& tangentialStepAssembler;

    NormalBlocks  normalStepBlocks;
    TangentialBlocks tangentialStepBlocks;

    bool isSymmetric;
    int nThreads;

    struct NormalStepApplyScaleAdd
    {
      NormalStepAssembler const& assembler;
      // Constructs the functor for matrix-vector multiplication and addition. If init is true, the first block
      // encountered in each row is called with apply (using alpha=1 and initializing the rhs to zero), otherwise
      // applyscaleadd is used. This prevents a double access of the right hand side, first for setting to zero,
      // then for adding the MV product.
      NormalStepApplyScaleAdd(NormalStepAssembler const& as, Scalar a_, Domain const& x_, Range& y_, bool init, Scalar* dp_=0): assembler(as), a(a_), x(x_), y(y_), dp(dp_? dp_: &dummy)
      {
        for (int i=0; i<lastRow; ++i)
          doInit[i] = init;
        *dp = 0;
      }

      template <class Block>
      void operator()(Block const& b) const {
        if( ( (Block::rowId==stateId && Block::colId==stateId) || (Block::rowId==stateId && Block::colId==controlId) || (Block::colId==stateId && Block::rowId==controlId) ) ) return;
//        if( (Block::rowId==stateId && Block::colId==stateId) ) return;
        using namespace boost::fusion;

        typename result_of::at_c<typename Domain::Functions const,Block::colId>::type xi = at_c<Block::colId>(x.data);
        typename result_of::at_c<typename Range::Functions,Block::rowId>::type  yi = at_c<Block::rowId>(y.data);
        if (!b.threadedMatrix.get())
          abort();

        if (doInit[Block::rowId])
        {
          *dp += b.threadedMatrix->mv(xi,yi);
          doInit[Block::rowId] = false;
        }
        else
          *dp += b.threadedMatrix->usmv(a,xi,yi);

        return;
      }

    private:
      Scalar a;
      Domain const& x;
      Range& y;
      Scalar* dp;
      Scalar dummy;                 // something to point to in case no target pointer has been provided
      mutable bool doInit[lastRow]; // ugly, but boost::fusion::for_each only works with const functors
    };


    struct TangentialStepApplyScaleAdd
    {
      TangentialStepAssembler const& assembler;
      // Constructs the functor for matrix-vector multiplication and addition. If init is true, the first block
      // encountered in each row is called with apply (using alpha=1 and initializing the rhs to zero), otherwise
      // applyscaleadd is used. This prevents a double access of the right hand side, first for setting to zero,
      // then for adding the MV product.
      TangentialStepApplyScaleAdd(TangentialStepAssembler const& as, Scalar a_, Domain const& x_, Range& y_, bool init, Scalar* dp_=0): assembler(as), a(a_), x(x_), y(y_), dp(dp_? dp_: &dummy)
      {
        for (int i=0; i<lastRow; ++i)
          doInit[i] = init;
        *dp = 0;
      }

      template <class Block>
      void operator()(Block const& b) const {
        if( !( (Block::rowId==stateId && Block::colId==stateId) || (Block::rowId==stateId && Block::colId==controlId) || (Block::colId==stateId && Block::rowId==controlId) ) ) return;
        using namespace boost::fusion;

        typename result_of::at_c<typename Domain::Functions const,Block::colId>::type xi = at_c<Block::colId>(x.data);
        typename result_of::at_c<typename Range::Functions,Block::rowId>::type  yi = at_c<Block::rowId>(y.data);
        if (!b.threadedMatrix.get())
          abort();

        if (doInit[Block::rowId])
        {
          *dp += b.threadedMatrix->mv(xi,yi);
          doInit[Block::rowId] = false;
        }
        else
          *dp += b.threadedMatrix->usmv(a,xi,yi);

        return;
      }

    private:
      Scalar a;
      Domain const& x;
      Range& y;
      Scalar* dp;
      Scalar dummy;                 // something to point to in case no target pointer has been provided
      mutable bool doInit[lastRow]; // ugly, but boost::fusion::for_each only works with const functors
    };
  };

  template <class Operator>
  class TransposedOperator : public Dune::LinearOperator<typename Operator::Range,typename Operator::Domain>
  {
  public:
    typedef typename Operator::Scalar Scalar;
    typedef typename Operator::Range Domain;
    typedef typename Operator::Domain Range;
    typedef typename Operator::Assembler Assembler;

    explicit TransposedOperator(Operator const& A_) : A(A_.getTriplet().transpose()){}

    virtual ~TransposedOperator(){}

    virtual void apply(Domain const& x, Range& b) const
    {
      A.apply(x,b);
    }

    virtual void applyscaleadd(Scalar alpha, Domain const& x, Range& b) const
    {
      A.applyscaleadd(alpha,x,b);
    }

    MatrixAsTriplet<Scalar> getTriplet() const { return A.getTriplet(); }

  private:
    MatrixRepresentedOperator<MatrixAsTriplet<double>,Domain,Range> A;
  };

} // end of namespace Kaskade

#endif
