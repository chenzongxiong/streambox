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

#ifndef HIEARARCHIC_ERROR_ESTIMATOR_HH
#define HIEARARCHIC_ERROR_ESTIMATOR_HH

#include "fem/variables.hh"
/**
 *
 * @file
 * @brief  Error estimation via hierachic FEM
 * @author Martin Weiser
 *
   This file provides components for error estimation and adaptive grid refinement. 
 */

namespace Kaskade
{
  /** \addtogroup adapt */
  /**@{*/

  /**
   * \cond internals
   */
  namespace HierarchicErrorEstimatorDetail {

    /**
     * A class that models the HierarchicErrorEstimatorD2Info concept and
     * can be used as the fourth template parameter to
     * HierarchicErrorEstimator.
     *
     * It relates the block info to its template parameter, but sets the
     * lumped flag to true and takes only those blocks of Functional::D2
     * as present which are symmetric.
     */
    template <class Functional>
    struct DefaultD2
    {
      template <int row, int col>
      class D2: public Functional::template D2<row,col>
      {
        typedef typename Functional::template D2<row,col> d2;

      public:
        static bool const present = d2::present && d2::symmetric;
        static bool const lumped = true;
      };
    };

    template <class Functional>
    struct LumpSymmetricD2
    {
      template <int row, int col>
      class D2: public Functional::template D2<row,col>
      {
      public:
        static constexpr bool lumped = Functional::template D2<row,col>::symmetric;
      };
    };

    /**
     * A class that models the HierarchicErrorEstimatorD2Info concept and
     * can be used as the fourth template parameter to
     * HierarchicErrorEstimator.
     *
     * It relates the block info to its template parameter, without
     * changing anything. This is less useful for actual error estimation,
     * but may help in analyzing algorithms and errors.
     */
    template <class Functional>
    struct TakeAllD2
    {
      template <int row, int col>
      class D2: public Functional::template D2<row,col> {};
    };


  } // End of namespace HierarchicErrorEstimatorDetail

  /**
   * \endcond
   */



  /**
   * \brief Defines assembly of hierarchically extended problems for defining DLY style error estimators.
   * 
   * This class defines the weak formulation of hierarchical error
   * estimation. It is a model of the LinearVariationalProblem concept.
   *
   * Assume you have discretized a linear(ized) variational functional
   * using a Galerkin ansatz space \f$ V_l \f$. The Hessian is \f$
   * A_{ll} \f$ and the gradient at \f$x_0=0 \f$ is \f$ b_l \f$, such
   * that the solution (minimizer) \f$ x_l \f$ satisfies \f$ A_{ll}x_l =
   * -b_l \f$.
   *
   * Now we extend the ansatz space by another (probably higher order)
   * ansatz space \f$ V_h \f$. The complete system to be solved would be
   *
   * \f[ \begin{bmatrix} A_{ll} & A_{lh} \\ A_{hl} & A_{hh} \end{bmatrix}
   *     \begin{bmatrix} x_l \\ x_h \end{bmatrix} =
   *     - \begin{bmatrix} b_l \\ b_h \end{bmatrix}.
   * \f]
   *
   * Since the assembly and solution of the large system is quite
   * expensive, we fix \f$ x_l \f$ obtained before and solve the defect equation
   *
   * \f[ A_{hh} x_h = -(b_h+A_{hl}x_l) \f]
   *
   * instead. And since the error \f$ x_h \f$ is structurally
   * oscillatory, one Jacobi iteration will do for error estimation
   * purposes, such that we only need to assemble the diagonal of \f$
   * A_{hh} \f$ (lumping).
   *
   * Important things to note:
   *
   * - The assembly produces \f$ A_{hh} \f$ as matrix and \f$
   *   b_h+A_{hl}x_l \f$ as right hand side. Make sure to provide \f$ x_l
   *   \f$ with the correct sign!
   *
   * - During the assembly, both original and extension ansatz functions need to be evaluated.
   *   Remember to provide all(!) required spaces. The easiest method to do so is to extend the
   *   space list of the original problem by the extension spaces, and to shift the space index
   *   of the extension variables accordingly.
   *
   *
   * \tparam LinearFunctional         the linear variational weak formulation for which the
   *                                  error will be estimated
   * \tparam ExtensionAnsatzVariables a description of the ansatz variables and hierarchic extension spaces
   * \tparam ExtensionTestVariables   a description of the test variables and hierarchic extension spaces
   * \tparam D2BlockInfo              proivdes static block information
   *
   * D2BlockInfo defines which blocks should be assembled and which ones
   * should be lumped (only diagonal assembled). The default is
   * HierarchicErrorEstimatorDetail::DefaultD2, which takes only the
   * symmetric blocks and lumps all of them.
   */
  template <class LinearFunctional,
            class ExtensionAnsatzVariables,
            class ExtensionTestVariables=ExtensionAnsatzVariables,
            class D2BlockInfo = HierarchicErrorEstimatorDetail::DefaultD2<LinearFunctional> >
  class HierarchicErrorEstimator
  {
    typedef HierarchicErrorEstimator<LinearFunctional,ExtensionAnsatzVariables,ExtensionTestVariables> Self;
    typedef typename LinearFunctional::AnsatzVars FEAnsatzVars;
    typedef typename ExtensionAnsatzVariables::Grid Grid;
    typedef typename Grid::template Codim<0>::Entity Cell;

  public:
    static ProblemType const type = LinearFunctional::type;
    typedef typename LinearFunctional::Scalar Scalar;
    typedef ExtensionAnsatzVariables AnsatzVars;
    typedef ExtensionTestVariables   TestVars;

    typedef typename FEAnsatzVars::VariableSet FESolution;

    /**
     * \cond internals
     */
    class DomainCache
    {
    public:
      DomainCache(LinearFunctional const& lf, FESolution const& ul_, int flags):
        ul(ul_), dc(lf.createDomainCache(flags))
      {}

      void moveTo(Cell const& e) { dc.moveTo(e); }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        dc.evaluateAt(x,evaluators);
        ulEval = evaluateVariables(ul,evaluators,valueMethod);
        ulGradEval = evaluateVariables(ul,evaluators,derivativeMethod);
      }

      Scalar d0() const { return dc.d0(); }

      template <int row, int dim>
      Dune::FieldVector<Scalar,TestVars::template Components<row>::m>
      d1(VariationalArg<Scalar,dim> const& arg) const
      {
        Dune::FieldVector<Scalar,TestVars::template Components<row>::m> result = dc.template d1<row>(arg);

        boost::fusion::for_each(typename FEAnsatzVars::Variables(),
            addLowerOrderD2<row,dim>(dc,ulEval,ulGradEval,arg,result));

        return  result;
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,
      AnsatzVars::template Components<col>::m>
      d2(VariationalArg<Scalar,dim> const& arg1, VariationalArg<Scalar,dim> const& arg2) const
      {
        return dc.template d2<row,col>(arg1,arg2);
      }

    private:
      FESolution const&                      ul;
      typename LinearFunctional::DomainCache dc;
      EvaluateVariables<FESolution,ValueMethod> ulEval;
      EvaluateVariables<FESolution,DerivativeMethod> ulGradEval;

      template <int row, int dim>
      struct addLowerOrderD2
      {
        addLowerOrderD2(typename LinearFunctional::DomainCache const& dc_,
                        EvaluateVariables<typename FEAnsatzVars::VariableSet,ValueMethod> const& ulEval_,
                        EvaluateVariables<typename FEAnsatzVars::VariableSet,DerivativeMethod> const& ulGradEval_,
                        VariationalArg<Scalar,dim> const& arg_,
                        Dune::FieldVector<Scalar,TestVars::template Components<row>::m>& result_):
              dc(dc_), ulEval(ulEval_), ulGradEval(ulGradEval_), arg(arg_), result(result_)
        {}

        template <class VarDescription>
        void operator()(VarDescription const& vd) const {
          int const col = VarDescription::id;
          bool const mirror = row<col && LinearFunctional::type==VariationalFunctional;
          bool const original = row>=col || LinearFunctional::type==WeakFormulation;


          // Check whether this block is present in the first
          // place. Note that it is important to use D2 from
          // LinearFunctional instead of our own, because for the error
          // estimation matrix we may deliberately drop blocks!
          //
          // If we have a symmetric problem, we only access subdiagonal
          // blocks and mirror them to the upper triangular part.
          if (! ( (original && LinearFunctional::template D2<row,col>::present) ||
              (mirror && LinearFunctional::template D2<col,row>::present) ) )
            return;


          typedef typename SpaceType<typename FEAnsatzVars::Spaces,
          VarDescription::spaceIndex>::type FEColSpace;

          VariationalArg<Scalar,dim> uArg;
          // Step through components of variable
          for (int c=0; c<VarDescription::m; ++c) {
            int off = c*FEColSpace::sfComponents;

            // Fill variational arg with values from that
            // component. Handle vectorial shape functions correctly.
            for (int i=0; i<FEColSpace::sfComponents; ++i) {
              uArg.value[i] = boost::fusion::at_c<col>(ulEval)[off+i];

              for (int j=0; j<dim; ++j)
              {
                uArg.derivative[i][j] = boost::fusion::at_c<col>(ulGradEval)[off+i][j];
              }
            }

            // Compute local matrix. Note that since we're working
            // currently on the c-th component, only the c-th column is
            // used.
            if (original && LinearFunctional::template D2<row,col>::present) {
              Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,
              AnsatzVars::template Components<col>::m>
              d2 = dc.template d2<row,col>(arg,uArg);

              for (int i=0; i<TestVars::template Components<row>::m; ++i)
              {
                result[i] += d2[i][c];
              }
            }

            // Apply transposed mirror-blocks if available.
            if (mirror && LinearFunctional::template D2<col,row>::present) {
              Dune::FieldMatrix<Scalar,AnsatzVars::template Components<col>::m,
              TestVars::template Components<row>::m>
              d2 = dc.template d2<col,row>(uArg,arg);

              for (int i=0; i<TestVars::template Components<row>::m; ++i)
              {
                result[i] += d2[c][i];
              }
            }
          }
        }

      private:
        typename LinearFunctional::DomainCache const&                               dc;
        EvaluateVariables<typename FEAnsatzVars::VariableSet,ValueMethod>           ulEval;
        EvaluateVariables<typename FEAnsatzVars::VariableSet,DerivativeMethod> const& ulGradEval;
        VariationalArg<Scalar,dim> const&                                           arg;
        Dune::FieldVector<Scalar,TestVars::template Components<row>::m>&            result;
      };
    };

    class BoundaryCache
    {
    public:
      // static const bool hasInteriorFaces = LinearFunctional::BoundaryCache::hasInteriorFaces;
      typedef typename Grid::LeafIntersectionIterator FaceIterator;

      BoundaryCache(LinearFunctional const& lf, FESolution const& ul_, int flags):
        ul(ul_), bc(lf.createBoundaryCache(flags))
      {}

      void moveTo(FaceIterator const& entity) { bc.moveTo(entity); }

      template <class Evaluators>
      void evaluateAt(Dune::FieldVector<typename Grid::ctype,Grid::dimension-1> const& x, Evaluators const& evaluators)
      {
        bc.evaluateAt(x,evaluators);
        ulEval = evaluateVariables(ul,evaluators,valueMethod);
      }

      Scalar d0() const { return bc.d0(); }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
      d1 (VariationalArg<Scalar,dim> const& arg) const
      {
        Dune::FieldVector<Scalar,TestVars::template Components<row>::m> result = bc.template d1<row>(arg);
        boost::fusion::for_each(typename FEAnsatzVars::Variables(),
            addLowerOrderD2<row,dim>(bc,ulEval,arg,result));
        return  result;
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const
      {
        return bc.template d2<row,col>(arg1,arg2);
      }

    private:
      FESolution const&                            ul;
      typename LinearFunctional::BoundaryCache     bc;
      EvaluateVariables<typename FEAnsatzVars::VariableSet,ValueMethod> ulEval;

      template <int row, int dim>
      struct addLowerOrderD2
      {
        addLowerOrderD2(typename LinearFunctional::BoundaryCache const& bc_,
                        EvaluateVariables<typename FEAnsatzVars::VariableSet,ValueMethod> const& ulEval_,
                        VariationalArg<Scalar,dim> const& arg_,
                        Dune::FieldVector<Scalar,TestVars::template Components<row>::m>& result_):
              bc(bc_), ulEval(ulEval_), arg(arg_), result(result_)
        {}

        template <class VarDescription>
        void operator()(VarDescription const& vd) const {
          int const col = VarDescription::id;
          bool const mirror = row<col && LinearFunctional::type==VariationalFunctional;
          bool const original = row>=col || LinearFunctional::type==WeakFormulation;

          // Check whether this block is present in the first
          // place. Note that it is important to use D2 from
          // LinearFunctional instead of our own, because for the error
          // estimation matrix we may deliberately drop blocks!
          //
          // If we have a symmetric problem, we only access subdiagonal
          // blocks and mirror them to the upper triangular part.
          if (! ( (original && LinearFunctional::template D2<row,col>::present) ||
              (mirror && LinearFunctional::template D2<col,row>::present) ) )
            return;

          typedef typename SpaceType<typename FEAnsatzVars::Spaces,
          VarDescription::spaceIndex>::type FEColSpace;

          VariationalArg<Scalar,dim> uArg;
          // Step through components of variable
          for (int c=0; c<VarDescription::m; ++c) {
            int off = c*FEColSpace::sfComponents;

            // Fill variational arg with values from that
            // component. Handle vectorial shape functions correctly.
            for (int i=0; i<FEColSpace::sfComponents; ++i)
              uArg.value[i] = boost::fusion::at_c<col>(ulEval)[off+i];

            // Compute local matrix. Note that since we're working
            // currently on the c-th component, only the c-th column is
            // used.
            if (original && LinearFunctional::template D2<row,col>::present) {
              Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,
              AnsatzVars::template Components<col>::m>
              d2 = bc.template d2<row,col>(arg,uArg);

              for (int i=0; i<TestVars::template Components<row>::m; ++i)
                result[i] += d2[i][c];
            }

            // Apply transposed mirror-blocks if available.
            if (mirror && LinearFunctional::template D2<col,row>::present) {
              Dune::FieldMatrix<Scalar,AnsatzVars::template Components<col>::m,
              TestVars::template Components<row>::m>
              d2 = bc.template d2<col,row>(uArg,arg);

              for (int i=0; i<TestVars::template Components<row>::m; ++i)
                result[i] += d2[c][i];
            }
          }
        }

      private:
        typename LinearFunctional::BoundaryCache const&                          bc;
        EvaluateVariables<typename FEAnsatzVars::VariableSet,ValueMethod> const& ulEval;
        VariationalArg<Scalar,dim> const&                                        arg;
        Dune::FieldVector<Scalar,TestVars::template Components<row>::m>&         result;
      };
    };

    template <int row>
    struct D1: public LinearFunctional::template D1<row>
    {};

    template <int row, int col>
    struct D2: public D2BlockInfo::template D2<row,col> {
      // Make sure that we only claim availability of blocks we can actually provide...
      static bool const present = D2BlockInfo::template D2<row,col>::present && LinearFunctional::template D2<row,col>::present;
    };

    /**
     * \endcond
     */

    /**
     * \param lf_ the linear weak formulation for which the error will
     * be estimated. This object can be temporary, since a copy is held
     * internally. Thus, it should not be too large.
     *
     * \param ul_ the approximate solution of lf_ on the lower order
     * ansatz space.
     */
    HierarchicErrorEstimator(LinearFunctional const& lf_, FESolution const& ul_):  lf(lf_), ul(ul_) {}


    DomainCache createDomainCache(int flags) const
    {
      return DomainCache(lf,ul,flags);
    }

    BoundaryCache createBoundaryCache(int flags) const
    {
      return BoundaryCache(lf,ul,flags);
    }

    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return lf.integrationOrder(cell,shapeFunctionOrder,boundary);
    }


    private:
    LinearFunctional  lf;
    FESolution const& ul;
  };

  //---------------------------------------------------------------------

  template <class ExtensionVariables, class LinearFunctional>
  HierarchicErrorEstimator<LinearFunctional,ExtensionVariables>
  hierarchicErrorEstimator(LinearFunctional const& lf, typename LinearFunctional::AnsatzVars::VariableSet const& u)
  {
    return HierarchicErrorEstimator<LinearFunctional,ExtensionVariables>(lf,u);
  }


  /**
@}
   */
} // namespace Kaskade
#endif
