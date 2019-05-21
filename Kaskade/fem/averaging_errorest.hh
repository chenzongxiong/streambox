/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AVERAGING_ERROREST_HH
#define AVERAGING_ERROREST_HH

#include "dune/geometry/quadraturerules.hh"

#include "fem/functionviews.hh"
#include "fem/assemble.hh"
#include "fem/fetransfer.hh"
#include "fem/integration"
#include "fem/celldata.hh"
//#include "nedelecspace.hh"
//#include "fem/gridscanner.hh"
#include "linalg/factorization.hh"
#include "linalg/pardiso_solve.hh"
#include "linalg/linearsystem.hh"

namespace Kaskade
{
  /**
   *
   * @file
   * @brief  Error estimation via gradient averaging and helper routines
   * @author Anton Schiela
   *
   This file provides components for error estimation and adaptive grid refinement. 

   In particular CellData \<Grid\> is a data structure that stores error indicators 
   and provides service routine for marking and refining. To be used together with 
   GridManager \<Grid\>.

   With gradientAveraging there is a simple and rather universal error estimator 
   provided. Independently of an equation it measures the difference of a
   discontinuous gradient to a continuous interpolation, which can be interpreted
   as the distance of the discrete solution to a smoother function space. 

   * <b> Example for usage </b>
\verbatim
 do
 {
    //Assemble and solve discrete Problem...

    CellData<Grid> indicator(gradientAveraging(currentsolution));
    double totalH1Error=std::sqrt(indicator.sum());
    gm.mark(markByBulkCriterion(indicator,theta));
    gm.adaptAtOnce();
 } while(totalH1Error > desiredAccuracy)
\endverbatim
   */

  /** \addtogroup adapt */
  /**@{*/
  namespace FunctionViews
  {
    /** View on FunctionSpaceElements. If .value() is called, then the
     * result is the difference of the two FunctionSpaceElements from
     * which this view was constructed. It is assumed that both functions locally
     * live on the same space. This means that they have the same shapefunctions on
     * each cell, but they may differ in their continuity properties
     *
     * For more information on the construction of such views cf. FunctionViews
     */
    template <typename Function1, typename Function2>
    class Difference
    {
    public:
      typedef typename Function1::Space Space;
      typedef typename Function1::ValueType ValueType;
      typedef typename Function1::GradientType GradientType;
      typedef typename Function1::Scalar Scalar;

      template<class SFS>
      int order(SFS const& sfs) const {return std::max(f1.order(sfs),f2.order(sfs)); }

      Space const& space() const {return f1.space();}

      Difference(Function1 const& f1_, Function2 const& f2_) : f1(f1_), f2(f2_) {}

      ValueType value(typename Function1::Space::Evaluator eval) const
      {
        return f1.value(eval)-f2.value(eval);
      }

      template<class EPtr, class V>
      ValueType value(EPtr const& ci, V const& v) const
      {
        return f1.value(ci,v)-f2.value(ci,v);
      }

      GradientType gradient(typename Function1::Space::Evaluator eval) const
      {
        return f1.gradient(eval)-f2.gradient(eval);
      }

    private:
      Function1 const& f1;
      Function2 const& f2;
    };

  }


  /// Create a CellData by computing local integrals over each Cell
  /** Class for creating a vector of pair<ValueType value, EntityPointer ptr>
   * from a FunctionSpaceElement with the following semantics:
   *
   * value = Integral of the FunctionSpaceElement on the entity *p.
   *
   * The integral is computed approximately by quadrature formulas.
   * For each leaf entity this vector will have exactly one entry. Possible
   * application: For a given pointwise error estimate we can compute a
   * cell-wise error indicator.
   */
  template<typename Space>
  class LocalIntegral
  {
    typedef typename Space::Grid Grid;
    typedef typename Space::template Element<1>::type Function;
    typedef typename Function::ValueType ValueType;

    class LocalIntegrator
    {
      typedef typename Grid::template Codim<0>::Entity Cell;
      typedef Dune::QuadraturePoint<typename Grid::ctype,Grid::dimension> QuadPoint;
    
    public:
      void operator()(typename Cell::EntityPointer const& ci, QuadPoint const& q, ValueType v)
      {
        v *= q.weight()*ci->geometry().integrationElement(q.position());
        if(localIntegral.empty() || !(localIntegral[localIntegral.size()-1].second == ci))
        {
          localIntegral.push_back(std::make_pair(v,ci));
        }
        else
        {
          localIntegral[localIntegral.size()-1].first += v;
        }
      }
      
      typename CellData<Grid, ValueType>::CellDataVector& getLocalIntegral() { assert(&localIntegral); return localIntegral; }
    
    private:
      typename CellData<Grid, ValueType>::CellDataVector localIntegral;
    };

  public:
    /// Perform local integration for WeakFunctionViews (less efficient)
    template<typename Function>
    typename CellData<Grid, ValueType>::CellDataVector operator()(Function f, Space s)
    {
      LocalIntegrator integrator;
      forEachQuadPoint(f,integrator,s);
      return integrator.getLocalIntegral();
    }

    /// Perform local integration for FunctionSpaceElement or FunctionViews (more efficient)
    template<typename Function>
    typename CellData<Grid, ValueType>::CellDataVector operator()(Function f)
    {
      LocalIntegrator integrator;
      forEachQuadPoint(f,integrator);
      return integrator.getLocalIntegral();
    }
  };

  /// Perform gradient averaging in order to construct an error estimator
  /** The result is a CellDataVector, from which an object of the class
    CellData can be constructed.  With the functionality of this class
    an adaptive refinement procedure can be driven easily */
  template<class Function>
  typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector gradientAveraging(Function const& f)
  {
    using namespace FunctionViews;
    typedef typename Function::Space Space;
    const int dim(Space::dim);
    typename Space::template Element<dim>::type smoothGradient(f.space());
    interpolateGloballyWeak<InverseVolume>(smoothGradient,makeView<Gradient>(f));
    LocalIntegral<Space> localIntegral;
    typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector
    errorIndicator(localIntegral(
        makeView<AbsSquare>(
            makeView<Difference>(smoothGradient,makeView<Gradient>(f))),f.space()));
    return errorIndicator;
  }

  /// Perform gradient averaging in order to construct a representation of the gradient error
  /** The result is a CellDataVector, from which an object of the class
    CellData can be constructed.  With the functionality of this class
    an adaptive refinement procedure can be driven easily */
  template<class GradientError, class Function>
  void gradientAveraging(GradientError & dx, Function const& x)
  {
    using namespace FunctionViews;
    const int dim(Function::Space::dim);
    typename Function::Space::template Element<dim>::type smoothGradient(x.space());
    interpolateGlobally<InverseVolume>(smoothGradient,makeView<Gradient>(x));
    interpolateGlobally<PlainAverage>(dx,makeView<Difference>(smoothGradient,makeView<Gradient>(x)));
  }

  /// Perform lower order interpolation in order to construct an error estimator
  /** The result is a CellDataVector, from which an object of the class
    CellData can be constructed.  With the functionality of this class
    an adaptive refinement procedure can be driven easily */
  template<class Function>
  typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector lowOrderInterpolation(Function const& f)
  {
    using namespace FunctionViews;
    typedef typename Function::Space Space;
    Space lowOrderSpace(f.space().gridManager(), f.space().indexSet(), f.space().mapper().getOrder()-1);
    typename Space::template Element<1>::type lowOrderInterpolant(lowOrderSpace);
    interpolateGlobally<PlainAverage>(lowOrderInterpolant,f);
    LocalIntegral<Space> localIntegral;
    typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector
    errorIndicator(localIntegral(
        makeView<AbsSquare>(
            makeView<Difference>(lowOrderInterpolant,f)
        ),f.space()
    )
    );
    return errorIndicator;
  }

  /// Construct a higher order estimate for the finite element solution from gradient information
  /** Using the gradient information from an averaging estimator
   * a higher order finite element function is constructed via a least squares fit
   *
   */
  template<class Grid, class Space>
  class HigherOrderRecovery
  {

    typedef FEFunctionSpace<DiscontinuousLagrangeMapper<double,Grid> > DiscontSpace;
    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid> > ContSpace;

  public:

    typedef DiscontSpace ErrorGradientFunctionSpace;
    typedef Space ErrorFunctionSpace;

    typedef typename ErrorGradientFunctionSpace::template Element<Grid::dimension>::type ErrorGradientFunction;
    typedef typename ErrorFunctionSpace::template Element<1>::type ErrorFunction;

  private:

    typedef boost::fusion::vector<const ErrorFunctionSpace*,const ErrorGradientFunctionSpace*> SpacesForAssembly;
    typedef boost::fusion::vector<VariableDescription<0,1,0> > ErrorFunctionVariables;

    typedef VariableSetDescription<SpacesForAssembly,ErrorFunctionVariables> VariableDefinition;

    template <class RType, class Variables>
    class HigherOrderRecoveryFunctional
    {
    public:
      typedef RType RT;
      typedef Variables TestVars;
      typedef Variables AnsatzVars;
      static ProblemType const type = VariationalFunctional;

    public:

      template<class Arg>
      bool inDomain(Arg const& a) const { return true; }

      HigherOrderRecoveryFunctional(ErrorGradientFunction const& x_) : x(x_){}

      struct DomainCache : public EvalCacheBase
      {
        DomainCache(ErrorGradientFunction const& u_) :  gradientError(u_)
        {
        };

        template<class Position, class Evaluators>
        void evaluateAt(Position const& x, Evaluators const& evaluators)
        {
          using namespace boost::fusion;
          du = gradientError.value(at_c<1>(evaluators));
        }

        RT d0() const { return 0; }

        template <int row, int dim> Dune::FieldVector<RT,1> d1(VariationalArg<RT,dim> const& arg) const { return du*arg.gradient[0]; }

        template <int row, int col, int dim>
        Dune::FieldMatrix<RT,1,1> d2(VariationalArg<RT,dim> const& arg1, VariationalArg<RT,dim> const& arg2) const
        {
          return arg1.gradient[0]*arg2.gradient[0] + arg1.value*arg2.value;
        }

      private:
        ErrorGradientFunction gradientError;
        Dune::FieldVector<RT,Grid::dimension> du;
      };

      struct BoundaryCache : public HomNeumannBoundaryCache<RT>
      {
      };

      DomainCache createDomainCache(int flags=7) const {return DomainCache(x);}
      BoundaryCache createBoundaryCache(int flags=7) const {return BoundaryCache();}

      ErrorGradientFunction const& x;

      public:
      template <int row, int col>
      struct D2
      {
        static bool const present = true;
        static bool const lumped = false;
        static bool const symmetric = true;
      };

      template <class Cell>
      int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
      {
        int matrixOrder =  boundary? 2*shapeFunctionOrder: 2*(shapeFunctionOrder-1);
        int lastIterateOrder = 0;

        return matrixOrder+lastIterateOrder;
      }
    };

  public:

    void getErrorFunction(ErrorFunction& ex, ErrorGradientFunction const& gex)
    {

      SpacesForAssembly spaces(&ex.space(),
          &gex.space());
      // construct variable list.
      VariableDefinition variableSet(spaces);

      // construct variational functional.
      typedef HigherOrderRecoveryFunctional<double,VariableDefinition> Functional;

      Functional HORF(gex);

      // construct Galerkin representation
      VariationalFunctionalAssembler<Functional> gop(spaces);

      size_t nnz = gop.nnz(0,1,0,1,false);
      size_t size = variableSet.degreesOfFreedom(0,1);

      std::vector<int> ridx(nnz), cidx(nnz);
      std::vector<double> data(nnz), rhs(size), sol(size);

      gop.assemble(HORF,6);
      gop.toTriplet(0,1,0,1,
          ridx.begin(), cidx.begin(),
          data.begin(),false);
      gop.toSequence(0,1,rhs.begin());

      MatrixAsTriplet<double> m;

      m.ridx=ridx; m.cidx=cidx; m.data=data;

      m.deleteLowerTriangle();

      PARDISOFactorization<double> pf(size,2,m.ridx,m.cidx,m.data,MatrixProperties::GENERAL);

      pf.solve(rhs,sol);
      assert((*ex).size()==size);
      for(int i=0; i<size;++i) (*ex)[i]=sol[i];
    }
  };
} /* end of namespace Kaskade */
/**
@}
 */

#endif


