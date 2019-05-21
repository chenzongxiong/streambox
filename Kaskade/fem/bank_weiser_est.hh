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

#ifndef BANK_WEISER_EST_HH
#define BANK_WEISER_EST_HH

#include "dune/grid/common/quadraturerules.hh"

#include "fem/functionviews.hh"
#include "fem/assemble.hh"
#include "fem/fetransfer.hh"
#include "fem/celldata.hh"
#include "fem/gridscanner.hh"
#include "linalg/pardiso_solve.hh"
#include "linalg/linearsystem.hh"

namespace Kaskade
{
  /**
   *
   * @file
   * @brief  Error estimation via bank-weiser type approach, for optimal control problems only
   * @author Anton Schiela
   *
   This file provides components for error estimation and adaptive grid refinement. 
   */
  /** \addtogroup adapt */
  /**@{*/
  template <int Id, int PDEId, int SpaceId, int components>
  struct EstimatorDescription
  {
    static int const id = Id;
    static int const pdeid = PDEId;
    static int const spaceIndex = SpaceId;
    static int const m = components;
  };


  // Construction of an error representation function...

  static int const indQ[10] = {0,1,2,3,4,5,6,7,8,9};
  static int const indQFull[10] = {0,1,2,3,4,5,6,7,8,9};

  static bool normalize = true;

  //static int const indQ[6] = {1,3,4,6,7,8};

  //template<int dim>
  //struct SystDim
  //{
  //  static const int value = 3*(dim-1);
  //};

  template<int dim>
  struct SystDim
  {
    static const int value = (dim+1)+3*(dim-1);
  };

  template<int dim>
  struct SystDimFull
  {
    static const int value = (dim+1)+3*(dim-1);
  };

  template <class LocalInd, class RT, class Cache, int dim,class Evaluators>
  struct HalfGradientJump
  {

    HalfGradientJump(LocalInd& localInd_,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        Cache const& cacheNeigh_,
        Evaluators const& evaluators_):
          localInd(localInd_),  cache(cache_), cacheNeigh(cacheNeigh_), evaluators(evaluators_)
    {
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      for(int row=0; row< localInd[Variable::id].size; ++row)
        localInd[Variable::id][row] += 0.5*(cacheNeigh.template d1<Variable::pdeid>(modarg)-cache.template d1<Variable::pdeid>(modarg))
        *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].value[0];
    }

    VariationalArg<RT,dim> modarg;
    LocalInd&         localInd;
    Cache const&      cache;
    Cache const&      cacheNeigh;
    Evaluators const & evaluators;
  };


  template <class LocalInd, class RT, class Cache, class BoundaryCache, int dim, class Evaluators>
  struct StrongBoundaryValues
  {
    StrongBoundaryValues(LocalInd& localInd_,
        RT& integrationElement,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        BoundaryCache const& cacheBoundary_,
        Evaluators const& evaluators_):
          localInd(localInd_), cache(cache_), boundaryCache(cacheBoundary_),evaluators(evaluators_)
    {
      modargB.value=integrationElement;
      modargB.gradient[0]=0;
      modargTest.value=1;
      modargTest.gradient[0]=0;
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      if(boundaryCache.template d2<Variable::pdeid,Variable::id>(modargTest,modargTest) < 1e5)
        for(int row=0; row< localInd[Variable::id].size; ++row)
        {
          localInd[Variable::id][row] -=
              (cache.template d1<Variable::pdeid>(modarg)+boundaryCache.template d1<Variable::pdeid>(modargB))
              *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].value[0];
        }
    }

    VariationalArg<RT,dim> modarg,modargB,modargTest;
    LocalInd&         localInd;
    Cache const& cache;
    BoundaryCache const&   boundaryCache;
    Evaluators const& evaluators;
  };

  template <class LocalRHS, class RT, class Cache, class Evaluators>
  struct WeakResiduum
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;
    WeakResiduum(LocalRHS& localRHS_, RT intEl_,Cache const& cache_, Evaluators const& evaluators_):
      localRHS(localRHS_), intEl(intEl_),  evaluators(evaluators_), cache(cache_)
    {
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      size_t n = (localRHS[Variable::id]).N();
      for(int row=0; row<n; ++row)
      {
        (localRHS[Variable::id])[row]
                                 -= intEl * cache.template d1<Variable::pdeid>(boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[row]);
      }
    }

    LocalRHS&   localRHS;
    RT intEl;
    Evaluators const& evaluators;
    Cache const& cache;
  };


  template <class LocalRHS, class RT, class Cache, class Evaluators>
  struct WeakResiduumMainPart
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;
    WeakResiduumMainPart(LocalRHS& localRHS_, RT intEl_,Cache const& cache_, Evaluators const& evaluators_):
      localRHS(localRHS_), intEl(intEl_),  evaluators(evaluators_), cache(cache_)
    {
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      VariationalArg<RT,dim> modarg1;
      modarg1.value[0]=0;
      size_t n = (localRHS[Variable::id]).N();
      for(int row=0; row<n; ++row)
      {
        modarg1.gradient[0]=boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].gradient[0];
        (localRHS[Variable::id])[row] -= intEl * cache.template d1<Variable::pdeid>(modarg1);
      }
    }

    LocalRHS&   localRHS;
    RT intEl;
    Evaluators const& evaluators;
    Cache const& cache;
  };


  template <class LocalInd, class RT, class Cache, int dim,class Evaluators>
  struct GradientAverage
  {

    GradientAverage(LocalInd& localInd_,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        Cache const& cacheNeigh_,
        Evaluators const& evaluators_):
          localInd(localInd_),  cache(cache_), cacheNeigh(cacheNeigh_), evaluators(evaluators_)
    {
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      for(int row=0; row< localInd[Variable::id].size; ++row)
        localInd[Variable::id][row] += 0.5*(cache.template d1<Variable::pdeid>(modarg)+cacheNeigh.template d1<Variable::pdeid>(modarg))
        *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].value[0];
    }

    VariationalArg<RT,dim> modarg;
    LocalInd&         localInd;
    Cache const&      cache;
    Cache const&      cacheNeigh;
    Evaluators const & evaluators;
  };



  template <class LocalInd, class RT, class Cache, class BoundaryCache, int dim, class Evaluators>
  struct WeakBoundaryValues
  {
    WeakBoundaryValues(LocalInd& localInd_,
        RT& integrationElement,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        BoundaryCache const& cacheBoundary_,
        Evaluators const& evaluators_):
          localInd(localInd_), cache(cache_), boundaryCache(cacheBoundary_),evaluators(evaluators_)
    {
      modargB.value=integrationElement;
      modargB.gradient[0]=0;
      modargTest.value=1;
      modargTest.gradient[0]=0;
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      if(boundaryCache.template d2<Variable::pdeid,Variable::id>(modargTest,modargTest) < 1e5)
        for(int row=0; row< localInd[Variable::id].size; ++row)
        {
          localInd[Variable::id][row] +=
              boundaryCache.template d1<Variable::pdeid>(modargB)
              *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].value[0];
        }
      else
        for(int row=0; row< localInd[Variable::id].size; ++row)
        {
          localInd[Variable::id][row] +=
              cache.template d1<Variable::pdeid>(modarg)
              *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQ[row]].value[0];
        }
    }

    VariationalArg<RT,dim> modarg,modargB,modargTest;
    LocalInd&         localInd;
    Cache const& cache;
    BoundaryCache const&   boundaryCache;
    Evaluators const& evaluators;
  };

  template <class LocalVector, class LocalMatrix>
  struct SolveLocalSystem
  {
    SolveLocalSystem(LocalVector& c_,
        LocalMatrix & A_,
        LocalVector& b_):
          c(c_), A(A_), b(b_)
    {
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      //Since Neumann problems are solved, achieve sum_i b_i = 0
      if(normalize)
      {
        int n(b[Variable::id].N());
        double sum(0.0);
        for(int i=0; i<n; ++i) sum+=b[Variable::id][i];
        sum /= n;
        for(int i=0; i<n; ++i) b[Variable::id][i]-=sum;
      }
      try
      {
        A.solve(c[Variable::id],b[Variable::id]);
      }
      catch(...)
      {
        std::cout << "Warning: Factorization failed: setting indicator to zero!" << std::endl;
        for(int i =0; i<c[Variable::id].N();++i)
          c[Variable::id][i]=0.0;
      }
    }

    LocalVector & c;
    LocalMatrix & A;
    LocalVector & b;
  };


  template <class LocalMatrix, class Evaluators, class RT, class Cache>
  void localNormalMatrix(LocalMatrix& localmatrix, RT intEl, Evaluators const& evaluators,Cache const& cache)
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;
    VariationalArg<RT,dim> modarg1,modarg2;
    modarg1.value[0]=0;
    modarg2.value[0]=0;
    size_t n = (localmatrix[0]).N();
    for(int row=0; row<n; ++row)
      for (size_t i=0; i<n; ++i)
      {
        modarg1.gradient[0]=boost::fusion::at_c<0>(evaluators).evalData[indQ[row]].gradient[0];
        modarg2.gradient[0]=boost::fusion::at_c<0>(evaluators).evalData[indQ[i]].gradient[0];

        localmatrix[row][i]
                         += 1e-3*intEl*boost::fusion::at_c<0>(evaluators).evalData[indQ[row]].value[0]*
                         boost::fusion::at_c<0>(evaluators).evalData[indQ[i]].value[0];
        localmatrix[row][i]+= intEl*cache.template d2<0,1>(modarg1,modarg2);
      }
  }




  template <class LocalVector, class Evaluators, class RT>
  void rank1KernelVector(LocalVector& localvector, RT intEl, Evaluators const& evaluators)
  {
    size_t n = (localvector).N();
    for (size_t i=0; i<n; ++i) 
    {
      (localvector)[i]
                    += intEl*boost::fusion::at_c<0>(evaluators).evalData[indQ[i]].value[0];
    }
  }

  template <class LocalSolution, class Evaluators, class GlobalSolution>
  struct ScatterLocalSolution
  {
    ScatterLocalSolution(LocalSolution const& localSolution_, Evaluators const& evaluators_, GlobalSolution& globalSolution_):
      localSolution(localSolution_), eval(evaluators_), globalSolution(globalSolution_)
    {}

    template <class Variable>
    void operator()(Variable const& ) const
    {
      using namespace boost::fusion;
      typename result_of::value_at_c<Evaluators,Variable::spaceIndex>::type const& reval = at_c<Variable::spaceIndex>(eval);
      assert(indQ[localSolution[Variable::id].N()-1]<reval.size());
      for (size_t i=0; i<localSolution[Variable::id].N(); i++)
        (*at_c<Variable::id>(globalSolution))[reval.globalIndices()[indQ[i]]] += (localSolution)[Variable::id][i];
    }

    LocalSolution const&   localSolution;
    Evaluators const& eval;
    GlobalSolution&        globalSolution;
  };

  /// Computes an error representation function, which is useful for goal oriented adaptivity
  /** Arguments
   *  - result : output: variable set of error functions
   *  - varDesc : description of the testfunctions
   *  - f  : functional, the residual of which is evaluated
   *  - spaces : spaces of the ansatz functions x, neede for the evaluation of f
   *  - extended spaces: spaces, where the error functions live
   *
   *  \attention: this routine is not yet fully general
   */
  template <class VariableDescriptions, int sysdim, class Result, class F, class Spaces, class ExtendedSpaces>
  void edgeAveraging(Result& result,
      F const& f,
      Spaces const& spaces,
      ExtendedSpaces const& extendedSpaces)
  {
    using namespace boost::fusion;
    using namespace Dune;

    VariableDescriptions const varDesc;

    typedef typename SpaceType<Spaces,0>::type::Grid Grid;
    typedef typename SpaceType<Spaces,0>::type::IndexSet IndexSet;
    typedef typename SpaceType<Spaces,0>::type::RT RT;
    typedef typename Grid::ctype CoordType;
    int const dim = Grid::dimension;

    IndexSet const& indexSet = at_c<0>(spaces)->grid().leafIndexSet();


    // Shape function cache. Remember that every thread has to use its own cache!
    typedef ShapeFunctionCache<Grid,RT> SfCache;
    SfCache sfCache,sfCacheN, sfCacheT;

    // Evaluators for shape functions of all FE spaces. Remember that
    // every thread has to use its own collection of evaluators!
    typedef typename result_of::as_vector<
    typename result_of::transform<Spaces, GetEvaluators<SfCache> >::type
    >::type Evaluators;

    typedef typename result_of::as_vector<typename result_of::transform<ExtendedSpaces, GetEvaluators<SfCache> >::type>::type ExtendedEvaluators;

    Evaluators evaluators(transform(spaces,GetEvaluators<SfCache>(&sfCache)));
    ExtendedEvaluators extendedEvaluators(transform(extendedSpaces,GetEvaluators<SfCache>(&sfCacheT)));
    Evaluators evaluatorsNeigh(transform(spaces,GetEvaluators<SfCache>(&sfCacheN)));

    // Define iterator and entity type for stepping through all cells of the grid.
    typedef typename IndexSet::template Codim<0>::template Partition<All_Partition>::Iterator CellIterator;
    typedef typename CellIterator::Entity::LeafIntersectionIterator LII;
    typedef typename CellIterator::Entity Entity;

    typename F::DomainCache domainCache(f.createDomainCache(6+8));         // 6 means RHS+Op
    typename F::DomainCache domainCacheNeigh(f.createDomainCache(6+8));
    typename F::BoundaryCache boundaryCache(f.createBoundaryCache(6+8));

    const int ssize = result_of::size<VariableDescriptions>::type::value;

    typedef Dune::FieldMatrix<RT, sysdim,sysdim> LocalMatrix;
    typedef Dune::FieldVector<Dune::FieldVector<RT, sysdim>,ssize> LocalVector;

    LocalMatrix A;
    LocalVector b,c;
    Dune::FieldVector<RT, sysdim> rk1;
    // Iterate over all cells.
    CellIterator end = indexSet.template end<0,All_Partition>();
    for (CellIterator ci=indexSet.template begin<0,All_Partition>(); ci!=end; ++ci) {

      domainCache.moveTo(*ci);

      moveEvaluatorsToCell(evaluators,*ci);
      moveEvaluatorsToCell(extendedEvaluators,*ci);
      int const shapeFunctionMaxOrder = maxOrder(evaluators);
      shapeFunctionMaxOrder = std::max(shapeFunctionMaxOrder,maxOrder(extendedEvaluators));
      int const p = f.integrationOrder(*ci,shapeFunctionMaxOrder,false);
      int const pb = f.integrationOrder(*ci,shapeFunctionMaxOrder,true);

      LII faceEnd = ci->ileafend();

      A=0;
      rk1=0;
      for(int i=0; i<ssize; ++i) { b[i]=0; c[i]=0; }

      GeometryType gti = ci->type();

      QuadratureRule<CoordType, dim> const& qri=QuadratureRules<CoordType,dim>::rule(gti,p);

      size_t nQuadPosi = qri.size();
      for (size_t g=0; g<nQuadPosi; ++g) {
        FieldVector<CoordType,dim> const& quadPos = qri[g].position();

        moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qri,g,0);
        moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPos,qri,g,0);
        domainCache.evaluateAt(quadPos,evaluators);
        CoordType weightTimesDetJac(ci->geometry().integrationElement(quadPos)); // determinant of jacobian
        weightTimesDetJac *= qri[g].weight();
        rank1KernelVector(rk1,weightTimesDetJac, extendedEvaluators);
        localNormalMatrix(A,weightTimesDetJac, extendedEvaluators, domainCache);

        for_each(varDesc,WeakResiduumMainPart<LocalVector, RT,typename F::DomainCache,ExtendedEvaluators>
        (b, weightTimesDetJac,domainCache, extendedEvaluators));
      }

      for (LII face=ci->ileafbegin(); face!=faceEnd; ++face) {
        if(face.neighbor()) {

          typename LII::EntityPointer o=face.outside();
          domainCacheNeigh.moveTo(*o);
          moveEvaluatorsToCell(evaluatorsNeigh,*o);

        } else if(face.boundary())
        {
          boundaryCache.moveTo(face);
        }
        GeometryType gt=face.intersectionSelfLocal().type();
        QuadratureRule<CoordType, dim-1> const& qr=QuadratureRules<CoordType,dim-1>::rule(gt,pb);
        size_t nQuadPos=qr.size();
        for(size_t g=0; g<nQuadPos; ++g){
          FieldVector<CoordType, dim-1> const& quadPos = qr[g].position();
          FieldVector<CoordType, dim> const& quadPosInSelf = face.intersectionSelfLocal().global(quadPos);
          moveEvaluatorsToIntegrationPoint(evaluators,quadPosInSelf,qr,g,face.numberInSelf());
          moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPosInSelf,qr,g,face.numberInSelf());
          domainCache.evaluateAt(quadPosInSelf,evaluators);
          FieldVector<RT, dim> outerNormal = face.integrationOuterNormal(quadPos);
          outerNormal *= qr[g].weight();
          if(face.neighbor())
          {
            Dune::FieldVector<typename Grid::ctype, Grid::dimension> const& quadPosInNeigh = face.intersectionNeighborLocal().global(quadPos);

            moveEvaluatorsToIntegrationPoint(evaluatorsNeigh,quadPosInNeigh,qr,g,face.numberInNeighbor());
            domainCacheNeigh.evaluateAt(quadPosInNeigh,evaluatorsNeigh);
            for_each(varDesc,GradientAverage<LocalVector,RT,typename F::DomainCache, dim,ExtendedEvaluators>
            (b,outerNormal,domainCache,domainCacheNeigh, extendedEvaluators));

          } else if(face.boundary())
          {
            boundaryCache.evaluateAt(quadPos,evaluators);

            CoordType weightTimesDetJac = face.intersectionGlobal().integrationElement(quadPos); // determinant of jacobian
            weightTimesDetJac *= qr[g].weight();
            for_each(varDesc,WeakBoundaryValues<LocalVector,RT, typename F::DomainCache, typename F::BoundaryCache, dim, ExtendedEvaluators>
            (b, weightTimesDetJac, outerNormal, domainCache, boundaryCache, extendedEvaluators));
          }
        }
      }
      for(int i=0; i<sysdim; ++i) for(int j=0; j<sysdim; ++j) A[i][j] += rk1[i]*rk1[j];
      for_each(varDesc, SolveLocalSystem<LocalVector, LocalMatrix>(c,A,b));
      for_each(varDesc, ScatterLocalSolution<LocalVector, ExtendedEvaluators,typename Result::Functions>(c,extendedEvaluators,result.vars));
    }
  };

  template <class LocalMatrix, class Evaluators, class RT, class Cache, class Row>
  struct LocalNormalMatrixBlock
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;
    LocalNormalMatrixBlock(LocalMatrix& localmatrix_, RT intEl_, Evaluators const& evaluators_,Cache const& cache_):
      localmatrix(localmatrix_), intEl(intEl_),  evaluators(evaluators_), cache(cache_)
    {
    }

    template <class Col>
    void operator()(Col const& ) const
    {
      size_t n = localmatrix.N()/2;
      for(int row=0; row<n; ++row)
        for (size_t col=0; col<n; ++col)
        {

          (localmatrix)[(Row::id)*n+row][(Col::id)*n+col]
                                         += intEl*cache.template d2<Row::id,Col::id>(boost::fusion::at_c<Row::spaceIndex>(evaluators).evalData[indQFull[row]],
                                             boost::fusion::at_c<Col::spaceIndex>(evaluators).evalData[indQFull[col]]);
        }
    }

    LocalMatrix&   localmatrix;
    RT intEl;
    Evaluators const& evaluators;
    Cache const& cache;
  };

  template <class LocalMatrix, class Evaluators, class RT, class Cache, class VarDesc>
  struct LocalNormalMatrixRow
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;
    LocalNormalMatrixRow(LocalMatrix& localmatrix_, RT intEl_, Evaluators const& evaluators_,Cache const& cache_):
      localmatrix(localmatrix_), intEl(intEl_),  evaluators(evaluators_), cache(cache_)
    {
    }

    template <class Row>
    void operator()(Row const&) const
    {
      using namespace boost::fusion;
      for_each(varDesc,LocalNormalMatrixBlock<LocalMatrix, Evaluators, RT,Cache, Row>
      (localmatrix,intEl, evaluators, cache));
    }

    VarDesc varDesc;
    LocalMatrix&   localmatrix;
    RT intEl;
    Evaluators const& evaluators;
    Cache const& cache;
  };

  template <class LocalVector, class LocalMatrix>
  void solveLocalFullSystem(LocalVector& c,
      LocalMatrix & A,
      LocalVector& b)
  {
    double sum;
    //Since Neumann problems are solved, achieve sum_i b_i = 0
    int n = b.N()/2;
    for(int row=0; row < 2; ++row)
    {
      sum=0.0;
      for(int i=0; i<n; ++i) sum+=b[i+n*row];
      sum /= n;
      for(int i=0; i<n; ++i) b[i+n*row]-=sum;
    }
    A.solve(c,b);
  }

  template <class LocalRHS, class RT, class Cache, class Evaluators>
  struct WeakResiduumFull
  {
    WeakResiduumFull(LocalRHS& localRHS_, RT intEl_,Cache const& cache_, Evaluators const& evaluators_):
      localRHS(localRHS_), intEl(intEl_),  evaluators(evaluators_), cache(cache_)
    {
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      size_t n = localRHS.N()/2;
      for(int row=0; row<n; ++row)
      {
        localRHS[row+(Variable::id)*n]
                 -= intEl * cache.template d1<Variable::id>(boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQFull[row]]);
      }
    }
    LocalRHS&   localRHS;
    RT intEl;
    Evaluators const& evaluators;
    Cache const& cache;
  };

  template <class LocalSolution, class Evaluators, class GlobalSolution>
  struct ScatterFullLocalSolution
  {
    ScatterFullLocalSolution(LocalSolution const& localSolution_, Evaluators const& evaluators_, GlobalSolution& globalSolution_):
      localSolution(localSolution_), eval(evaluators_), globalSolution(globalSolution_)
    {}

    template <class Variable>
    void operator()(Variable const& ) const
    {
      using namespace boost::fusion;
      typename result_of::value_at_c<Evaluators,Variable::spaceIndex>::type const& reval = at_c<Variable::spaceIndex>(eval);
      int n=localSolution.N()/2;
      for (size_t i=0; i<n; i++)
        (*at_c<Variable::id>(globalSolution))[reval.globalIndices()[indQFull[i]]] += (localSolution)[i+(Variable::id)*n];
    }

    LocalSolution const&   localSolution;
    Evaluators const& eval;
    GlobalSolution&        globalSolution;
  };


  template <class LocalInd, class RT, class Cache, int dim,class Evaluators>
  struct GradientAverageFull
  {

    GradientAverageFull(LocalInd& localInd_,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        Cache const& cacheNeigh_,
        Evaluators const& evaluators_):
          localInd(localInd_),  cache(cache_), cacheNeigh(cacheNeigh_), evaluators(evaluators_)
    {
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      int n=localInd.N()/2;

      for(int row=0; row<n; ++row)
        localInd[row+(Variable::id)*n] +=0.5*(cache.template d1<Variable::id>(modarg)+cacheNeigh.template d1<Variable::id>(modarg))
        *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQFull[row]].value[0];
    }

    VariationalArg<RT,dim> modarg;
    LocalInd&         localInd;
    Cache const&      cache;
    Cache const&      cacheNeigh;
    Evaluators const & evaluators;
  };

  template <class LocalInd, class RT, class Cache, class BoundaryCache, int dim, class Evaluators>
  struct WeakBoundaryValuesFull
  {
    WeakBoundaryValuesFull(LocalInd& localInd_,
        RT& integrationElement,
        Dune::FieldVector<RT,dim>const & outernormal,
        Cache const& cache_,
        BoundaryCache const& cacheBoundary_,
        Evaluators const& evaluators_):
          localInd(localInd_), cache(cache_), boundaryCache(cacheBoundary_),evaluators(evaluators_)
    {
      modargB.value=integrationElement;
      modargB.gradient[0]=0;
      modargTest.value=1;
      modargTest.gradient[0]=0;
      modarg.value=0;
      modarg.gradient[0]=outernormal;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      int n=localInd.size/2;
      if(boundaryCache.template d2<Variable::pdeid,Variable::id>(modargTest,modargTest) < 1e5)
        for(int row=0; row< n; ++row)
        {
          localInd[row+(Variable::id)*n] +=
              boundaryCache.template d1<Variable::id>(modargB)
              *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQFull[row]].value[0];
        }
      else
        for(int row=0; row< n; ++row)
        {
          localInd[row+(Variable::id)*n] +=
              cache.template d1<Variable::id>(modarg)
              *boost::fusion::at_c<Variable::spaceIndex>(evaluators).evalData[indQFull[row]].value[0];
        }

    }

    VariationalArg<RT,dim> modarg,modargB,modargTest;
    LocalInd&         localInd;
    Cache const& cache;
    BoundaryCache const&   boundaryCache;
    Evaluators const& evaluators;
  };


  /// Computes an error representation function, which is useful for goal oriented adaptivity
  /** Arguments
   *  - result : output: variable set of error functions
   *  - varDesc : description of the testfunctions
   *  - f  : functional, the residual of which is evaluated
   *  - spaces : spaces of the ansatz functions x, neede for the evaluation of f
   *  - extended spaces: spaces, where the error functions live
   *
   *  \attention: this routine is not yet fully general
   */
  template <class VariableDescriptions, int sysdim, class Result, class F, class FS, class Spaces, class ExtendedSpaces>
  void edgeAveragingFull(Result& result,
      F const& f,
      FS const& fs,
      Spaces const& spaces,
      ExtendedSpaces const& extendedSpaces)
  {
    using namespace boost::fusion;
    using namespace Dune;

    VariableDescriptions const varDesc;

    typedef typename SpaceType<Spaces,0>::type::Grid Grid;
    typedef typename SpaceType<Spaces,0>::type::IndexSet IndexSet;
    typedef typename SpaceType<Spaces,0>::type::RT RT;
    typedef typename Grid::ctype CoordType;
    int const dim = Grid::dimension;

    IndexSet const& indexSet = at_c<0>(spaces)->grid().leafIndexSet();


    // Shape function cache. Remember that every thread has to use its own cache!
    typedef ShapeFunctionCache<Grid,RT> SfCache;
    SfCache sfCache,sfCacheN, sfCacheT;

    // Evaluators for shape functions of all FE spaces. Remember that
    // every thread has to use its own collection of evaluators!
    typedef typename result_of::as_vector<
    typename result_of::transform<Spaces, GetEvaluators<SfCache> >::type
    >::type Evaluators;

    typedef typename result_of::as_vector<typename result_of::transform<ExtendedSpaces, GetEvaluators<SfCache> >::type>::type ExtendedEvaluators;

    Evaluators evaluators(transform(spaces,GetEvaluators<SfCache>(&sfCache)));
    ExtendedEvaluators extendedEvaluators(transform(extendedSpaces,GetEvaluators<SfCache>(&sfCacheT)));
    Evaluators evaluatorsNeigh(transform(spaces,GetEvaluators<SfCache>(&sfCacheN)));

    // Define iterator and entity type for stepping through all cells of the grid.
    typedef typename IndexSet::template Codim<0>::template Partition<All_Partition>::Iterator CellIterator;
    typedef typename CellIterator::Entity::LeafIntersectionIterator LII;
    typedef typename CellIterator::Entity Entity;

    typename F::DomainCache domainCache(f.createDomainCache(6+8));         // 6 means RHS+Op
    typename FS::DomainCache domainCacheNeigh(fs.createDomainCache(6+8));
    typename FS::DomainCache domainCacheSimpl(fs.createDomainCache(6+8));
    typename FS::BoundaryCache boundaryCache(fs.createBoundaryCache(6+8));

    typedef Dune::FieldMatrix<RT, 2*sysdim,2*sysdim> LocalMatrix;
    typedef Dune::FieldVector<RT, 2*sysdim> LocalVector;

    LocalMatrix A;
    LocalVector b,c;


    // Iterate over all cells.
    CellIterator end = indexSet.template end<0,All_Partition>();
    for (CellIterator ci=indexSet.template begin<0,All_Partition>(); ci!=end; ++ci) {

      A=0; b=0; c=0;
      domainCache.moveTo(*ci);
      domainCacheSimpl.moveTo(*ci);

      moveEvaluatorsToCell(evaluators,*ci);
      moveEvaluatorsToCell(extendedEvaluators,*ci);
      int const shapeFunctionMaxOrder = maxOrder(evaluators);
      shapeFunctionMaxOrder = std::max(shapeFunctionMaxOrder,maxOrder(extendedEvaluators));
      int const p = f.integrationOrder(*ci,shapeFunctionMaxOrder,false);
      int const pb = f.integrationOrder(*ci,shapeFunctionMaxOrder,true);

      LII faceEnd = ci->ileafend();

      GeometryType gti = ci->type();

      QuadratureRule<CoordType, dim> const& qri=QuadratureRules<CoordType,dim>::rule(gti,p);

      size_t nQuadPosi = qri.size();
      for (size_t g=0; g<nQuadPosi; ++g) {
        FieldVector<CoordType,dim> const& quadPos = qri[g].position();

        moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qri,g,0);
        moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPos,qri,g,0);
        domainCache.evaluateAt(quadPos,evaluators);
        CoordType weightTimesDetJac(ci->geometry().integrationElement(quadPos)); // determinant of jacobian
        weightTimesDetJac *= qri[g].weight();
        for_each(varDesc,LocalNormalMatrixRow<LocalMatrix, ExtendedEvaluators, RT,typename F::DomainCache, VariableDescriptions>
        (A,weightTimesDetJac, extendedEvaluators, domainCache));
        for_each(varDesc,WeakResiduumFull<LocalVector, RT,typename F::DomainCache,ExtendedEvaluators>
        (b, weightTimesDetJac,domainCache, extendedEvaluators));
      }

      for (LII face=ci->ileafbegin(); face!=faceEnd; ++face) {
        if(face.neighbor()) {

          typename LII::EntityPointer o=face.outside();
          domainCacheNeigh.moveTo(*o);
          moveEvaluatorsToCell(evaluatorsNeigh,*o);

        } else if(face.boundary())
        {
          boundaryCache.moveTo(face);
        }
        GeometryType gt=face.intersectionSelfLocal().type();
        QuadratureRule<CoordType, dim-1> const& qr=QuadratureRules<CoordType,dim-1>::rule(gt,pb);
        size_t nQuadPos=qr.size();
        for(size_t g=0; g<nQuadPos; ++g){
          FieldVector<CoordType, dim-1> const& quadPos = qr[g].position();
          FieldVector<CoordType, dim> const& quadPosInSelf = face.intersectionSelfLocal().global(quadPos);
          moveEvaluatorsToIntegrationPoint(evaluators,quadPosInSelf,qr,g,face.numberInSelf());
          moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPosInSelf,qr,g,face.numberInSelf());
          domainCacheSimpl.evaluateAt(quadPosInSelf,evaluators);
          FieldVector<RT, dim> outerNormal = face.integrationOuterNormal(quadPos);
          outerNormal *= qr[g].weight();
          if(face.neighbor())
          {
            Dune::FieldVector<typename Grid::ctype, Grid::dimension> const& quadPosInNeigh = face.intersectionNeighborLocal().global(quadPos);

            moveEvaluatorsToIntegrationPoint(evaluatorsNeigh,quadPosInNeigh,qr,g,face.numberInNeighbor());
            domainCacheNeigh.evaluateAt(quadPosInNeigh,evaluatorsNeigh);
            for_each(varDesc,GradientAverageFull<LocalVector,RT,typename FS::DomainCache, dim,ExtendedEvaluators>
            (b,outerNormal,domainCacheSimpl,domainCacheNeigh, extendedEvaluators));

          } else if(face.boundary())
          {
            boundaryCache.evaluateAt(quadPos,evaluators);

            CoordType weightTimesDetJac = face.intersectionGlobal().integrationElement(quadPos); // determinant of jacobian
            weightTimesDetJac *= qr[g].weight();
            for_each(varDesc,WeakBoundaryValuesFull<LocalVector,RT, typename FS::DomainCache, typename FS::BoundaryCache, dim, ExtendedEvaluators>
            (b, weightTimesDetJac, outerNormal, domainCacheSimpl, boundaryCache, extendedEvaluators));
          }
        }
      }
      solveLocalFullSystem(c,A,b);
      for_each(varDesc, ScatterFullLocalSolution<LocalVector, ExtendedEvaluators,typename Result::Functions>(c,extendedEvaluators,result.vars));
    }
  };

  /// Obtain an error function by solving local Neumann problems in the flavour of Bank/Weiser
  template< class Description, class Grid, class Functional, class FunctionalSimpl>
  class BWErrorFunction
  {

  public:

    typedef FEFunctionSpace<DiscontinuousLagrangeMapper<double,Grid> > RecoverySpace;
    typedef boost::fusion::vector<RecoverySpace const*> RecoverySpaces;
    typedef VariableSetDescription<RecoverySpaces,Description> RecoveryVariableSet;
    typedef typename RecoveryVariableSet::VariableSet RecoveryRepresentation;
    static int const syst =  SystDim<RecoverySpace::dim>::value;
    static int const systFull =  SystDimFull<RecoverySpace::dim>::value;

    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid> > ContRecoverySpace;
    typedef boost::fusion::vector<ContRecoverySpace const*> ContRecoverySpaces;
    typedef VariableSetDescription<ContRecoverySpaces,Description> ContRecoveryVariableSet;
    typedef typename ContRecoveryVariableSet::VariableSet ContRecoveryRepresentation;

    template<class Spaces>
    void getErrorFunction(RecoveryRepresentation & recovery, Spaces const& spaces, FunctionalSimpl const& vf)
    {
      edgeAveraging<Description,syst>(recovery, vf, spaces, recovery.descriptions.spaces);
    }

    template<class Spaces>
    void getFullErrorFunction(RecoveryRepresentation & recovery, Spaces const& spaces, Functional const& vf, FunctionalSimpl const& vfs)
    {
      edgeAveragingFull<Description,systFull>(recovery, vf, vfs, spaces, recovery.descriptions.spaces);
    }
  };

  /// Construct an error indicator in the flavour of Bank/Weiser
  template< class Description, class Function, class Functional>
  typename CellData<typename Function::Grid>::CellDataVector
  BWErrorIndicator(Function const& f, Functional const& vf)
  {
    typedef typename Function::Grid Grid;
    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid> > Space;

    using namespace FunctionViews;
    using namespace boost::fusion;
    typedef FEFunctionSpace<DiscontinuousLagrangeMapper<double, Grid> > RecoverySpace;
    RecoverySpace rspace(at_c<0>(f.vars).space().gridManager(),at_c<0>(f.vars).space().indexSet(),2);
    RecoverySpace r1space(at_c<0>(f.vars).space().gridManager(),at_c<0>(f.vars).space().indexSet(),1);
    typedef boost::fusion::vector<RecoverySpace const*> RecoverySpaces;
    RecoverySpaces rspaces(&rspace);

    typedef VariableSetDescription<RecoverySpaces,typename Functional::TestVars::Variables> VariableSet;

    VariableSet variableSet(rspaces);

    typedef typename VariableSet::VariableSet RecoveryRepresentation;

    RecoveryRepresentation recovery(variableSet);
    int const dim(Space::dim);
    static int const syst =  SystDim<RecoverySpace::dim>::value;

    edgeAveraging<Description,syst>(recovery, vf, f.descriptions.spaces, rspaces);


    typename RecoverySpace::template Element<dim>::type gradientFkt(r1space);
    interpolateGloballyWeak<PlainAverage>(gradientFkt,makeView<Gradient>(at_c<0>(recovery.vars)));

    LocalIntegral<Space> localIntegral;
    typename CellData<Grid>::CellDataVector
    errorIndicator(localIntegral(
        makeView<AbsSquare>(gradientFkt)
        ,at_c<0>(f.vars).space())
    );
    return errorIndicator;
  }
} /* end of namespace Kaskade */
/**
@}
 */

#endif


