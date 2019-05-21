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

#ifndef RESIDUAL_ERROREST_HH
#define RESIDUAL_ERROREST_HH

#include "fem/assemble.hh"

namespace Kaskade
{
  template <class LocalInd,  class Evaluators, class RT, class TestFunctions, class Cache>
  struct GradientJumpsTimesTestFunction
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;

    GradientJumpsTimesTestFunction(LocalInd& localInd_,
        Dune::FieldVector<RT,dim>const & outernormal_,
        Evaluators const& evaluators_,
        TestFunctions const& tf_,
        Cache const& cache_,
        Cache const& cacheNeigh_):
          localInd(localInd_), evaluators(evaluators_), tf(tf_), cache(cache_), cacheNeigh(cacheNeigh_)
    {
      modarg.value=0;
      modarg.gradient[0]=outernormal_;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      localInd[Variable::id] -= 0.5*(cache.template d1<Variable::id>(modarg)-cacheNeigh.template d1<Variable::id>(modarg))*
          boost::fusion::at_c<Variable::id>(tf.vars).value(boost::fusion::at_c<Variable::spaceIndex>(evaluators));
    }

    VariationalArg<RT,dim> modarg;
    LocalInd&         localInd;
    Evaluators const& evaluators;
    TestFunctions const& tf;
    Cache const&      cache;
    Cache const&      cacheNeigh;
  };

  template <class LocalInd,  class Evaluators, class RT, class TestFunctions, class Cache, class BoundaryCache>
  struct BoundaryJumpsTimesTestFunction
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;

    BoundaryJumpsTimesTestFunction(LocalInd& localInd_,
        Dune::FieldVector<RT,dim>const & outernormal_,
        RT const & integrationElement_,
        Evaluators const& evaluators_,
        TestFunctions const& tf_,
        Cache const& cache_,
        BoundaryCache const& cacheBoundary_):
          localInd(localInd_), evaluators(evaluators_), tf(tf_), cache(cache_), cacheBoundary(cacheBoundary_)
    {
      modarg.value=0;
      modarg.gradient[0]=outernormal_;
      modargB.value=integrationElement_;
      modargB.gradient[0]=0;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      if(cacheBoundary.template d2<Variable::pdeid,Variable::id>(modargB,modargB) < 1e5)
        localInd[Variable::id] -=
            (cache.template d1<Variable::id>(modarg)+cacheBoundary.template d1<Variable::id>(modargB))*
            boost::fusion::at_c<Variable::id>(tf.vars).value(boost::fusion::at_c<Variable::spaceIndex>(evaluators));
    }

    VariationalArg<RT,dim> modarg,modargB;
    LocalInd&         localInd;
    Evaluators const& evaluators;
    TestFunctions const& tf;
    Cache const&      cache;
    BoundaryCache const&      cacheBoundary;
  };

  template <class LocalInd,  class Evaluators, class RT, class TestFunctions, class Cache>
  struct StrongResidualsTimesTestFunction
  {
    static const int dim = boost::fusion::result_of::value_at_c<Evaluators,0>::type::Space::dim;

    StrongResidualsTimesTestFunction(LocalInd& localInd_,
        RT const & integrationElement,
        Evaluators const& evaluators_,
        TestFunctions const& tf_,
        Cache const& cache_):
          localInd(localInd_), evaluators(evaluators_), tf(tf_), cache(cache_)
    {
      modarg.value=integrationElement;
      modarg.gradient[0]=0;
    }

    template <class Variable>
    void operator()(Variable const& ) const
    {
      localInd[Variable::id] += cache.template d1<Variable::id>(modarg)*
          boost::fusion::at_c<Variable::id>(tf.vars).value(boost::fusion::at_c<Variable::spaceIndex>(evaluators));
    }

    VariationalArg<RT,dim> modarg;
    LocalInd&         localInd;
    Evaluators const& evaluators;
    TestFunctions const& tf;
    Cache const&      cache;
  };


  /** \ingroup adapt
   * \brief Evaluates strong residual, applied to some test-function <f(x), tf>
   *
   * Arguments
   *  - ind : output vector: for each equation of the system and each cell the result of the computation
   *  - varDesc : description of the testfunctions
   *  - f  : functional, the residual of which is evaluated
   *  - tf : testfunctions the residual is applied to
   *  - spaces : spaces of the ansatz functions x
   *  - testSpaces : spaces of the test functions tf
   */
  template <class ErrorDescriptions, class F, class FS, class TestFunctions, class Spaces, class TestSpaces>
  void residualTimesTestFunction(std::vector< Dune::FieldVector<typename SpaceType<Spaces,0>::type::RT,
      boost::fusion::result_of::size<ErrorDescriptions>::type::value> >& ind,
      ErrorDescriptions const& varDesc,
      F const& f,
      FS const& fs,
      TestFunctions const& tf,
      Spaces const& spaces,
      TestSpaces const& testSpaces)
  {
    using namespace boost::fusion;
    using namespace Dune;


    typedef typename SpaceType<Spaces,0>::type::Grid Grid;
    typedef typename SpaceType<Spaces,0>::type::IndexSet IndexSet;
    typedef typename SpaceType<Spaces,0>::type::RT RT;
    typedef typename Grid::ctype CoordType;
    int const dim = Grid::dimension;

    IndexSet const& indexSet = at_c<0>(spaces)->grid().leafIndexSet();

    ind.resize(indexSet.size(0));

    for(int i=0; i<ind.size(); ++i)
      for(int j=0; j<ind[i].N(); ++j)
        ind[i][j]=0.0;

    // Shape function cache. Remember that every thread has to use its own cache!
    typedef ShapeFunctionCache<Grid,RT> SfCache;
    SfCache sfCache,sfCacheN,sfCacheT;

    // Evaluators for shape functions of all FE spaces. Remember that
    // every thread has to use its own collection of evaluators!
    typedef typename result_of::as_vector<
    typename result_of::transform<Spaces, GetEvaluators<SfCache> >::type
    >::type Evaluators;

    typedef typename result_of::as_vector<
    typename result_of::transform<TestSpaces, GetEvaluators<SfCache> >::type
    >::type TestEvaluators;

    Evaluators evaluators(transform(spaces,GetEvaluators<SfCache>(&sfCache)));
    TestEvaluators testEvaluators(transform(testSpaces,GetEvaluators<SfCache>(&sfCacheT)));
    Evaluators evaluatorsNeigh(transform(spaces,GetEvaluators<SfCache>(&sfCacheN)));

    // Define iterator and entity type for stepping through all cells of the grid.
    typedef typename IndexSet::template Codim<0>::template Partition<All_Partition>::Iterator CellIterator;
    typedef typename CellIterator::Entity::LeafIntersectionIterator LII;
    typedef typename CellIterator::Entity Entity;

    typename F::DomainCache domainCache(f.createDomainCache(2));         // 2 means RHS
    typename FS::DomainCache domainCacheSimpl(fs.createDomainCache(2+8));
    typename FS::DomainCache domainCacheNeigh(fs.createDomainCache(2+8));
    typename FS::BoundaryCache boundaryCache(fs.createBoundaryCache(6));
    typedef Dune::FieldVector<RT,result_of::size<ErrorDescriptions>::type::value> LocalInd;

    // Iterate over all cells.
    CellIterator end = indexSet.template end<0,All_Partition>();
    for (CellIterator ci=indexSet.template begin<0,All_Partition>(); ci!=end; ++ci) {

      domainCache.moveTo(*ci);
      domainCacheSimpl.moveTo(*ci);

      moveEvaluatorsToCell(evaluators,*ci);
      moveEvaluatorsToCell(testEvaluators,*ci);
      int shapeFunctionMaxOrder = maxOrder(evaluators);
      shapeFunctionMaxOrder = std::max(shapeFunctionMaxOrder,maxOrder(testEvaluators));
      int p = f.integrationOrder(*ci,shapeFunctionMaxOrder,false);
      int pb = f.integrationOrder(*ci,shapeFunctionMaxOrder,true);
      pb=shapeFunctionMaxOrder+2;
      p = pb+1;

      GeometryType gti = ci->type();

      QuadratureRule<CoordType, dim> const& qri=QuadratureRules<CoordType,dim>::rule(gti,p);

      size_t nQuadPosi = qri.size();
      for (size_t g=0; g<nQuadPosi; ++g) {
        FieldVector<CoordType,dim> const& quadPos = qri[g].position();

        moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qri,g,0);
        moveEvaluatorsToIntegrationPoint(testEvaluators,quadPos,qri,g,0);
        domainCache.evaluateAt(quadPos,evaluators);
        CoordType weightTimesDetJac(ci->geometry().integrationElement(quadPos)); // determinant of jacobian
        weightTimesDetJac *= qri[g].weight();
        for_each(varDesc, StrongResidualsTimesTestFunction<LocalInd,TestEvaluators,RT, TestFunctions, typename F::DomainCache>
        (ind[indexSet.index(*ci)],weightTimesDetJac,testEvaluators, tf,domainCache));
      }

      LII faceEnd = ci->ileafend();
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
          moveEvaluatorsToIntegrationPoint(testEvaluators,quadPosInSelf,qr,g,face.numberInSelf());
          domainCacheSimpl.evaluateAt(quadPosInSelf,evaluators);
          FieldVector<RT, dim> outerNormal = face.integrationOuterNormal(quadPos);
          outerNormal *= qr[g].weight();
          if(face.neighbor())
          {
            Dune::FieldVector<typename Grid::ctype, Grid::dimension> const& quadPosInNeigh = face.intersectionNeighborLocal().global(quadPos);

            moveEvaluatorsToIntegrationPoint(evaluators,Neigh,quadPosInNeigh,qr,g,face.numberInNeighbor());
            domainCacheNeigh.evaluateAt(quadPosInNeigh,evaluatorsNeigh);

            for_each(varDesc,GradientJumpsTimesTestFunction<LocalInd,TestEvaluators,RT, TestFunctions,
                typename FS::DomainCache>(ind[indexSet.index(*ci)],outerNormal,testEvaluators, tf,domainCacheSimpl,domainCacheNeigh));

          } else if(face.boundary())
          {
            boundaryCache.evaluateAt(quadPos,evaluators);

            CoordType weightTimesDetJac = face.intersectionGlobal().integrationElement(quadPos); // determinant of jacobian
            weightTimesDetJac *= qr[g].weight();

            for_each(varDesc,BoundaryJumpsTimesTestFunction<LocalInd,TestEvaluators,RT, TestFunctions,
                typename FS::DomainCache, typename FS::BoundaryCache>(ind[indexSet.index(*ci)],outerNormal, weightTimesDetJac,
                    testEvaluators, tf,domainCacheSimpl,boundaryCache));
          }
        }
      }
    }
  };
} // end of namespace Kaskade
#endif
