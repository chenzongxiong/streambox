/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ERROR_DISTRIBUTION_HH
#define ERROR_DISTRIBUTION_HH

#include <string>

#include <boost/utility.hpp>
#include <boost/fusion/include/as_vector.hpp>

#include "fem/assemble.hh"
#include "fem/functionspace.hh"
#include "fem/lagrangespace.hh"
#include "fem/shapefunctioncache.hh"

// forward declarations
namespace Dune
{
  template <class,int> class FieldVector;
  template <class,int,int> class FieldMatrix;
  template <class, int> class QuadratureRule;
}

namespace Kaskade
{
  using namespace boost::fusion;
  using namespace AssemblyDetail;

  enum class ErrorNorm { Energy, L2, H1, H1_half };

  template <class OriginalEvaluators, class ExtensionEvaluators, class Functional,
            class ArgYH, class ArgYE, class ArgUH, class ArgUE, class ArgPH, class ArgPE>
  void evaluateData(OriginalEvaluators& originalEvaluators, ExtensionEvaluators& extendedEvaluators, Functional const& f,
                    ArgYH& yl, ArgYE& yh, ArgUH& ul, ArgUE& uh, ArgPH& pl, ArgPE& ph, double w)
  {
      if( f.considerControlVariable_ )
      {
        uh.value += w * at_c<Functional::uIdx>(f.errorEstimateH.data).value(at_c<Functional::uSHIdx>(extendedEvaluators));
        ul.value += w * at_c<Functional::uIdx>(f.errorEstimateL.data).value(at_c<Functional::uSLIdx>(originalEvaluators));
        uh.gradient += w * at_c<Functional::uIdx>(f.errorEstimateH.data).gradient(at_c<Functional::uSHIdx>(extendedEvaluators));
        ul.gradient += w * at_c<Functional::uIdx>(f.errorEstimateL.data).gradient(at_c<Functional::uSLIdx>(originalEvaluators));
      }
      if( f.considerStateVariable_ )
      {
        yh.value += w * at_c<Functional::yIdx>(f.errorEstimateH.data).value(at_c<Functional::ySHIdx>(extendedEvaluators));
        yl.value += w * at_c<Functional::yIdx>(f.errorEstimateL.data).value(at_c<Functional::ySLIdx>(originalEvaluators));
        yh.gradient += w * at_c<Functional::yIdx>(f.errorEstimateH.data).gradient(at_c<Functional::ySHIdx>(extendedEvaluators));
        yl.gradient += w * at_c<Functional::yIdx>(f.errorEstimateL.data).gradient(at_c<Functional::ySLIdx>(originalEvaluators));
      }
      if( f.considerAdjointVariable_ )
      {
        ph.value += w * at_c<Functional::pIdx>(f.errorEstimateH.data).value(at_c<Functional::pSHIdx>(extendedEvaluators));
        pl.value += w * at_c<Functional::pIdx>(f.errorEstimateL.data).value(at_c<Functional::pSLIdx>(originalEvaluators));
        ph.gradient += w * at_c<Functional::pIdx>(f.errorEstimateH.data).gradient(at_c<Functional::pSHIdx>(extendedEvaluators));
        pl.gradient += w * at_c<Functional::pIdx>(f.errorEstimateL.data).gradient(at_c<Functional::pSLIdx>(extendedEvaluators));
      }
  }

  template <class ArgYL, class ArgYH, class ArgUL, class ArgUH, class ArgPL, class ArgPH>
  void clearVarArgs(ArgYL& yl, ArgYH& yh, ArgUL& ul, ArgUH& uh, ArgPL& pl, ArgPH& ph)
  {
    ul.value = 0;  ul.gradient = 0;
    yl.value = 0;  yl.gradient = 0;
    pl.value = 0;  pl.gradient = 0;
    uh.value = 0;  uh.gradient = 0;
    yh.value = 0;  yh.gradient = 0;
    ph.value = 0;  ph.gradient = 0;
  }

  template <class ArgU, class ArgY, class ArgP>
  double computeL2Error(ArgU const& u, ArgY const& y, ArgP const& p)
  {
    return u.value*u.value + y.value*y.value + p.value*p.value;
  }

  template <class ArgU, class ArgY, class ArgP>
  double computeH1HalfError(ArgU const& u, ArgY const& y, ArgP const& p)
  {
    LinAlg::EuclideanScalarProduct sp;
    return sp(u.gradient,u.gradient) + sp(y.gradient,y.gradient) + sp(p.gradient,p.gradient);
  }


  template <class Functional, class ExtendedAnsatzVars >
  class ErrorDistribution
  {
    typedef typename Functional::AnsatzVars OriginalAnsatzVars;
    typedef typename OriginalAnsatzVars::Grid Grid;
    typedef ShapeFunctionCache<Grid, typename Functional::Scalar> SfCache;
    typedef ShapeFunctionCache<Grid, typename Functional::Scalar> SfCache2;
    typedef typename OriginalAnsatzVars::Spaces OriginalSpaces;
    typedef typename ExtendedAnsatzVars::Spaces ExtendedSpaces;
    typedef typename result_of::as_vector<typename result_of::transform<OriginalSpaces, GetEvaluators<SfCache> >::type>::type OriginalEvaluators;
    typedef typename result_of::as_vector<typename result_of::transform<ExtendedSpaces, GetEvaluators<SfCache2> >::type>::type ExtendedEvaluators;
    typedef typename Grid::ctype CoordType;
    typedef Dune::QuadratureRule<typename Functional::AnsatzVars::Grid::ctype, Functional::AnsatzVars::Grid::dimension> QuadRule;
    typedef Dune::QuadratureRule<typename Functional::AnsatzVars::Grid::ctype, Functional::AnsatzVars::Grid::dimension-1> FaceQuadRule;
  public:
    typedef typename Functional::Scalar Scalar;
    typedef FEFunctionSpace<DiscontinuousLagrangeMapper<Scalar,typename Grid::LeafGridView> > AnsatzSpace;
    typedef boost::fusion::vector<AnsatzSpace const*> AnsatzSpaces;
    typedef Variable<SpaceIndex<0>,Components<1>,VariableId<0> > AnsatzVariableInformation;
    typedef boost::fusion::vector<AnsatzVariableInformation> VariableDescriptions;
    typedef VariableSetDescription<AnsatzSpaces,VariableDescriptions> AnsatzVars;
    typedef typename AnsatzVars::template CoefficientVectorRepresentation<0,1>::type ErrorVector;
    typedef AnsatzVars TestVars;
    typedef OriginalAnsatzVars OriginVars;
    //typedef AnsatzVars OriginVars;
    static int const dim = Grid::dimension;
    static ProblemType const type = Functional::type;

    static constexpr int yIdx = getStateId<Functional>();
    static constexpr int uIdx = getControlId<Functional>();
    static constexpr int pIdx = getAdjointId<Functional>();
    static constexpr int uSLIdx = result_of::value_at_c<typename OriginalAnsatzVars::Variables, uIdx>::type::spaceIndex;
    static constexpr int ySLIdx = result_of::value_at_c<typename OriginalAnsatzVars::Variables, yIdx>::type::spaceIndex;
    static constexpr int pSLIdx = result_of::value_at_c<typename OriginalAnsatzVars::Variables, pIdx>::type::spaceIndex;
    static constexpr int uSHIdx = result_of::value_at_c<typename ExtendedAnsatzVars::Variables, uIdx>::type::spaceIndex;
    static constexpr int ySHIdx = result_of::value_at_c<typename ExtendedAnsatzVars::Variables, yIdx>::type::spaceIndex;
    static constexpr int pSHIdx = result_of::value_at_c<typename ExtendedAnsatzVars::Variables, pIdx>::type::spaceIndex;

    class DomainCache
    {
    public:
      DomainCache(ErrorDistribution const& f_, typename OriginVars::VariableSet const& /*vars*/, int flags=7)
        : f(f_), domainCache(f.functional,f.iterate,flags)
      {
      }

      template <class Entity>
      void moveTo(Entity const &entity)
      {
        e  = &entity;
        domainCache.moveTo(entity);
      }

      template <class Position, class Evaluators>
      void evaluateAt(Position const&, Evaluators const& evaluators)
      {
        clearVarArgs(yl,yh,ul,uh,pl,ph);
        OriginalEvaluators originalEvaluators(transform(f.iterate.descriptions.spaces,GetEvaluators<SfCache>(&sfCache)));
        ExtendedEvaluators extendedEvaluators(transform(f.errorEstimateH.descriptions.spaces,GetEvaluators<SfCache>(&extendedSFCache)));
        QuadRule qr = QuadratureTraits<QuadRule>().rule(e->type(),f.qOrder);
        moveEvaluatorsToCell(originalEvaluators,*e);
        moveEvaluatorsToCell(extendedEvaluators,*e);
        useQuadratureRuleInEvaluators(originalEvaluators,qr,0);
        useQuadratureRuleInEvaluators(extendedEvaluators,qr,0);

        size_t nQuadPos = qr.size();
        for (size_t g=0; g<nQuadPos; ++g)
        {
          // pos of integration point
          Dune::FieldVector<CoordType,dim> quadPos = qr[g].position();
          // for all spaces involved, update the evaluators associated
          // to this quadrature point
          moveEvaluatorsToIntegrationPoint(originalEvaluators,quadPos);
          moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPos);
          // prepare evaluation of functional
          domainCache.evaluateAt(qr[g].position(),originalEvaluators);

          evaluateData(originalEvaluators, extendedEvaluators, f, yl, yh, ul, uh, pl, ph, qr[g].weight());
        }
      }

      Scalar d0() const { return 0; }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
      d1 (VariationalArg<Scalar,dim> const& arg) const
      {
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<yIdx>::m> y(yh);
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<uIdx>::m> u(uh);
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<pIdx>::m> p(ph);

        if( !f.onlyH_ )
        {
          u.value += ul.value;
          u.gradient += ul.gradient;
          y.value += yl.value;
          y.gradient += yl.gradient;
          p.value += pl.value;
          p.gradient += pl.gradient;
        }

        if( f.errorNorm == ErrorNorm::L2 ) return computeL2Error(u,y,p);
        if( f.errorNorm == ErrorNorm::H1_half ) return computeH1HalfError(u,y,p);
        if( f.errorNorm == ErrorNorm::H1 ) return computeL2Error(u,y,p) + computeH1HalfError(u,y,p);

        Scalar result = 0;
        result += Functional::template D2<yIdx,yIdx>::present ? domainCache.template d2_impl<yIdx,yIdx>(y,y) : 0.;
        result += Functional::template D2<yIdx,uIdx>::present ? domainCache.template d2_impl<yIdx,uIdx>(y,u) : 0.;
        result += Functional::template D2<uIdx,yIdx>::present ? domainCache.template d2_impl<uIdx,yIdx>(u,y) : 0.;
        result += Functional::template D2<uIdx,uIdx>::present ? domainCache.template d2_impl<uIdx,uIdx>(u,u) : 0.;
        if( !f.useStateNormForAdjoint_ )
        {
          result += Functional::template D2<yIdx,pIdx>::present ? domainCache.template d2_impl<yIdx,pIdx>(y,p) : 0.;
          result += Functional::template D2<pIdx,yIdx>::present ? domainCache.template d2_impl<pIdx,yIdx>(p,y) : 0.;
          result += Functional::template D2<pIdx,uIdx>::present ? domainCache.template d2_impl<pIdx,uIdx>(p,u) : 0.;
          result += Functional::template D2<uIdx,pIdx>::present ? domainCache.template d2_impl<uIdx,pIdx>(u,p) : 0.;
          result += Functional::template D2<pIdx,pIdx>::present ? domainCache.template d2_impl<pIdx,pIdx>(p,p) : 0.;
        }
        else
        {
          result += Functional::template D2<pIdx,pIdx>::present ? domainCache.template d2_impl<yIdx,yIdx>(p,p) : 0.;
        }
        return result*arg.value[0];
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar,1,1>	d2 (VariationalArg<Scalar,dim> const&, VariationalArg<Scalar,dim> const&) const { return Dune::FieldMatrix<Scalar,1,1>(0); }

    private:
      ErrorDistribution const& f;
      typename Functional::DomainCache domainCache;
      typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<yIdx>::m> yl;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<yIdx>::m> yh;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<uIdx>::m> ul;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<uIdx>::m> uh;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<pIdx>::m> pl;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<pIdx>::m> ph;
      SfCache sfCache;
      SfCache2 extendedSFCache;
    };

    class BoundaryCache
    {
    public:
      BoundaryCache(ErrorDistribution const& f_, typename OriginVars::VariableSet const& vars, int flags=7)
        : f(f_), boundaryCache(f.functional,f.iterate,flags)
      {}

      template <class FaceIterator>
      void moveTo(FaceIterator const& entity)
      {
        e = &entity;
        boundaryCache.moveTo(entity);
      }

      template <class Evaluators>
      void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype, AnsatzVars::Grid::dimension-1> const& xlocal, Evaluators const& evaluators)
      {
        clearVarArgs(yl,yh,ul,uh,pl,ph);
        OriginalEvaluators originalEvaluators(transform(f.iterate.descriptions.spaces,GetEvaluators<SfCache>(&sfCache)));
        ExtendedEvaluators extendedEvaluators(transform(f.errorEstimateH.descriptions.spaces,GetEvaluators<SfCache>(&extendedSFCache)));
        FaceQuadRule qr = QuadratureTraits<FaceQuadRule>().rule((*e)->geometryInInside().type(),f.qOrder);
        moveEvaluatorsToCell(originalEvaluators,*at_c<ySLIdx>(evaluators).cell_);
        moveEvaluatorsToCell(extendedEvaluators,*at_c<ySLIdx>(evaluators).cell_);
        moveEvaluatorsToCell(originalEvaluators,*e);
        moveEvaluatorsToCell(extendedEvaluators,*e);
        useQuadratureRuleInEvaluators(originalEvaluators,qr,(*e)->indexInInside());
        useQuadratureRuleInEvaluators(extendedEvaluators,qr,(*e)->indexInInside());

        size_t nQuadPos = qr.size();
        for (size_t g=0; g<nQuadPos; ++g)
        {
          // pos of integration point
          Dune::FieldVector<CoordType,dim> quadPos = (*e)->geometryInInside().global(qr[g].position());
          // for all spaces involved, update the evaluators associated
          // to this quadrature point
          moveEvaluatorsToIntegrationPoint(originalEvaluators,quadPos);
          moveEvaluatorsToIntegrationPoint(extendedEvaluators,quadPos);
          // prepare evaluation of functional
          boundaryCache.evaluateAt(qr[g].position(),originalEvaluators);

          evaluateData(originalEvaluators, extendedEvaluators, f, yl, yh, ul, uh, pl, ph, qr[g].weight());
        }
      }

      Scalar d0() const { return 0; }

      template<int row, int dim>
      Dune::FieldVector<Scalar,1> d1 (VariationalArg<Scalar,dim> const& arg) const
      {
        return Dune::FieldVector<Scalar,1>(0);
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<yIdx>::m> y(yh);
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<uIdx>::m> u(uh);
        VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<pIdx>::m> p(ph);

        if( !f.onlyH_ )
        {
          u.value += ul.value;
          u.gradient += ul.gradient;
          y.value += yl.value;
          y.gradient += yl.gradient;
          p.value += pl.value;
          p.gradient += pl.gradient;
        }

        if( f.errorNorm != ErrorNorm::Energy ) return computeL2Error(u,y,p);

        Scalar result = Functional::template D2<yIdx,yIdx>::present ? boundaryCache.template d2_impl<yIdx,yIdx>(y,y) : 0.;
        result += Functional::template D2<yIdx,uIdx>::present ? boundaryCache.template d2_impl<yIdx,uIdx>(y,u) : 0.;
        result += Functional::template D2<uIdx,yIdx>::present ? boundaryCache.template d2_impl<uIdx,yIdx>(u,y) : 0.;
        result += Functional::template D2<uIdx,uIdx>::present ? boundaryCache.template d2_impl<uIdx,uIdx>(u,u) : 0.;
        if( !f.useStateNormForAdjoint_ )
        {
          result += Functional::template D2<yIdx,pIdx>::present ? boundaryCache.template d2_impl<yIdx,pIdx>(y,p) : 0.;
          result += Functional::template D2<pIdx,yIdx>::present ? boundaryCache.template d2_impl<pIdx,yIdx>(p,y) : 0.;
          result += Functional::template D2<pIdx,uIdx>::present ? boundaryCache.template d2_impl<pIdx,uIdx>(p,u) : 0.;
          result += Functional::template D2<uIdx,pIdx>::present ? boundaryCache.template d2_impl<uIdx,pIdx>(u,p) : 0.;
          result += Functional::template D2<pIdx,pIdx>::present ? boundaryCache.template d2_impl<pIdx,pIdx>(p,p) : 0.;
        }
        else
        {
          result += Functional::template D2<pIdx,pIdx>::present ? boundaryCache.template d2_impl<yIdx,yIdx>(p,p) : 0.;
        }
        return result*arg.value[0];
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar,1,1> d2 (VariationalArg<Scalar,dim> const&, VariationalArg<Scalar,dim> const&) const { return Dune::FieldMatrix<Scalar,1,1>(0); }

    private:
      ErrorDistribution const& f;
      typename Functional::BoundaryCache boundaryCache;
      typename AnsatzVars::Grid::LeafGridView::IntersectionIterator const* e;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<yIdx>::m> yl;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<yIdx>::m> yh;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<uIdx>::m> ul;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<uIdx>::m> uh;
      VariationalArg<Scalar,dim,OriginalAnsatzVars::template Components<pIdx>::m> pl;
      VariationalArg<Scalar,dim,ExtendedAnsatzVars::template Components<pIdx>::m> ph;
      SfCache sfCache;
      SfCache2 extendedSFCache;
    };

    ErrorDistribution(Functional const& functional_, typename OriginalAnsatzVars::VariableSet const& iterate_,
                      typename OriginalAnsatzVars::VariableSet const& errorEstimateL_, typename ExtendedAnsatzVars::VariableSet const& errorEstimateH_)
      : functional(functional_), iterate(iterate_), errorEstimateL(errorEstimateL_), errorEstimateH(errorEstimateH_), ansatzSpace(boost::fusion::at_c<0>(iterate.descriptions.spaces)->gridManager(), iterate.descriptions.gridView,0),
        ansatzSpaces(&ansatzSpace), varName( {"error"} ), ansatzVars(ansatzSpaces,varName)

    {}

    template <int row>
    struct D1
    {
      static bool const present   = true;
      static bool const constant  = false;
    };

    template <int row, int col>
    struct D2
    {
      static bool const present = false;
      static bool const symmetric = false;
      static bool const lumped = false;
    };

    template <class Cell>
    int integrationOrder(Cell const&, int, bool) const { return 0; }

    AnsatzSpaces const& getSpaces() const { return ansatzSpaces; }

    AnsatzVars const& getVariableSetDescription() const { return ansatzVars; }

    void ignoreLowerOrderError(bool ignore) { onlyH_ = ignore; }
    void considerStateVariable(bool consider) { considerStateVariable_ = consider; }
    void considerControlVariable(bool consider) { considerControlVariable_ = consider; }
    void considerAdjointVariable(bool consider) { considerAdjointVariable_ = consider; }
    void setErrorNorm(ErrorNorm errNorm) { errorNorm = errNorm; }
    void useStateNormForAdjoint(bool useStateNorm) { useStateNormForAdjoint_ = useStateNorm; }

    friend class DomainCache;
    friend class BoundaryCache;

    Functional const& functional;
    typename OriginalAnsatzVars::VariableSet const& iterate;
    typename OriginalAnsatzVars::VariableSet const& errorEstimateL;
    typename ExtendedAnsatzVars::VariableSet const& errorEstimateH;
    AnsatzSpace ansatzSpace;
    AnsatzSpaces ansatzSpaces;
    std::string varName[1];
    AnsatzVars ansatzVars;
    ErrorNorm errorNorm = ErrorNorm::Energy;
    bool onlyH_ = false, considerStateVariable_ = true, considerControlVariable_ = true, considerAdjointVariable_ = false, useStateNormForAdjoint_ = false;
    int qOrder = 6;
  };
}
#endif
