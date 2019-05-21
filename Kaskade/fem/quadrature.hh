/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef QUADRATURE_HH
#define QUADRATURE_HH

#include <map>

#include "dune/geometry/quadraturerules.hh"

#include "utilities/threading.hh"

namespace Kaskade
{
  /**
   * \todo docme
   */
  template<class QuadRule>
  struct QuadratureTraits
  {
    typedef QuadRule Rule;

    static Rule const& rule(const Dune::GeometryType& t, int p)
    {
      return instance().rl;
    }

  private:
    static QuadratureTraits& instance()
    {
      static QuadratureTraits instance;
      return instance;
    }

    QuadRule const rl;
  };

  /**
   * \brief A cache that stores suitable quadrature rules for quick retrieval.
   * 
   * Essentially this replicates the functionality of Dune::QuadratureRules<GType, dim>::rule, 
   * but with higher performance (compared to Dune 2.2).
   */
  template<class GType, int dim>
  struct QuadratureTraits<Dune::QuadratureRule<GType, dim> >
  {
    typedef Dune::QuadratureRule<GType, dim> Rule;

    QuadratureTraits(): cacheP(-1), cacheRule(nullptr)
    {
    }

    /**
     * \brief Returns a quadrature rule for the given geometry type and integration order.
     *
     * This method is *not* thread-safe. Take care that each thread has its own QuadratureTraits object.
     */
    Rule const& rule(const Dune::GeometryType& t, int p)
    {
      // The higher performance is based
      // (i)   caching the previous rule
      // (ii)  using smaller std::maps by one map for each order (results in fewer comparisons if multiple orders are used)
      // (iii) avoiding the dumb double look-up done in Dune
      // (iv)  using one cache per thread in order to avoid locking (which is necessary for the Dune::QuadratureRules singleton)

      assert(p>=0);

      // Check if the previously returned rule satisfies the needs (this is a frequent case)
      if (p==cacheP && t==cacheT)
      {
        return *cacheRule;
      }

      // Not immediately available - look for the wanted rule
      if (p<=pmax)
      {
        // Low order rule - look in the specialized maps
        typename DirectMap::const_iterator i = rules[p].find(t);
        if (i == rules[p].end())
        {
          // Rule is not yet in our own map - fetch from Dune and remember here.
          // Remember to acquire a lock, as the Dune factory is not thread-safe.
#ifndef KASKADE_SEQUENTIAL
          boost::lock_guard<boost::mutex> lock(Kaskade::DuneQuadratureRulesMutex);
#endif
          cacheRule = &Dune::QuadratureRules<GType, dim>::rule(t,p);
          rules[p][t] = cacheRule;
        }
        else
          cacheRule = i->second;
      }
      else
      {
        // High order rule - look in the fallback map
        typename OverspillMap::const_iterator i = ruleOverspill.find(std::make_pair(t,p));
        if (i==ruleOverspill.end())
        {
          // Rule is not yet in our own map - fetch from Dune and remember here
          // Remember to acquire a lock, as the Dune factory is not thread-safe.
#ifndef KASKADE_SEQUENTIAL
          boost::lock_guard<boost::mutex> lock(Kaskade::DuneQuadratureRulesMutex);
#endif
          cacheRule = &Dune::QuadratureRules<GType, dim>::rule(t,p);
          ruleOverspill[std::make_pair(t,p)] = cacheRule;
        }
        else
          cacheRule = i->second;
      }

      // cache the result
      cacheT = t;
      cacheP = p;

      return *cacheRule;
    }

  private:
    static int const pmax = 8;

    // One would assume that due to frequent lookup and rare insertion, a boost::container::flat_map would
    // provide better performance. However, this appears not to be the case.
    //typedef boost::container::flat_map<Dune::GeometryType,Rule const*> DirectMap;
    //typedef boost::container::flat_map<std::pair<Dune::GeometryType,int>,Rule const*> OverspillMap;
    typedef std::map<Dune::GeometryType,Rule const*> DirectMap;
    typedef std::map<std::pair<Dune::GeometryType,int>,Rule const*> OverspillMap;
    DirectMap rules[pmax+1];
    OverspillMap ruleOverspill;

    Dune::GeometryType cacheT;
    int                cacheP;
    Rule const*        cacheRule;
  };


}

#endif