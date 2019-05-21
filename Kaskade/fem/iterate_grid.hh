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

#ifndef ITERATE_GRID_HH
#define ITERATE_GRID_HH

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

#include "dune/grid/common/grid.hh"
#include "dune/geometry/quadraturerules.hh"

#include "fem/functionspace.hh"
#include "fem/quadrature.hh"
#include "fem/variables.hh"

#include "utilities/threading.hh"

namespace Kaskade
{
  namespace GridIterateDetail
  {
    template <class VariableDescriptions, class Functor, class Functions, class Spaces, class CellRange>
    void gridIterateRange(VariableDescriptions const& varDesc, Functions const& functions, Spaces const& spaces, Functor& f, CellRange cells)
    {
      using namespace boost::fusion;
      using namespace Dune;
      
      typedef typename SpaceType<Spaces,0>::type::Grid Grid;
      typedef typename SpaceType<Spaces,0>::type::GridView GridView;
      typedef typename SpaceType<Spaces,0>::type::Scalar Scalar;
      typedef typename Grid::ctype CoordType;
      int const dim = Grid::dimension;
      
      GridView const& gridView = at_c<0>(spaces)->gridView();
      
      // Shape function cache. Remember that every thread has to use its own cache!
      typedef ShapeFunctionCache<Grid,Scalar> SfCache;
      SfCache sfCache;
      
      // Quadrature rule cache
      typedef Dune::QuadratureRule<CoordType,dim> QuadRule;
      QuadratureTraits<QuadRule> quadratureCache;
      
      // Evaluators for shape functions of all FE spaces. Remember that
      // every thread has to use its own collection of evaluators!
      auto evaluators = getEvaluators(spaces,&sfCache);
      
      // Iterate over all cells.
      for (auto ci=cells.begin(); ci!=cells.end(); ++ci) {
        
        // for all spaces involved, compute the shape functions and
        // their global indices, which are needed for evaluating the
        // functional's derivative.
        auto idx = gridView.indexSet().index(*ci);
        moveEvaluatorsToCell(evaluators,*ci,idx);
        
        // loop over all quadrature points and evaluate variables
        int const p = f.integrationOrder(*ci,maxOrder(evaluators));
        
        QuadRule const& qr = quadratureCache.rule(ci->type(),p);
        useQuadratureRuleInEvaluators(evaluators,qr,0);
        
        size_t nQuadPos = qr.size();
        for (size_t g=0; g<nQuadPos; ++g) {
          // pos of integration point
          auto const& quadPos = qr[g].position();
          
          // for all spaces involved, update the evaluators associated
          // to this quadrature point
          moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qr,g,0);
          
          // evaluate functions and call functor
          f(ci,idx,quadPos,evaluateFunctions<VariableDescriptions>(functions,evaluators,valueMethod),
            ci->geometry().integrationElement(quadPos)*qr[g].weight());
        }
      }
    };
    
  } // end of namespace GridIterateDetail

  /**
   * \brief A function that supports general iterations over a spatial domain and
   * supports efficient evaluation of FE functions at the iterated
   * quadrature points.
   *
   * \param varDesc a heterogeneous container of VariableDescription structures
   *
   * \param functions a heterogeneous container of FE functions or function
   * views that can be evaluated given an Evaluator of the associated space.
   *
   * \param spaces a heterogeneous container of pointers to FE function spaces.
   * VariableDescription classes that gives the mapping of functions in
   * data to their associated spaces.
   *
   * \param f the functor that does the work. On each quadrature point,
   * is operator() is called with a heterogeneous container of the
   * values of the functions in data.
   */
  template <class VariableDescriptions, class Functor, class Functions, class Spaces>
  void gridIterate(VariableDescriptions const& varDesc, Functions const& functions, Spaces const& spaces, Functor& f)
  {
    using namespace boost::fusion;
    auto const& cellRanges = deref(begin(functions)).space().gridManager().leafCellRanges(); // WARNING: what if the grid view is not leaf?
    
    int n = std::min(NumaThreadPool::instance().cpus(),cellRanges.maxRanges());
    std::vector<Functor> fs(n,f);

    parallelFor([&varDesc,&functions,&spaces,&cellRanges,&fs](int i, int n) 
    {
      GridIterateDetail::gridIterateRange(varDesc,functions,spaces,fs[i],cellRanges.range(n,i));
    },n);
    
    for (auto& g: fs)
      f.join(g);
  }
  
  //---------------------------------------------------------------------

  /**
   * \brief A trivial Scaling.
   *
   * The identity scaling which leaves its operand simply as it is. This
   * is a convenience class to evaluate plain, unscaled norms.
   */
  struct IdentityScaling
  {
    template <class Cell, class F>
    void scale(Cell const&,
               Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const&,
               F&) const {}
  };


  /**
   * \brief A \ref Collector that sums up the weighted contributions.
   * 
   * This realizes a plain integration.
   */
  struct SummationCollector
  {
    template <class Cell>
    int integrationOrder(Cell const& /* ci */, int shapeFunctionOrder) const
    {
      return shapeFunctionOrder;
    }

    template <class CellIterator, class Index, class Sequence>
    void operator()(CellIterator const&, Index, double weight, Sequence const& x) 
    {
      if (sum.empty())
        sum.resize(boost::fusion::size(x),0);
      auto i = sum.begin();                           // provide a reference to the iterator
      boost::fusion::for_each(x,Add(weight,i));       // because that is incremented on each call of Add
    }                                                 // (hack: for_each requires an immutable functor)
    
    void join(SummationCollector const& c)
    {
      if (sum.empty())
        sum = c.sum;
      else
      {
        assert(sum.size()==c.sum.size());
        for (int i=0; i<sum.size(); ++i)
          sum[i] += c.sum[i];
      }
    }

    std::vector<double> sum;

  private:
    struct Add
    {
      Add(double w_, std::vector<double>::iterator& i_): w(w_), i(i_) {}
      template <class T>
      void operator()(T const& t) const { *i += w*t; ++i; }
      
    private:
      double w;
      std::vector<double>::iterator& i;
    };
  };




  /**
   * \brief A functor for computing scaled \f$L^2\f$-norms of a set of (possibly vector-valued) functions. 
   * 
   * This is intended to be used with gridIterate. The plain \f$L^2\f$-norm is computed if the
   * IdentityScaling is provided.
   */
  template <class Functions, class Scaling, class Collector>
  class ScaledTwoNorm2Collector
  {
  public:
    
    /**
     * \brief Constructor.
     * 
     * \param scaling
     * \param collector 
     */
    ScaledTwoNorm2Collector(Scaling const& scaling_, Collector& collector_): scaling(scaling_), collector(collector_) {    }

    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder) const
    {
      return 2*collector.integrationOrder(cell,shapeFunctionOrder);
    }

    template <class CellIterator, class Index, class Sequence>
    void operator()(CellIterator const& ci, Index idx,
                    Dune::FieldVector<typename CellIterator::Entity::Geometry::ctype,CellIterator::Entity::dimension> const& pos,
                    Sequence const& seq, double weight)
    {
      using namespace boost::fusion;
      typename result_of::as_vector<Sequence>::type values(seq);
      scaling.scale(*ci,pos,values);
      collector(ci,idx,weight,transform(values,TwoNorm2()));
    }

    void join(ScaledTwoNorm2Collector const& c)
    {
      collector.join(c.collector);
    }
    
    Collector const& data() const { return collector; }

  private:
    // a boost::fusion functor computing the squared euclidean norm of given vectors
    struct TwoNorm2
    {
      template <class T> struct result {};

      template <class Vec>
      struct result<TwoNorm2(Vec)> { typedef  Dune::FieldVector<double,1> type; };

      template <class Vec>
      typename result<TwoNorm2(Vec)>::type operator()(Vec const& v) const { return v.two_norm2(); }
    };

    Scaling const& scaling;
    Collector collector; // keep collector by value as we write concurrently in different objects.
  };

  /**
   * \brief Evaluates the square of the scaled \f$L^2\f$-norms of a set of functions.
   * \tparam Variables a boost::fusion sequence of VariableDescription types
   * \tparam Functions a boost::fusion sequence of FE function types (or types providing the required interface subset)
   * \tparam Spaces
   * \tparam Scaling
   * \tparam Collector
   * 
   * \param varDesc a boost::fusion sequence of VariableDescription entries, e.g. obtained from VariableSetDescription::Variables()
   * 
   * The required interface subset of FE functions includes the value() and the space() methods.
   */
  template <class Variables, class Functions, class Spaces, class Scaling, class Collector>
  void scaledTwoNormSquared(Variables const& varDesc, Functions const& f, Spaces const& spaces,
                            Scaling const& scaling, Collector& sum)
  {
    ScaledTwoNorm2Collector<Functions,Scaling,Collector> collector(scaling,sum);
    gridIterate(varDesc,f,spaces,collector);
    sum = collector.data();
  }


  //---------------------------------------------------------------------

  /**
   * \cond internals
   */
  namespace relativeErrorDetail {

    template <class F>
    struct Difference
    {
      typedef typename F::Scalar RT;
      typedef RT Scalar;

      static int const Components = F::components;

      typedef typename F::Space     Space;
      typedef typename F::ValueType ValueType;

      Difference(F const& f1_, F const& f2_): f1(f1_), f2(f2_) { }

      ValueType value(typename Space::Evaluator const& evaluator) const
      {
        return f1.value(evaluator)-f2.value(evaluator);
      }
      
      Space const& space() const { return f1.space(); }

    private:
      F const& f1;
      F const& f2;
    };


    struct MakeDifference
    {
      template <class T> struct result {};

      template <class Pair>
      struct result<MakeDifference(Pair)>
      {
        typedef typename boost::fusion::result_of::value_at_c<Pair,0>::type T;

        typedef Difference<typename boost::remove_const<typename boost::remove_reference<T>::type>::type> type;
      };

      template <class Pair>
      typename result<MakeDifference(Pair)>::type operator()(Pair const& pair) const
      {
        return typename result<MakeDifference(Pair)>::type(boost::fusion::at_c<0>(pair),boost::fusion::at_c<1>(pair));
      }
    };


  } // End of namespace relativeErrorDetail
  /**
   * \endcond
   */


  /**
   * For each variable, this function computes the following pair of values:
   * \f[ (\|f_1-f_2\|,\|f_3\|) \f]
   * \tparam Variables a boost::fusion sequence of variable descriptions
   * \tparam OutIter an output iterator with value type std::pair(double,double)
   * \tparam Functions a boost::fusion sequence of finite element functions, referenced by the variable descriptions
   * \tparam Spaces a boost::fusion sequence of pointers to spaces, referenced by the variable descriptions
   * \tparam Scaling
   */
  template <class Variables, class OutIter, class Functions, class Spaces, class Scaling>
  void relativeError(Variables const& varDesc, Functions const& f1, Functions const& f2, Functions const& f3,
                     Spaces const& spaces, Scaling const& scaling, OutIter out)
  {
    using namespace boost::fusion;

    int const s = size(f1);
    SummationCollector sum;
    //   scaledTwoNormSquared(join(varDesc,varDesc),
    //                        join( f3, as_vector(transform(zip(f1,f2),relativeErrorDetail::MakeDifference())) ),
    //                        spaces,scaling,sum);
    scaledTwoNormSquared(join(varDesc,varDesc),
                         join(as_vector(transform(zip(f1,f2),relativeErrorDetail::MakeDifference())),f3),
                         spaces,scaling,sum);
    // To my understanding, the as_vector call in the statement above
    // should not be necessary. However, it prevents a (rare but
    // reproducible) segfault with gcc 4.3.4 -O2. Funny thing is, some
    // output statements to cerr have the same effect. The cause may be
    // either a wild pointer bug in the code exposed by -O2 or a code
    // generation bug in gcc.  Sigh. At least, since Difference objects
    // are lightweight, the explicit vector construction should not be a
    // performance penalty. ws 2009-12-03

    for (int i=0; i<s; ++i) {
      *out = std::make_pair(std::sqrt(sum.sum[i]),std::sqrt(sum.sum[i+s]));
      ++out;
    }
  }
} // end of namespace Kaskade


#endif
