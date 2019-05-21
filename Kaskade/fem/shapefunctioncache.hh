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

#ifndef SHAPEFUNCTIONCACHE_HH
#define SHAPEFUNCTIONCACHE_HH
#include <map>
#include <tuple>

#include "boost/multi_array.hpp"
#include "boost/signals2.hpp"

#ifndef KASKADE_SEQUENTIAL
#include "boost/thread.hpp"
#include "utilities/threading.hh"
#endif

#include <boost/mpl/int.hpp>
#include <boost/mpl/range_c.hpp>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

#include "dune/geometry/quadraturerules.hh"
#include "dune/istl/bvector.hh"

#include "fem/fetransfer.hh"
#include "fem/pshapefunctions.hh"


//---------------------------------------------------------------------
namespace Kaskade
{
  /**
   * \ingroup fem
   * \brief A class that stores values, gradients, and Hessians of evaluated
   * FE functions / test functions.
   *
   * This is the argument type used during assembly to provide test
   * function values and derivatives to variational functionals.
   */
  template <typename Scalar, int dim, int components=1>
  struct VariationalArg
  {
    VariationalArg(): gradient(derivative) {}
    
    /**
     * \brief Constructor
     * 
     * This initializes the value as provied and initializes derivative and hessian to zero.
     */
    VariationalArg(Dune::FieldVector<Scalar,components> const& value_)
    : value(value_), derivative(0), gradient(derivative), hessian(0) {}
    
    
    VariationalArg(Dune::FieldVector<Scalar,components> const& value_,
                   Dune::FieldMatrix<Scalar,components,dim> const& derivative_)
    : value(value_), derivative(derivative_), gradient(derivative), hessian(0) {}
    
    VariationalArg(Dune::FieldVector<Scalar,components> const& value_,
                   Dune::FieldMatrix<Scalar,components,dim> const& derivative_,
                   Tensor3<Scalar,components,dim,dim> const& hessian_)
    : value(value_), derivative(derivative_), gradient(derivative), hessian(hessian_) {}
    
    // required due to gradient backwards compatibility reference that has always to refer to the *own* derivatives vector
    VariationalArg(VariationalArg<Scalar,dim,components> const& other)
    : value(other.value), derivative(other.derivative), gradient(derivative), hessian(other.hessian) {}
    VariationalArg& operator=(VariationalArg<Scalar,dim,components> const& other) { value = other.value; derivative = other.derivative; hessian = other.hessian; return *this; }
    
    Dune::FieldVector<Scalar,components>     value;
    Dune::FieldMatrix<Scalar,components,dim> derivative;
    Dune::FieldMatrix<Scalar,components,dim>& gradient; // deprecated, reference for backwards compatibility only
    Tensor3<Scalar,components,dim,dim>       hessian;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \ingroup fem
   * \brief This class caches values and derivatives of shape functions at quadrature points.
   *
   * Limitations are:
   *
   * - The number of components of shape functions is limited (to one
   *   below the template parameter ComponentsEnd).
   *
   * - Since evaluating the shape functions may modify the internal cache
   *   data structure, this is NOT thread safe. Every thread should use
   *   its own ShapeFunctionCache.
   *
   * \tparam G the grid type 
   * \tparam T the scalar field type
   * 
   * \todo make shape function cache independent of grid type - only coordinate type and dimension are required
   */
  template <class G, class T, int ComponentsEnd=4>
  class ShapeFunctionCache
  {

  private:
    typedef typename boost::mpl::range_c<int,1,ComponentsEnd>::type ComponentRange;

    // This is a functor that, given a number of
    // components (a boost::mpl integer), defines an associative
    // container type with a tuple (shape function set address,
    // quadrature rule address) as key
    // and a twodimensional boost::multi_array of pairs of vector and
    // matrix with given components as value type.
    struct MakeMapType {
      template <class Int>
      auto operator()(Int i) 
      {
        static int const nComp = Int::value;

        typedef VariationalArg<T,G::dimension,nComp> EntryType;
        typedef boost::multi_array<EntryType,2> DataType;
        typedef ShapeFunctionSet<typename G::ctype,G::dimension,T,Int::value> Sfs;
        typedef std::tuple<Sfs const*,void const*> KeyType;
        return std::map<KeyType,DataType>();
      }
    };

    typedef typename boost::fusion::result_of::as_vector<
                typename boost::fusion::result_of::transform<ComponentRange,MakeMapType>::type
        >::type MapsType;

  public:

    /**
     * \brief Defines the type returned by evaluate(sfs,qr,subId).
     *
     * This is a two-dimensional array type to be indexed by integration point and shape function number. The entries
     * are VariationalArg structures containing shape function values, derivatives, and hessians at the integration points.
     */
    template <int nComp>
    struct DataType 
    {
      typedef typename std::result_of<MakeMapType(boost::mpl::int_<nComp>)>::type::mapped_type type;
    };
    
    /**
     * \brief Defines the type returned by evaluate(sfs,qr,ip,subId).
     *
     * This is a one-dimensional array type to be indexed by shape function number. The entries
     * are pairs of shape function values and derivatives at the integration points.
     */
    template <int nComp>
    struct LocalDataType 
    {
      typedef typename DataType<nComp>::type::template const_array_view<1>::type type;
    };

    
    /**
     * \brief Evaluate comp-th component of isf-th shape function of the shape
     * function set sfs ath the ip-th integration point of quadrature
     * rule qr.
     *
     * This method is explicitly NOT thread safe. Use one cache object per thread.
     */
    template <int nComp, int subDim>
    Dune::FieldVector<T,nComp> evaluateFunction(ShapeFunctionSet<typename G::ctype,G::dimension,T,nComp> const& sfs, int isf,
                                                Dune::QuadratureRule<typename G::ctype,subDim> const& qr, int ip, int subId)
    {
      assert(0<=isf && isf<sfs.size());
      assert(0<=ip && ip<qr.size());
      return get(sfs,qr,subId)[ip][isf].value;
    }

    /**
     * \brief Evaluate derivative of comp-th component of isf-th shape function
     * of the shape function set sfs ath the ip-th integration point of
     * quadrature rule qr.
     *
     * This method is explicitly NOT thread safe. Use one cache object per thread.
     */
    template <int nComp, int subDim>
    Dune::FieldMatrix<T,nComp,G::dimension> const& evaluateDerivative(ShapeFunctionSet<typename G::ctype,G::dimension,T,nComp> const& sfs, int isf,
                                                                      Dune::QuadratureRule<typename G::ctype,subDim> const& qr, int ip,
                                                                      int subId)
    {
      assert(0<=isf && isf<sfs.size());
      assert(0<=ip && ip<qr.size());
      return get(sfs,qr,subId)[ip][isf].derivative;
    }

    /**
     * \brief Returns the values of all shape functions at given integration point. 
     *
     * The return value can be indexed by shape function number,
     * yielding a pair of value and derivative of the shape function's
     * component.
     *
     * E.g., the value can be obtained by evaluate(sfs,qr,ip,subId)[isf].value
     *
     * This method is explicitly NOT thread safe. Use one cache object per thread.
     *
     * \param subId the index of the subentity in the cell (has to be between 0 and the maximum supported subentity number)
     */
    template <int nComp, int subDim>
    typename LocalDataType<nComp>::type
    evaluate(ShapeFunctionSet<typename G::ctype,G::dimension,T,nComp> const& sfs,
             Dune::QuadratureRule<typename G::ctype,subDim> const& qr, int ip, int subId)
    {
      typedef typename DataType<nComp>::type DType;
      typedef typename DType::index_range range;
      typename DType::index_gen indices;
      assert(sfs.size()>=0);
      return get(sfs,qr,subId)[indices[ip][range(0,sfs.size())]];
    }

  
   /**
    * \brief Returns the values of all shape functions. 
    *
    * The return value can be indexed by shape function number,
    * yielding a pair of value and derivative of the shape function's
    * component.
    *
    * E.g., the value can be obtained by evaluate(sfs,qr,subId)[ip][isf].value
    * 
    * This method is explicitly NOT thread safe. Use one cache object per thread.
    *
    * \param subId the index of the subentity in the cell (has to be between 0 and the maximum supported subentity number)
    */
   template <int nComp, int subDim>
   typename DataType<nComp>::type const& evaluate(ShapeFunctionSet<typename G::ctype,G::dimension,T,nComp> const& sfs,
                                                  Dune::QuadratureRule<typename G::ctype,subDim> const& qr, int subId) 
   {
     return get(sfs,qr,subId);
   }
   
   
 
  /**
   * \brief Reports the maximum subentity number that is allowed.
   */
  int maximumSubentityIndex() const 
  {
    return maxSubId;
  }

  private:
    // TODO: use boost::container::flat_map from boost 1.49 on
    static int const maxSubId = 11; // 12 edges in a cube - that's the upper limit (up to 3d, of course...)
    
    // For each subentity number, we maintain a different cache for direct access
    MapsType caches[maxSubId+1];
    
    template <int nComp, int subDim>
    typename DataType<nComp>::type const&
    get(ShapeFunctionSet<typename G::ctype,G::dimension,T,nComp> const& sfs, Dune::QuadratureRule<typename G::ctype,subDim> const& qr, int subId)
    {
      assert(0<=subId && subId<=maxSubId);
      using MapType = typename std::result_of<MakeMapType(boost::mpl::int_<nComp>)>::type;
      typedef typename MapType::key_type KeyType;
      typedef typename DataType<nComp>::type DType;
      MapType& cache = boost::fusion::at_c<nComp-1>(caches[subId]);

      typename MapType::const_iterator i = cache.find(KeyType(&sfs,&qr));
      
      // check whether the shape functions' values are already available
      if (i==cache.end()) {
        // No, they are not. Evaluate shape functions and store the result.
        Dune::FieldVector<T,G::dimension> pos;
        auto const& refElem = sfs.referenceElement();
        assert(sfs.size()>=0);
        DType data(boost::extents[qr.size()][sfs.size()]);
#ifndef KASKADE_SEQUENTIAL     
        if(G::dimension != subDim)
          for (int isf=0; isf<sfs.size(); ++isf)
            for (int ip=0; ip<qr.size(); ++ip) {
              {
                boost::mutex::scoped_lock lock(refElementMutex);
                pos = refElem.template geometry<G::dimension-subDim>(subId).global(qr[ip].position());
              }
              data[ip][isf].value = sfs[isf].evaluateFunction(pos);
              data[ip][isf].derivative = sfs[isf].evaluateDerivative(pos);
              data[ip][isf].hessian = sfs[isf].evaluate2ndDerivative(pos);
            }
        else
#endif
          for (int isf=0; isf<sfs.size(); ++isf)
            for (int ip=0; ip<qr.size(); ++ip) {
              pos = refElem.template geometry<G::dimension-subDim>(subId).global(qr[ip].position());
              data[ip][isf].value = sfs[isf].evaluateFunction(pos);
              data[ip][isf].derivative = sfs[isf].evaluateDerivative(pos);
              data[ip][isf].hessian = sfs[isf].evaluate2ndDerivative(pos);
            }

        // store the values, obtaining an iterator pointing to the newly inserted values
        i = cache.insert(cache.begin(),typename MapType::value_type(KeyType(&sfs,&qr),data));
      }
      
      // From here on, the values are in the map (at the position pointed to by i).
      return i->second;
    }
  };
} // end of namespace Kaskade
//---------------------------------------------------------------------
//---------------------------------------------------------------------


#endif
