/*
 * scalarproducts.hh
 *
 *  Created on: Sep 16, 2011
 *      Author: bzflubko
 */

#ifndef SCALARPRODUCTS_HH_
#define SCALARPRODUCTS_HH_

#include <cmath>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Kaskade
{
  namespace {

    /// Norm traits.
    /**
     * Provides the underlying scalar type of Type in type value_type and
     * a const iterator in type iterator
     */
    template <class Type>
    struct NormTraits{
      typedef typename Type::value_type return_type;
      typedef typename Type::const_iterator iterator;
    };

    template <class value_type, int rows>
    struct NormTraits<Dune::FieldVector<value_type,rows> >{
      typedef value_type return_type;
      typedef typename Dune::FieldVector<value_type,rows>::const_iterator iterator;
    };

    template <class value_type, int rows, int cols>
    struct NormTraits<Dune::FieldMatrix<value_type,rows,cols> >{
      typedef value_type return_type;
      typedef typename Dune::FieldMatrix<value_type,rows,cols>::const_iterator iterator;
    };

  } // end of anonymous namespace
  namespace LinAlg{
    /// Euclidean scalar product.
    /**
     * Implemented for stl-compliant classes providing an const iterator under const_iterator
     * and the underlying value type under value_type. For type identification the class NormTraits
     * may be specialized.
     */
    struct EuclideanScalarProduct{
      template <class Type>
      typename NormTraits<Type>::return_type operator()(Type const& v1, Type const& v2) const
      {
        typename NormTraits<Type>::return_type result(0);
        typename NormTraits<Type>::iterator iend = v1.end(), it1 = v1.begin(), it2 = v2.begin();
        for(; it1 != iend; ++it1, ++it2)
          result += (*it1) * (*it2);
        return result;
      }
    };

    // forward declaration
    template <class> struct Norm;

    /// Norm squared defined via ScalarProduct
    template <class ScalarProduct_>
    struct NormSquared{
      typedef ScalarProduct_ ScalarProduct;
      typedef Norm<ScalarProduct> AssociatedNorm;

      NormSquared(ScalarProduct const sp=ScalarProduct()) : scalarProduct(sp){}

      template <class Type>
      typename NormTraits<Type>::return_type operator()(Type const& type) const
      {
        return scalarProduct(type,type);
      }

    private:
      ScalarProduct const scalarProduct;
    };

    /// Norm defined via ScalarProdcut
    template <class ScalarProduct_>
    struct Norm{
      typedef ScalarProduct_ ScalarProduct;
      typedef NormSquared<ScalarProduct> AssociatedNormSquared;

      Norm(ScalarProduct const sp=ScalarProduct()) : normSquared(sp){}

      template <class Type>
      typename NormTraits<Type>::return_type operator()(Type const& type) const
      {
        return sqrt(normSquared(type));
      }

    private:
      NormSquared<ScalarProduct> normSquared;
    };

    /// Euclidean norm
    typedef Norm<EuclideanScalarProduct> EuclideanNorm;

    /// Euclidean norm squared
    typedef NormSquared<EuclideanScalarProduct> EuclideanNormSquared;

    /// Infinity norm for stl-container, Dune::FieldVector and Dune::FieldMatrix
    /**
     * May be used for other types as well. Therefore provide a specialization of
     * NormTraits<class> with type value_type for the underlying scalar type and
     * iterator as a const_iterator
     */
    struct InfinityNorm{
      template<class Type>
      typename NormTraits<Type>::return_type operator()(Type const& type) const
      {
        typename NormTraits<Type>::return_type result(0);
        typename NormTraits<Type>::iterator iend = type.end(), iter = type.begin();
        for(; iter != iend; ++iter) if(result < fabs(*iter)) result = fabs(*iter);
        return result;
      }
    };

    /// One norm for stl-container, Dune::FieldVector and Dune::FieldMatrix
    /**
     * May be used for other types as well. Therefore provide a specialization of
     * NormTraits<class> with type value_type for the underlying scalar type and
     * iterator as const_iterator.
     */
    struct OneNorm{
      template<class Type>
      typename NormTraits<Type>::return_type operator()(Type const& type) const
      {
        typename NormTraits<Type>::return_type result(0);
        typename NormTraits<Type>::iterator iend = type.end(), iter = type.begin();
        for(; iter != iend; ++iter)
        	result += fabs(*iter);
        return result;
      }
    };

  }
}
#endif /* SCALARPRODUCTS_HH_ */
