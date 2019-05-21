/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SCALAR_HH_
#define SCALAR_HH_

#include <complex>
#include <type_traits>

/// \internal 
// forward declarations
namespace Dune
{
 template <class Scalar, int n>        class FieldVector; 
 template <class Scalar, int n, int m> class FieldMatrix; 
}
/// \endinternal

namespace Kaskade {

  /**
   * \brief Helper class for working with scalar field types
   */
  template <class Scalar>
  class ScalarTraits {
  public:
    /**
     * \brief The real type on which the scalar field is based.
     * 
     * For real types such as \a float or \a double, this is just \a Scalar. For complex types
     * such as \a std::complex<T>, this is \a T.
     */
    typedef Scalar Real;
    
    /**
     * \brief Conversion to the real type, ignoring the imaginary part if nonzero.
     */
    static Real real(Scalar const& s) { return s; }
  };

  // Specialization for standard complex types
  template <class T>
  class ScalarTraits<std::complex<T> > {
  public:
    typedef T Real;
    
    static Real real(std::complex<T> const& z) { return z.real(); }
  };
  
  //-----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup utilities
   * \brief Reports the converted type.
   */
  template <class T, class Real>
  struct ConvertTo 
  {
    using type = std::conditional_t<std::is_convertible<T,Real>::value,Real,T>;
  };
  
  template <class T, int n, int m, class Real>
  struct ConvertTo<Dune::FieldMatrix<T,n,m>,Real>
  {
    using type = Dune::FieldMatrix<typename ConvertTo<T,Real>::type,n,m>;
  };
}

#endif
