/*
 * norm.hh
 *
 *  Created on: May 2, 2015
 *      Author: sunayana_ghosh
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef KASKADE_TIMESTEPPING_EULERSDC_NORM_HH_
#define KASKADE_TIMESTEPPING_EULERSDC_NORM_HH_

#include <numeric>
#include <vector>

//includes from Kaskade7.3
#include "utilities/linalg/scalarproducts.hh"

namespace Kaskade {

  //=======================================================================================================
  //  Implementation for enum class NormType
  //=======================================================================================================

  /**
   * \ingroup eulerSDC
   * \brief Enum class to define the norm type in class Norm.
   */

  //---C++11 strongly typed enums. Would be used in the work model to call different functions.
  enum class NormType {
    ONE_NORM,         //!< ONE_NORM
    MAX_NORM          //!< MAX_NORM
  };

  //=======================================================================================================
  //  Implementation abstract base class Norm
  //=======================================================================================================
  /**
   * \ingroup eulerSDC
   * \brief Abstract base class of various norms
   *
   * This class represents a generic base class for an arbitrary norm of a vector of X's.
   * It also implements certain functions necessary for integrate function in EulerSDC and are
   * dependent on the norm in use.
   *
   * Vector is usually a vector of type Dune::FieldVector and Matrix is of type Dune::FieldMatrix respectively.
   */
  template <class Vector>
  class Norm
  {
  public:
    using field_type = typename Vector::value_type;



    //pure virtual functions
    /**
     *
     * @param vecs : the member function value takes in a vector of Vectors and computes its norm.
     * @return It returns the norm, which is usually of type double.
     */
    virtual field_type value(std::vector<Vector> const& vecs) = 0;

    //every abstract base class has a virtual destructor
    virtual ~Norm (){}

  private:
    NormType norm_t;

  };

  //=================================================================================================
  //  Implementation for derive class one norm
  //=================================================================================================
  /**
   * \ingroup eulerSDC
   * \brief OneNorm class derived from the abstract base class Norm.
   *
   * Here the one norm represents the one norm of a vector of Vectors or a vector of Matrices.
   * For a vector of Vectors \f$ \mathbf{V} := \{\mathbf{v_1},\ldots, \mathbf{v_m} \} \f$, where
   * \f$ \mathbf{v_i} = \{v_{i,1},\ldots, v_{i,n}\} \f$. Such that
   * \f$ \| V \|_1 := \sum_{i=1}^m \|\mathbf{v_i}\|_1 \f$ where,
   * \f$ \|\mathbf{v_i}\|_1:= \sum_{j=1}^n |v_{ij}|  \f$.
   */

  template <class Vector>
  class OneNorm: public Norm<Vector>
  {
  public:
    using field_type = typename Vector::value_type;



    virtual field_type value(std::vector<Vector> const& vecs)
    {

      //Using Kaskade::LinAlg::OneNorm instead of one_norm for Dune::FieldVector since the one norm defined
      //in Kaskade works for stl-container, Dune::FieldVector and Dune::FieldMatrix making it more generic.
      return std::accumulate(vecs.begin(),vecs.end(),0.0,[](field_type init, Vector const& v)
          {
        Kaskade::LinAlg::OneNorm on;
        return init + on(v);
          });
    }

    virtual ~OneNorm () {}
  };

  //=======================================================================================================
  //  Implementation for derived class max norm
  //=======================================================================================================
  /**
   * \ingroup eulerSDC
   * \brief MaxNorm class derived from the abstract base class Norm.
   * Here the max norm represents the max norm of a vector of Vectors or a vector of Matrices.
   * For a vector of Vectors \f$ \mathbf{V} := \{\mathbf{v_1},\ldots, \mathbf{v_m} \}\f$, where
   * \f$ \mathbf{v_i} = \{v_{i1},\ldots, v_{i,n}\}\f$. Such that \f$ \|V\|_{\infty} := \max_{1 \leq i \leq m} \|\mathbf{v_i}\|_{\infty}\f$ such that
   * \f$\|\mathbf{v_i}\|_{\infty}:= \max_{1 \leq j \leq n} |v_{ij}| \f$.
   */
  template <class Vector>
  class MaxNorm : public Norm<Vector>
  {
  public:
    using field_type = typename Vector::value_type;

    virtual field_type value(std::vector<Vector> const& vecs)
    {
      //Using Kaskade::LinAlg::OneNorm instead of one_norm for Dune::FieldVector since the one norm defined
      //in Kaskade works for stl-container, Dune::FieldVector and Dune::FieldMatrix making it more generic.
      Kaskade::LinAlg::InfinityNorm in;
      field_type max = in(vecs[0]);
      for (auto i = 1u; i < vecs.size(); ++i)
        max = std::max(max, in(vecs[i]));

      return max;
    }

    virtual ~MaxNorm () {}
  };


} //end namespace Kaskade


#endif /* KASKADE_TIMESTEPPING_EULERSDC_NORM_HH_ */
