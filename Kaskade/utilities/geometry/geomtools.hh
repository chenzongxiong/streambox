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

#ifndef GEOMTOOLS_HH
#define GEOMTOOLS_HH

#include <dune/common/fvector.hh>

/**
 * @file
 * @brief  Some simple tools for geometric calculations. Please extend.
 * @author Lars Lubkoll
 */

namespace GeomTools{

    /// Compute cross product for Dune::FieldVector<Scalar,3>.
    /*
     * \param Scalar floating point type (= double/float)
     * \param v1 input vector
     * \param v2 input vector
     * \return v1 x v2
     */
    template<class Scalar>
    Dune::FieldVector<Scalar,3> crossProduct(Dune::FieldVector<Scalar,3> const& v1, Dune::FieldVector<Scalar,3> const& v2){
    	Dune::FieldVector<Scalar,3> result;
    	result[0] = v1[1]*v2[2] - v1[2]*v2[1];
      result[1] = v1[2]*v2[0] - v1[0]*v2[2];
      result[2] = v1[0]*v2[1] - v1[1]*v2[0];
      return result;
    }

    /// Compute cross product for Dune::FieldVector<Scalar,3>.
    /*
     * \param Scalar floating point type (= double/float)
     * \param v1 input vector
     * \param v2 input vector
     * \param result = v1 x v2
     */
    template<class Scalar>
    void crossProduct(Dune::FieldVector<Scalar,3> const& v1, Dune::FieldVector<Scalar,3> const& v2, Dune::FieldVector<Scalar,3> &result){
      result[0] = v1[1]*v2[2] - v1[2]*v2[1];
      result[1] = v1[2]*v2[0] - v1[0]*v2[2];
      result[2] = v1[0]*v2[1] - v1[1]*v2[0];
    }

    /// Normalize vector.
    /*
     * \param Vector Type of the input vector
     * \param vector vector to be normalized
     */
    template<class Vector>
    inline Vector normalize(Vector &vector){
      if(vector.two_norm() < 1e-16) return vector;
      vector /= vector.two_norm();
      return vector;
    }

    /// Normalize vector.
    /*
     * \param Vector Type of the input vector, must provide a method called "two_norm()" for normalization,
     * a copy constructor and /= operator.
     * \param vector vector to be normalized
     * \result normalized copy of vector
     */
    template<class Vector>
    inline Vector getNormalized(Vector const& vector){
    	if(vector.two_norm() < 1e-16) return vector;
    	Vector result(vector);
    	result /= result.two_norm();
    	return result;
    }

    /// Project Vector vec on plane given through planeNormal. No translation is performed.
    /*
     * It is assumed that the origin of the coordinate system lies in the plane (no translation of vec).
     *
     * \param Vector Type of the input vector, must provide standard vector space operations.
     * \param vec vector to project, will be replaced by the projection
     * \param planeNormal unit(!) normal vector of the plane
     */
    template<class Vector>
    inline void projectOnPlane(Vector &vec, Vector const& planeNormal){
      vec -= (vec*planeNormal)*planeNormal;
    }

    /// Project Vector vec on plane given through planeNormal. No translation is performed.
    /*
     * It is assumed that the origin of the coordinate system lies in the plane (no translation of vec).
     *
     * \param Vector Type of the input vector, must provide standard vector space operations.
     * \param vec vector to project, will be replaced by the projection
     * \param projectionUnitNormal unit(!) normal vector (projection direction)
     */
    template <class Vector>
    inline void project(Vector &vec, Vector const& projectionUnitNormal)
    {
      vec -=(vec*projectionUnitNormal)*projectionUnitNormal;
    }
}


#endif // GEOMTOOLS_HH

