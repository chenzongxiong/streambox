/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HERMITE_INTERPOLATION_BASE
#define HERMITE_INTERPOLATION_BASE

#include <vector>
#include <dune/common/fvector.hh>

namespace Kaskade
{
  /**
   * @file hermiteinterpolationbase.hh
   * @brief Base class for hermite interpolation assuming zero values on the interpolation nodes and using
   * 			prescribed surface normals.
   * @author Lars Lubkoll
   */

  /// Base class for hermite interpolation.
  template<class Scalar_, int dimension>
  class HIPBase
  {
  public:
    typedef Scalar_ Scalar;
    static int const dim = dimension;
    typedef Dune::FieldVector<Scalar,dim> Vector;

    ~HIPBase(){}

    /// Evaluate interpolation polynomial at position x (in local coordinates).
    /*
     * \param entity entity containing x
     * \param x evaluation point in local coordinates
     * \param shapeFunctionSet shape function set of the reference entity
     * \result interpolation polynomial evaluated at x
     */
    template <class ShapeFunctionSet>
    Vector evaluate(Vector const &x, ShapeFunctionSet const& shapeFunctionSet) const
    {
      Vector result(0);
      for(int sfId=0; sfId<shapeFunctionSet.size(); ++sfId) result += values[sfId] * shapeFunctionSet[sfId].evaluateFunction( x )[0] * directions[sfId];
      return result;
    }

  protected:
    /// Constructor.
    explicit HIPBase(int numberOfShapeFunctions)
    {
      values = std::vector<Scalar>(numberOfShapeFunctions);
      directions = std::vector<Vector>(numberOfShapeFunctions,Vector(0));
    }

    /// Store value and direction associated with the id's shape function.
    void insertEntry(Scalar const value, Vector const& dir, int const id){
      values[id] = value;
      directions[id] = dir;
    }

    std::vector<Scalar> values;
    std::vector<Vector> directions;
  };


  // class declaration
  /// Implementation of interpolation polynomials from given vertex and/or edge directions ("normals").
  /*
   * Only 2D and 3D implemented.
   * In 2D and 3D for every vertex a "normal" must be given.
   * In 3D also "normals" for the vertices are necessary.
   * The given directions are interpreted as desired normals of the surface.
   * The polynomials calculate the necessary deformations of the vertices. In 2D this is done exactly.
   * In 3D the displacement is exact at the edges. On the faces it is approximated as deformation in direction
   * of the actual normal of the face.
   */
  template<class Scalar, class GridView, class InterpolationPolicy, int dimension> class HIPImpl;
} /* end of namespace Kaskade */
#endif
