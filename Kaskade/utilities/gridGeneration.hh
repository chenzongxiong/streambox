/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  Created on: 25.06.2013
 *      Author: bzflubko
 */

#ifndef KASKADE_GRID_GENERATION_HH_
#define KASKADE_GRID_GENERATION_HH_

#include <algorithm>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include "GridGeneration/cube.hh"

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace GridGeneration_Detail
  {
    /// Extract vector of vertices and tetrahedra's vertex ids.
    /**
     * Taking a list of cubes, which is each represented as union of a set of tetrahedra this function collects all the local vertices,
     * removes duplicates and stores them in the first entry of the result-pair. The second contains the connectivity of the tetrahedra
     * wrt. the returned global vertex list.
     *
     * \param cubes list of cubes
     * \return pair of vectors, vertices in the first, connectivity in the second entry
     */
//    template <class Scalar,int dim>
    template <class Type,
               int dim = std::conditional<std::is_same<Type,Cube<typename Type::Real> >::value,std::integral_constant<int,3>,std::integral_constant<int,2> >::type::value>
    std::pair<std::vector<Dune::FieldVector<typename Type::Real,dim>>, std::vector<std::vector<unsigned int>>> extractSimplices(std::vector<Type> const& cubes)
    {
      std::vector<Dune::FieldVector<typename Type::Real,dim> > vertices;
      std::vector<std::vector<unsigned int> > tetrahedra;

      for(Type const& cube : cubes)
      {
        auto const& localVertices = cube.getVertices();
        auto const& localTetrahedra = cube.getSimplices();
        size_t const nVertices = localVertices.size();

        // maps local vertex ids to global ones (i.e. the entry in 'vertices')
        std::vector<size_t> vertexIdMap(nVertices);
        for(size_t i=0; i<nVertices; ++i) 
          vertexIdMap[i] = addToVector(localVertices[i], vertices);

        // map local connectivity onto global one (wrt. numbering in vector 'vertices')
        for(std::vector<unsigned int> const& tet : localTetrahedra)
        {
          std::vector<unsigned int> newTet(tet.size());

          for(size_t i=0; i<newTet.size(); ++i) newTet[i] = vertexIdMap[tet[i]];
          tetrahedra.push_back(newTet);
        }
      }

      return std::make_pair(vertices,tetrahedra);
    }


    template <class Grid, bool isUGGrid>
    struct SetDefaultHeapSizeImpl
    {
      static void apply(size_t){}
    };

    template <class Grid>
    struct SetDefaultHeapSizeImpl<Grid,true>
    {
      static void apply(size_t heapSize)
      {
        Grid::setDefaultHeapSize(heapSize);
      }
    };


    /**
     * \brief Fill rectangle with squares.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Square<Scalar> > fillRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Dune::FieldVector<Scalar,2> const& dh, bool symmetric=true)
    {
      std::vector<Square<Scalar> > squares;
      size_t xi_end = (size_t)round(fabs(dx[0]/dh[0])),
          yi_end = (size_t)round(fabs(dx[1]/dh[1]));
      Dune::FieldVector<Scalar,2> dsquare(dh), x(x0);

      for(size_t xi=0; xi<xi_end; ++xi)
      {
        x[0] = x0[0] + xi*dh[0];
        for(size_t yi=0; yi<yi_end; ++yi)
        {
          x[1] = x0[1] + yi*dh[1];
          squares.push_back(Square<Scalar>(x,dsquare,symmetric));
        }
      }

      return squares;
    }

    /**
     * \brief Fill rectangle with squares.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Square<Scalar> > fillRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Scalar dh, bool symmetric=true)
    {
      return fillRectangle(x0,dx,Dune::FieldVector<Scalar,2>(dh),symmetric);
    }

    /**
     * \brief Fill cuboid with smaller cuboids.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side lengths of cuboids that are used to fill the cuboid
     * \param symmetric true: each cuboid provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Cube<Scalar> > fillCuboid(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx, Dune::FieldVector<Scalar,3> const& dh, bool symmetric=true)
    {
      std::vector<Cube<Scalar> > cubes;
      size_t xi_end = (size_t)round(fabs(dx[0]/dh[0])),
          yi_end = (size_t)round(fabs(dx[1]/dh[1])),
          zi_end = (size_t)round(fabs(dx[2]/dh[2]));
      Dune::FieldVector<Scalar,3> x(x0);

      for(size_t xi=0; xi<xi_end; ++xi)
      {
        x[0] = x0[0] + xi*dh[0];
        for(size_t yi=0; yi<yi_end; ++yi)
        {
          x[1] = x0[1] + yi*dh[1];
          for(size_t zi=0; zi<zi_end; ++zi)
          {
            x[2] = x0[2] + zi*dh[2];
            cubes.push_back(Cube<Scalar>(x,dh,symmetric));
          }
        }
      }

      return cubes;
    }


    /**
     * \brief Fill cuboid with cubes
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Cube<Scalar>> fillCuboid(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx, Scalar dh, bool symmetric=true)
    {
      return fillCuboid(x0,dx,Dune::FieldVector<Scalar,3>(dh),symmetric);
    }



    template<class Scalar, int dim>
    std::vector<Cube<Scalar>> fillTwoLayerCuboid(Dune::FieldVector<Scalar,dim> x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                                    Dune::FieldVector<Scalar,dim> const& dh1, Scalar dh2, bool symmetric=true)
    {
      auto cubes = fillCuboid(x0,dx1,dh1,symmetric);

      x0[dim-1] += dx1[dim-1];
      auto cubes2 = fillCuboid(x0,dx2,dh2,symmetric);
      for( Cube<Scalar> const& cube : cubes2 ) cubes.push_back(cube);

      return cubes;
    }

    template<class Scalar, int dim>
    std::vector<Cube<Scalar> > fillTwoLayerCuboid(Dune::FieldVector<Scalar,dim> x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                                    Dune::FieldVector<Scalar,dim> const& dh1, Dune::FieldVector<Scalar,dim> const& dh2, bool symmetric=true)
    {
      auto cubes = fillCuboid(x0,dx1,dh1,symmetric);

      x0[dim-1] += dx1[dim-1];
      auto cubes2 = fillCuboid(x0,dx2,dh2,symmetric);
      for( Cube<Scalar> const& cube : cubes2 ) cubes.push_back(cube);

      return cubes;
    }

    template <int dim>
    struct FillCuboid;

    template <>
    struct FillCuboid<3>
    {
      template <class Scalar, class EdgeLength>
      static std::vector<Cube<Scalar> > apply(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx, EdgeLength const& dh, bool symmetric=true)
      {
        return fillCuboid(x0, dx, dh, symmetric);
      }
    };

    template <>
    struct FillCuboid<2>
    {
      template <class Scalar, class EdgeLength>
      static std::vector<Square<Scalar> > apply(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, EdgeLength const& dh, bool symmetric=true)
      {
        return fillRectangle(x0, dx, dh, symmetric);
      }
    };


    /**
     * \brief Fill L-shaped domain with cubes
     *
     * \param x0 left, lower, front corner of the L-shaped domain
     * \param dx length of horizontal line
     * \param dy length of vertical line
     * \param thickness
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar, int dim>
    std::vector< typename std::conditional<dim==3,Cube<Scalar>,Square<Scalar> >::type > fillLShape(Dune::FieldVector<Scalar,dim> x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true)
    {
      // fill horizontal line
      Dune::FieldVector<Scalar,dim> d0(thickness);    d0[0] = dx;
      auto cubes = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // fill vertical line
      d0 = thickness;   d0[1] = dy;
      x0[1] += thickness[1];
      auto cubes2 = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // merge
      for(auto const& cube : cubes2) cubes.push_back(cube);

      return cubes;
    }

    /**
     * \brief Fill L-shaped domain with cubes
     *
     * \param x0 center on lower boundary of the T-shaped domain
     * \param dx length of horizontal line
     * \param dy length of vertical line
     * \param thickness
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar, int dim>
    std::vector< typename std::conditional<dim==3,Cube<Scalar>,Square<Scalar>>::type> fillTShape(Dune::FieldVector<Scalar,dim> x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true)
    {
      // fill vertical line
      x0[0] -= 0.5*thickness[0];
      Dune::FieldVector<Scalar,dim> d0(thickness);    d0[1] = dy - thickness[1];
      auto cubes = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // fill vertical line
      x0[0] += 0.5*thickness[0] - 0.5*dx;
      x0[1] += d0[1];
      d0[1] = thickness[1];   d0[0] = dx;
      auto cubes2 = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // merge
      for(auto const& cube : cubes2) cubes.push_back(cube);

      return cubes;
    }


//#ifdef HAVE_UG
    template <class Grid>
    using SetDefaultHeapSize = SetDefaultHeapSizeImpl<Grid,std::is_same<Grid,Dune::UGGrid<Grid::dimension> >::value>;
/*#else
    template <class Grid>
    using SetDefaultHeapSize = SetDefaultHeapSizeImpl<Grid,false>;
#endif*/
  } // end of namespace GridGeneration_Detail
  /**
   * \endcond
   */


  /// Extract simplicial grid from list of cubes.
  /**
   * From a given list of (conforming) cubes or rectangles, this creates a simplicial grid.
   * 
   * \ingroup gridInput
   * \tparam Grid the type of grid to create 
   * \tparam Type
   * 
   * Usage: 
   * \code
   * grid = createGrid<Grid>(GridGeneration_Detail::fillCuboid(x0,dx,dh,true))
   * \endcode
   */
  template <class Grid, class Type>
  std::unique_ptr<Grid> createGrid(std::vector<Type> const& cubes, size_t heapSize=0)
  {
    auto simplices = GridGeneration_Detail::extractSimplices(cubes);
    if(heapSize > 0) GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,Type::dim);

    for(Dune::FieldVector<double,Type::dim> const& vertex : simplices.first) factory.insertVertex(vertex);
    for(std::vector<unsigned int> const& elem : simplices.second) factory.insertElement(gt,elem);

    return std::unique_ptr<Grid>(factory.createGrid());
  }

  /**
   * \ingroup gridInput
   * \brief Creates a uniform simplicial grid on a rectangular or cuboid domain.
   * 
   * \tparam Grid type of grid to create
   * \tparam Scalar
   * \tparam dim spatial dimension (usually 2 or 3)
   *
   * \param x0 left, lower, front corner of the cuboid
   * \param dx side lengths
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar, int dim, class EdgeLength>
  std::unique_ptr<Grid> createCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx, EdgeLength const& dh, bool symmetric=true, size_t heapSize=0)
  {
    return createGrid<Grid>(GridGeneration_Detail::FillCuboid<dim>::apply(x0,dx,dh,symmetric),heapSize);
  }

  template <class Grid>
  std::unique_ptr<Grid> createUnitSquare(double dh = 1., bool symmetric=true, size_t heapSize=0)
  {
    typedef Dune::FieldVector<double,2> Vec;
    return createCuboid<Grid>(Vec(0.), Vec(1.), dh, symmetric, heapSize);
  }

  template <class Grid>
  std::unique_ptr<Grid> createUnitCube(double dh = 1., bool symmetric=true, size_t heapSize=0)
  {
    typedef Dune::FieldVector<double,3> Vec;
    return createCuboid<Grid>(Vec(0.), Vec(1.), dh, symmetric, heapSize);
  }


  /**
   * \ingroup gridInput
   * \brief Fill rectangle with squares.
   *
   * \param x0 left, lower, front corner of the cuboid
   * \param dx side lengths
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar>
  std::unique_ptr<Grid> createRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Scalar dh, bool symmetric=true, size_t heapSize=0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillRectangle(x0,dx,dh,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Fill L-shaped domain with cubes
   *
   * \param x0 left, lower, front corner of the L-shaped domain
   * \param dx length of horizontal line
   * \param dy length of vertical line
   * \param thickness
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createLShape(Dune::FieldVector<Scalar,dim> const& x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillLShape(x0,dx,dy,thickness,dh,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Fill L-shaped domain with cubes
   *
   * \param x0 center on lower boundary of the T-shaped domain
   * \param dx length of horizontal line
   * \param dy length of vertical line
   * \param thickness
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTShape(Dune::FieldVector<Scalar,dim> const& x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTShape(x0,dx,dy,thickness,dh,symmetric),heapSize);
  }

  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTwoLayerCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                             Dune::FieldVector<Scalar,dim> const& dh1, Scalar dh2, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTwoLayerCuboid(x0,dx1,dx2,dh1,dh2,symmetric),heapSize);
  }

  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTwoLayerCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                             Dune::FieldVector<Scalar,dim> const& dh1, Dune::FieldVector<Scalar,dim> const& dh2, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTwoLayerCuboid(x0,dx1,dx2,dh1,dh2,symmetric),heapSize);
  }

} // end of namespace Kaskade

#endif /* GRIDGENERATION_HH_ */
