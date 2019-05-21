/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * foreach.hh
 *
 * @brief: for_each implementations for easy iteration over grids using functors and lambdas.
 * TODO: Maybe some methods should be removed since not really needed?
 *
 *  Created on: 28.04.2013
 *      Author: Lars Lubkoll
 */

#ifndef KASKADE_FOREACH_HH_
#define KASKADE_FOREACH_HH_

#include <algorithm>

namespace Kaskade
{
  /// iterates over each cell and applies functor to cell. Each cell is visited exactly once.
  template <class GridView, class Functor>
  void forEachCell(GridView const& gridView, Functor functor)
  {
    for(auto const& element : elements(gridView)) {
      functor(element);
    }
  }

  /// iterates over each vertex and applies functor to vertex. Each vertex is visited exactly once.
  template <class GridView, class Functor>
  void forEachVertex(GridView const& gridView, Functor functor)
  {
    for(auto const& vertex : vertices(gridView)) {
      functor(vertex);
    }
  }

  /// iterates over each face and applies functor to face. Each boundary face is visited exactly once
  /// and each inner face is visited exactly twice (from each neighbouring cell) and so the functor is applied twice to inner faces.
  template <class GridView, class Functor>
  void forEachFace(GridView const& gridView, Functor functor)
  {
    for(auto const& element : elements(gridView)) {
      for(auto const& intersection : intersections(gridView,element)) {
        functor(intersection);
      }
    }
  }

  /// same as forEachFace, but the functor is applied to the face and the cell from which the face is visited
  template <class GridView, class Functor>
  void forEachFaceOfCell(GridView const& gridView, Functor functor)
  {
    for(auto const& element : elements(gridView)) {
      for(auto const& intersection : intersections(gridView,element)) {
        functor(element, intersection);
      }
    }
  }

//  template <class GridView, class Functor>
//  void forEachBoundaryFace(GridView const& gridView, Functor functor)
//  {
//    forEachFace(gridView, [&](typename GridView::Intersection const& face)
//    {
//      if(face.boundary()) functor(face);
//    });
//  }

  // seems to be faster, at least for UG grid (tested in experiment); for the old version see above
  /// iterates over each boundary face and applies functor to face. Each boundary face is visited exactly once.
  template <class GridView, class Functor>
  void forEachBoundaryFace(GridView const& gridView, Functor functor)
  {
    for(auto const& element : elements(gridView)) {
      if(element.hasBoundaryIntersections()) {// this if clause accelerates iteration over boundary faces
        for(auto const& intersection : intersections(gridView,element)) {
          if(intersection.boundary()) functor(intersection);
        }
      }
    }
  }

  /// iterates over each inner face and applies functor to face. Each inner face is visited exactly twice (and so functor is applied twice).
  template <class GridView, class Functor>
  void forEachInnerFace(GridView const& gridView, Functor functor)
  {
    forEachFace(gridView, [&](typename GridView::Intersection const& face)
    {
      if(!face.boundary()) functor(face);
    });
  }

  /// iterates over each cell and applies functor to cell and then iterates over each face of cell and applies functor to face
  template <class GridView, class Functor>
  void forEachCellAndEachFace(GridView const& gridView, Functor functor)
  {
    for(auto const& element : elements(gridView)) {
      functor(element);
      for(auto const& intersection : intersections(gridView,element)) {
        functor(intersection);
      }
    }
  }

  template <class GridManager, class Functor>
  void markCells(GridManager& gridManager, Functor functor)
  {
    forEachCell(gridManager.grid().leafView(), [&](typename GridManager::Grid::LeafGridView::template Codim<0>::Entity const& cell)
    {
#ifdef TESTOUTPUT
      size_t marked = 0, notMarked = 0;
#endif

      bool mark = false;
      mark = mark || functor(cell);
      std::for_each(gridManager.grid().leafView().ibegin(cell), gridManager.grid().leafView().iend(cell), [&](typename GridManager::Grid::LeafGridView::Intersection const& face)
      {
        mark = mark || functor(face);
      });
      if(mark) gridManager.mark(1,cell);

#ifdef TESTOUTPUT
      if(mark) ++marked;
      else ++notMarked;
      std::cout << "Marked " << marked << " cells of " << (marked+notMarked) << " cells." << std::endl;
#endif
    });
  }
} // end of namespace Kaskade


#endif /* KASKADE_FOREACH_HH_ */
