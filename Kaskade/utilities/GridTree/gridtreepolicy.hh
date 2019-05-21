/*
 * gridtreepolicy.h
 *
 *  Created on: 17.12.2011
 *      Author: Lars
 */

#ifndef GRIDTREEPOLICY_HH_
#define GRIDTREEPOLICY_HH_

#include <cmath>
#include <utility>
#include <vector>
#include <dune/common/fvector.hh>
#include "utilities/get_protected_constructor.hh"
#include "utilities/geometry/geometric_objects.hh"

namespace {
  using namespace GeometricObject;

  /// helper class storing information for current search direction (i.e. for bisection)
  template <class Scalar, int dim>
  struct SearchDetail{
    SearchDetail(Scalar const f, Scalar const l, Scalar const s, Direction const dir) : first(f), last(l), split(s), delta(last-split), direction(dir)
    {}

    template <class BoundingBoxType>
    SearchDetail(BoundingBoxType const& boundingBox, Direction dir) :
      first(boundingBox.coord[dir].first), last(boundingBox.coord[dir].second),
      split(0.5*(boundingBox.coord[dir].first+boundingBox.coord[dir].second)), delta(split-first), direction(dir)
    {}

    Scalar const first, last;
    Scalar split, delta;
    Direction const direction;
  };

  /**
   * Create a new bounding box using the existing BoundingBox boundingBox and the additional information of searchDetail.
   * The two bounding boxes created using lower=true and lower=false cover the same domain as boundingBox does.
   */
  template <class Scalar, int dim, class BoundingBoxType>
  BoundingBoxType createBoundingBox(BoundingBoxType const& boundingBox, SearchDetail<Scalar,dim> const& searchDetail, bool const lower)
  {
    BoundingBoxType newBoundingBox;
    if(lower) newBoundingBox.coord[searchDetail.direction].first = boundingBox.coord[searchDetail.direction].first,
              newBoundingBox.coord[searchDetail.direction].second = searchDetail.split;
    else newBoundingBox.coord[searchDetail.direction].first = searchDetail.split,
         newBoundingBox.coord[searchDetail.direction].second = boundingBox.coord[searchDetail.direction].second;

    size_t dir = (1+searchDetail.direction)%dim;
    newBoundingBox.coord[dir].first = boundingBox.coord[dir].first, newBoundingBox.coord[dir].second = boundingBox.coord[dir].second;
    dir = (2+searchDetail.direction)%dim;
    newBoundingBox.coord[dir].first = boundingBox.coord[dir].first, newBoundingBox.coord[dir].second = boundingBox.coord[dir].second;
    return newBoundingBox;
  }

  /// Helper function for adding a cell to the left or right part of a partition
  template <class Cell, class Partition>
  void addCell(Cell const& cell, SearchDetail<typename Cell::ctype, Cell::dimension> const searchDetail, Partition &partition)
  {
    bool left = false, right = false;

    for(int cornerId=0; cornerId<cell.geometry().corners(); ++cornerId)
      (cell.geometry().corner(cornerId)[searchDetail.direction] < (searchDetail.split)) ? left = true : right = true;

    if(left) partition.first.push_back(&cell);
    if(right) partition.second.push_back(&cell);
  }

  template <class Partition>
  void swapPartition(Partition &p1, Partition &p2)
  {
    std::swap(p1.first,p2.first);
    std::swap(p1.second,p2.second);
  }

  template <class Partition>
  void clear(Partition &partition)
  {
    partition.first.clear();
    partition.second.clear();
  }

  /// Get difference in the partitions sizes
  template <class Partition>
  size_t getDifference(Partition const& partition)
  {
    return (size_t) std::abs((int)partition.first.size()-(int)partition.second.size());
  }

  /// check if difference is admissible according to an absolute threshold value.
  template <class Partition, class Scalar>
  bool checkPartition(Partition const& partition, Scalar const difference)
  {
     return (getDifference(partition) < difference);
  }

  Direction changeDirection(Direction dir)
  {
    if(dir==X) return Y;
    if(dir==Y) return Z;
    return X;
  }
}

/// The default splitting policy
template <class Cell, class BoundingBox>
struct BisectionPolicy {
  static int const dim = Cell::dimension;
  typedef typename Cell::ctype Scalar;
  typedef GetProtectedConstructor<Cell> GPC;
  typedef std::pair<std::vector<GPC const*>, std::vector<GPC const*> > Partition;

  // create tree node's children
  template <class TreeNode>
  void createChildren(BoundingBox const& boundingBox, std::vector<GPC const*> const& cells_, Direction dir, std::pair<TreeNode*, TreeNode*> &children, Scalar const reltol, Scalar const height, Scalar const maxHeight)
  {

    // split cell vector
    Scalar const MINIMAL_DELTA = 0.01;
    size_t const numberOfCells = cells_.size();
    size_t const tolerance = std::max((size_t)std::fabs(numberOfCells*reltol), (size_t)25);
    SearchDetail<Scalar,dim> searchDetail(boundingBox, dir);
    Partition partition, tmpPartition;

    splitCells(searchDetail, cells_, partition);
    size_t difference = getDifference(partition);
    while(difference > tolerance && searchDetail.delta > MINIMAL_DELTA)
    {
      searchDetail.delta *= 0.5;
      searchDetail.split += searchDetail.delta;
      splitCells(searchDetail, cells_, tmpPartition);
      if(checkPartition(tmpPartition, difference))
      {
        swapPartition(partition,tmpPartition);
        difference = getDifference(partition);
        continue;
      }
      else
      {
        clear(tmpPartition);
        searchDetail.split -= 2.*searchDetail.delta;
        splitCells(searchDetail, cells_, tmpPartition);
        if(checkPartition(tmpPartition, difference))
        {
          swapPartition(partition,tmpPartition);
          difference = getDifference(partition);
        }
        else searchDetail.split += searchDetail.delta;
      }
      clear(tmpPartition);
    }

    // create children
    children.first = new TreeNode(createBoundingBox(boundingBox,searchDetail,true), partition.first, changeDirection(searchDetail.direction), reltol, height+1, maxHeight);
    children.second = new TreeNode(createBoundingBox(boundingBox,searchDetail,false), partition.second, changeDirection(searchDetail.direction), reltol, height+1, maxHeight);
  }

  void splitCells(SearchDetail<Scalar,dim> const& searchDetail, std::vector<GPC const*> const& cells_, Partition &partition) const
  {
    typename std::vector<GPC const*>::const_iterator ci=cells_.begin(), cend=cells_.end();
    for(;ci!=cend; ++ci) addCell(**ci, searchDetail, partition);
    /*std::for_each(cells_.begin(), cells_.end(), [](Cell const* cell)
        {
          addCell(*cell, searchDetail, partition);
        });*/
  }
};

#endif /* GRIDTREEPOLICY_HH_ */
