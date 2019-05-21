/*
 * gridtree.hh
 *
 *  Created on: 16.12.2011
 *      Author: Lars
 */

#ifndef GRIDTREE_HH_
#define GRIDTREE_HH_

#include <algorithm>
#include <utility>
#include <vector>
#include "tools/linalg/scalarproducts.hh"
#include "gridtreepolicy.hh"



/// Stores a grid in a tree structure, allowing fast searches.
/**
 * \param GRID
 * \param BoundingBoxType type of the bounding box. Provided are BoundingBox and FastBoundingBox. The last inherits BoundingBox and additionally stores its corners, edges and faces
 *        for fast access.
 * \param lowerSplittingBorder if a node gets less cells it will not be split in grid tree creation process
 * \param SplittingPolicy policy that determines the way new tree nodes are created (adjust this only if you have to much time and nothing else to do)
 */
template <class GRID, template <class,int> class BoundingBoxType=GeometricObject::BoundingBox, int lowerSplittingBorder=1000, template <class,class> class SplittingPolicy=BisectionPolicy>
class GridTree {
private:
  /// grid treee node
  template <class Cell, template <class,int> class LocalBoundingBoxType, int localLowerSplittingBorder, template <class,class> class LocalSplittingPolicy=BisectionPolicy>
  class GridTreeNode : public LocalSplittingPolicy<Cell,LocalBoundingBoxType<typename Cell::ctype, Cell::dimension> >
  {
  public:
    static int const dim = Cell::dimension;
    typedef typename Cell::ctype Scalar;
    typedef typename Cell::Geometry::GlobalCoordinate GlobalCoordinate;
    typedef GetProtectedConstructor<Cell> GPC;
    typedef GridTreeNode<Cell,LocalBoundingBoxType,localLowerSplittingBorder,LocalSplittingPolicy> This;
    typedef LocalBoundingBoxType<Scalar,dim> BoundingBox;

    GridTreeNode(BoundingBox boundingBox_, std::vector<GPC const*> const& cells_, Direction dir, Scalar const reltol, Scalar const height, Scalar const maxHeight, bool isRootNode=false)
    : boundingBox(boundingBox_)
    {
      initialize(cells_, dir, reltol, height, maxHeight);
      if(isRootNode) cells = cells_;
      BoundingBoxWrapper<BoundingBox>::init(boundingBox);
      isRoot=isRootNode;
    }

    GridTreeNode(This const& node)
    : cells(node.cells), boundingBox(node.boundingBox), children(node.children)
    {}

    ~GridTreeNode()
    {
      if(children.first) delete children.first;
      if(children.second) delete children.second;
      if(isRoot){
        typename std::vector<GPC const*>::iterator iter = cells.begin(), iend = cells.end();
        for(;iter!=iend; ++iter) delete *iter;
      }
    }

    template <class Grid, class Scalar>
    static This createRootNode(Grid const& grid, Scalar const maxH, Scalar const reltol)
    {
      std::vector<GPC const*> cells_(grid.size(0));
      typedef typename Grid::LeafGridView::template Codim<0>::Iterator CellIterator;
      CellIterator ci = grid.leafView().template begin<0>(), cend = grid.leafView().template end<0>();
      BoundingBox boundingBox(ci->geometry().corner(0));

      int cellIndex = 0, corners;
      for(;ci!=cend; ++ci)
      {
        (cells_[cellIndex]) = new GPC(*ci);
        ++cellIndex;
        corners = ci->geometry().corners();
        for(int cornerId=0; cornerId<corners; ++cornerId) boundingBox.update(ci->geometry().corner(cornerId));
      }
      /*std::for_each(grid.leafView().begin<0>(), grid.leafView.end<0>(), [](Cell const &cell)
          {
            cells_[cellIndex] = &cell;
            ++cellIndex;
            int corners = cell.geometry().corners();
            for(int cornerId=0; cornerId<corners; +cornerId) boundingBox.update(cell.geometry().corner(cornerId));
          });*/
      BoundingBoxWrapper<BoundingBox>::init(boundingBox);
      return This(boundingBox, cells_, X, reltol, 0, maxH, true);
    }

    bool checkInside(GlobalCoordinate const& x) const
    {
      if(!boundingBox.contains(x)) return false;

      if(children.first && children.second) return ( children.first->checkInside(x) || children.second->checkInside(x) );

      typename std::vector<GPC const*>::const_iterator iter=cells.begin(), iend=cells.end();
      const Dune::GenericReferenceElement<double,3> *simplex_ptr = &Dune::GenericReferenceElements<double,3>::simplex();
      for(;iter!=iend; ++iter)
      {
        if(simplex_ptr->checkInside((*iter)->geometry().local(x))) return true;
      }
      return false;
    }

    Cell const* getCell(GlobalCoordinate const& x) const
    {
      if(!boundingBox.contains(x)) return 0;
      if(children.first && children.second) return (children.first->getCell(x)==0) ? children.second->getCell(x) : children.first->getCell(x);
      typename std::vector<GPC const*>::const_iterator iter=cells.begin(), iend=cells.end();

      for(;iter!=iend; ++iter) if((*iter)->checkInside((*iter)->geometry().local(x))) return *iter;
      //else std::for_each(cells.begin(), cells.end(), [](Cell const* cell){ if(cell->checkInside(cell->geometry().local(x))) return cell; });
      return 0;
    }

    void print(std::ostream &stream, int const level) const
    {
      using std::endl;
      stream << endl;
      stream << "Level: " << level << endl;
      stream << boundingBox;
      stream << "Number of cells: " << cells.size() << endl;
      stream << "Has first child: " << (children.first!=0) << endl;
      stream << "Has second child: " << (children.second!=0) << endl;
      stream << endl;
    }

    void print(std::ostream &stream, int const level, int const printLevel ) const
    {
      if(level==printLevel){
        print(stream, level);
        return;
      }
      if(children.first) children.first->print(stream,level+1,printLevel);
      if(children.second) children.second->print(stream,level+1,printLevel);
    }

    int getTreeHeight() const
    {
      return (children.first && children.second) ? 1+std::max(children.first->getTreeHeight(),children.second->getTreeHeight()) : 0;
    }

    std::vector<GPC const*> getCellsInBall(Ball<Scalar,dim> const& ball) const
    {
      if(GeometricObject::intersects<LinAlg::EuclideanNorm>(ball, boundingBox))
      {
        // if children exist delegate to them
        if(children.first && children.second)
        {
          std::vector<GPC const*> result = children.first->getCellsInBall(ball),
                                     tmp = children.second->getCellsInBall(ball);
          result.insert(result.begin(), tmp.begin(), tmp.end());
          return result;
        }
        else // else check list of cells
        {
          typename std::vector<GPC const*>::const_iterator iter=cells.begin(), iend=cells.end();

        }
      } // end if -> ball lies completely outside the bounding box
      // else
      return std::vector<GPC const*>();
    }

    std::vector<GPC const*> cells;
    BoundingBox boundingBox;
    std::pair<This*,This*> children ;
  private:
    bool isRoot;
    /// not implemented
    GridTreeNode();
    /// not implemented
    This operator=(This const&);

    void initialize(std::vector<GPC const*> const& cells_, Direction dir, Scalar const reltol, Scalar const height, Scalar const maxHeight)
    {
      if(cells_.size() > localLowerSplittingBorder && (height<maxHeight)) createChildren(boundingBox, cells_, changeDirection(dir), children, reltol, height, maxHeight);
      else
      {
        children.first=0;
        children.second=0;
        cells = cells_;
      }
    }
  };

  typedef GridTreeNode<typename GRID::LeafGridView::template Codim<0>::Entity,BoundingBoxType,lowerSplittingBorder,SplittingPolicy> NodeType;

public:
  typedef GRID Grid;
  typedef typename Grid::LeafGridView::template Codim<0>::Entity Cell;
  typedef typename Cell::ctype Scalar;
  typedef typename Cell::Geometry::GlobalCoordinate GlobalCoordinate;
  typedef BoundingBoxType<Scalar,Grid::dimension> BoundingBox;

  /// Constructor
  /**
   * \param grid
   * \param maxHeight maximal grid tree height
   * \param reltol determine allowed relative difference in the number of cells of a nodes children (optional)
   */
  GridTree(Grid const& grid, Scalar const maxHeight, Scalar const reltol = 0.2 ) : rootNode(NodeType::createRootNode(grid, maxHeight, reltol))
  {}

  /// Copy constructor
  GridTree(GridTree<Grid,BoundingBoxType,lowerSplittingBorder,SplittingPolicy> const& gt) : rootNode(gt.rootNode)
  {}

  ~GridTree(){}

  /// true if tree contains a cell containing x, else false
  bool contains(GlobalCoordinate const& x) const{ return rootNode.checkInside(x); }

  /**
   * \param x coordinate for which a cell containing it is searched
   * \result cell containing x if tree contains x, else a null pointer
   */
  Cell const* getCell(GlobalCoordinate const& x) const{ return rootNode.getCell(x); }

  int height() const
  {
    return rootNode.getTreeHeight();
  }

  void print() const
  {
    int h = height();
    for(int level=0; level<h; ++level)
      rootNode.print(std::cout,0,level);
  }

  BoundingBox const& getBoundingBox() const
  {
    return rootNode.boundingBox;
  }

// @TODO
//  std::vector<Cell const*> getCellsInBall(Ball<Scalar,Grid::dimension> &ball) const
//  {
//    return rootNode.getCellsInBall(ball);
//  }

private:
   NodeType const rootNode;
};

#endif /* GRIDTREE_HH_ */
