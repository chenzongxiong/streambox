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

#ifndef MLLGEOMETRY_HH
#define MLLGEOMETRY_HH

/**
 * @file
 * @brief  Class for transformations between ancestors and their children
 * @author Anton Schiela
 */

#include <limits>
#include <stack>

#include "dune/common/fmatrix.hh"
#include <dune/geometry/referenceelements.hh>

#include "fem/fixdune.hh"

namespace Kaskade
{
  /** 
   * \brief Transformation of coordinates between ancestor and child.
   * 
   * Provides a convenient transfer between a cell and its (possibly higher degree) ancestor.
   * 
   * This defines a mapping \f$ g:L\to G \f$, where \f$ L \f$ is the "local" grid cell and 
   * \f$ G \f$ the "global" grid cell. The evaluation of the Jacobian \f$ J = g'(x) \f$ is also supported.
   * 
   * It is possible to choose, which of child and ancestor are to be considered as "global".
   * Moreover, it is possible to provide a flag that decides if (chained) transformations
   * should be performed via geometry.global/local, or via localGeometry.global/local. The
   * first option is more efficient for high degree ancestors, but it is only sensible,
   * if the child is actually geometrically a subset of the father. This requirement might be
   * violated in some special cases. In general the transformation from child to father is more
   * efficient than the other way round. This has to do with the data structures in the grid.
   * @todo: implement pull-back and pull-forward operation for gradients
   */
  template <class Grid>
  class MultiLevelLocalGeometry
  {
    // Changed in DUNE 2.1
    //typedef typename Grid::template Codim<0>::GlobalIdSet::IdType Id;
    typedef typename Grid::GlobalIdSet::IdType Id;
    typedef typename Grid::template Codim<0>::EntityPointer       CellPointer;

  protected:
    CellPointer child, ancestor;

  private:

    // Changed in DUNE 2.1
    // typename Grid::template Codim<0>::GlobalIdSet const& idSet;
    typename Grid::GlobalIdSet const& idSet;
    int spanwidth;
    bool geometricallyNested, transformViaGlobal;

  public:

    enum Direction{
      ChildIsGlobal,
      FatherIsGlobal
    };

    /// The entity type between which this geometry object maps.
    typedef typename Grid::template Codim<0>::Entity Entity;


  private:

    Direction dir;

    typedef typename Grid::ctype ctype;
    static int const dim = Grid::dimension;
    static int const dimw = Grid::dimensionworld;

  public:

    /** Constructor Arguments:
     * 
     *  @param grid: Grid in which everthing takes place
     *  @param child_: Pointer to the child entity
     *  @param ancestor_: Pointer to the ancestor entity
     *  @param direction_: determine, which of both entities is seen as the global entitiy
     *  @param geometricallyNested_: true: assume that father and child are
     *  geometrically nested. This allows sometimes more efficient
     *  transformation, but is not alwasy feasible to assume this (for example
     *  in boundary adaptation).
     *  @param tol_: geometric tolerance for checking points to be inside cells or not. Use this to make sure
     *  that points inside are detected as such, even if by rounding errors they appear to be outside (by some
     *  negligible margin). Defaults to \f( 10 \epsilon \f).
     * 
     *  @todo: (i) Simplify code by using entity pointer's comparison
     *  operator instead of ids (caveat: can ids lead to shorter
     *  hierarchies? Remember that entities copied to higher grid levels
     *  have the same id. If the comparison operators between entity
     *  pointers evaluate to "not equal" in this case, comparing ids may be
     *  more efficient. (ii) Check whether spanwidth can be computed as
     *  difference of the entity levels instead of using the loop. The same
     *  consideration as for (i) applies.
     */
    MultiLevelLocalGeometry(Grid const &grid,
                            CellPointer const& child_,
                            CellPointer const& ancestor_,
                            MultiLevelLocalGeometry<Grid>::Direction direction_,
                            bool geometricallyNested_= (dim==dimw),
                            ctype tol_=-100):
      child(child_),
      ancestor(ancestor_),
      idSet(grid.globalIdSet()),
      geometricallyNested(geometricallyNested_),
      dir(direction_),
      tol(tol_<=-100? 10*std::numeric_limits<ctype>::epsilon(): tol_)
    {
      CellPointer ci(child);
      spanwidth = 0;
      while (idSet.id(*ci) != idSet.id(*ancestor))
      {
        ci = ci->father();
        spanwidth++;
      }
      transformViaGlobal = geometricallyNested && (spanwidth >= 3);
    }

    /// A reference to the local entity
    Entity const& localEntity() const { return dir==ChildIsGlobal? *ancestor: *child; }


    /// Transformation from local to global
    void global(std::vector<Dune::FieldVector<ctype,dim> > & x,
                std::vector<int> &componentsInside) const
    {
      if (dir==ChildIsGlobal) {
        if(transformViaGlobal) childToFatherViaGlobal(x);
        else                   childToFather(x);
      }
      else {
        if(transformViaGlobal) fatherToChildViaGlobal(x,componentsInside);
        else                   fatherToChild(x,componentsInside);
      }
    }

    ///Transformation from global to local
    void local(std::vector< Dune::FieldVector<ctype,dim> >& x,
               std::vector<int>& componentsInside) const
    {
      if (dir==FatherIsGlobal) {
        if(transformViaGlobal) fatherToChildViaGlobal(x,componentsInside);
        else                   fatherToChild(x,componentsInside);
      }
      else {
        if(transformViaGlobal) childToFatherViaGlobal(x);
        else                   childToFather(x);
      }
    }

    /**
     * @brief inverse of transposed Jacobian.
     *
     * The Jacobian is just the derivative of the local to global
     * map \f$ g \f$. This methhod returns the transpose of its inverse.
     */
    Dune::FieldMatrix<ctype,dim,dim> jacobianInverseTransposed (Dune::FieldVector<ctype,dim> const& xLocal) const
    {
      if (transformViaGlobal)
        if (dir==FatherIsGlobal) {
          Dune::FieldVector<ctype,dim> xGlobal = ancestor->geometry().local(child->geometry().global(xLocal));
          return ancestor->geometry().jacobianTransposed(xGlobal) * child->geometry().jacobianInverseTransposed(xLocal);
        } else {
          Dune::FieldVector<ctype,dim> xGlobal = child->geometry().local(ancestor->geometry().global(xLocal));
          Dune::FieldMatrix<ctype,dim,dim> Jct = child->geometry().jacobianInverseTransposed(xGlobal);
          Jct.invert();
          return Jct * ancestor->geometry().jacobianInverseTransposed(xLocal);
        }
      else {
        Dune::FieldMatrix<ctype,dim,dim> invt = unitMatrix<ctype,dim>();
        Dune::FieldVector<ctype,dim> x = xLocal;
        CellPointer ci(child);
        if (dir==FatherIsGlobal) {
          while(idSet.id(*ci) != idSet.id(*ancestor)) {
            invt.leftmultiply(ci->geometryInFather().jacobianInverseTransposed(x));
            x = ci->geometryInFather().global(x);
          }
        } else {
          std::stack<CellPointer> stack;
          while(idSet.id(*ci) != idSet.id(*ancestor)) {
            stack.push(ci);
            ci = ci->father();
          }
          while (!stack.empty()) {
            x = stack.top()->geometryInFather().local(x);
            // TODO: throw exception if x not in *ci
            invt.rightmultiply(stack.top()->geometryInFather().jacobianInverseTransposed(x));
            stack.pop();
          }
          invt.invert();
        }
        return invt;
      }
    }

  private:

    /** 
     * \brief In-place transformation from a father entity to a child entity.
     */
    void childToFather(std::vector< Dune::FieldVector<ctype,dim> >& x) const
    {
      CellPointer ci(child);
      while(idSet.id(*ci) != idSet.id(*ancestor))
      {
        auto const& geo = ci->geometryInFather();
        for(auto& xi: x)
          xi = geo.global(xi);
        ci = ci->father();
      }
    }

    void childToFatherViaGlobal(std::vector<Dune::FieldVector<ctype,dim>>& x) const
    {
      auto const& geoChild = child->geometry();
      auto const& geoAnces = ancestor->geometry();
      for(auto& xi: x) 
        xi = geoAnces.local(geoChild.global(xi));
    }

    /** 
     * \brief In-place transformation from a father entity to a child entity. 
     * 
     * Only those points are returned that lie inside the child
     * entity. All others are deleted during the process. Note that "insideness" can be
     * controlled by setting the tol_ parameter on construction.
     */
    void fatherToChild(std::vector<Dune::FieldVector<ctype,dim> >& x,
        std::vector<int> &componentsInside) const
    {
      assert(x.size()==componentsInside.size());

      // We need to go from father cells up to child cells. However, we can only obtain the father
      // directly. Thus we have to traverse from child to ancestor and store the whole line
      // of fathers to be traced back in reverse direction.
      CellPointer ci(child);
      std::stack<CellPointer> ciStack;
      while(idSet.id(*ci) != idSet.id(*ancestor)) {
        ciStack.push(ci);
        ci = ci->father();
      }
      // Now ciStack.top()->father is ancestor.

      // Trace back the line from ancestor to child and transform the coordinates given in x
      // to the childrens' local coordinates on the way.
      while(!ciStack.empty()) {
        ci = ciStack.top();
        typename std::vector<Dune::FieldVector<ctype,dim> >::iterator vi(x.begin());
        std::vector<int>::iterator compi(componentsInside.begin());
        while (vi!=x.end()) {
          *vi = ci->geometryInFather().local(*vi);
          if (checkInside(ci->type(),*vi)<tol) { // (probably) inside
            ++vi;
            ++compi;
          } else { // (almost sure) outside
            vi = x.erase(vi);
            compi = componentsInside.erase(compi);
          }
        }
        ciStack.pop();
      }
    }

    void fatherToChildViaGlobal(std::vector<Dune::FieldVector<ctype,dim>>& x,
                                std::vector<int>& componentsInside) const
    {
      for(int i=0; i<x.size(); i++)
        x[i] = child->geometry().local(ancestor->geometry().global(x[i]));

      auto vi = x.begin();
      auto compi = componentsInside.begin(); 
      while (vi!=x.end()) {
        if (checkInside(child->type(),*vi)<tol) { // (probably) inside
          ++vi;
          ++compi;
        } else { // (almost sure) outside
          vi = x.erase(vi);
          compi = componentsInside.erase(compi);
        }
      }
    }

  private:
    ctype tol;
  };

  /**
   * \todo doc me
   */
  template<class Grid>
  class LocalGeometryInCoarseGridAncestor : public MultiLevelLocalGeometry<Grid>
  {
    typedef typename Grid::template Codim<0>::EntityPointer CellPointer;
    typedef typename Grid::ctype ctype;
    static int const dim = Grid::dimension;
  public:

    LocalGeometryInCoarseGridAncestor(Grid const &grid,
                                      CellPointer const& child,
                                      bool geometricallyNested=false) :
    MultiLevelLocalGeometry<Grid>(grid,child,coarseGridAncestor(child),MultiLevelLocalGeometry<Grid>::ChildIsGlobal,geometricallyNested)
    {
    }

    /// Transformation from local to global
    Dune::FieldVector<ctype,dim> global(Dune::FieldVector<ctype,dim>const& x) const
    {
      std::vector<Dune::FieldVector<ctype,dim> > xi(1,x);
      std::vector<int> isi(1,0);
      MultiLevelLocalGeometry<Grid>::global(xi,isi);
      return xi[0];
    }

    CellPointer const& getFather() const  { return this->ancestor; }

  private:

    static CellPointer coarseGridAncestor(CellPointer const& child)
    {
      CellPointer ancestor(child);

      while (ancestor->level()>0) 
        ancestor = ancestor->father();
      return ancestor;
    }
  };
} // end of namespace Kaskade

#endif
