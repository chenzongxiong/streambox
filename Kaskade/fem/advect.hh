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

#ifndef ADVECT_HH
#define ADVECT_HH

#include <cassert>
#include <utility>

#include "fem/fixdune.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  
  /// \cond internals
  namespace AdvectDetail
  {
    /**
     * \ingroup fetransfer
     * \brief Finds the face that is intersected by a ray starting from an interior point of a reference cell.
     * 
     * \param gt the type of cell
     * \param xi the point in the reference cell
     * \param dxi the direction in reference coordinates (nonzero)
     * 
     * \return a pair (face subentity index,step size)
     * 
     * If a zero direction dxi is provided, there is no face that could be encountered. In that case, 
     * a DetailedException is thrown.
     */
    template <class LocalCoordinate>
    std::pair<int,double> intersectedFace(Dune::GeometryType const& gt, LocalCoordinate const& xi, LocalCoordinate const& dxi)
    {
      using Scalar = typename LocalCoordinate::field_type;
      
      if (gt.isSimplex())
      {
        // Check minimal step size until faces are encountered
        int imin = -1;
        double minStep = std::numeric_limits<Scalar>::infinity();
        
        for (int i=0; i<xi.N(); ++i)       // first step through coordinate planes
          if (dxi[i]<0)                    // only encountered if direction has negative value in that coordinate
          {
            Scalar stepi = -xi[i]/dxi[i];  // compute the step size until intersection
            if (stepi < minStep)           // select the one that is encountered first
            {
              imin = i;
              minStep = stepi;
            }
          }
          
        LocalCoordinate one(1.0);                       // test the face normal to [1,1,1]
        Scalar dxione = one*dxi;                        // normal component of direction
        if (dxione > 0 && (1-one*xi)/dxione < minStep)  // if positive, compute intersection step size
        {                                               // and check whether this is encountered first
          minStep = (1-one*xi)/dxione;
          imin = xi.N();
        }
        
        // Make sure we don't return nonsense (like negative entity indices) just becaus we couldn't find
        // an intersection.
        if (imin < 0)
          throw DetailedException("Vanishing direction provided for face-ray intersection.",__FILE__,__LINE__);
        
        // Find the face subentity index in the reference simplex of Dune
        if (LocalCoordinate::dimension==1)              
          return std::make_pair(imin,minStep);
        if (LocalCoordinate::dimension==2)
        {
          constexpr int sube[] = { 1, 0, 2 };
          return std::make_pair(sube[imin],minStep);
        }
        if (LocalCoordinate::dimension==3)
        {
          constexpr int sube[] = { 2, 1, 0, 3 };
          return std::make_pair(sube[imin],minStep);
        }
        // never get here
      }
      else if (gt.isCube())
      {
        // not yet implemented
        abort();
      }
      
      // Never get here
      abort();
      return std::make_pair(-1,-1.0);
    }
    
  }
  /// \endcond
  
  /**
   * \ingroup fetransfer
   * \brief A weak function view that defines an advected function.
   *
   * The function view's value is defined by advection along the velocity field \f$ v\f$ as
   * \f$ \hat f(x) = f(y(\tau)) \f$, where \f$ y(t) \f$ satisfies
   * \f[ \dot y = -v(y), \quad y(0) = x. \f]
   * 
   * The integration is performed with \f$ n \f$ equidistant explicit Euler steps. The integration is stopped 
   * prematurely if the domain boundary is encountered.
   * 
   * \tparam Function a weak function view (e.g. a finite element function) type
   * \tparam Velocity a vector-valued WeakFunctionView. The vector size has to be the spatial dimension of the grid.
   */
  template <class Function, class Velocity>
  class AdvectedFunctionView
  {
  public:
    using GridView = typename Function::Space::GridView;
    using Cell = typename GridView::template Codim<0>::Entity;
    using Position = Dune::FieldVector<typename GridView::ctype,GridView::dimension>;
    
  private:
    using GlobalPosition = std::pair<Cell,Position>;
    
  public:

    /**
     * \brief Constructor.
     * \param f the function to be advected
     * \param v the velocity field 
     * \param tau the integration length
     * \param n the number of Euler integration steps (>= 1)
     */
    AdvectedFunctionView(Function const& f_, Velocity const& v_, double tau_, int n_=1): f(f_), v(v_), tau(tau_), n(n_)
    {
      assert(n>=1);
    }
    
    typename Function::ValueType value(Cell const& cell, Position const& xi) const 
    {
      // perform n explicit Euler steps in backwards direction for finding the source point y(-tau)
      // TODO: consider using a higher order integrator. Is this more efficient? Or do 
      //       discontinuities (of higher derivatives) of v at faces prevent this to 
      //       achieve higher order? Is a restriction of higher order integration to 
      //       single cells more appropriate?
      GlobalPosition y(cell,xi);
      for (int i=0; i<n; ++i)
        y = advectPosition(y,-tau/n*v.value(y.first,y.second));

      // evaluate function
      return f.value(y.first,y.second);
    }
    
  private:
    Function const& f;
    Velocity const& v;    // velocity field
    double tau;           // integration length
    int n;                // number of Euler steps
    
    // A single Euler step. As the Euler step can leave the current cell, we have to 
    // search for the cell containing the target point. Currently we follow the 
    // Euler step line through the grid as it passes through the cells.
    // TODO: For large step sizes relative to cell size this is inefficient, as many cells 
    //       have to be crossed. Then, going down to coarser levels with larger cells should
    //       speed up the process, finally going up again with a hierarchic search. Where is
    //       the trade-off?
    GlobalPosition advectPosition(GlobalPosition const& y, Position const& dy) const
    {
      auto const& geo = y.first.geometry();
      
      auto xi = y.second;                                           // current local point
      auto dxi = transpose(geo.jacobianInverseTransposed(xi))*dy;   // step in reference coordinates
      
      // TODO: xinext = xi+dxi could be obtained by geo.local(geo.global(xi)+dy). But what happens
      // if this is outside the cell?

      if (checkInside(geo.type(),xi+dxi) <= 1e-8)             // if (almost) inside current cell...
        return GlobalPosition(y.first,xi+dxi);                // ... we've found the target point
        
      // Now comes the hard part. We've to find the cell containing the target point. First we 
      // identify the face and intersection point.
      auto encounter = AdvectDetail::intersectedFace(geo.type(),xi,dxi);
      auto xinext = xi+encounter.second*dxi;
      
      // Then we find the neighboring cell (if any) and proceed with the 
      for (auto const& face: intersections(f.space().gridView(),y.first))
        if (face.indexInInside()==encounter.first)           // found the face
        {
          if (face.neighbor())                               // if there is a neighboring cell, step over and go ahead
            return advectPosition(GlobalPosition(face.outside(),face.geometryInOutside().global(face.geometryInInside().local(xinext))),
                                  (1-encounter.second)*dy);  // (tail recursion should allow compiler optimization)
          else                                               // no neighboring cell - we blundered into the domain boundary
            return GlobalPosition(y.first,xinext);           // ... return the boundary point
        }
      // never get here
      abort(); return y;
    }
  };

  // -------------------------------------------------------------------------------------------
  
  /**
   * \ingroup fetransfer
   * \brief Creates a weak function view for advecting FE functions.
   * \relates AdvectedFunctionView
   */
  template <class Function, class Velocity>
  auto makeAdvectedFunctionView(Function const& f, Velocity const& v, double tau, int n=1)
  {
    return AdvectedFunctionView<Function,Velocity>(f,v,tau,n);
  }
  
}

#endif