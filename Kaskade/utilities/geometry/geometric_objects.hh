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
/*
 * geometric_objects.h
 *
 *  Created on: 17.12.2011
 *      Author: Lars Lubkoll
 */

#ifndef GEOMETRIC_OBJECTS_HH_
#define GEOMETRIC_OBJECTS_HH_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <dune/common/fvector.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include "utilities/geometry/geomtools.hh"
#include "tools/linalg/scalarproducts.hh"

namespace GeometricObject
{
  enum Direction{ X, Y, Z };

  /*******************************************************************************************************************************/

  template <class Scalar, int dim>
  class Point : public Dune::FieldVector<Scalar,dim>{
  public:
    typedef Dune::FieldVector<Scalar,dim> Base;
    Point() : Base(){}
    explicit Point(Dune::FieldVector<Scalar,dim> const& p) : Base(p) {}
  };

  /*******************************************************************************************************************************/

  template <class Scalar, int dim>
  struct Line{
    Line() : start(), end()
    {}

    Line(Point<Scalar,dim> const& s, Point<Scalar,dim> const& e) : start(s), end(e)
    {}

    Point<Scalar,dim> start, end;
  };

  /*******************************************************************************************************************************/

  // TODO implement 2D-Polygon
  template <class Scalar, int dim>
  struct Triangle{
    typedef Point<Scalar,dim> Vertex;
    Triangle()
    {
      init(Vertex(), Vertex(), Vertex());
    }

    Triangle(Vertex const& v1, Vertex const& v2, Vertex const& v3)
    {
      init(v1,v2,v3);
    }

    boost::array<Vertex,3> corners;

  private:
    void init(Vertex const& v1, Vertex const& v2, Vertex const& v3)
    {
      corners[0] = v1; corners[1] = v2; corners[2] = v3;
    }
  };

  /*******************************************************************************************************************************/

  template <class Scalar, int dim>
  struct Rectangle{
    Rectangle()
    {
      init(Point<Scalar,dim>(), Point<Scalar,dim>(), Point<Scalar,dim>(), Point<Scalar,dim>());
    }

    explicit Rectangle(boost::array<Point<Scalar,dim>, 4> const& e)
    {
      init(e[0], e[1], e[2], e[3]);
    }

    Rectangle(Point<Scalar,dim> const& p0, Point<Scalar,dim> const& p1, Point<Scalar,dim> const& p2, Point<Scalar,dim> const& p3)
    {
      init(p0, p1, p2, p3);
    }

    boost::array<Point<Scalar,dim>,4> corners;

  private:
    void init(Point<Scalar,dim> const& p0, Point<Scalar,dim> const& p1, Point<Scalar,dim> const& p2, Point<Scalar,dim> const& p3)
    {
      corners[0] = p0;  corners[1] = p1;  corners[2] = p2;  corners[3] = p3;
    }
  };

  /*******************************************************************************************************************************/

  namespace {
    template <int dim>
    struct BoundingBoxContainsImpl{
      template <class Array, class Coordinate>
      static bool apply(Array const& array, Coordinate const& x)
      {
        return (BoundingBoxContainsImpl<dim-1>::apply(array, x) && array[dim-1].first<=x[dim-1] && array[dim-1].second>=x[dim-1]);
      }
    };

    template <>
    struct BoundingBoxContainsImpl<1>{
      template <class Array, class Coordinate>
      static bool apply(Array const& array, Coordinate const& x)
      {
        return (array[0].first <= x[0] && array[0].second >= x[0]);
      }
    };

    template <class Scalar, int dim, int id>
    struct BoundingBoxGetCornersImpl
    {
      template <class Array>
      static void apply(std::vector<Point<Scalar,dim> > &corners, Array const& array)
      {
        size_t offset = (size_t)(corners.size()/pow(2,id+1));
        bool first = true;
        size_t changes = 0;
        for(size_t i=0; i<corners.size(); ++i)
        {
          if( (i-changes*offset)==offset )
          {
            ++changes;
            first = !first;
          }

          if(first) corners[i][id] = array[id].first;
          else corners[i][id] = array[id].second;
        }
        BoundingBoxGetCornersImpl<Scalar,dim,id+1>::apply(corners, array);
      }
    };

    template <class Scalar, int dim>
    struct BoundingBoxGetCornersImpl<Scalar,dim,dim>
    {
      template <class Array>
      static void apply(std::vector<Point<Scalar,dim> > &corners, Array const& array){}
    };


    template <class Scalar, int dim, int id>
    struct BoundingBoxGetEdgesImpl
    {
      template <class Array>
      static void apply(std::vector<Line<Scalar,dim> > &edges, Array const& array, int indexOffset)
      {
        int localNumberOfEdges = pow(2,dim-1);
        for(int edge=0; edge<localNumberOfEdges; ++edge)
        {
          Point<Scalar,dim> start, end;
          start[id] = array[id].first;
          end[id] = array[id].second;

          for(int i=0; i<dim; ++i)
          {
            if(i==id) continue;
            if(edge==0){
              insert(start, end, array[i].first, i);
              continue;
            }
            if(edge==1)
            {
              if(i==1) insert(start, end, array[i].second, i);
              if(i==2) insert(start, end, array[i].first, i);
              continue;
            }
            if(edge==2)
            {
              if(i==1) insert(start, end, array[i].first, i);
              if(i==2) insert(start, end, array[i].second, i);
              continue;
            }
            if(edge==3)
            {
              if(i==1) insert(start, end, array[i].second, i);
              if(i==2) insert(start, end, array[i].second, i);
              continue;
            }
          }
          edges[edge] = Line<Scalar,dim>(start,end);
        }
        BoundingBoxGetEdgesImpl<Scalar,dim,id+1>::apply(edges, array, localNumberOfEdges);
      } // end apply

    private:
      static void insert(Point<Scalar,dim>& start, Point<Scalar,dim>& end, Scalar value, int index)
      {
        start[index] = value, end[index] = value;
      }
    };

    template <class Scalar, int dim>
    struct BoundingBoxGetEdgesImpl<Scalar,dim,dim>
    {
      template <class Array>
      static void apply(std::vector<Line<Scalar,dim> >&, Array const&, int){};
    };

    template <int v1, int v2>
    struct MyLess{
      static bool const value=v1<v2;
    };

    template <int v1, int v2>
    struct Minus{
      static int const value = v1-v2;
    };

    template <int val>
    struct Decrement{
      static int const value = val-1;
    };

    template <class Scalar, int dim, int dimMinusCodim>
    struct GetCodimTypeImpl;

    template <class Scalar, int dim>
    struct GetCodimTypeImpl<Scalar,dim,0>{
      typedef typename boost::enable_if<MyLess<0,dim>, Point<Scalar,dim> >::type type;
    };

    template <class Scalar, int dim>
    struct GetCodimTypeImpl<Scalar,dim,1>{
      typedef typename boost::enable_if<MyLess<1,dim>, Line<Scalar,dim> >::type type;
    };

    template <class Scalar, int dim>
    struct GetCodimTypeImpl<Scalar,dim,2>{
      typedef typename boost::enable_if<MyLess<2,dim>, Rectangle<Scalar,dim> >::type type;
    };

    template <class Scalar, int dim, int codim>
    struct GetCodimType{
      typedef typename GetCodimTypeImpl<Scalar,dim,dim-codim>::type type;
    };

    template <class Scalar, int dim, int dimMinusCodim>
    struct CodimBaseImpl{};

    template <class Scalar, int dim>
    struct CodimBaseImpl<Scalar,dim,2>
    {
      template <class Array>
      std::vector<typename GetCodimType<Scalar,dim,dim-2>::type > getFaces(Array const& array) const
      {
        std::vector<Rectangle<Scalar,3> > faces(6);
        std::vector<Point<Scalar,3> > corners=CodimBaseImpl<Scalar,dim,0>::getCorners(array);
        faces[0] = Rectangle<Scalar,3>(corners[0], corners[2], corners[3], corners[1]);
        faces[1] = Rectangle<Scalar,3>(corners[4], corners[6], corners[7], corners[5]);
        faces[2] = Rectangle<Scalar,3>(corners[0], corners[4], corners[6], corners[2]);
        faces[3] = Rectangle<Scalar,3>(corners[1], corners[5], corners[7], corners[3]);
        faces[4] = Rectangle<Scalar,3>(corners[0], corners[4], corners[5], corners[1]);
        faces[5] = Rectangle<Scalar,3>(corners[2], corners[6], corners[7], corners[3]);
        return faces;
      }
    };

    template <class Scalar, int dim>
    struct CodimBaseImpl<Scalar,dim,1>
    {
      template <class Array>
      std::vector<Line<Scalar,dim> > getEdges(Array const& array) const
      {
        if(dim==0) return std::vector<Line<Scalar,dim> >();

        size_t numEdges = 1;
        for(int i=1; i<dim; ++i) numEdges = 2*numEdges + pow(2,i);

        std::vector<Line<Scalar,dim> > edges(numEdges);
        BoundingBoxGetEdgesImpl<Scalar,dim, 0>::apply(edges, array, 0);
        return edges;
      }
    };

    template <class Scalar, int dim>
    struct CodimBaseImpl<Scalar, dim, 0>
    {
      template <class Array>
      static std::vector<Point<Scalar,dim> > getCorners(Array const& array)
      {
        std::vector<Point<Scalar,dim> > corners(pow(2,dim));
        BoundingBoxGetCornersImpl<Scalar,dim,0>::apply(corners, array);
        return corners;
      }
    };

    template <class Scalar, int dim, int codim>
    struct CodimBase : public CodimBaseImpl<Scalar,dim,dim-codim>, public CodimBase<Scalar,dim,codim-1>{};

    template <class Scalar, int dim>
    struct CodimBase<Scalar,dim,0>{};
  }

  /// A bounding box
  template <class Scalar, int dimension>
  struct BoundingBox : public CodimBase<Scalar,dimension,dimension>{
    static int const dim=dimension;
    BoundingBox()
    {
      for(int i=0; i<dim; ++i) coord[i].first = coord[i].second = 0;
    }

    BoundingBox(BoundingBox const& boundingBox)
    {
      coord = boundingBox.coord;
    }

    explicit BoundingBox(typename Point<Scalar,dim>::Base const& x)
    {
      for(int i=0; i<dim; ++i) coord[i].first = coord[i].second = x[i];
    }

    template <class Coordinate>
    void update(Coordinate const& x)
    {
      for(int i=0; i<dim; ++i)
      {
        if(coord[i].first > x[i]) coord[i].first = x[i];
        if(coord[i].second < x[i]) coord[i].second = x[i];
      }
    }

    std::ostream& print(std::ostream &stream) const
    {
      stream << "BoundingBox: ";
      for(int i=0; i<dim-1; ++i)
        stream << coord[i].first << " - " << coord[i].second << " | ";
      stream << coord[dim-1].first << " - " << coord[dim-1].second << std::endl;
      return stream;
    }

    template <class Coordinate>
    bool contains(Coordinate const& x) const
    {
      return BoundingBoxContainsImpl<dim>::apply(coord, x);
    }

    friend std::ostream& operator<<(std::ostream& stream, BoundingBox<Scalar,dim> const& boundingBox)
    {
      return boundingBox.print(stream);
    }

    Dune::FieldVector<std::pair<Scalar,Scalar>,dim> coord;
  };

  template <class Scalar, int dim>
  struct FastBoundingBox : public BoundingBox<Scalar,dim>
  {
    typedef BoundingBox<Scalar,dim> Base;

    FastBoundingBox() : Base(){}

    FastBoundingBox(FastBoundingBox const& boundingBox) : Base(boundingBox),
        corners(boundingBox.corners), edges(boundingBox.edges), faces(boundingBox.faces)
    {}

    template <class Coordinate>
    explicit FastBoundingBox(Coordinate const& x) : Base(x){}

    void initCorners()
    {
      corners = getCorners(this->coord);
    }

    void initEdges()
    {
      edges = getEdges(this->coord);
    }

    void initFaces()
    {
      faces = getFaces(this->coord);
    }

    void initAll()
    {
      initCorners();    initEdges();    initFaces();
    }

    std::vector<Point<Scalar,dim> > const& getStoredCorners() const
                    {
      return corners;
                    }

    std::vector<Line<Scalar,dim> > const& getStoredEdges() const
                    {
      return edges;
                    }

    std::vector<Rectangle<Scalar,dim> > const& getStoredFaces() const
                    {
      return faces;
                    }

    std::vector<Point<Scalar,dim> > corners;
    std::vector<Line<Scalar,dim> > edges;
    std::vector<Rectangle<Scalar,dim> > faces;
  };


  template <class BB>
  struct BoundingBoxWrapper;

  template <class Scalar, int dim>
  struct BoundingBoxWrapper<BoundingBox<Scalar,dim> >
  {
    typedef BoundingBox<Scalar,dim> BB;

    static void init(BB& boundingBox){}
    static void initCorners(BB& boundingBox){}
    static void initEdges(BB& boundingBox){}
    static void initFaces(BB& boundingBox){}
  };

  template <class Scalar, int dim>
  struct BoundingBoxWrapper<FastBoundingBox<Scalar,dim> >
  {
    typedef FastBoundingBox<Scalar,dim> BB;

    static void init(BB& boundingBox){ boundingBox.initAll(); }
    static void initCorners(BB& boundingBox){ boundingBox.initCorners(); }
    static void initEdges(BB& boundingBox){ boundingBox.initEdges(); }
    static void initFaces(BB& boundingBox){ boundingBox.initFaces(); }
  };

  /*******************************************************************************************************************************/

  /// A ball 
  template <class Scalar, int dim>
  struct Ball{
    Ball(Dune::FieldVector<Scalar,dim> &c, Scalar r) : center(c), radius(r)
    {}

    template <class Position>
    bool contains(Position const& x) const
    {
      return (x-center).two_norm() <= radius;
    }

    template <class Position>
    Scalar distanceFromCenter(Position const& x) const
    {
      return (x-center).two_norm();
    }

    Dune::FieldVector<Scalar,dim> center;
    Scalar radius;
  };
  /*******************************************************************************************************************************/
  // function declarations, definitions follow after namespace ImplementationDetail
  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& first, Point<Scalar,dim> const& second);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Line<Scalar,dim> const& line);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Line<Scalar,dim> const& line, Point<Scalar,dim> const& point);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Point<Scalar,dim> const& point);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Ball<Scalar,dim> const& ball);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Line<Scalar,dim> const& line);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Line<Scalar,dim> const& line, Ball<Scalar,dim> const& ball);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Rectangle<Scalar,dim> const& rectangle);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Rectangle<Scalar,dim> const& rectangle, Ball<Scalar,dim> const& ball);

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, BoundingBox<Scalar,dim> const& boundingBox);

  namespace ImplementationDetail{

    template <class FirstObject, class SecondObject, class Scalar, int dim, class ScalarProduct>
    struct ProjectionImpl;

    template <class Scalar, int dim, class ScalarProduct>
    struct ProjectionImpl<Point<Scalar,dim>, Line<Scalar,dim>, Scalar, dim, ScalarProduct>
    {
      static Point<Scalar,dim> apply(Point<Scalar,dim> const& point, Line<Scalar,dim> const& line, Scalar& a, ScalarProduct const scalarProduct)
      {
        Point<Scalar,dim> const direction = GeomTools::normalize(line.end-line.start);
        return line.start + scalarProduct(point-line.start, direction)*direction;
      }
    };

    template <class Scalar, int dim, class ScalarProduct>
    struct ProjectionImpl<Point<Scalar,dim>, Rectangle<Scalar,dim>, Scalar, dim, ScalarProduct>
    {
      static Point<Scalar,dim> apply(Point<Scalar,dim> const& point, Rectangle<Scalar,dim> const& rectangle, Scalar& a, Scalar& b, ScalarProduct const scalarProduct)
      {
        Point<Scalar,dim> const dir1 = GeomTools::normalize(rectangle.corners[1]-rectangle.corners[0]),
            dir2 = GeomTools::normalize(rectangle.corners[3]-rectangle.corners[0]);
        Point<Scalar,dim> const vec = point-rectangle.corners[0];
        return rectangle.corners[0] + scalarProduct(vec, dir1)*dir1 + scalarProduct(vec, dir2)*dir2;
      }
    };

    template <class Metric, class Scalar, int dim, class FirstObject, class SecondObject>
    struct DistanceImpl;

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Point<Scalar,dim> >
    {
      static Scalar apply(Point<Scalar,dim> const& firstPoint, Point<Scalar,dim> const& secondPoint, Metric metric)
      {
        return metric(firstPoint-secondPoint);
      }
    };

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Line<Scalar,dim> >
    {
      static Scalar apply(Point<Scalar,dim> const& point, Line<Scalar,dim> const& line, Metric metric)
      {
        Scalar a;
        Point<Scalar,dim> projectedPoint = ProjectionImpl<Point<Scalar,dim>, Line<Scalar,dim>, Scalar,dim,typename Metric::ScalarProduct>::apply(point,line,a);
        Scalar length = distance<Metric>(line.start, line.end);
        if(0 < a &&  a < length) return distance<Metric>(point,projectedPoint,metric);
        if(a <= 0) return distance<Metric>(point, line.start,metric);
        return distance<Metric>(point, line.end,metric);
      }
    };

    // template <class Metric

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Rectangle<Scalar,dim> >
    {
      typedef Line<Scalar,dim> Edge;
      static Scalar apply(Point<Scalar,dim> const& point, Rectangle<Scalar,dim> const& rectangle, Metric metric)
      {
        Scalar a,b;
        Point<Scalar,dim> projectedPoint = ProjectionImpl<Point<Scalar,dim>, Rectangle<Scalar,dim>, Scalar,dim,typename Metric::ScalarProduct>::apply(point, rectangle, a, b);
        Scalar l0 = distance<Metric>(rectangle.corners[0], rectangle.corners[1], metric),
            l1 = distance<Metric>(rectangle.corners[0], rectangle.corners[3], metric);
        if(0 < a && a < l0 && 0 < b && b < l1) return distance<Metric>(point, projectedPoint, metric);

        if(a <= 0)
        {
          if(b <= 0) return distance<Metric>(point, rectangle.corners[0], metric);
          if(b >= l1) return distance<Metric>(point, rectangle.corners[3], metric);
          return distance<Metric>(point, Edge(rectangle.corners[0], rectangle.corners[3]), metric);
        }
        if(a >= l0)
        {
          if(b <= 0) return distance<Metric>(point, rectangle.corners[1], metric);
          if(b >= l1) return distance<Metric>(point, rectangle.corners[2], metric);
          return distance<Metric>(point, Edge(rectangle.corners[1], rectangle.corners[2]), metric);
        }
        if(b <= 0) return distance<Metric>(point, Edge(rectangle.corners[0], rectangle.corners[1]), metric);
        return distance<Metric>(point, Edge(rectangle.corners[2], rectangle.corners[3]), metric);
      }
    };

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Point<Scalar,dim> >
    {
      static Scalar apply(Ball<Scalar,dim> const& ball, Point<Scalar,dim> const& point, Metric metric)
      {
        Scalar distanceToCenter = distance<Metric>(point, ball.center, metric);
        return (distanceToCenter < ball.radius) ? 0 : distanceToCenter-ball.radius;
      }
    };

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Line<Scalar,dim> >
    {
      static Scalar apply(Ball<Scalar,dim> const& ball, Line<Scalar,dim> const& line, Metric metric)
      {

        int distanceToCenter = distance<Metric>(ball.center, line, metric);
        return (distanceToCenter < ball.radius) ? 0 : distanceToCenter - ball.radius;
      }
    };

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Rectangle<Scalar,dim> >
    {
      static Scalar apply(Ball<Scalar,dim> const& ball, Rectangle<Scalar,dim> const& rectangle, Metric metric)
      {
        Scalar distanceToCenter = distance<Metric>(ball.center, rectangle, metric);
        return (distanceToCenter < ball.radius) ? 0 : distanceToCenter - ball.radius;
      }
    };

    template <class Metric, class Scalar, int dim>
    struct DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, BoundingBox<Scalar,dim> >
    {
      static Scalar apply(Ball<Scalar,dim> const& ball, BoundingBox<Scalar,dim> const& boundingBox, Metric metric)
      {
        if(boundingBox.contains(ball.center)) return 0;
        Scalar result=9999;
        // check corners
        std::vector<Point<Scalar,dim> > const& corners = boundingBox.getCorners();
        for(size_t cornerId=0; cornerId<corners.size(); ++cornerId)
        {
          Scalar tmp = distance<Metric>(ball, corners[cornerId], metric);
          if(tmp > 0) result = std::min(tmp,result); else return 0;
        }
        // check edges
        std::vector<Line<Scalar,dim> > const edges = boundingBox.getEdges();
        for(size_t edgeId=0; edgeId<edges.size(); ++edgeId)
        {
          Scalar tmp = distance<Metric>(ball, edges[edgeId], metric);
          if(tmp > 0) result = std::min(tmp, result); else return 0;
        }
        // check faces
        std::vector<Rectangle<Scalar,dim> > const faces = boundingBox.getFaces();
        for(size_t faceId=0; faceId<faces.size(); ++faceId)
        {
          Scalar tmp = distance<Metric>(ball, faces[faceId], metric);
          if(tmp > 0) result = std::min(tmp, result); else return 0;
        }

        return result;
      } // end distance
    };

//    template <class Metric, class Scalar, int dim, class DuneSimplex>
//    struct DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, DuneSimplex>
//    {
//      static Scalar apply(Ball<Scalar,dim> const& ball, DuneSimplex const& duneEntity)
//      {
//        if(Dune::GenericReferenceElements<double,3>::simplex().checkInside(duneEntity.geometry().local(ball.center))) return 0;
//        Scalar result = 9999;
//        int const numberOfCorners = duneEntity.geometry().corners();
//        // check corners
//        for(int cornerId=0; cornerId<numberOfCorners; ++cornerId)
//        {
//          Scalar tmp = distance<Metric>(ball, duneEntity.geometry().corner(cornerId));
//          if(tmp > 0) result = std::min(tmp, result);
//          else return 0;
//        }
//        // check edges
//        for(int startCornerId=0; startCornerId<numberOfCorners; ++startCornerId)
//        {
//          for(int endCornerId=startCornerId+1; endCornerId<numberOfCorners; ++endCornerId)
//          {
//            Scalar tmp = distance<Metric>(ball, Line<Scalar,dim>(duneEntity.geometry().corner(startCornerId), duneEntity.geometry(),corner(endCornerId)));
//            if(tmp > 0) result = std::min(tmp, result);
//            else return 0;
//          }
//        }
//        // check faces
//        if(dim>2)
//        {
//
//        }
//
//        return result;
//      }
//    };

    /*----------------------------------------------------------------------------------*/

    template <class Metric, class Scalar, int dim, class FirstObject, class SecondObject>
    struct IntersectionCheckImpl;

    template <class Metric, class Scalar, int dim>
    struct IntersectionCheckImpl<Metric, Scalar, dim, Ball<Scalar,dim>, BoundingBox<Scalar,dim> >
    {
      static bool apply(Ball<Scalar,dim> const& ball, BoundingBox<Scalar,dim> const& boundingBox, Metric metric)
      {
        return DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, BoundingBox<Scalar,dim> >::apply(ball,boundingBox,metric) <= 0;
      }
    };

//    template <class Metric, class Scalar, int dim, class DuneSimplex>
//    struct IntersectionCheckImpl<Metric, Ball<Scalar,dim>, DuneSimplex, Scalar, dim>
//    {
//      static bool apply(Ball<Scalar,dim> const& ball, DuneSimplex const& duneEntity)
//      {
//        return DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, DuneSimplex >::apply(ball,duneEntity) <= 0;
//      }
//    };

  }

  /*******************************************************************************************************************************/

  // project first object on second
  template <class Scalar, int dim, class FirstObject, class SecondObject, class ScalarProduct>
  FirstObject projectFirstOnSecond(FirstObject const& first, SecondObject const& second, ScalarProduct const& scalarProduct = ScalarProduct())
  {
    return ImplementationDetail::ProjectionImpl<FirstObject, SecondObject, Scalar,dim,ScalarProduct>::apply(first,second,scalarProduct);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& first, Point<Scalar,dim> const& second, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Point<Scalar,dim> >::apply(first,second,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Line<Scalar,dim> const& line, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Line<Scalar,dim> >::apply(point,line,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Line<Scalar,dim> const& line, Point<Scalar,dim> const& point, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Line<Scalar,dim> >::apply(point,line,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Triangle<Scalar,dim> const& triangle, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Triangle<Scalar,dim> >::apply(point,triangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Triangle<Scalar,dim> const& triangle, Point<Scalar,dim> const& point, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Triangle<Scalar,dim> >::apply(point,triangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Rectangle<Scalar,dim> const& rectangle, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Rectangle<Scalar,dim> >::apply(point,rectangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Rectangle<Scalar,dim> const& rectangle, Point<Scalar,dim> const& point, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Point<Scalar,dim>, Rectangle<Scalar,dim> >::apply(point,rectangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Point<Scalar,dim> const& point, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Point<Scalar,dim> >::apply(ball,point,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Point<Scalar,dim> const& point, Ball<Scalar,dim> const& ball, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Point<Scalar,dim> >::apply(ball,point,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Line<Scalar,dim> const& line, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Line<Scalar,dim> >::apply(ball,line,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Line<Scalar,dim> const& line, Ball<Scalar,dim> const& ball, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Line<Scalar,dim> >::apply(ball,line,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Triangle<Scalar,dim> const& triangle, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Triangle<Scalar,dim> >::apply(ball,triangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Triangle<Scalar,dim> const& triangle, Ball<Scalar,dim> const& ball, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Triangle<Scalar,dim> >::apply(ball,triangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, Rectangle<Scalar,dim> const& rectangle, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Rectangle<Scalar,dim> >::apply(ball,rectangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Rectangle<Scalar,dim> const& rectangle, Ball<Scalar,dim> const& ball, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, Rectangle<Scalar,dim> >::apply(ball,rectangle,metric);
  }

  template <class Scalar, int dim, class Metric>
  Scalar distance(Ball<Scalar,dim> const& ball, BoundingBox<Scalar,dim> const& boundingBox, Metric const& metric = Metric())
  {
    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, BoundingBox<Scalar,dim> >::apply(ball,boundingBox,metric);
  }

//  template <class Metric, class Scalar, int dim, class DuneSimplex>
//  Scalar distance(Ball<Scalar,dim> const& ball, DuneSimplex const& duneEntity)
//  {
//    typedef typename DuneSimplex::Geometry Geom;
//    return ImplementationDetail::DistanceImpl<Metric, Scalar, dim, Ball<Scalar,dim>, DuneSimplex >::apply(ball,duneEntity);
//  }

  template <class Scalar, int dim, class Metric>
  Scalar intersects(Ball<Scalar,dim> const& ball, BoundingBox<Scalar,dim> const& boundingBox, Metric const& metric = Metric())
  {
    return ImplementationDetail::IntersectionCheckImpl<Metric, Scalar, dim, Ball<Scalar,dim>, BoundingBox<Scalar,dim> >::apply(ball,boundingBox,metric);
  }

//  template <class Metric, class Scalar, int dim, class DuneSimplex>
//  Scalar intersects(Ball<Scalar,dim> const& ball, DuneSimplex const& boundingBox)
//  {
//    typedef typename DuneSimplex::Geometry Geom;
//    return ImplementationDetail::IntersectionCheckImpl<Metric, Scalar, dim, Ball<Scalar,dim>, DuneSimplex>::apply(ball,duneEntity);
//  }
}

#endif /* GEOMETRIC_OBJECTS_HH_ */
