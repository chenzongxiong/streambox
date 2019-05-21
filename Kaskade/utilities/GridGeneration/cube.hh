#ifndef KASKADE_CUBE_HH
#define KASKADE_CUBE_HH

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include "utilities/vectorTools.hh"

namespace Kaskade
{
  template <class Real, int dim_>
  class BasicGridElement
  {
  public:
    static constexpr int dim = dim_;

    BasicGridElement() : vertices(), simplices(), pi(3.141592653589793){}

    BasicGridElement(BasicGridElement const& other) = default;

    BasicGridElement(BasicGridElement&& other)
    {
      *this = std::move(other);
    }

    BasicGridElement& operator=(BasicGridElement const& other) = default;

    BasicGridElement& operator=(BasicGridElement&& other)
    {
      vertices = std::move(other.vertices);
      simplices = std::move(other.simplices);
      pi = other.pi;
      return *this;
    }

    std::vector<Dune::FieldVector<Real,dim> > const& getVertices() const { return vertices; }
    std::vector<std::vector<unsigned int> > const& getSimplices() const { return simplices; }

  protected:
    std::vector<Dune::FieldVector<Real,dim> > vertices;
    std::vector<std::vector<unsigned int> > simplices;
    Real pi;
  };


  template <class Real_=double>
  class Square : public BasicGridElement<Real_,2>
  {
  public:
    typedef Real_ Real;

    explicit Square(bool symmetric=false) : Square(Dune::FieldVector<Real,2>(0.0), Dune::FieldVector<Real,2>(1.0),symmetric)
    {}

    Square(Dune::FieldVector<Real,2> const& x0, Dune::FieldVector<Real,2> const& dx, bool symmetric=false)
    {
      if(!symmetric)
      {
        initVertices(x0[0],x0[1],dx[0],dx[1]);
        initTriangles();
      }
      else
      {
        initSymmetricVertices(x0[0],x0[1],dx[0],dx[1]);
        initSymmetricTriangles();
      }
    }

  private:
    using BasicGridElement<Real,2>::vertices;
    using BasicGridElement<Real,2>::simplices;

    void initVertices(Real x, Real y, Real dx, Real dy)
    {
      Dune::FieldVector<Real,2> v;
      v[0] = x;       v[1] = y;       vertices.push_back(v);
      v[0] = x+dx;    v[1] = y;       vertices.push_back(v);
      v[0] = x;       v[1] = y+dy;    vertices.push_back(v);
      v[0] = x+dx;    v[1] = y+dy;    vertices.push_back(v);
    }

    void initSymmetricVertices(Real x, Real y, Real dx, Real dy)
    {
      Dune::FieldVector<Real,2> v;
      v[0] = x;         v[1] = y;         vertices.push_back(v);
      v[0] = x+dx;      v[1] = y;         vertices.push_back(v);
      v[0] = x+dx;      v[1] = y+dy;      vertices.push_back(v);
      v[0] = x;         v[1] = y+dy;      vertices.push_back(v);
      v[0] = x+0.5*dx;  v[1] = y+0.5*dx;  vertices.push_back(v);
    }

    void initTriangles()
    {
      std::vector<unsigned int> p(3,0);

      p[0] = 0;     p[1] = 1;     p[2] = 3;   simplices.push_back(p);
      p[0] = 0;     p[1] = 3;     p[2] = 2;   simplices.push_back(p);
    }

    void initSymmetricTriangles()
    {
      std::vector<unsigned int> p(3,0);

      p[0] = 0;     p[1] = 1;     p[2] = 4;   simplices.push_back(p);
      p[0] = 1;     p[1] = 2;     p[2] = 4;   simplices.push_back(p);
      p[0] = 2;     p[1] = 3;     p[2] = 4;   simplices.push_back(p);
      p[0] = 3;     p[1] = 0;     p[2] = 4;   simplices.push_back(p);
    }
  };

  template <class Real_=double>
  class Cube : public BasicGridElement<Real_,3>
  {
  public:
    typedef Real_ Real;

    // creates a standard cube
    explicit Cube(bool symmetric=false) : Cube(Dune::FieldVector<Real,3>(0.0), Dune::FieldVector<Real,3>(1.0), symmetric)
    {}

    // creates a standard cube with the left, down, front corner in the point(x,y,z)
    Cube(Dune::FieldVector<Real,3> const& x0, Dune::FieldVector<Real,3> const& dx, bool symmetric)
    {
      if(!symmetric)
      {
        initVertices(x0[0],x0[1],x0[2],dx[0],dx[1],dx[2]);
        initTetrahedra();
      }
      else
      {
        initSymmetricVertices(x0[0],x0[1],x0[2],dx[0],dx[1],dx[2]);
        initSymmetricTetrahedra();
      }
    }

  private:
    using BasicGridElement<Real_,3>::vertices;
    using BasicGridElement<Real_,3>::simplices;

    void initVertices(Real x, Real y, Real z, Real dx, Real dy, Real dz)
    {
      Dune::FieldVector<Real,3> v;
      //points in the frontside
      v[0]=x;     v[1]=y;     v[2]=z;      vertices.push_back(v);
      v[0]=x+dx;  v[1]=y;     v[2]=z;      vertices.push_back(v);
      v[0]=x+dx;  v[1]=y;     v[2]=z+dz;   vertices.push_back(v);
      v[0]=x;     v[1]=y;     v[2]=z+dz;   vertices.push_back(v);

      //points in the backside
      v[0]=x;     v[1]=y+dy;  v[2]=z;      vertices.push_back(v);
      v[0]=x+dx;  v[1]=y+dy;  v[2]=z;      vertices.push_back(v);
      v[0]=x+dx;  v[1]=y+dy;  v[2]=z+dz;   vertices.push_back(v);
      v[0]=x;     v[1]=y+dy;  v[2]=z+dz;   vertices.push_back(v);

      //point in the center
      v[0]=x+0.5*dx; v[1]=y+0.5*dy ;v[2]=z+0.5*dz; vertices.push_back(v);
    }

    void initSymmetricVertices(Real x, Real y, Real z, Real dx, Real dy, Real dz)
    {
      Dune::FieldVector<Real,3> v;

      // points on the frontside
      v[0]=x;     v[1]=y;     v[2]=z;      vertices.push_back(v);           // 0
      v[0]=x+dx;  v[1]=y;     v[2]=z;      vertices.push_back(v);           // 1
      v[0]=x+dx;  v[1]=y;     v[2]=z+dz;   vertices.push_back(v);           // 2
      v[0]=x;     v[1]=y;     v[2]=z+dz;   vertices.push_back(v);           // 3

      // points on the backside
      v[0]=x;     v[1]=y+dy;  v[2]=z;      vertices.push_back(v);           // 4
      v[0]=x+dx;  v[1]=y+dy;  v[2]=z;      vertices.push_back(v);           // 5
      v[0]=x+dx;  v[1]=y+dy;  v[2]=z+dz;   vertices.push_back(v);           // 6
      v[0]=x;     v[1]=y+dy;  v[2]=z+dz;   vertices.push_back(v);           // 7

      // point in the center
      v[0]=x+0.5*dx; v[1]=y+0.5*dy; v[2]=z+0.5*dz; vertices.push_back(v);   // 8

      // points in the center of face
      v[0]=x;        v[1]=y+0.5*dy; v[2]=z+0.5*dz; vertices.push_back(v);   // 9
      v[0]=x+dx;     v[1]=y+0.5*dy; v[2]=z+0.5*dz; vertices.push_back(v);   // 10
      v[0]=x+0.5*dx; v[1]=y;        v[2]=z+0.5*dz; vertices.push_back(v);   // 11
      v[0]=x+0.5*dx; v[1]=y+dy;     v[2]=z+0.5*dz; vertices.push_back(v);   // 12
      v[0]=x+0.5*dx; v[1]=y+0.5*dy; v[2]=z;        vertices.push_back(v);   // 13
      v[0]=x+0.5*dx; v[1]=y+0.5*dy; v[2]=z+dz;     vertices.push_back(v);   // 14
    }

    void initSymmetricTetrahedra()
    {
      std::vector<unsigned int> p(4,0);

      // bottom
      p[0]=13;    p[1]=1;   p[2]=5;   p[3]=8; simplices.push_back(p);
      p[0]=13;    p[1]=5;   p[2]=4;   p[3]=8; simplices.push_back(p);
      p[0]=13;    p[1]=4;   p[2]=0;   p[3]=8; simplices.push_back(p);
      p[0]=13;    p[1]=0;   p[2]=1;   p[3]=8; simplices.push_back(p);

      // top
      p[0]=14;    p[1]=3;   p[2]=2;   p[3]=8; simplices.push_back(p);
      p[0]=14;    p[1]=2;   p[2]=6;   p[3]=8; simplices.push_back(p);
      p[0]=14;    p[1]=6;   p[2]=7;   p[3]=8; simplices.push_back(p);
      p[0]=14;    p[1]=7;   p[2]=3;   p[3]=8; simplices.push_back(p);

      // right
      p[0]=10;    p[1]=6;   p[2]=5;   p[3]=8; simplices.push_back(p);
      p[0]=10;    p[1]=5;   p[2]=1;   p[3]=8; simplices.push_back(p);
      p[0]=10;    p[1]=1;   p[2]=2;   p[3]=8; simplices.push_back(p);
      p[0]=10;    p[1]=2;   p[2]=6;   p[3]=8; simplices.push_back(p);

      // left
      p[0]=9;     p[1]=3;   p[2]=0;   p[3]=8; simplices.push_back(p);
      p[0]=9;     p[1]=0;   p[2]=4;   p[3]=8; simplices.push_back(p);
      p[0]=9;     p[1]=4;   p[2]=7;   p[3]=8; simplices.push_back(p);
      p[0]=9;     p[1]=7;   p[2]=3;   p[3]=8; simplices.push_back(p);

      // front
      p[0]=11;    p[1]=1;   p[2]=0;   p[3]=8; simplices.push_back(p);
      p[0]=11;    p[1]=0;   p[2]=3;   p[3]=8; simplices.push_back(p);
      p[0]=11;    p[1]=3;   p[2]=2;   p[3]=8; simplices.push_back(p);
      p[0]=11;    p[1]=2;   p[2]=1;   p[3]=8; simplices.push_back(p);

      // back
      p[0]=12;    p[1]=4;   p[2]=5;   p[3]=8; simplices.push_back(p);
      p[0]=12;    p[1]=5;   p[2]=6;   p[3]=8; simplices.push_back(p);
      p[0]=12;    p[1]=6;   p[2]=7;   p[3]=8; simplices.push_back(p);
      p[0]=12;    p[1]=7;   p[2]=4;   p[3]=8; simplices.push_back(p);
    }

    void initTetrahedra()
    {
      std::vector<unsigned int> p(4,0);

      // front side
      p[0]=0; p[1]=1; p[2]=2; p[3]=8; simplices.push_back(p);
      p[0]=2; p[1]=3; p[2]=0; p[3]=8; simplices.push_back(p);

      // back side
      p[0]=6; p[1]=5; p[2]=4; p[3]=8; simplices.push_back(p);
      p[0]=4; p[1]=6; p[2]=7; p[3]=8; simplices.push_back(p);

      // right side
      p[0]=6; p[1]=1; p[2]=5; p[3]=8; simplices.push_back(p);
      p[0]=2; p[1]=1; p[2]=6; p[3]=8; simplices.push_back(p);

      // left side
      p[0]=0; p[1]=4; p[2]=7; p[3]=8; simplices.push_back(p);
      p[0]=0; p[1]=7; p[2]=3; p[3]=8; simplices.push_back(p);

      // bottom
      p[0]=0; p[1]=1; p[2]=5; p[3]=8; simplices.push_back(p);
      p[0]=0; p[1]=5; p[2]=4; p[3]=8; simplices.push_back(p);

      // top
      p[0]=3; p[1]=2; p[2]=6; p[3]=8; simplices.push_back(p);
      p[0]=3; p[1]=6; p[2]=7; p[3]=8; simplices.push_back(p);
    }
  };

  template <class Real_=double>
  class Circle : public BasicGridElement<Real_,2>
  {
  public:
    typedef Real_ Real;

    // creates a circle
    explicit Circle(bool symmetric=false) : Circle(Dune::FieldVector<Real,2>(0.0), 1., 8, symmetric)
    {}

    // creates a standard cube with the left, down, front corner in the point(x,y,z)
    Circle(Dune::FieldVector<Real,2> const& center, Real radius, size_t nCorners, bool symmetric)
    {
      assert(nCorners>0);
      Real localAngle = (2*pi)/nCorners;
      Real ratio = 1./localAngle;

      if(!symmetric)
      {
        initVertices(center, radius, localAngle, nCorners, symmetric);
        initSimplices(nCorners, localAngle);
      }
      else
      {
        initSymmetricVertices(center, radius, localAngle, nCorners, symmetric);
        initSymmetricTetrahedra();
      }
    }

  private:
    using BasicGridElement<Real,2>::vertices;
    using BasicGridElement<Real,2>::simplices;
    using BasicGridElement<Real,2>::pi;

    Real log2(Real val)
    {
      return log(val)/log(2);
    }

    void initVertices(Dune::FieldVector<Real,2> const& center, Real radius, Real localAngle, size_t nCorners, bool symmetric)
    {
      size_t regFactor = (nCorners > 8) ? (size_t)round(log2(1/localAngle)) : 1;
      Real div = 1./regFactor;
      if(div < 1 ) div=0.5;

      vertices.push_back(center);

      Dune::FieldVector<Real,2> v;
      size_t innerCorners = 8;
      Real angle = pi*0.25;
      for(size_t i=0; i<innerCorners; ++i)
      {
        v[0] = div * radius * sin(i*angle);
        v[0] = div * radius * cos(i*angle);
        vertices.push_back(v);
      }

//      if(nCorners <= innerCorners) return;
//
//      innerCorners *= 2;
//      angle *= 0.5;
//      div = 1;
//      for(size_t i=0; i<innerCorners; ++i)
//      {
//        v[0] = div * radius * sin(i*angle);
//        v[0] = div * radius * cos(i*angle);
//        vertices.push_back(v);
//      }

//      for(size_t i=0; i<regFactor; ++i)
//      {
//
//      }
//
//      if(regFactor > 1)
//      {
//        for(size_t i=0; i<nCorners; ++i)
//        {
//          v[0] = 0.5*radius * sin(i*localAngle);
//          v[1] = 0.5*radius * cos(i*localAngle);
//          vertices.push_back(v);
//        }
//      }
//
//      for(size_t i=0; i<nCorners; ++i)
//      {
//        v[0] = radius * sin(i*localAngle);
//        v[1] = radius * cos(i*localAngle);
//        vertices.push_back(v);
//      }


    }

    void initSimplices(size_t nCorners, Real localAngle)
    {
      size_t regFactor = (nCorners > 8) ? (size_t)round(log2(1/localAngle)) : 1;
      Real div = 1./regFactor;
      if(div < 1 ) div=0.5;

      // inner ring
      std::vector<unsigned int> p(3,0);
      for(size_t i=0; i<8; ++i)
      {
        p[1] = i+1;
        if(i+2 > 8) p[2] = i+2-8;
        else p[2] = i+2;
        simplices.push_back(p);
      }

//      // outer ring
//      size_t lastRingNoNodes = 8;
//      unsigned int a, b, c;
//      for(size_t i=1; i<=lastRingNoNodes; ++i)
//      {
//        a = lastRingNoNodes + 2*i - 1;
//        b = lastRingNoNodes + 2*i;
//        c = i;
//        addSimplex(a,b,c);
//
//        a = i+1;
//        addSimplex(a,b,c);
//
//        c = (b+1)==vertices.size() ? lastRingNoNodes+1 : b+1;
//        addSimplex(a,b,c);
//      }
    }

    void addSimplex(unsigned int a, unsigned int b, unsigned int c)
    {
      std::vector<unsigned int> p = { a, b, c };
      simplices.push_back(p);
    }

    void initSymmetricVertices(Dune::FieldVector<Real,2> const& center, Real radius, Real localAngle, size_t nCorners, bool symmetric)
    {
      vertices.push_back(center);

      Dune::FieldVector<Real,2> v;
      for(size_t i=0; i<nCorners; ++i)
      {
        v[0] = radius * sin(i*localAngle);
        v[1] = radius * cos(i*localAngle);
        vertices.push_back(v);
      }
    }

    void initSymmetricTetrahedra()
    {

    }
  };

} // end namespace Kaskade

#endif
