// -*- tab-width: 4; indent-tabs-mode: nil -*-
// $Id: vtkwriter.hh 4915 2009-03-10 10:50:41Z mnolte $

#ifndef QUAD_GNUPLOTWRITER_HH
#define QUAD_GNUPLOTWRITER_HH

#include <memory> // std::unique_ptr
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <list>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#include "io/base64.hh"
#include "io/check_endianness.hh"

/** @file
    @author Peter Bastian, modified by Irina Gossmann (irina.gossmann@gmail.com)
    @brief Provides file i/o for the visualization toolkit
 */

/** 
    This is a modified version of Dune's vtkwriter.hh that writes out the grid cells
    using Gnuplot data. Also, support for binary output was added.

    Implementation details: added new iterators over corner and midpoint vertices of a reference 
    element. 

    Base64 support for binary output is provided by a separate class (in base64.hh)
    allowed for free re-distribution by its author. See base64.hh for copyright notice.
 */


/** \brief options for Gnuplot output
        \ingroup Gnuplot */


using namespace Dune;

struct GnuplotOptions
{
  enum OutputType {
    /** @brief Output to the file is in ascii. */
    ascii,
    /** @brief Output to the file is binary. */
    binary,
    /** @brief Ouput is appended to the binary file. */
    binaryappended
  };
  enum DataMode {
    /** @brief Output conforming data. */
    conforming,
    /** @brief Output non conforming data. */
    nonconforming
  };
  enum Order {
    /** @brief Linear output. */
    linear,
    /** @brief Quadratic output. */
    quadratic
  };
  enum Info { 
    /** @brief No info messages. */
    none, 
    /** @brief Statistical summary. */
    summary, 
    /** @brief Statistical summary and listing of all points which exceed the standard deviation. */
    detail
  };
};


// map type to name in data array
template<class T>
struct GnuplotTypeNameTraits {
  std::string operator () (){
    return "";
  }
};

template<>
struct GnuplotTypeNameTraits<char> {
  std::string operator () () {
    return "Int8";
  }
  typedef int PrintType;
};

template<>
struct GnuplotTypeNameTraits<unsigned char> {
  std::string operator () () {
    return "UInt8";
  }
  typedef int PrintType;
};

template<>
struct GnuplotTypeNameTraits<short> {
  std::string operator () () {
    return "Int16";
  }
  typedef short PrintType;
};

template<>
struct GnuplotTypeNameTraits<unsigned short> {
  std::string operator () () {
    return "UInt16";
  }
  typedef unsigned short PrintType;
};

template<>
struct GnuplotTypeNameTraits<int> {
  std::string operator () () {
    return "Int32";
  }
  typedef int PrintType;
};

template<>
struct GnuplotTypeNameTraits<unsigned int> {
  std::string operator () () {
    return "UInt32";
  }
  typedef unsigned int PrintType;
};

template<>
struct GnuplotTypeNameTraits<float> {
  std::string operator () () {
    return "Float32";
  }
  typedef float PrintType;
};

template<>
struct GnuplotTypeNameTraits<double> {
  std::string operator () () {
    return "Float64";
  }
  typedef double PrintType;
};


/**
 * @brief Writer for the output of grid functions in the vtk format.
 * @ingroup Gnuplot
 *
 * Writes arbitrary grid functions (living on cells or vertices of a grid)
 * to a file suitable for easy visualization with
 * <a href="http://public.kitware.com/Gnuplot/">The Visualization Toolkit (Gnuplot)</a>.
 */
template< class GridView >
class GnuplotWriter {
  template<int dim>
  struct P0Layout
  {
    bool contains (Dune::GeometryType gt)
    {
      return gt.dim()==dim;
    }
  };

  template<int dim>
  struct P1Layout
  {
    bool contains (Dune::GeometryType gt)
    {
      return gt.dim()==0;
    }
  };

  template<int dim>
  struct P3Layout
  {
    bool contains (Dune::GeometryType gt)
    {
      return (gt.dim()==1 || gt.dim()==0);
    }
  };

  // extract types
  typedef typename GridView::Grid Grid;
  typedef typename Grid::ctype DT;
  enum { n = GridView::dimension };
  enum { w = GridView::dimensionworld };

  typedef typename GridView::template Codim< 0 >::Entity Cell;
  typedef typename GridView::template Codim< n >::Entity Vertex;
  typedef Cell Entity;

  typedef typename GridView::IndexSet IndexSet;

  static const PartitionIteratorType Gnuplot_Partition = InteriorBorder_Partition;

  typedef typename GridView::template Codim< 0 >
  ::template Partition< Gnuplot_Partition >::Iterator
   GridCellIterator;
  typedef typename GridView::template Codim< n >
  ::template Partition< Gnuplot_Partition >::Iterator
   GridVertexIterator;

  typedef MultipleCodimMultipleGeomTypeMapper< GridView, P1Layout > VertexMapper;
  typedef MultipleCodimMultipleGeomTypeMapper< GridView, P3Layout > QuadVertexMapper;

public:
  /** \brief A base class for grid functions with any return type and dimension
          \ingroup Gnuplot

          Trick : use double as return type
   */
  class GnuplotFunction
  {
  public:
    //! return number of components
    virtual int ncomps () const = 0;

    //! evaluate single component comp in the entity e at local coordinates xi
    /*! Evaluate the function in an entity at local coordinates.
        @param[in]  comp   number of component to be evaluated
        @param[in]  e      reference to grid entity of codimension 0
        @param[in]  xi     point in local coordinates of the reference element of e
        \return            value of the component
     */
    virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const = 0;

    //! get name
    virtual std::string name () const = 0;

    //! virtual destructor
    virtual ~GnuplotFunction () {}
  };

protected:
  typedef typename std::list<GnuplotFunction*>::iterator FunctionIterator;

class CellIterator :
public ForwardIteratorFacade<CellIterator, Entity, Entity&, int>
{
  GridCellIterator git;
  GridCellIterator gend;
public:
  CellIterator(const GridCellIterator & x, const GridCellIterator & end) : git(x), gend(end) {};
  void increment ()
  {
    ++git;
    while (git!=gend && git->partitionType()!=InteriorEntity) ++git;
  }
  bool equals (const CellIterator & cit) const
  {
    return git == cit.git;
  }
  Entity& dereference() const
  {
    return (Entity&) *git;
  }
  const FieldVector<DT,n> position() const
            {
    return ReferenceElements<DT,n>::general(git->type()).position(0,0);
}
};

CellIterator cellBegin() const
{
  return CellIterator( gridView_.template begin< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >() );
}

CellIterator cellEnd() const
{
  return CellIterator( gridView_.template end< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >() );
}

class VertexIterator :
public ForwardIteratorFacade<VertexIterator, Entity, Entity&, int>
{
  GridCellIterator git;
  GridCellIterator gend;
  GnuplotOptions::DataMode datamode;
  int index;
  const VertexMapper & vertexmapper;
  std::vector<bool> visited;
  const std::vector<int> & number;
  int offset;
protected:
  void basicIncrement ()
  {
    if( git == gend )
      return;
    ++index;
    const int numCorners = git->template count< n >();
    if( index == numCorners )
    {
      offset += numCorners;
      index = 0;

      ++git;
      while( (git != gend) && (git->partitionType() != InteriorEntity) )
        ++git;
    }
  }
public:
  VertexIterator(const GridCellIterator & x,
      const GridCellIterator & end,
      const GnuplotOptions::DataMode & dm,
      const VertexMapper & vm,
      const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), index(0),
        vertexmapper(vm), visited(vm.size(), false),
        number(num), offset(0)
  {
    if (datamode == GnuplotOptions::conforming && git != gend)
      visited[vertexmapper.map(*git,index,n)] = true;
  };
  void increment ()
  {
    switch (datamode)
    {
    case GnuplotOptions::conforming:
      while(visited[vertexmapper.map(*git,index,n)])
      {
        basicIncrement();
        if (git == gend) return;
      }
      visited[vertexmapper.map(*git,index,n)] = true;
      break;
    case GnuplotOptions::nonconforming:
      basicIncrement();
      break;
    }
  }
  bool equals (const VertexIterator & cit) const
  {
    return git == cit.git
        && index == cit.index && datamode == cit.datamode;
  }
  Entity& dereference() const
  {
    return (Entity&) *git;
  }
  int id () const
  {
    switch (datamode)
    {
    case GnuplotOptions::conforming:
      return
      number[vertexmapper.map(*git,renumber(*git,index),n)];
    case GnuplotOptions::nonconforming:
      return offset + renumber(*git,index);
    default:
      DUNE_THROW(IOError,"GnuplotWriter: unsupported DataMode" << datamode);
    }
  }
  int localindex () const
  {
    return index;
  }
  const FieldVector<DT,n> & position () const
            {
    return ReferenceElements<DT,n>::general(git->type()).position(index,n);
            }
};

VertexIterator vertexBegin () const
{
  return VertexIterator( gridView_.template begin< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >(),
      datamode, *vertexmapper, number );
}

VertexIterator vertexEnd () const
{
  return VertexIterator( gridView_.template end< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >(),
      datamode, *vertexmapper, number );
}


//same as VertexIterator, but iterates over midpoints vertices in addition to corner vertices
//iteration over midpoints is implemented as iteration over edges of the element
class QuadVertexIterator :
public ForwardIteratorFacade<QuadVertexIterator, Entity, Entity&, int>
{
  GridCellIterator git;
  GridCellIterator gend;
  GnuplotOptions::DataMode datamode;
  int index;
  const QuadVertexMapper & quadvertexmapper;
  std::vector<bool> visited;
  const std::vector<int> & number;
  int offset;
  int codim; //codim of the current element in quadvertexmapper (\n-1\ for midpoint vertex or \n\ for corner vertex)
  int subindex; //index of the current element in quadvertexmapper with current codim
protected:
  void basicIncrement ()
  {
    if( git == gend )
      return;
    ++index;
    const int numCorners = git->template count<n>();
    const int numEdges = git-> template count<n-1>();
    const int numVertices = numCorners + numEdges;
    if( index == numVertices )
    {
      offset += numVertices;
      index = 0;

      ++git;
      while( (git != gend) && (git->partitionType() != InteriorEntity) )
        ++git;
    }
    if (index < numCorners) {
      codim = n;
      subindex = index;
    }
    else {
      codim = n-1;
      subindex = index-numCorners;
    }
  }

  int qvmIndex() const{
    if (codim==n)
      return quadvertexmapper.map(*git,subindex,n);
    else //codim == n-1
        return quadvertexmapper.map(*git,subindex,n-1);
  }

  int renumberedqvmIndex() const{
    if (codim==n)
      return quadvertexmapper.map(*git,renumber(*git,subindex,codim),n);
    else
      return quadvertexmapper.map(*git,renumber(*git,subindex,codim),n-1);
  }


public:
  QuadVertexIterator(const GridCellIterator & x,
      const GridCellIterator & end,
      const GnuplotOptions::DataMode & dm,
      const QuadVertexMapper & qvm,
      const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), index(0), 
        quadvertexmapper(qvm), visited(qvm.size(), false),
        number(num), offset(0), codim(n), subindex(0)
  {
    if (datamode == GnuplotOptions::conforming && git != gend)
      visited[quadvertexmapper.map(*git,index,n)] = true;
  };
  void increment ()
  {
    switch (datamode)
    {
    case GnuplotOptions::conforming:
      while(visited[qvmIndex()])
      {
        basicIncrement();
        if (git == gend) return;
      }
      visited[qvmIndex()] = true;
      break;
    case GnuplotOptions::nonconforming:
      basicIncrement();
      break;
    }
  }
  bool equals (const QuadVertexIterator & cit) const
  {
    return git == cit.git
        && index == cit.index && datamode == cit.datamode;
  }
  Entity& dereference() const
  {
    return (Entity&) *git;
  }
  int id () const
  {
    switch (datamode)
    {
    case GnuplotOptions::conforming:
      return number[renumberedqvmIndex()];
    case GnuplotOptions::nonconforming:
    {if (codim==n)
      return offset + renumber(*git,subindex,codim);
    else
      return offset + renumber(*git,subindex,codim) + index - subindex;
    }
    default:
      DUNE_THROW(IOError,"GnuplotWriter: unsupported DataMode" << datamode);
    }
  }
  int localindex () const
  {
    return index;
  }
  const FieldVector<DT,n> & position () const
            {
    return ReferenceElements<DT,n>::general(git->type()).position(subindex,codim);
            }
};

QuadVertexIterator quadvertexBegin () const
{
  return QuadVertexIterator( gridView_.template begin< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >(),
      datamode, *quadvertexmapper, number );
}

QuadVertexIterator quadvertexEnd () const
{
  return QuadVertexIterator( gridView_.template end< 0, Gnuplot_Partition >(),
      gridView_.template end< 0, Gnuplot_Partition >(),
      datamode, *quadvertexmapper, number );
}

private:
/** \brief take a vector and interpret it as cell data
          \ingroup Gnuplot
 */
template<class V>
class P0VectorWrapper : public GnuplotFunction
{
  typedef MultipleCodimMultipleGeomTypeMapper< GridView, P0Layout > VM0;
public:
  //! return number of components
  virtual int ncomps () const
  {
    return 1;
  }

  //! evaluate
  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    return v[mapper.map(e)];
  }

  //! get name
  virtual std::string name () const
  {
    return s;
  }

  //! construct from a vector and a name
  P0VectorWrapper ( const Grid &g_, const IndexSet &is_, const V &v_, std::string s_)
  : g(g_), is(is_), v(v_), s(s_), mapper(g_,is_)
  {
    if (v.size()!=(unsigned int)mapper.size())
      DUNE_THROW(IOError,"GnuplotWriter::P0VectorWrapper: size mismatch");
  }

  virtual ~P0VectorWrapper() {}

private:
  const Grid& g;
  const IndexSet &is;
  const V& v;
  std::string s;
  VM0 mapper;
};

/** \brief take a vector and interpret it as vertex data
          \ingroup Gnuplot
 */
template<class V>
class P1VectorWrapper : public GnuplotFunction
{
  typedef MultipleCodimMultipleGeomTypeMapper< GridView, P1Layout > VM1;

public:
  //! return number of components
  virtual int ncomps () const
  {
    return 1;
  }

  //! evaluate
  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    double min=1E100;
    int imin=-1;
    Dune::GeometryType gt = e.type();
    for (int i=0; i<e.template count<n>(); ++i)
    {
      Dune::FieldVector<DT,n>
      local = Dune::ReferenceElements<DT,n>::general(gt).position(i,n);
      local -= xi;
      if (local.infinity_norm()<min)
      {
        min = local.infinity_norm();
        imin = i;
      }
    }
    return v[mapper.template map<n>(e,imin)];
  }

  //! get name
  virtual std::string name () const
  {
    return s;
  }

  //! construct from a vector and a name
  P1VectorWrapper ( const Grid &g_, const IndexSet &is_, const V &v_, std::string s_ )
  : g(g_), is(is_), v(v_), s(s_), mapper(g_,is_)
  {
    if (v.size()!=(unsigned int)mapper.size())
      DUNE_THROW(IOError,"GnuplotWriter::P1VectorWrapper: size mismatch");
  }

  virtual ~P1VectorWrapper() {}

private:
  const Grid& g;
  const IndexSet &is;
  const V& v;
  std::string s;
  VM1 mapper;
};

public:
/**
 * @brief Construct a GnuplotWriter working on a specific GridView.
 *
 *
 * @param gridView The gridView the grid functions live on. (E. g. a LevelGridView.)
 * @param dm The data mode.
 */
explicit GnuplotWriter ( const GridView &gridView,
    GnuplotOptions::Order od,
    GnuplotOptions::Info in,
    GnuplotOptions::DataMode dm = GnuplotOptions::conforming )
: gridView_( gridView ),
  grid( gridView.grid() ),
  is( gridView_.indexSet() ),
  order( od ),
  inf( in ),
  datamode( dm )
{
  indentCount = 0;
  numPerLine = 4*3; //should be a multiple of 3 !
}

/**
 * @brief Add a grid function that lives on the cells of the grid to the visualization.
 * @param p The function to visualize.  The GnuplotWriter object will take
 *          ownership of the GnuplotFunction *p and delete it when it's done.
 */
void addCellData (GnuplotFunction* p)
{
  celldata.push_back(p);
}

/**
 * @brief Add a grid function (represented by container) that lives on the cells of
 * the grid to the visualization.
 *
 * The container has to have random access via operator[] (e. g. std::vector). The
 * value of the grid function for an arbitrary element
 * will be accessed by calling operator[] with the id of the element.
 *
 * @param v The container with the values of the grid function for each cell.
 * @param name A name to identify the grid function.
 */
template<class V>
void addCellData (const V& v, std::string name)
{
  GnuplotFunction* p = new P0VectorWrapper<V>(grid,is,v,name);
  celldata.push_back(p);
}

/**
 * @brief Add a grid function that lives on the vertices of the grid to the visualization.
 * @param p The function to visualize.  The GnuplotWriter object will take
 *          ownership of the GnuplotFunction *p and delete it when it's done.
 */
void addVertexData (GnuplotFunction* p)
{
  vertexdata.push_back(p);
}

/**
 * @brief Add a grid function (represented by container) that lives on the cells of the
 * grid to the visualization output.
 *
 * The container has to have random access via operator[] (e. g. std::vector). The value
 * of the grid function for an arbitrary element
 * will be accessed by calling operator[] with the id of the element.
 *
 * @param v The container with the values of the grid function for each cell.
 * @param name A name to identify the grid function.
 */
template<class V>
void addVertexData (const V& v, std::string name)
{
  GnuplotFunction* p = new P1VectorWrapper<V>(grid,is,v,name);
  vertexdata.push_back(p);
}

//! clear list of registered functions
void clear ()
{
  for (FunctionIterator it=celldata.begin();
      it!=celldata.end(); ++it)
    delete *it;
  celldata.clear();
  for (FunctionIterator it=vertexdata.begin();
      it!=vertexdata.end(); ++it)
    delete *it;
  vertexdata.clear();
}

//! destructor
virtual ~GnuplotWriter ()
{
  this->clear();
}

/** \brief write output (interface might change later)
 *
 *  \param[in]  name  basic name to write (may not contain a path)
 *  \param[in]  type  type of output (e.g,, ASCII) (optional)
 */
std::string write ( const std::string &name,
    GnuplotOptions::Info inf,
    GnuplotOptions::OutputType type = GnuplotOptions::ascii )
{
  return write( name, inf, type, gridView_.comm().rank(), gridView_.comm().size() );
}

/** \brief write output (interface might change later)
 *
 *  \param[in]  name        basic name to write (may not contain a path)
 *  \param[in]  path        path to data output
 *  \param[in]  extendpath  path keyword for each process
 *  \param[in]  type        type of output (e.g,, ASCII) (optional)
 */
std::string pwrite ( const char* name,  const char* path, const char* extendpath,
    GnuplotOptions::OutputType type = GnuplotOptions::ascii )
{
  return pwrite( name, path, extendpath, type, gridView_.comm().rank(), gridView_.comm().size() );
}

protected:
std::string write ( const std::string &name,
    GnuplotOptions::Info in,
    GnuplotOptions::OutputType type,
    const int commRank,
    const int commSize )
{
  // make data mode visible to private functions
  outputtype = type;
  inf = in;

  // generate filename for process data
  std::ostringstream pieceName, pieceName2;
  if( commSize > 1 )
  {
    pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << commSize << ":";
    pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << commRank << ":";
  }
  pieceName << name << ".data";

  // write process data
  std::ofstream file;
  file.open( pieceName.str().c_str() );
  pieceName2 << name << ".gnu";

  // write process data
  std::ofstream filecmd;
//  filecmd.open( pieceName2.str().c_str() );
  writeDataFile( file, filecmd, name.c_str(),inf );
  file.close();  // filecmd.close();
  char gcmd[512];
  const char* gnuplot;
  gnuplot = getenv("GNUPLOT");
  if ( gnuplot == nullptr ) { gnuplot="gnuplot"; };
  sprintf(gcmd,"%s << //\nload \"%s.gnu\"\n//",gnuplot,name.c_str());
//  std::cout << "executing " << std::string(gcmd) << std::endl;
  system(gcmd);

  // for serial jobs we're done here
  if( commSize == 1 )
    return pieceName.str();

  // synchronize processes
  gridView_.comm().barrier();

  // generate name of parallel header
  std::ostringstream parallelName;
  parallelName << "s" << std::setfill( '0' ) << std::setw( 4 ) << commSize << ":";
  parallelName << name << (GridView::dimension > 1 ? ".pvtu" : ".pvtp");

  // on process 0: write out parallel header

  // synchronize processes
  gridView_.comm().barrier();
  return parallelName.str();
}

private:
//quadratic cell types
enum GnuplotGeometryType
{
  vtkLine = 3,
      vtkTriangle = 5,
      vtkQuadrilateral = 9,
      vtkTetrahedron = 10,
      vtkHexahedron = 12,
      vtkPrism = 13,
      vtkPyramid = 14,
      vtkquadLine = 21,
      vtkquadTriangle = 22,
      vtkquadQuadrilateral = 23,
      vtkquadTetrahedron = 24,
      vtkquadHexahedron = 25,
};

protected:
//! mapping from GeometryType to GnuplotGeometryType
static GnuplotGeometryType vtkType(const Dune::GeometryType & t, GnuplotOptions::Order od)
{if (od == GnuplotOptions::quadratic) {
  if (t.isLine())
    return vtkquadLine;
  if (t.isTriangle())
    return vtkquadTriangle;
  if (t.isQuadrilateral())
    return vtkquadQuadrilateral;
  if (t.isTetrahedron())
    return vtkquadTetrahedron;
  if (t.isHexahedron())
    return vtkquadHexahedron;
  DUNE_THROW(IOError,"GnuplotWriter: unsupported GeometryType " << t);
}
else {
  if (t.isLine())
    return vtkLine;
  if (t.isTriangle())
    return vtkTriangle;
  if (t.isQuadrilateral())
    return vtkQuadrilateral;
  if (t.isTetrahedron())
    return vtkTetrahedron;
  if (t.isPyramid())
    return vtkPyramid;
  if (t.isPrism())
    return vtkPrism;
  if (t.isHexahedron())
    return vtkHexahedron;
  DUNE_THROW(IOError,"GnuplotWriter: unsupported GeometryType " << t);
}
}

private:
//! write data file to stream
void writeDataFile (std::ostream& s,std::ostream& c, std::string name, GnuplotOptions::Info inf)
 {
  // Grid characteristics
  vertexmapper = new VertexMapper(gridView_);
  quadvertexmapper = new QuadVertexMapper(gridView_);
  if (datamode == GnuplotOptions::conforming)
  {
    if (order == GnuplotOptions::quadratic)
      number.resize(quadvertexmapper->size());
    else
      number.resize(vertexmapper->size());
    for (std::vector<int>::size_type i=0; i<number.size(); i++)
      number[i] = -1;
  }
  countEntities(nvertices, ncells, ncorners);
//   // PointData and Values
   writeGnuplotData(s,c,name,inf);
 }

protected:
//! count the vertices, cells and corners
//! with nonconforming datamode the number of vertices and corners is the same
virtual void countEntities(int &nvertices, int &ncells, int &ncorners)
{
  nvertices = 0;
  ncells = 0;
  ncorners = 0;
  if (order == GnuplotOptions::quadratic) {
    for (CellIterator it=cellBegin(); it!=cellEnd(); ++it) {
      ncells++;
      for (int i=0; i<it->template count<n>(); ++i) {
        ncorners++;
        if (datamode == GnuplotOptions::conforming) {
          int alpha = quadvertexmapper->map(*it,i,n);
          if (number[alpha]<0)
            number[alpha] = nvertices++;
        } else {
          nvertices++;
        }
      }
      for (int i=0; i<it->template count<n-1>(); ++i){
        ncorners++;
        if (datamode == GnuplotOptions::conforming) {
          int alpha = quadvertexmapper->map(*it,i,n-1);
          if (number[alpha]<0)
            number[alpha] = nvertices++;
        }
        else
          nvertices++;
      }
    }
  } else /* linear */ {
    for (CellIterator it=cellBegin(); it!=cellEnd(); ++it) {
      ncells++;
      for (int i=0; i<it->template count<n>(); ++i) {
        ncorners++;
        if (datamode == GnuplotOptions::conforming) {
          int alpha = vertexmapper->map(*it,i,n);
          if (number[alpha]<0)
            number[alpha] = nvertices++;
        } else
          nvertices++;
      }
    }
  }
}

template <int dim>
class SortNode
{
  public:
  SortNode* next;
  SortNode(Dune::FieldVector<double,dim> _coord, double* _values): coord(_coord), values(_values) {};
  ~SortNode() {free(values);};
  Dune::FieldVector<double,dim> getCoord() {return coord;};
  double* getValues() {return values;};
  private:
  Dune::FieldVector<double,dim> coord;
  double* values;
};

//! write vertex data
virtual void writeGnuplotData (std::ostream& s, std::ostream& c, std::string name, GnuplotOptions::Info inf)
{
  using std::max;
  const int dim=GridView::dimension;
  const int ptdim=max(2,dim);
  double mean, stddev, deviation, uh, pt[ptdim];
  int nver; 
  SortNode<dim> *cbegin, *clatest, *csearch, *cprev, *citer;
  double* vnew;
  cbegin=nullptr;  clatest=nullptr;
  for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it) {
    int ncomps = (*it)->ncomps();
    double *m1=(double *)(malloc(ncomps*sizeof(double))),
           *m2=(double *)(malloc(ncomps*sizeof(double))),
           *umin=(double *)(malloc(ncomps*sizeof(double))),
           *umax=(double *)(malloc(ncomps*sizeof(double))),
           *ptxmin=(double *)(malloc(ncomps*sizeof(double))),
           *ptymin=(double *)(malloc(ncomps*sizeof(double))),
           *ptxmax=(double *)(malloc(ncomps*sizeof(double))),
           *ptymax=(double *)(malloc(ncomps*sizeof(double)));
    for(int i=0; i<ncomps; i++) { m1[i]=0.0; m2[i]=0.0; umin[i]=DBL_MAX; umax[i]=-DBL_MAX; };
    if (order == GnuplotOptions::quadratic) {
      QuadVertexIterator qvEnd = quadvertexEnd(), qvBeg=quadvertexBegin();
      for (QuadVertexIterator qvit=qvBeg; qvit!=qvEnd; ++qvit) {
        /* sort points by first point component */
        vnew=(double *)(malloc(ncomps*sizeof(double)));
        for (int j=0; j<ncomps; j++) {
          vnew[j]=(*it)->evaluate(j,*qvit,qvit.position());
        };
        SortNode<dim> *cnew = new SortNode<dim>(qvit->geometry().global(qvit.position()),vnew);
        if ( clatest==nullptr )
        {
          cbegin=cnew; clatest=cnew; cnew->next=nullptr;
        }
        else
        {
          if ( cnew->getCoord()[0] > clatest->getCoord()[0] )
            {csearch=clatest->next; cprev=clatest;}
          else
            {csearch=cbegin; cprev=nullptr;};
          if ( csearch == nullptr )
            { cnew->next=nullptr; clatest->next=cnew;}
          else
          {
            while ( csearch != nullptr )
              if ( cnew->getCoord()[0] > csearch->getCoord()[0] )
                {cprev=csearch; csearch=csearch->next;}
              else
                { csearch=nullptr;};
            if ( cprev==nullptr )
              {cnew->next=cbegin; cbegin=cnew;}
            else
              {cnew->next=cprev->next; cprev->next=cnew;};
          };
          clatest=cnew;
        };
      };
    } else {
      VertexIterator vEnd = vertexEnd();
      for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit) {
        /* sort points by first point component */
        vnew=(double *)(malloc(ncomps*sizeof(double)));
        for (int j=0; j<ncomps; j++) {
          vnew[j]=(*it)->evaluate(j,*vit,vit.position());
        };
        SortNode<dim> *cnew = new SortNode<dim>(vit->geometry().global(vit.position()),vnew);
        if ( clatest==nullptr )
        {
          cbegin=cnew; clatest=cnew; cnew->next=nullptr;
        }
        else
        {
          if ( cnew->getCoord()[0] > clatest->getCoord()[0] )
            {csearch=clatest->next; cprev=clatest;}
          else
            {csearch=cbegin; cprev=nullptr;};
          if ( csearch == nullptr )
            { cnew->next=nullptr; clatest->next=cnew; }
          else
          {
            while ( csearch != nullptr )
              if ( cnew->getCoord()[0] > csearch->getCoord()[0] )
                {cprev=csearch; csearch=csearch->next;}
              else
                { csearch=nullptr;};
            if ( cprev==nullptr )
              {cnew->next=cbegin; cbegin=cnew;}
            else
              {cnew->next=cprev->next; cprev->next=cnew;};
          };
          clatest=cnew;
        };
      };
    };
    /* output sorted points */
    citer=cbegin;  nver=0;
    while ( citer != nullptr )
    {
      nver++;
      Dune::FieldVector<double,dim>  pt=citer->getCoord();
      for (int j=0; j<dim; j++) s << pt[j] << " " ;
      double* u=citer->getValues();
      for (int j=0; j<ncomps; j++) {
        uh = u[j];
        s << uh << " " ; m1[j] += uh; m2[j] += uh*uh;
        if ( uh<umin[j] ) { umin[j]=uh; ptxmin[j]=pt[0]; if (dim>1) ptymin[0]=pt[1]; };
        if ( uh>umax[j] ) { umax[j]=uh; ptxmax[j]=pt[0]; if (dim>1) ptymax[0]=pt[1]; };
      };
      s << std::endl;
      cprev=citer;
      citer=citer->next;
      delete cprev;
    };
    /* infomational output follows */
    if (inf != GnuplotOptions::none) {
      std::cout << "*** GnuplotWriter diagnostics ***" << std::endl;
      for (int j=0; j<ncomps; j++) {
        mean = m1[j]/nver;  stddev=sqrt(m2[j]/nver-mean*mean);
        std::cout << "mean value    of u component " << j << ": " << mean << std::endl;
        std::cout << "std-deviation of u component " << j << ": " << stddev << std::endl;
        if ( dim==2 )
        {
          std::cout << "minimum u of component " << j << ": " << umin[j] << " at point (" << ptxmin[j] << "," << ptymin[j] << ")" << std::endl;
          std::cout << "  (deviation is " << umin[j]-mean << ")" << std::endl;
          std::cout << "maximum u of component " << j << ": " << umax[j] << " at point (" << ptxmax[j] << "," << ptymax[j] << ")" << std::endl;
        } 
        else if ( dim==1 )
        {
          std::cout << "minimum u of component " << j << ": " << umin[j] << " at abzissa " << ptxmin[j] << std::endl;
          std::cout << "  (deviation is " << umin[j]-mean << ")" << std::endl;
          std::cout << "maximum u of component " << j << ": " << umax[j] << " at abzissa " << ptxmax[j] << std::endl;
        };
        std::cout << "  (deviation is " << umax[j]-mean << ")" << std::endl;
        if (inf == GnuplotOptions::detail) {
          if (order == GnuplotOptions::quadratic) {
            QuadVertexIterator qvEnd = quadvertexEnd(), qvBeg=quadvertexBegin();
            for (QuadVertexIterator qvit=qvBeg; qvit!=qvEnd; ++qvit) {
              uh = (*it)->evaluate(j,*qvit,qvit.position());
              deviation = uh-mean;
              if ( fabs(deviation)>1.5*stddev ) {
                for (int k=0; k<dim; k++) {
                  pt[k] = qvit->geometry().global(qvit.position())[k];
                };
                if ( dim==2 )
                {
                  std::cout << "deviation at (" << pt[0] << "," << pt[1] << ") exceeds 1.5*std-deviation: u=" << uh <<
                    ",deviation=" << deviation << std::endl;
                }
                else if ( dim==1 )
                {
                  std::cout << "deviation at " << pt[0] <<  " exceeds 1.5*std-deviation: u=" << uh <<
                    ",deviation=" << deviation << std::endl;
                }
              }
            }
          } else {
            VertexIterator vEnd = vertexEnd();
            for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit) {
              uh = (*it)->evaluate(j,*vit,vit.position());
              deviation = uh-mean;
              if ( fabs(deviation)>1.5*stddev ) {
                  for (int k=0; k<dim; k++) {
                    pt[k] = vit->geometry().global(vit.position())[k];
                  };
                  if ( dim==2 )
                  {
                    std::cout << "deviation at (" << pt[0] << "," << pt[1] << ") exceeds 1.5*std-deviation: u=" << uh <<
                      ",deviation=" << deviation << std::endl;
                  }
                  else if ( dim==1 )
                  {
                    std::cout << "deviation at " << pt[0] <<  " exceeds 1.5*std-deviation: u=" << uh <<
                      ",deviation=" << deviation << std::endl;
                  }
                }
              }
            }
          }
        };
      std::cout << "*********************************" << std::endl;
    };
    free(m1); free(m2); free(umin); free(umax); free(ptxmin); free(ptymin); free(ptxmax); free(ptymax);
  };
  
  
// c << "set terminal x11" << std::endl <<
//   "set dgrid3d 50,50 splines" << std::endl <<
//   "set dummy u,v" << std::endl <<
//   "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000" << std::endl <<
//   "set parametric" << std::endl <<
//   "set style fill solid" << std::endl <<
//   "set title \"" << name << "\"" << std::endl <<
//   "set pm3d map" << std::endl <<
//   "splot \"" << name << ".data\"" << std::endl <<
//   "pause 15" << std::endl;
}

// renumber Gnuplot -> Dune
static int renumber (const GeometryType &t, int i, int codim)
{
  static const int triangleRenumbering[2][3] = { {0,1,2}, {0,2,1} };
  static const int quadRenumbering[2][4] = { {0,1,3,2}, {2,1,3,0} };
  static const int hexaRenumbering[2][12] = { {0,1,3,2,4,5,7,6}, {8,5,9,4,10,7,11,6,0,1,3,2} };
  static const int prismRenumbering[6] = {0,2,1,3,5,4};
  static const int tetrahedronRenumbering[3][6] = { {0,1,2,3}, {0,2,1,3,4,5}, {3,2,1,0} } ;
  switch (vtkType(t,GnuplotOptions::linear))
  {
  case vtkQuadrilateral:
    return quadRenumbering[n-codim][i];
  case vtkHexahedron:
    return hexaRenumbering[n-codim][i];
  case vtkPrism:
    return prismRenumbering[i];
  case vtkTriangle:
    return triangleRenumbering[n-codim][i];
  case vtkTetrahedron:
    return tetrahedronRenumbering[n-codim][i] ;
  default:
    return i;
  }
}
static int renumber (const Entity& e, int i, int codim)
{
  return renumber(e.type(), i, codim);
}
static int renumber (const Entity& e, int i) {
  return renumber(e.type(), i, n);
}
static int renumber (const GeometryType& t, int i) {
  return renumber(t, i, n);
}



// the list of registered functions
std::list<GnuplotFunction*> celldata;
std::list<GnuplotFunction*> vertexdata;

private:
// the grid
GridView gridView_;
const Grid& grid;

// the indexset
const IndexSet& is;

// indend counter
int indentCount;
int numPerLine;

// temporary grid information
protected:
int ncells;
int nvertices;
int ncorners;
private:
VertexMapper* vertexmapper;
QuadVertexMapper* quadvertexmapper;
std::vector<int> number;
protected:
GnuplotOptions::Order order;
protected:
GnuplotOptions::OutputType outputtype;
GnuplotOptions::Info inf;
private:
GnuplotOptions::DataMode datamode;
};

/** \brief GnuplotWriter on the leaf grid
      \ingroup Gnuplot
 */
template< class Grid >
class LeafGnuplotWriter
: public GnuplotWriter< typename Grid::LeafGridView >
{
  typedef GnuplotWriter< typename Grid::LeafGridView > Base;

public:
  /** \brief Construct a Gnuplot writer for the leaf level of a given grid */
  explicit LeafGnuplotWriter ( const Grid &grid,
      GnuplotOptions::DataMode dm = GnuplotOptions::conforming ) DUNE_DEPRECATED
      : Base( grid.leafView(), dm )
        {}
};

/** \brief GnuplotWriter on a given level grid
      \ingroup Gnuplot
 */
template< class Grid >
class LevelGnuplotWriter
: public GnuplotWriter< typename Grid::LevelGridView >
{
  typedef GnuplotWriter< typename Grid::LevelGridView > Base;

public:
  /** \brief Construct a Gnuplot writer for a certain level of a given grid */
  LevelGnuplotWriter ( const Grid &grid, int level,
      GnuplotOptions::DataMode dm = GnuplotOptions::conforming ) DUNE_DEPRECATED
      : Base( grid.levelView( level ), dm )
        {}
};

#endif
