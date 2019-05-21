#ifndef MYAMIRAGRID_HH
#define MYAMIRAGRID_HH

#include <amiramesh/AmiraMesh.h>

/**
 * \cond internals
 */
namespace  {

template <class value_type> struct TypeId {};

template <> struct TypeId<unsigned char>        { static int const id = HxBYTE; };
template <> struct TypeId<int>                  { static int const id = HxINT32; };
template <> struct TypeId<float>                { static int const id = HxFLOAT; };
template <> struct TypeId<std::complex<float> > { static int const id = HxCOMPLEX; };


template <class value_type, int components>
void extractData(AmiraMesh::Data* dataObj, std::vector<Dune::FieldVector<value_type,components> >& data) 
{
  size_t const n= dataObj->location()->nNodes();
  data.clear();
  data.reserve(n);
  data.insert(data.begin(),n,Dune::FieldVector<value_type,components>(0));
  
  value_type const* dataPtr = static_cast<value_type const*>(dataObj->dataPtr());
  for (size_t i=0; i<n; ++i)
    for (int j=0; j<components; ++j)
      data[i][j] = dataPtr[i*components+j];
}

template <class value_type, int components>
void readData(AmiraMesh& am, std::string const& fname,
              char const* location, char const* name, std::vector<Dune::FieldVector<value_type,components> >& data) 
{
  int const id = TypeId<value_type>::id;
  AmiraMesh::Data* dataObj = am.findData(location,id,components,name);
  if (!dataObj) {
    
    std::cerr << "Amira mesh file " << fname << " does not contain data'" << name << "' of '";
    switch (id) {
    case HxBYTE: std::cerr << "unsigned char"; break;
    case HxINT32: std::cerr << "int"; break;
    case HxFLOAT: std::cerr << "float"; break;
    case HxCOMPLEX: std::cerr << "complex<float>"; break;
    default:  std::cerr << "unknown(" << id << ")"; break;
    }
    std::cerr << "[" << components << "]' at '" << location << "'!\n";
    abort();
  }
  extractData(dataObj,data);
}

std::pair<std::pair<int,int>,std::complex<float> > combine(Dune::FieldVector<int,2> const& nodes, std::complex<float> const& value) 
{
  return std::make_pair(std::make_pair(nodes[0],nodes[1]),value);
}

struct VertexData 
{
  VertexData(bool boundary_, Dune::FieldVector<double,3> pos_, int oldNr_):
    boundary(boundary_), pos(pos_), oldNr(oldNr_) 
  {}

  bool boundary;
  Dune::FieldVector<double,3> pos;
  int oldNr;
};

bool boundaryVertex(VertexData const&  v) 
{
  return v.boundary; // belongs to both body and non-body
}

} // End of namespace 
/**
 * \endcond 
 */




template <class FACTORY>
void ReadAmiraGrid(FACTORY &factory, std::string const& gridfile) {
  std::vector<Dune::FieldVector<float,3> >               vertices;
  std::vector<Dune::FieldVector<int,4> >                 tetrahedra;
  std::vector<Dune::FieldVector<int,3> >                 boundary;
  std::vector<Dune::FieldVector<unsigned char,1> >       tissue;
//  std::vector<Dune::FieldVector<int,2> >                 edges;
  
  // First open grid file.
  std::cout << "Reading amira mesh\n";
  std::cout << "opening " << gridfile << '\n';
  AmiraMesh* am = AmiraMesh::read(gridfile.c_str());

  // Obtain vertices, tetrahedra, tissues, and edges.
  readData(*am,gridfile,"Nodes","Coordinates",vertices);
  readData(*am,gridfile,"Tetrahedra","Nodes",tetrahedra);
  readData(*am,gridfile,"Tetrahedra","Materials",tissue);
//  readData(*am,gridfile,"Edges","fromTo",edges);
  readData(*am,gridfile,"BoundaryTriangles","Nodes",boundary);

  
  // close grid file, open sar file
  std::free(am);
  assert(tissue.size()==tetrahedra.size());

  // We use start indices 0. Amira uses 1. Correct indices.
  for (size_t i=0; i<tetrahedra.size(); ++i) {
    for (int j=0; j<4; ++j) {
      tetrahedra[i][j] -= 1;
      if (tetrahedra[i][j]>=vertices.size())
        std::cerr << "tetrahedron node " << i << "[" << j << "]=" << tetrahedra[i][j] << " exceeds vertex numbers!\n";
    }
  }
//   for (size_t i=0; i<edges.size(); ++i) {
//     edges[i][0] -= 1;
//     edges[i][1] -= 1;
//     if (edges[i][0]>=vertices.size() || edges[i][1]>=vertices.size())
//       std::cerr << "edge " << i << "=[" << edges[i][0] << ',' << edges[i][1] << "] exceeds edge numbers!\n";
//   }
  for (size_t i=0; i<boundary.size(); ++i)
    for (int j=0; j<3; ++j)
      boundary[i][j] -= 1;
  

  // Occasionally tissue types start by 0. Correct this to one.
  int min = 10000;
  for (size_t i=0; i<tissue.size(); ++i) {  
    min = std::min(min,static_cast<int>(tissue[i][0]));
  }
  if (min==0)
    for (size_t i=0; i<tissue.size(); ++i) 
      tissue[i] += 1;
  for (size_t i=0; i<tissue.size(); ++i) 
      if (tissue[i]<=0 || tissue[i]>=8)
        std::cerr << "Suspect tissue[" << i << "] = " << tissue[i] << '\n';
  
  
  // Spatial coordinates in amira meshes are usually given in
  // cm. Correct this to SI units.  This does also affect the edge
  // element shape functions, such that the shape functions have to be
  // scaled along.
//   for (size_t i=0; i<vertices.size(); ++i)
//     vertices[i] /= 100;

  // Create grid. Note that we take only those tetrahedra that are
  // contained in the body, which is detected by their tissue
  // type. Boundary vertices are moved to the front, since this is the
  // UG sort order.
  std::cout << "preparing grid creation\n";
  
  std::vector<VertexData> vertexData;
  for (size_t i=0; i<vertices.size(); ++i) {
    Dune::FieldVector<double,3> x;
    for (int j=0; j<3; ++j) x[j] = vertices[i][j];
    vertexData.push_back(VertexData(false,x,i));
  }
  for (size_t i=0; i<boundary.size(); ++i)
    for (int j=0; j<3; ++j)
      vertexData[boundary[i][j]].boundary = true;
  
  // move boundary vertices to the front, adjust tetrahedra and edge indices
  std::partition(vertexData.begin(),vertexData.end(),&boundaryVertex);
  std::vector<int> newVertexNumber(vertices.size(),-1);
  for (size_t i=0; i<vertexData.size(); ++i)
    newVertexNumber[vertexData[i].oldNr] = i;
  assert(*std::min_element(newVertexNumber.begin(),newVertexNumber.end())==0);

  for (size_t i=0; i<tetrahedra.size(); ++i)
    for (int j=0; j<4; ++j)
      tetrahedra[i][j] = newVertexNumber[tetrahedra[i][j]];
//   for (size_t i=0; i<edges.size(); ++i)
//     for (int j=0; j<2; ++j)
//       edges[i][j] = newVertexNumber[edges[i][j]];
  

  std::cout << "inserting vertices and elements\n";

  for (size_t i=0; i<vertexData.size(); ++i) {
    factory.insertVertex(vertexData[i].pos);
  }
  
  for (size_t i=0; i<tetrahedra.size(); ++i) {
    std::vector<unsigned int> tmp(4);
    for (int j=0; j<4; ++j) {
      assert(tetrahedra[i][j]>=0 && tetrahedra[i][j]<vertexData.size());
      tmp[j] = tetrahedra[i][j];
    }
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3),tmp);
  }
}

#endif
