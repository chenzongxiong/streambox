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

#include <cctype>
#include <fstream>
#include <sstream>

#include "io/readPoly.hh"
#include "utilities/detailed_exception.hh"

namespace  {

  // reads the next nonempty line from the given stream into the given string. 
  // Comments starting with '#' extend up to the next newline and are considered as "empty".
  void readLine(std::istream& in, std::string& line) 
  {
    do {
      std::getline(in,line);
      line.erase(std::min(line.find('#'),line.size()),line.size());
      
      size_t i=0;
      for (; i<line.size(); ++i)
        if (!std::isspace(line[i])) break;
        line.erase(0,i);
    } while(line.size()==0);
  }
  
  std::vector<Dune::FieldVector<double,2>> readNodeData(std::string const& name) 
  {
    using namespace std;
    
    //open file with nodes
    ifstream in(name.c_str());
    if (!in) // check for successful opening
      throw Kaskade::FileIOException("Cannot open for reading.",name,__FILE__,__LINE__);
    
    string line;
    readLine(in,line);
    istringstream is(line);
    
    int nVertices, dimension, nAttributes, nBoundaryMarkers;
    is >> nVertices >> dimension >> nAttributes >> nBoundaryMarkers;
    assert(dimension==2); // TODO: throw exception
    assert(0<=nBoundaryMarkers && nBoundaryMarkers<=1);
    
    std::vector<Dune::FieldVector<double,2>> nodes(nVertices);
    for (int i=0; i<nVertices; ++i) {
      readLine(in,line);
      istringstream is(line);
      int id; // the node number
      is >> id;
      
      if (id<0 || id>nVertices)
        throw Kaskade::FileIOException("Node index out of range.",name,__FILE__,__LINE__); 
      
      // Note that .poly files may be 0-based or 1-based. We just wrap around to cover both cases.
      id = id % nVertices;
      // Read in the coordinates
      is >> nodes[id][0] >> nodes[id][1];
    }
    
    return nodes;
  }
} // End of namespace 

std::unique_ptr<Dune::UGGrid<2>> readPolyData(std::string const& name, int heapSize, int envSize) 
{
  using namespace std;
  
  Dune::UGGrid<2>::setDefaultHeapSize(heapSize);
  Dune::UGGrid<2>* grid(new Dune::UGGrid<2>);
  Dune::GridFactory<Dune::UGGrid<2>> factory(grid);
  
  
   // read node file
  auto nodes = readNodeData(name+".node"); 
  for (auto node: nodes) 
    factory.insertVertex(node);

  // read elements file
  string eleName = name+".ele" ;
  ifstream in(eleName.c_str());
  if (!in)
    throw Kaskade::FileIOException("Cannot open for reading.",name,__FILE__,__LINE__);
  
  string line;
  readLine(in,line);
  istringstream is(line);

  // read header
  int nEle, nNode, nAttributes;
  is >> nEle >> nNode >> nAttributes;

  assert(nNode>=3); // TODO throw

  std::vector<unsigned int> elemNodes(3);
  for (int i=0; i<nEle; ++i) {
    readLine(in,line);
    istringstream is(line);
    int id;
    is >> id >> elemNodes[0] >> elemNodes[1] >> elemNodes[2];
    
    if (!(0<=elemNodes[0] && elemNodes[0]<=nodes.size() && 
          0<=elemNodes[1] && elemNodes[1]<=nodes.size() && 
          0<=elemNodes[2] && elemNodes[2]<=nodes.size()) )
      throw Kaskade::FileIOException("Node index of element out of range.",name,__FILE__,__LINE__);
    
    // TODO: read attributes
    
    // Remember that node numbers can be 0-based or 1-based. As in readNodeData, we just wrap around.
    for (auto& n: elemNodes) 
      n = n % nodes.size();
    
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),elemNodes);
  }
   
  return std::unique_ptr<Dune::UGGrid<2>>( factory.createGrid() );
}

  
