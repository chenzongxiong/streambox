/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cctype>
#include <sstream>
#include <cmath>

#include "io/readPoly.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/rotate.hh" //rotates 2D objects (so far) and allows to translate them in y-direction; to be used
//after reading in a grid, before handing over to grid manager

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

void readNodeData(std::string const& name, std::vector<Dune::FieldVector<double,2> >& nodes, int mortarCheck, Dune::BitSetVector<1> * obsField,
                  std::array<double,3>* pelvis)
{
    using namespace std;

    //open file with nodes
    ifstream in(name.c_str());
    if (!in) // check for successful opening
        throw Kaskade::FileIOException("Cannot open for reading.",name,__FILE__,__LINE__);

    string line;
    readLine(in,line);
    istringstream is(line);

    //Read out the number of vertices, the dimension, number of attributes and boundary markers from the first line of the file.
    int nVertices, dimension, nAttributes, nBoundaryMarkers;
    is >> nVertices >> dimension >> nAttributes >> nBoundaryMarkers;
    assert(dimension==2); // TODO: throw exception
    assert(0<=nBoundaryMarkers && nBoundaryMarkers<=1);

    if(nBoundaryMarkers==0 && obsField!=nullptr)
        std::cerr << "No boundary markers!\n" << "The BitSetVector obsField isn't properly initialised" << std::endl;

    //Check whether obsField is to be written
    nodes.resize(nVertices);
    if(obsField!=nullptr){
        obsField->clear();
        obsField->resize(nVertices);
    }

    for (int i=0; i<nVertices; ++i) {
        readLine(in,line);
        istringstream is(line);
        int id; // the node number
        is >> id;
        assert(0<=id && id<=nVertices); // TODO: throw
        // Note that .poly files may be 0-based or 1-based. We just wrap around to cover both cases.
        id = id % nVertices;
        // Read in the coordinates
        is >> nodes[id][0] >> nodes[id][1];

        //stores obstical field data (pelvis to be deleted; for testing included)
        if(obsField!=nullptr){
            //Boundary markers need to exist in order to initialise obsField
            if(nBoundaryMarkers==0)
                continue;

            if(nAttributes!=0){
                //Discard attribute
                float dummy;
                is >> dummy;
            }

            //Read boundary marker
            int n;
            is >> n;

            //Check for boundary (manually)
            if(n==1 && pelvis==nullptr){
                if(mortarCheck){  //lower body, gets mortarField
                    if(nodes[id][1]>-0.4)
                        (*obsField)[id][0] = n;
                }
                else{ //upper body, gets obsField
                    if(nodes[id][1]<0.2)
                        (*obsField)[id][0] = n;
                }
            }

            //If we look at the pelvis, the condition for setting field to 1 is to look in a certain radius of the barycenter in the joint region
            else{
                double radius = sqrt(pow(nodes[id][0]-std::get<0>(*pelvis),2)+pow(nodes[id][1]-std::get<1>(*pelvis),2));
                if(radius<=std::get<2>(*pelvis))
                    (*obsField)[id][0] = n;
            }
        }
    }
}
}
 // End of namespace

std::unique_ptr<Dune::UGGrid<2> > readPolyData(std::string const& name, int heapSize, double alpha, std::array<double,2>* barycenter, double * deviation,
                                               int mortarCheck, Dune::BitSetVector<1> * obsField, std::array<double,3>* pelvis)
{
  using namespace std;
  using namespace Kaskade;
  using Grid = Dune::UGGrid<2>;
  
  Grid::setDefaultHeapSize(heapSize);
  Grid * grid(new Grid);
  Dune::GridFactory<Grid > factory(grid);
  
  
//   grid->createBegin();
//     
   // read node file
  vector<Dune::FieldVector<double,2> > nodes;
  readNodeData(name+".node",nodes,mortarCheck,obsField,pelvis);
  if(alpha!=0){

      //Determine the biggest y coordinate to restore this position after rotation
      double y_max = std::numeric_limits<double>::lowest();

      for (int i=0; i<nodes.size(); ++i)
      {
          y_max = std::max(y_max, (nodes[i])[1]);
      }

      try{
          find2DBarycenter(nodes, barycenter);
          // Coordinates for 2D hip joint model
          //        barycenter[0] = 88.33;
          //        barycenter[1] = -104.47;
          rotate(alpha, barycenter, nodes);
      }
      catch(std::exception& e) {
          std::cerr << "error: " << e.what() << "\n";
          return nullptr;
      }

      *deviation = std::numeric_limits<double>::lowest();

      for (int i=0; i<nodes.size(); ++i)
      {
          *deviation = std::max(*deviation, (nodes[i])[1]);
      }

      *deviation -= y_max;

      if(*deviation != 0)
          translate(*deviation, nodes);
  }

  for (int i=0; i<nodes.size(); ++i) 
    factory.insertVertex(nodes[i]);

  string eleName = name+".ele" ;
  ifstream in(eleName.c_str());
  if (!in)
    throw Kaskade::FileIOException("Cannot open for reading.",name,__FILE__,__LINE__);
  
  string line;
  readLine(in,line);
  istringstream is(line);

  int nEle, nNode, nAttributes;
  is >> nEle >> nNode >> nAttributes;

  assert(nNode>=3); // TODO throw

  std::vector<unsigned int> elemNodes(3);
  for (int i=0; i<nEle; ++i) {
    readLine(in,line);
    istringstream is(line);
    int id;
    is >> id >> elemNodes[0] >> elemNodes[1] >> elemNodes[2];
    assert(0<=elemNodes[0] && elemNodes[0]<=nodes.size() && 
           0<=elemNodes[1] && elemNodes[1]<=nodes.size() && 
           0<=elemNodes[2] && elemNodes[2]<=nodes.size()); // TODO: throw
    // Remember that node numbers can be 0-based or 1-based. As in readNodeData, we just wrap around.
    for (int j=0; j<3; ++j) 
      elemNodes[j] = elemNodes[j]%nodes.size();
    
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),elemNodes);
  }
   
  return std::unique_ptr<Grid>( factory.createGrid() );
}
