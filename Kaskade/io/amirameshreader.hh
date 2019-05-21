/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AMIRAMESHREADER_HH
#define AMIRAMESHREADER_HH

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <dune/common/fvector.hh>

#include <amiramesh/AmiraMesh.h>
#include "fem/lagrangespace.hh"
#include "io/ioTools.hh"
#include "io/vtk.hh"

// forward declaration
template <> class std::complex<float>;

namespace Kaskade
{
  /// Reader for Amira meshes.
  /**
   * Note that the Amira mesh file must be in ASCII format. Currently UG and ALUSimplex grids are supported.
   *
   */
  namespace AmiraMeshReader
  {
    namespace ImplementationDetails{

      /**
       * \cond internals
       */
      template <class value_type> struct TypeId{};

      template <> struct TypeId<unsigned char>        { static int const id = HxBYTE; };
      template <> struct TypeId<int>                  { static int const id = HxINT32; };
      template <> struct TypeId<float>                { static int const id = HxFLOAT; };
      template <> struct TypeId<double>               { static int const id = HxDOUBLE; };
      template <> struct TypeId<std::complex<float> > { static int const id = HxCOMPLEX; };
      /**
       * \endcond
       */


      std::string exceptionMessage(std::string const& function, std::string const& gridfile, int line)
      {
        return std::string("In AmiraMeshReader::") + function + ", line " + boost::lexical_cast<std::string>(line) + ": Error reading " + gridfile + "\n";
      }

      std::string exceptionMessage(std::string const& function, std::string const& location, std::string const& component, int line)
      {
        return std::string("In AmiraMeshReader::") + function + ", line " + boost::lexical_cast<std::string>(line)  + std::string(": Error reading ") + component + " of " + location + "\n";
      }


      template <class value_type, int components>
      void extractData(AmiraMesh::Data* dataObj, std::vector<Dune::FieldVector<value_type,components> >& data)
      {
        size_t const n= dataObj->location()->nNodes();
        std::cout << n << " entries...";
        data.clear();
        data.reserve(n);
        data.insert(data.begin(),n,Dune::FieldVector<value_type,components>(0));

        value_type const* dataPtr = static_cast<value_type const*>(dataObj->dataPtr());
        for (size_t i=0; i<n; ++i)
          for (int j=0; j<components; ++j)
            data[i][j] = dataPtr[i*components+j];

      }

      template <class value_type, int components>
      bool readData(AmiraMesh& am, char const* location, char const* name, std::vector<Dune::FieldVector<value_type,components> >& data)
      {
        std::string pname=name;
        std::cout << "trying to read " << pname << " of " << location << ": " << std::flush;
        int const id = ImplementationDetails::TypeId<value_type>::id;
        AmiraMesh::Data* dataObj = am.findData(location,id,components,name);
        if (!dataObj){
          std::cout << "empty." << std::endl << std::flush;
          return false;
        }

        extractData(dataObj,data);
        std::cout << "done." << std::endl << std::flush;
        return true;
      }


    } // End of ImplementationDetails


    /// function reading an Amira mesh file in ASCII format.
    /**
     * This function reads vertices, tetrahedra and, if defined, boundary triangles
     * of an Amira mesh file into a grid and returns this grid.
     *
     * \param Scalar scalar type of the vertex coordinates
     * \param gridfile the name of the Amira mesh file
     * \param initialGridSize optional parameter for reserving memory for the grid
     * \return pointer on the created grid
     */
    template<class Grid, class Scalar>
    std::unique_ptr<Grid> readGrid(std::string const& gridfile, int initialGridSize = 0, bool measureTime = false){
      int const dim = Grid::dimension;
      std::vector<Dune::FieldVector<Scalar,dim> >    vertices;
      std::vector<Dune::FieldVector<int,dim+1> >     cells;
      std::vector<Dune::FieldVector<int,dim> >       boundary;

      // open grid file.
      std::cout << "Reading amira mesh\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am){
        std::string const mes = ImplementationDetails::exceptionMessage("readGrid(std::string const&, int)", gridfile, __LINE__);
        throw std::runtime_error(mes);
      }


      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, int)", "Nodes", "Coordinates", __LINE__));
      }
      // read cells
      std::string cellName;
      if(dim == 2) cellName = "Triangles";
      if(dim == 3) cellName = "Tetrahedra";
      if(dim != 2 && dim != 3){
        delete(am);
        throw std::runtime_error("Sorry!\nOnly 2D and 3D grids are supported.");
      }
      exists = ImplementationDetails::readData(*am,cellName.c_str(),"Nodes",cells);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, int)", cellName, "Nodes", __LINE__));
      }
      // read boundary segments if defined
      bool boundaryExists = ImplementationDetails::readData(*am,"BoundaryTriangles","Nodes",boundary);

      // close grid file
      delete(am);

      // We use start indices 0. Amira uses 1. Correct indices.
      for (size_t i=0; i<cells.size(); ++i) {
        for (int j=0; j<dim+1; ++j) {
          cells[i][j] -= 1;
          if (cells[i][j]>=vertices.size()){
            std::string message = std::string("In ") + __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__) + ": cell node " + boost::lexical_cast<std::string>(i) + "[" + boost::lexical_cast<std::string>(j) + "]=" + boost::lexical_cast<std::string>(cells[i][j]) + " exceeds vertex numbers!\n";
            throw std::runtime_error(message);
          }
        }
      }

      if(boundaryExists)
        for (size_t i=0; i<boundary.size(); ++i)
          for (int j=0; j<dim; ++j)
            boundary[i][j] -= 1;


      // Create grid.
      std::cout << std::endl;
      std::cout << "preparing grid creation\n";

      Dune::GridFactory<Grid> factory = IOTools::FactoryGenerator<dim,Grid>::createFactory(initialGridSize);
      std::cout << "inserting vertices...";
      for (size_t i=0; i<vertices.size(); ++i) {
        Dune::FieldVector<double,dim> pos;
        for(int j=0; j<dim; ++j)
          pos[j] = (double)vertices[i][j];
        factory.insertVertex(pos);
      }
      std::cout << "done." << std::endl;

      if(boundaryExists){
        std::cout << "inserting boundary segments...";
        for(size_t i=0; i<boundary.size(); ++i){
          std::vector<unsigned int> tmp(dim);
          for(size_t j=0; j<dim; ++j){
            assert(boundary[i][j]>=0 && boundary[i][j]<vertices.size());
            tmp[j] = boundary[i][j];
          }
          factory.insertBoundarySegment(tmp);
        }
        std::cout << "done." << std::endl;
      }

      std::cout << "inserting elements...";
      for (size_t i=0; i<cells.size(); ++i) {
        std::vector<unsigned int> tmp(dim+1);
        for (int j=0; j<dim+1; ++j) {
          assert(cells[i][j]>=0 && cells[i][j]<vertices.size());
          tmp[j] = cells[i][j];
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),tmp);
      }
      std::cout << "done." << std::endl;

      boost::timer::cpu_timer timer;
      std::unique_ptr<Grid> result(factory.createGrid());
      std::cout << "grid creation finished in " << boost::timer::format(timer.elapsed()) << std::endl;
      std::cout << std::endl;

      return result;
    }


    /// function reading an Amira mesh file in ASCII format.
    /**
     * This function reads vertices, tetrahedra and, if defined, boundary triangles
     * of an Amira mesh file into a grid and returns this grid.
     *
     * \param Scalar scalar type of the vertex coordinates
     * \param gridfile the name of the Amira mesh file
     * \param initialGridSize optional parameter for reserving memory for the grid
     * \return pointer on the created grid
     */
    template<class Grid, class Scalar>
    std::unique_ptr<Grid> readGrid(std::string const& gridfile, std::vector<int>& boundaryIds, int initialGridSize = 0){
      int const dim = Grid::dimension;
      std::vector<Dune::FieldVector<Scalar,dim> >    vertices;
      std::vector<Dune::FieldVector<int,dim+1> >     cells;
      std::vector<Dune::FieldVector<int,dim> >       boundary;

      // open grid file.
      std::cout << "Reading amira mesh\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", gridfile, __LINE__));

      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", "Nodes", "Coordinates", __LINE__));
      }
      // read cells
      std::string cellName;
      if(dim == 2) cellName = "Triangles";
      if(dim == 3) cellName = "Tetrahedra";
      if(dim != 2 && dim != 3){
        delete(am);
        throw std::runtime_error("Sorry!\nOnly 2D and 3D grids are supported.");
      }
      exists = ImplementationDetails::readData(*am,cellName.c_str(),"Nodes",cells);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", cellName, "Nodes", __LINE__));
      }
      // read boundary segments if defined
      exists = ImplementationDetails::readData(*am,"BoundaryTriangles","Nodes",boundary);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", "BoundaryTriangles", "Nodes", __LINE__));
      }

      typedef std::vector<Dune::FieldVector<unsigned char, 1> > TmpIdVector;
      TmpIdVector boundaryID;
      ImplementationDetails::readData(*am,"BoundaryTriangles","Id",boundaryID);

      // close grid file
      delete(am);

      // We use start indices 0. Amira uses 1. Correct indices.
      for (size_t i=0; i<cells.size(); ++i) {
        for (int j=0; j<dim+1; ++j) {
          cells[i][j] -= 1;
          if (cells[i][j]>=vertices.size()){
            std::string message = std::string("In ") + __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__) + ": cell node " + boost::lexical_cast<std::string>(i) + "[" + boost::lexical_cast<std::string>(j) + "]=" + boost::lexical_cast<std::string>(cells[i][j]) + " exceeds vertex numbers!\n";
            throw std::runtime_error(message);
          }
        }
      }

      for (size_t i=0; i<boundary.size(); ++i)
        for (int j=0; j<dim; ++j)
          boundary[i][j] -= 1;


      // Create grid.
      std::cout << std::endl;
      std::cout << "preparing grid creation\n";

      Dune::GridFactory<Grid> factory = IOTools::FactoryGenerator<dim,Grid>::createFactory(initialGridSize);
      std::cout << "inserting vertices...";
      for (size_t i=0; i<vertices.size(); ++i) {
        Dune::FieldVector<double,dim> pos;
        for(int j=0; j<dim; ++j)
          pos[j] = (double)vertices[i][j];
        factory.insertVertex(pos);
      }
      std::cout << "done." << std::endl;

      std::cout << "inserting boundary segments...";
      for(size_t i=0; i<boundary.size(); ++i){
        std::vector<unsigned int> tmp(dim);
        for(size_t j=0; j<dim; ++j){
          assert(boundary[i][j]>=0 && boundary[i][j]<vertices.size());
          tmp[j] = boundary[i][j];
        }
        factory.insertBoundarySegment(tmp);
      }
      std::cout << "done." << std::endl;

      std::cout << "inserting elements...";
      for (size_t i=0; i<cells.size(); ++i) {
        std::vector<unsigned int> tmp(dim+1);
        for (int j=0; j<dim+1; ++j) {
          assert(cells[i][j]>=0 && cells[i][j]<vertices.size());
          tmp[j] = cells[i][j];
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),tmp);
      }
      std::cout << "done." << std::endl;

      boost::timer::cpu_timer timer;
      std::unique_ptr<Grid> grid(factory.createGrid());
      std::cout << "grid creation finished in " << boost::timer::format(timer.elapsed()) << std::endl;

      // get boundary segment indices
      std::cout << "sorting boundary ids...";
      IOTools::BoundarySegmentIndexWrapper<dim, Grid>::readBoundarySegmentIndices(*grid, factory, vertices, boundary, boundaryID, boundaryIds);
      std::cout << "done." << std::endl;

      return grid;
    }

    /// function reading additional data from an Amira mesh file in ASCII format.
    /**
     * \param Scalar scalar type of the data in the Amira mesh file
     * \param FunctionSpaceElement type of FE-Function to store the additional data in
     * \param datafile the name of the Amira mesh file containing additional data
     * \param dataname the name of the data in the Amira mesh file
     * \param componentname the name of the data's components in the Amira mesh file
     * \param data fe-function object for storing the additional data
     */
    template<class Scalar, class FunctionSpaceElement>
    void readData(std::string const& datafile, std::string const& dataname, std::string const& componentname, FunctionSpaceElement &data){
      constexpr int numberOfComponents = FunctionSpaceElement::StorageValueType::dimension;
      // FunctionSpaceElement::StorageValueType in general uses floating point precision.
      // Thus, in order to read integer values or characters, the DataVector is defined manually.
      // When reading the information into the function space element the data is
      // casted to the desired type.
      typedef std::vector<Dune::FieldVector<Scalar, numberOfComponents> > DataVector;
      DataVector dataVector;

      // open data file
      std::cout << std::endl;
      std::cout << "Reading additional data" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      AmiraMesh* am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, FunctionSpaceElement&)", datafile, __LINE__));


      bool exists = ImplementationDetails::readData(*am, dataname.c_str(), componentname.c_str(), dataVector);
      if(!exists){
        delete(am);
        throw std::runtime_error( ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, FunctionSpaceElement&)", dataname, componentname, __LINE__) );
      }

      // close data file
      delete(am);

      std::cout << "number of entries: " << dataVector.size() << std::endl;
      std::cout << "number of components: " << numberOfComponents << std::endl;
      std::cout << "data object: " << data.coefficients().size() << "x" << data.coefficients()[0].size() << std::endl;

      // reading prescribed data into function space element
      for(size_t i=0; i<dataVector.size(); ++i)
        for(size_t j=0; j<numberOfComponents; ++j)
          data.coefficients()[i][j] = (typename FunctionSpaceElement::Scalar) dataVector[i][j];

    }

    /// function reading additional data for each vertex from an Amira mesh file in ASCII format.
    /**
     * \param Scalar scalar type of the data in the Amira mesh file
     * \param numberOfComponents number of components of each entry of the data
     * \param datafile the name of the Amira mesh file containing additional node data
     * \param dataname the name of the data in the Amira mesh file
     * \param componentname the name of the data's components in the Amira mesh file
     * \param data vector containing the data
     */
    template<class Scalar, int numberOfComponents>
    void readData(std::string const& datafile, std::string const& dataname, std::string const& componentname,
        std::vector<Dune::FieldVector<Scalar,numberOfComponents> > &data){
      // open data file
      std::cout << std::endl;
      std::cout << "Reading additional vertex data" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      AmiraMesh* am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, std::vector<Dune::FieldVector>&)", datafile, __LINE__));

      bool exists = ImplementationDetails::readData(*am, dataname.c_str(), componentname.c_str(), data);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, std::vector<Dune::FieldVector>&)", dataname, componentname, __LINE__));
      }

      // close data file
      delete(am);
    }

    /// function reading return a FE-function describing the deformation.
    /**
     * Function reading a deformed and an undeformed geometry from Amira mesh files
     * in ASCII format. These geometries must possess the same number of nodes. For each
     * vertex the difference between the geometries is calculated. This difference is stored
     * into the FE-function at index 0 of the variable set data.
     *
     * \param datafile the name of the Amira mesh file containing the deformed geometry
     * \param gridfile the name of the Amira mesh file containing the undeformed geometry
     * \param data fe-function object for storing the deformation of type VariableSet::VariableSet
     */
    template <class VarSetDesc>
    void readDeformationIntoVariableSetRepresentation(std::string const& gridfile, std::string const& datafile, 
						      typename VarSetDesc::VariableSet& data){
      using namespace boost::fusion;
      int const numberOfComponents = result_of::at_c<typename VarSetDesc::VariableSet::Functions, 0>::type::Components;
      typedef std::vector<Dune::FieldVector<float,numberOfComponents> > DataVector;
      std::vector<Dune::FieldVector<float,numberOfComponents> > vertices;
      DataVector nodeData;

      // open grid file.
      std::cout << "Reading undeformed geometry\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", gridfile, __LINE__));

      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)","Nodes","Coordinates",__LINE__));
      }

      // close grid file
      delete(am);

      // open data file
      std::cout << std::endl;
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", datafile, __LINE__));

      exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", "Nodes", "Coordinates", __LINE__);
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      for(size_t i=0; i<nodeData.size(); ++i)
        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<double> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<numberOfComponents; ++j)
          tmpVec[numberOfComponents*i+j] = (double)nodeData[i][j];

      data.read(tmpVec.begin());
    }

    /// function reading return a fe-function describing the deformation.
    /**
     * Function reading a deformed and a undeformed geometry from Amira mesh files
     * in ASCII format. These geometry must possess the same number of nodes. For each
     * vertex the difference between the geometries is calculated. This difference is read
     * into the fe-function object data.
     *
     * \param datafile the name of the Amira mesh file containing the deformed geometry
     * \param gridfile the name of the Amira mesh file containing the undeformed geometry
     * \param data fe-function object for storing the deformation
     */
    template<class FunctionSpaceElement, class FileScalar=float>
    void readDeformation(std::string const& gridfile, std::string const& datafile, FunctionSpaceElement &data) {
      int const numberOfComponents = FunctionSpaceElement::Components;
      typedef std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > DataVector;
      typedef typename FunctionSpaceElement::Scalar Scalar;
      std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > vertices;
      DataVector nodeData;

      // open grid file.
      std::cout << "Reading undeformed geometry\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", gridfile, __LINE__));

      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close grid file
      delete(am);

      // open data file
      std::cout << std::endl;
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", datafile, __LINE__));

      exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      for(size_t i=0; i<nodeData.size(); ++i)
        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<Scalar> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<FunctionSpaceElement::Components; ++j)
          data.coefficients()[i][j] = (Scalar) nodeData[i][j];
    }


    template<class GridView, class FunctionSpaceElement, class FileScalar=float>
    void readDeformation2(GridView const& gridView, std::string const& datafile, FunctionSpaceElement &data){
      int const numberOfComponents = FunctionSpaceElement::Components;
      typedef std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > DataVector;
      typedef typename FunctionSpaceElement::Scalar Scalar;
      std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > vertices;
      DataVector nodeData;

      // open data file
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      AmiraMesh *am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", datafile, __LINE__));

      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      auto vIter = gridView.template begin<GridView::dimension>();
      auto vend = gridView.template end<GridView::dimension>();
      size_t i=0;
      for(;vIter!=vend;++vIter)
      {
        nodeData[i] -= vIter->geometry().corner(0);
        ++i;
      }
      //      for(size_t i=0; i<nodeData.size(); ++i)
//        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<Scalar> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<FunctionSpaceElement::Components; ++j)
          data.coefficients()[i][j] = (Scalar) nodeData[i][j];
    }


    /// Read boundary indices and save as scalar field in .vtu-file
    /**
     * \param gridfile name of the Amira mesh file
     * \param savefilename name of the output file (optional), default: "boundaryConditions.vtu"
     * \param initialGridSize initial grid size in mb(only for UGGrid)
     */
    template <class Grid>
    void boundaryConditionsToScalarField(std::string const& gridfile, std::string const& savefilename=std::string("boundaryConditions"), int initialGridSize = 0) {

      int const dim = Grid::dimension;
      // read grid
      typedef typename Grid::LeafGridView LeafView;
      std::vector<int> boundaryIndices;
      std::unique_ptr<Grid> grid = readGrid<Grid,float>(gridfile, boundaryIndices);
      GridManager<Grid> gridManager(grid);

      // create variable set
      typedef FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> > Space;
      Space space(gridManager, gridManager.grid().leafView(), 1);
      typedef boost::fusion::vector<Space const*> Spaces;
      Spaces spaces(&space);
      typedef boost::fusion::vector< VariableDescription<0,1,0> > VariableDescriptions;
      typedef VariableSetDescription<Spaces,VariableDescriptions> VarSetDesc;
      std::string name[1] = { "boundary ids" };
      VarSetDesc varSetDesc(spaces, name);
      typename VarSetDesc::VariableSet indices(varSetDesc);

      // read boundary indices
      std::vector<Dune::FieldVector<int,dim> > boundaryVertices;

      std::string const comp_name("BoundaryTriangles");
      std::string const node_name("Nodes");

      ImplementationDetails::readData<int,dim>(gridfile, comp_name, node_name, boundaryVertices);
      std::vector<double> tmp(gridManager.grid().size(dim),0);
      for(int i=0; i<boundaryIndices.size(); ++i){
        tmp[boundaryVertices[i][0]-1] = boundaryIndices[i];
        tmp[boundaryVertices[i][1]-1] = boundaryIndices[i];
        tmp[boundaryVertices[i][2]-1] = boundaryIndices[i];
      }
      indices.read(tmp.begin());

      // write file
      IoOptions options;
      options.outputType = IoOptions::ascii;
      writeVTKFile(gridManager.grid().leafView(), varSetDesc, indices, savefilename,options);
    }
  }
} // namespace Kaskade
#endif
