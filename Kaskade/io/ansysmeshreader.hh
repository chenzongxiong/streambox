/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ANSYS_MESH_READER_HH
#define ANSYS_MESH_READER_HH

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <limits>

#include <boost/lexical_cast.hpp>

#include <dune/common/fvector.hh>
#include "dune/grid/config.h"


#if HAVE_UG
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

#if HAVE_ALBERTA
#include <dune/grid/albertagrid/gridfactory.hh>
#endif

#if ENABLE_ALUGRID
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif

#include "fem/lagrangespace.hh"
#include "io/ioTools.hh"
// #include "io/vtk.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"

#ifndef Cygwin
using std::to_string;
#else
std::string to_string(size_t n)
{
  std::stringstream s; s << n; return s.str();
}
#endif

namespace Kaskade
{
  /**
   * \ingroup gridInput
   * \brief Reader for Ansys *.inp grid files.
   */
  class AnsysMeshReader
  {
  public:
    /**
     * \brief Constructor.
     * \param scale Scale geometry by the given factor. Use this to convert the geometry to SI units.
     * \param verbose Output level to standard out (0=no output (default), 1=status messages, 2=debug messages).
     */
    explicit AnsysMeshReader(double scale_=1.0, int verbose_=0) : scale(scale_), verbose(verbose_)
    {}

    /**
     * \brief Reads an Ansys file and creates a corresponding grid.
     * \param fileName
     * \param initialGridSize
     */
    template <class Grid, class Scalar>
    std::unique_ptr<Grid> readGrid(std::string const& fileName, int initialGridSize = 1000)
    {
      constexpr int dim = Grid::dimension;
      std::ifstream file(fileName);
      currentLine = 0;

      if(!(bool)file) throw IOException(std::string("Could not read file ") + fileName, __FILE__, __LINE__);
      else if(verbose>=1) std::cout << "Reading " << fileName << std::endl;

      BoundingBox<Scalar,dim> boundingBox;
      std::vector<Dune::FieldVector<Scalar,dim> > vertices;
      std::vector<std::vector<unsigned int> > cubes;
      std::vector<std::pair<unsigned int,unsigned int> > offsets(1,std::make_pair(0,0));

      auto NOT_FOUND = std::string::npos;
      std::string buffer;

      // First line should contain "*HEADING". Check if this is the case.
      getline(file,buffer);
      if(buffer.find("*HEADING") == NOT_FOUND) 
        throw IOException(std::string("Error reading file ") + fileName + "at line " + to_string(currentLine),__FILE__,__LINE__);
      getline(file,buffer);
      short elementSetId = 0;

      while(file.good())
      {

        if(buffer.empty())
        {
          getline(file,buffer);
          continue;
        }

        if(verbose >= 1) std::cout << buffer << std::endl;
        if(buffer[0] != '*') throw IOException(std::string("Expected comment at line " + to_string(currentLine) + "."),__FILE__,__LINE__);
        if(buffer.find("**") != NOT_FOUND)
        {
          if(buffer.find("Materialdaten") != NOT_FOUND) break;
          if(verbose >= 1) std::cout << "End of section. Starting new one." << std::endl;
          getline(file,buffer);
          continue;
        }

        // read vertices
        if(buffer.find("*NODE") != NOT_FOUND)
        {
          if(verbose >= 1) std::cout << "Reading vertices...";
          readVertices<double>(file, buffer, vertices, offsets);
          if(verbose >= 1) std::cout << "done: " << vertices.size() << " entries.\n" << std::endl;
          continue;
        }

        // read elements
        if(buffer.find("*ELEMENT") != NOT_FOUND)
        {
          if(buffer.find("TYPE=C3D8") != NOT_FOUND) // hexahedral element section found
          {
            std::string::size_type pos = buffer.find("ELSET=");
            if(pos == NOT_FOUND) throw IOException(std::string("Could not determine element set."),__FILE__,__LINE__);


            if (verbose >= 1) std::cout << "Reading elements..." << std::flush;
            readElements(file, buffer, cubes, offsets);
            if (verbose >= 1) std::cout << "done: " << cubes.size() << " entries.\n" << std::endl;
            continue;
          }

          if(buffer.find("TYPE=C3D4") != NOT_FOUND) // tetrahedral element section found
          {
            std::string::size_type pos = buffer.find("ELSET=");
            if(pos == NOT_FOUND) throw IOException(std::string("Could not determine element set."),__FILE__,__LINE__);

            if (verbose >= 1) std::cout << "Reading tetrahedra..." << std::flush;
            readTetrahedra(file, buffer, cubes, offsets, elementSetId++);
            if (verbose >= 1) std::cout << "done: " << cubes.size() << " entries.\n" << std::endl;
            continue;
          }

          throw IOException(std::string("For Ansys meshes currently only supported for cubic and tetrahedral geometries."),__FILE__,__LINE__);
        }
      }

      // create factory
      Dune::GridFactory<Grid> factory = IOTools::FactoryGenerator<dim,Grid>::createFactory(initialGridSize);
      /****************************************/
      if (verbose >= 1) std::cout << "inserting vertices...";
      for(auto const& vertex : vertices) factory.insertVertex(vertex);
      if (verbose >= 1) std::cout << "done." << std::endl;

      /****************************************/
      if (verbose >= 1) std::cout << "inserting elements...";

      for (auto const& vertexIds: cubes) {
        // include each vertex into the bounding box
        if (verbose >= 2)
          for (auto vertexId: vertexIds)
            boundingBox.update(vertices[vertexId]);

        // create the cell. TODO: is this restricted to simplicial cells? what about hexahedral ones?
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertexIds);
      }
      if (verbose >= 1) std::cout << "done." << std::endl;

      if (verbose >= 2) boundingBox.print();

      if (verbose >= 1) std::cout << "grid creation finished\n\n";

      return std::unique_ptr<Grid>(factory.createGrid());
    }


    /**
     * \brief Extracts the element set id of each cell and stores it in the coefficient of a discontinuous P0 FE function.
     */
    template <class FSElement>
    void readElementSetIds(FSElement& element) const
    {
      assert(element.coefficients().N() == elementSetIds.size());
      for(size_t i=0; i<elementSetIds.size(); ++i)
        element.coefficients()[i][0] = static_cast<typename FSElement::Scalar>(elementSetIds[i]);
    }


  private:
      template <class Scalar, int dim>
      struct BoundingBox
      {
        BoundingBox() : max(-1.0*std::numeric_limits<Scalar>::max()), min(std::numeric_limits<Scalar>::max())
        {}

        explicit BoundingBox(Dune::FieldVector<Scalar,dim> const& v) : max(v), min(v)
        {}

        void update(Dune::FieldVector<Scalar,dim> const& v)
        {
//           if(v[1] > 0) return; // ???
          for(int i=0; i<dim; ++i)
          {
            min[i] = std::min(min[i],v[i]);
            max[i] = std::max(max[i],v[i]);
          }
        }

        bool isInside(Dune::FieldVector<Scalar,dim> const& v) const
        {
          for(int i=0; i<dim; ++i) 
            if ((min[i] > v[i]) || (max[i] < v[i])) // outside
              return false;
          return true;
        }

        void print() const
        {
          std::cout << "Bounding Box" << std::endl;
          for(int i=0; i<dim; ++i) std::cout << i << ": [" << min[i] << ", " << max[i] << "]" << std::endl;
          std::cout << std::endl;
        }

        Dune::FieldVector<Scalar,dim> max, min;
      };

      static const char delimiter = ',';

      /// Correction of local vertex order such that the requirements of the Dune::GridFactory are satisfied.
      /**
       * The local vertex order will be changed for cubic elements only.
       * For pyramids and simplices it is assumed, that the ordering is already correct.
       * (See method AnsysMeshReader::readGrid).
       */
      void sortElements(std::vector<std::vector<unsigned int> >& elements) const;

      /// Adjust vertex indices stored in element using the information of offset.
      /**
       * \param element global vertex indices of the element
       * \param offsets vector with offset information. The offset is accessed via offsets[i].second.
       *                offsets[i].first is the first index from which on the offset has to be considered
       */
      void adjustElement(std::vector<unsigned int>& element, std::vector<std::pair<unsigned int,unsigned int> > const& offsets);

      inline void removeLeadingWhiteSpaces(std::string& buffer)
      {
        buffer.erase(0,buffer.find_first_not_of(' '));
      }

      template <class Type>
      inline Type readEntry(std::ifstream& file, std::string& buffer)
      {
        getline(file,buffer,delimiter);
        removeLeadingWhiteSpaces(buffer);
        return boost::lexical_cast<Type>(buffer);
      }

      template <class Scalar, int dim>
      void readVertices(std::ifstream& file, std::string& buffer, std::vector<Dune::FieldVector<Scalar,dim> >& vertices, 
                        std::vector<std::pair<unsigned int,unsigned int> >& offsets)
      {
        Dune::FieldVector<Scalar,dim> vertex;

        unsigned int currentIndex = 1, givenIndex;
        char tmp;

        getline(file,buffer,delimiter);
        removeLeadingWhiteSpaces(buffer);

        while(file.good())
        {
          givenIndex = boost::lexical_cast<unsigned int>(buffer);

          // throw if no ordering exists in file
          if(currentIndex+offsets.back().second > givenIndex)
          {
            std::cerr << vertices.size() << ", " << currentIndex << ", " << offsets.back().second << ", " << givenIndex << std::endl;
            throw IOException(std::string("currentIndex > givenIndex in line ") + to_string(currentLine),__FILE__,__LINE__);
          }
          // store offsets
          if(currentIndex + offsets.back().second < givenIndex) 
            offsets.push_back(std::make_pair(givenIndex-1, givenIndex - currentIndex));
          ++currentIndex;

          // read data
          vertex[0] = readEntry<Scalar>(file,buffer);
          vertex[1] = readEntry<Scalar>(file,buffer);
          file >> vertex[2];
          
          // scale the position as requested
          vertex *= scale;

          vertices.push_back(vertex);

          getline(file,buffer);
          file >> tmp;
          file.unget();
          if(tmp=='*')
          {
            getline(file,buffer);
            return;
          }
          getline(file,buffer,delimiter);
          removeLeadingWhiteSpaces(buffer);
        }
      }

      void readElements(std::ifstream& file, std::string& buffer, std::vector<std::vector<unsigned int> >& elements, 
                        std::vector<std::pair<unsigned int,unsigned int> > const& offsets);
      

      void readTetrahedra(std::ifstream& file, std::string& buffer, std::vector<std::vector<unsigned int> >& elements, 
                          const std::vector<std::pair<unsigned int,unsigned int> >& offsets, short elementSetId);

      template <class Scalar, int dim>
      inline bool equal(Dune::FieldVector<Scalar,dim> const& v0, Dune::FieldVector<Scalar,dim> const& v1)
      {
        return (v0-v1).two_norm2() <= power<2>(100*std::numeric_limits<Scalar>::epsilon()) * (v0.two_norm2() + v1.two_norm2());
      }

      template <class Scalar, int dim>
      inline unsigned int getIndex(Dune::FieldVector<Scalar,dim> const& vertex, std::vector<Dune::FieldVector<Scalar,dim> > const& vertices)
      {
        for(size_t i=0; i<vertices.size(); ++i) 
          if(equal(vertex,vertices[i])) 
            return i;
        return vertices.size();
      }


      template <class Scalar, int dim>
      void getElementsVertices(std::vector<Dune::FieldVector<Scalar,dim> > const& vertices, std::vector<Dune::FieldVector<Scalar,dim> >& elementVertices,
                               std::vector<unsigned int> const& element, std::vector<unsigned int>& newElement, unsigned int index)
      {
        for(size_t i=0; i<element.size(); ++i)
        {
          unsigned int index = getIndex(vertices[element[i]], elementVertices);
          if(index==elementVertices.size())
          {
            newElement[i] = elementVertices.size();
            elementVertices.push_back(vertices[element[i]]);
          }
          else
            newElement[i] = index;
        }
      }
      
      // Read a line from the file, advancing the currentLine counter appropriately.
      void getline(std::ifstream& file, std::string& line, char const delim='\n');

      std::vector<short> elementSetIds;
      double scale;
      int verbose;
      size_t currentLine; // storing the line number on reading
  };
} // namespace Kaskade
#endif
