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


#ifndef VTKWRITER_HH
#define VTKWRITER_HH

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory> 
#include <sstream>
#include <type_traits>
#include <vector>

#include <boost/fusion/include/zip.hpp>
#include <boost/range/iterator_range.hpp>

#include <dune/grid/common/grid.hh>
#include <dune/geometry/referenceelements.hh>

#include "io/iobase.hh"
#include "utilities/detailed_exception.hh"
#undef min





namespace Kaskade
{
  /**
   * \cond internals
   */

  namespace FunctionSpace_Detail {
    // forward declaration
    template <class> class ToDomainRepresentation;
  }

  namespace VTKWriterDetail
  {
    // returns the textual VTK notation of this system's endianness.
    std::string vtkEndianness();
    
    enum VTKGeometryType
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
    
    //! mapping from GeometryType to VTKGeometryType
    VTKGeometryType vtkType(Dune::GeometryType const& t, int order);
    
    /**
     * \brief Returns the subentity number and codimension for any VTK point in a cell.
     */
    std::pair<int,int> vtkPointToDune(Dune::GeometryType const& type, int vtkPoint);
    

    // forward declaration
    class Base64Writer;
    
    //! base class for VTK data array writers
    // T  is a scalar type (float,double,int,...)
    template <class T>
    class VTKDataArrayWriter
    {
    public:
    /** 
     * @brief Constructor.
     * 
     * @param s           The stream to write to
     * @param outputType  ASCII or binary
     * @param name        The name of the vtk array
     * \param ncomps      number of vectorial components (1 or 3)
     * @param entries     the number of data entries to write
     */
      VTKDataArrayWriter(std::ostream &s, IoOptions::OutputType outputType, std::string const& name, int ncomps,
                         size_t entries, int indentCount, int precision);
      
      //! write one data element
      template <int m>
      void write (Dune::FieldVector<T,m> const&);
      
      void write(T);
      
      //!  destructor
      ~VTKDataArrayWriter();
      
    private:
      std::ostream& s;
      IoOptions::OutputType outputType;
      int ncomps;
      int counter;
      int numPerLine;
      int indentCount;
      int precision;
      std::unique_ptr<Base64Writer> base64;
    };
 
    // ---------------------------------------------------------------------------

    //! write indentation to stream
    void indent(std::ostream& s, int indentCount);
    
    /// This gathers some information about the grid and supports the output of the grid
    template <class GridView> 
    struct VTKGridInfo
    {
      using Cell         = typename GridView::template Codim<0>::Entity;
      using IndexSet     = typename GridView::IndexSet;
      using CellIterator = typename GridView::template Codim<0>::Iterator;
      
      static int constexpr dim = GridView::dimension;
      
      
      VTKGridInfo(GridView const& gridView_, IoOptions const& options_)
      : gridView(gridView_),
        cells(gridView.template begin<0>(),gridView.template end<0>()),
        ncells(gridView.size(0)),
        options(options_)
      {        
        // For counting the corners (i.e. all corners of all cells) we step through the 
        // cells and count their respective corners (and edges for order two). This is 
        // general enough to treat mixed grids containing different cell types.
        ncorners = 0;
        for (auto const& cell: cells) 
        {
          ncorners += cell.subEntities(dim);
          if (options.order==2)
            ncorners += cell.subEntities(dim-1);
        }
        
        // In conforming mode, the number of points is just the number of vertices 
        // (plus edges for oder two). Otherwise, points and corners coincide.
        if (options.dataMode==IoOptions::conforming)
        {
          npoints = gridView.size(dim);
          if (options.order==2)
            npoints += gridView.size(dim-1);
        }
        else
          npoints = ncorners;
      }
      
      // returns entity center coordinates
      template <class Cell>
      auto center(Cell const& cell, int subindex, int codim) const
      {
        return Dune::ReferenceElements<typename GridView::ctype,dim>::general(cell.type()).position(subindex,codim);
      }
      
      // returns a globally unique index for vertices and edges in conforming mode. 
      // Vertices come first, then edges.
      template <class Cell>
      int conformingIndex(Cell const& cell, int subindex, int codim) const
      {
        // due to the sorting vertices first then edges there is no need to switch between order 1 and 2
        if (codim==dim)
          return gridView.indexSet().subIndex(cell,subindex,codim);
        else
          return gridView.size(dim) + gridView.indexSet().subIndex(cell,subindex,codim);
      }
      
      // Writes the point coordinates to the VTK XML stream
      void writeGridPoints (int indentCount, std::ostream& s) const
      {
        VTKDataArrayWriter<typename GridView::ctype> p(s, options.outputType, "Coordinates", 3, npoints, indentCount, options.precision); 
        
        if (options.dataMode==IoOptions::conforming)
        {
          std::vector<Dune::FieldVector<typename GridView::ctype,GridView::dimensionworld>> xs(npoints);
          for (auto const& cell: cells)
          {
            for (int corner=0; corner<cell.subEntities(dim); ++corner)   // global coordinates of all corners
              xs[conformingIndex(cell,corner,dim)] = cell.geometry().corner(corner); 
            
            if (options.order == 2)
              for (int edge=0; edge<cell.subEntities(dim-1); ++edge)     // global coordinates of all edges
                xs[conformingIndex(cell,edge,dim-1)] = cell.geometry().global(center(cell,edge,dim-1));
          }
          
          for (auto const& x: xs)
            p.write(x);
        }
        else // nonconforming
        {
          for (auto const& cell: cells)
          {
            for (int corner=0; corner<cell.subEntities(dim); ++corner)
              p.write(cell.geometry().corner(corner));
            
            if (options.order == 2)
              for (int edge=0; edge<cell.subEntities(dim-1); ++edge)
                p.write(cell.geometry().global(center(cell,edge,dim-1)));
          }
        }
      }
      
      void writeGridCells (int indentCount, std::ostream& s) const
      {
        // connectivity
        auto p1 = std::make_unique<VTKDataArrayWriter<int>>(s, options.outputType, "connectivity", 1, ncorners, indentCount, options.precision);
        
        int offset = 0;
        std::vector<int> offsets; offsets.reserve(ncells);
        for (auto const& cell: cells)                                             // step through all cells
        {
          int const numCorners = cell.subEntities(dim);
          int const numEdges = cell.subEntities(dim-1);
          int const numPoints = numCorners + (options.order==2? numEdges: 0);
          
          for (int vtkPoint=0; vtkPoint<numPoints; ++vtkPoint)                   // visit each point in that cell
          {
            auto subentity = vtkPointToDune(cell.type(),vtkPoint);               // and find its Dune coordinates
            p1->write( options.dataMode==IoOptions::conforming?
                          conformingIndex(cell,subentity.first,subentity.second):
                          offset + (subentity.second==dim? 0: numCorners) + subentity.first );
          }
          offset += numPoints;                                                   // in nonconforming mode, the points are stored blockwise,
          offsets.push_back(offset);
        }                                                                        // cell by cell, just advance the offset
        p1.reset(); // write end tag on destruction
        
        // write offsets
        auto p2 = std::make_unique<VTKDataArrayWriter<int>>(s, options.outputType, "offsets", 1, ncells, indentCount, options.precision);
        for (int off: offsets)
          p2->write(off);
        p2.reset(); // write end tag on destruction
        
        // types only if dimension greater than 1 - there is only one 1D cell type: line
        if (dim>1)
        {
          VTKDataArrayWriter<unsigned char> p(s, options.outputType, "types", 1, ncells, indentCount, options.precision);
          for (auto const& cell: cells)
            p.write(static_cast<unsigned char>(vtkType(cell.type(),options.order)));
        }
      }
    
    
      GridView const& gridView;
      boost::iterator_range<CellIterator> cells;
      size_t ncorners;
      size_t npoints;
      size_t ncells;
      IoOptions const& options;
    };
    
    
    // ----------------------------------------------------------------------------------------

    
    template <class GridView, class Function, class Names>
    void writeCellData (VTKGridInfo<GridView> const& gridInfo, Function const& pair, Names const& names, int indentCount, std::ostream& s)
    {
      auto const& f = boost::fusion::at_c<0>(pair);
      
      // Check whether this is best represented as cell data: only if it is (piecewise) constant.
      if (f.space().mapper().maxOrder() > 0)
        return;
      
      auto varDesc = boost::fusion::at_c<1>(pair);
      VTKDataArrayWriter<double> p(s, gridInfo.options.outputType, names[varDesc.id], varDesc.m, gridInfo.ncells, indentCount, gridInfo.options.precision);
      for (auto const& cell: gridInfo.cells)
        p.write(f.value(cell,gridInfo.center(cell,0,0)));
    }
    
    
    
    template <class GridView, class Function, class Names>
    void writeVertexData (VTKGridInfo<GridView> const& gridInfo, Function const& pair, Names const& names, int indentCount, std::ostream& s)
    {
      auto const& f = boost::fusion::at_c<0>(pair);
      using ValueType = typename std::remove_reference_t<typename boost::fusion::result_of::value_at_c<Function,0>::type>::ValueType;
      constexpr int dim = GridView::dimension;
     
     // Check whether this is best represented as vertex data: only if it is not (piecewise) constant.
      if (f.space().mapper().maxOrder() == 0)
        return;

      auto varDesc = boost::fusion::at_c<1>(pair);
      VTKDataArrayWriter<double> p(s, gridInfo.options.outputType, names[varDesc.id], varDesc.m, gridInfo.npoints, indentCount, gridInfo.options.precision);
      
      if (gridInfo.options.dataMode==Kaskade::IoOptions::conforming)
      {
        std::vector<ValueType> fs(gridInfo.npoints,ValueType(0.0));   // function values, in conforming mode average over multiple points
        std::vector<short> count(gridInfo.npoints,0);                  // count how many points contribute to an entity

        // TODO: Currently boundary FE functions are zero on cells which only touch the boundary by one vertex (or edge in 3D). This will lead to unwanted results when values are averaged.
        // Therefore, we use the following two lines to transform boundary FE functions to continuous FE functions over the whole domain, then averaging works as expected.
        // Remove this when not necessary any longer. (And also delete ToDomainRepresetation from functionspace.hh)
        FunctionSpace_Detail::ToDomainRepresentation<std::remove_const_t<std::remove_reference_t<typename boost::fusion::result_of::value_at_c<Function,0>::type>>> tdr(f);
        auto const& df = tdr.get();

        for (auto const& cell: gridInfo.cells)
        {
          int const numPoints = cell.subEntities(dim) + (gridInfo.options.order==2? cell.subEntities(dim-1): 0);
          for (int point=0; point<numPoints; ++point)
          {
            auto subentity = vtkPointToDune(cell.type(),point);
            auto xi = gridInfo.center(cell,subentity.first,subentity.second);          // local coordinate of point in current cell
            int idx = gridInfo.conformingIndex(cell,subentity.first,subentity.second); // global index of point
            ++count[idx];
            fs[idx] += df.value(cell,xi);
          }          
        }
        
        for (int i=0; i<gridInfo.npoints; ++i)          // now write out the values in vertex order
          p.write(fs[i]/static_cast<double>(count[i])); // and average them on writing
      }
      else // nonconforming
      {
        for (auto const& cell: gridInfo.cells)
        {
          // First write corner values
          for (int corner=0; corner<cell.subEntities(dim); ++corner)
            p.write(f.value(cell,gridInfo.center(cell,corner,dim)));

          // If required, write edge midpoint values as well.
          if (gridInfo.options.order==2)
            for (int edge=0; edge<cell.subEntities(dim-1); ++edge)
              p.write(f.value(cell,gridInfo.center(cell,edge,dim-1)));
        }
      }      
    }
    
    // The following piece of code defines the "active" scalar and vectorial variables (those used to immediately color the output in Paraview). 
    // According to VTK file formats, there need not be any active scalar or vectorial variable. Taking simply the first one as is done here makes little sense. 
    // TODO: make this configurable via IoOptions.
    template <class Functions, class Names>
    std::string guessActiveFields(Functions const& functions, Names const& names, bool cellData) 
    {
      std::string scalar = "";
      std::string vector = "";
      
      boost::fusion::for_each(functions,[&](auto const& pair) 
      {
        auto const& f = boost::fusion::at_c<0>(pair);
        auto varDesc = boost::fusion::at_c<1>(pair);
      
        // Check whether this is best represented as cell or point data.
        if ( (f.space().mapper().maxOrder()==0 && cellData)
          || (f.space().mapper().maxOrder()>0 && !cellData) )
        {
          if (f.components>1 && vector.empty())
            vector = " Vectors=\"" + names[varDesc.id] + "\"";
          if (f.components==1 && scalar.empty())
            scalar = " Scalars=\"" + names[varDesc.id] + "\"" ;
        }
      });
      return scalar+vector;
    }
  
  }
  /**
   * \endcond
   */
  
  // ---------------------------------------------------------------------------------------------
  
  /**
   * \ingroup IO
   * \brief Writes a set of finite element functions in VTK XML format to a stream.
   * \param vars a set of finite element functions to write
   * \param options specifies how the data shall be written
   * \param s the output stream
   */
  template <class VariableSet>
  void writeVTK(VariableSet const& vars, IoOptions options, std::ostream& s)
  {
    using namespace VTKWriterDetail;
    options.order = std::max(1,std::min(options.order,2));
    
    // If data mode "inferred" is requested, try to select either conforming or nonconforming
    // based on rules.
        
    
    
    int indentCount = 0;
    constexpr int dim = VariableSet::Descriptions::Grid::dimension;
    
    VTKGridInfo<typename VariableSet::Descriptions::GridView> gridInfo(vars.descriptions.gridView,options);
    
    // xml header
    s << "<?xml version=\"1.0\"?>" << std::endl;
    
    // VTKFile
    std::string const gridTag = dim>1? "UnstructuredGrid": "PolyData";
    s << "<VTKFile type=\"" << gridTag << "\" version=\"0.1\" byte_order=\"" << vtkEndianness() << "\">" << std::endl;
    ++indentCount;
    
    // UnstructuredGrid 
    indent(s,indentCount);
    s << "<" << gridTag << ">" << std::endl;
    ++indentCount;
    
    // Piece
    indent(s,indentCount);
    if (dim>1)
      s << "<Piece NumberOfPoints=\"" << gridInfo.npoints << "\" NumberOfCells=\"" << gridInfo.ncells << "\">" << std::endl;
    else
      s << "<Piece NumberOfPoints=\"" << gridInfo.npoints << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << gridInfo.ncells << "\" NumberOfPolys=\"0\">" << std::endl;
    ++indentCount;
    
    // Write grid points
    indent(s,indentCount); s << "<Points>" << std::endl;    
    gridInfo.writeGridPoints(indentCount+1,s);
    indent(s,indentCount); s << "</Points>" << std::endl;
   
    // Write grid cells
    std::string const cellTagName = dim>1? "Cells": "Lines";        
    indent(s,indentCount);  s << '<' << cellTagName << '>' << std::endl;
    gridInfo.writeGridCells(indentCount+1,s);
    indent(s,indentCount); s << "</" << cellTagName << '>' << std::endl;
    
    // associate the FE functions with their variable descriptions
    auto varDesc = typename VariableSet::Descriptions::Variables();
    auto functions = boost::fusion::zip(vars.data,varDesc);
    
    indent(s,indentCount); s << "<CellData" << guessActiveFields(functions,vars.descriptions.names,true) << ">" << std::endl;
    boost::fusion::for_each(functions,[&](auto const& f) { writeCellData(gridInfo,f,vars.descriptions.names,indentCount+1,s); });
    indent(s,indentCount); s << "</CellData>" << std::endl;
    
    indent(s,indentCount); s << "<PointData" << guessActiveFields(functions,vars.descriptions.names,false) << ">" << std::endl;
    boost::fusion::for_each(functions,[&](auto const& f) { writeVertexData(gridInfo,f,vars.descriptions.names,indentCount+1,s); });
    indent(s,indentCount); s << "</PointData>" << std::endl;
    
    // /Piece
    --indentCount;
    indent(s,indentCount); s << "</Piece>" << std::endl;
    
    // /UnstructuredGrid
    --indentCount;
    indent(s,indentCount); s << "</" << gridTag << ">" << std::endl;
    
    // /VTKFile
    s << "</VTKFile>" << std::endl;
  }
  
  /**
   * \ingroup IO
   * \brief Writes a set of finite element functions in VTK XML format to a stream.
   * \param vars a set of finite element functions
   * \param name file name (without trailing .vtu)
   * 
   * Note that Paraview can create a time line implicitly from numbered files. use paddedString for easy creation 
   * of numbered file names.
   */
  template <class VariableSet>
  void writeVTK(VariableSet const& vars, std::string name, IoOptions const& options)
  {
    // generate filename for process data
    name += (VariableSet::Descriptions::GridView::dimension>1? ".vtu" : ".vtp");
    
    // write process data
    std::ofstream file(name);
    if (!file)
      throw Kaskade::FileIOException("Opening  failed.\n",name,__FILE__,__LINE__);
    writeVTK(vars,options,file);
    if (!file)
      throw Kaskade::FileIOException("Writing  failed.\n",name,__FILE__,__LINE__);
    file.close();
  }
 
}

// --------------------------------------------------------------------------------------------



#endif
