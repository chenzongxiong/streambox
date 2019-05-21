/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <memory>
#include <type_traits>

#include <boost/fusion/include/as_map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/fusion/include/map.hpp>


#include "dune/grid/config.h"

#include "io/check_endianness.hh"
#include "io/vtkwriter.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  namespace VTKWriterDetail
  {
    std::string vtkEndianness()
    {
      if(bigEndian())
        return "BigEndian";
      else if(littleEndian())
        return "LittleEndian";
      else 
        throw std::runtime_error("Could not determine endianness. No VTK-file will be created.");
      return "";
    }
    
    // ----------------------------------------------------------
    
    void indent(std::ostream& s, int indentCount)
    {
      for (int i=0; i<indentCount; i++)
        s << "  ";
    }
    
    // ----------------------------------------------------------
    
    // map scalar types to their VTK names
    auto vtkTypeNames = boost::fusion::as_map(boost::fusion::make_map<char,unsigned char, short, unsigned short, int, unsigned int, float, double>
                                                                     ("Int8","UInt8",     "Int16","UInt16",      "Int32","UInt32",  "Float32","Float64"));

    // type to use in ASCII output. This is needed to prevent Int8 and UInt8 to be printed as characters.
    template <class T>
    using PrintType = typename std::conditional<std::is_same<T,char>::value || std::is_same<T,unsigned char>::value, int, T>::type;
        

    //! mapping from GeometryType to VTKGeometryType
    VTKGeometryType vtkType(Dune::GeometryType const& t, int order)
    {
      if (order == 2) 
      {
        if (t.isLine())          return vtkquadLine;
        if (t.isTriangle())      return vtkquadTriangle;
        if (t.isQuadrilateral()) return vtkquadQuadrilateral;
        if (t.isTetrahedron())   return vtkquadTetrahedron;
        if (t.isHexahedron())    return vtkquadHexahedron;
      }
      else // order == 1
      {
        if (t.isLine())          return vtkLine;
        if (t.isTriangle())      return vtkTriangle;
        if (t.isQuadrilateral()) return vtkQuadrilateral;
        if (t.isTetrahedron())   return vtkTetrahedron;
        if (t.isPyramid())       return vtkPyramid;
        if (t.isPrism())         return vtkPrism;
        if (t.isHexahedron())    return vtkHexahedron;
      }
      throw DetailedException("VTKWriter: unsupported GeometryType",__FILE__,__LINE__);
      return vtkLine; // never get here
    }
    
    std::pair<int,int> vtkPointToDune(Dune::GeometryType const& type, int vtkPoint)
    {
      static std::pair<int,int> triangle[] = { {0,2}, {1,2}, {2,2},        {0,1}, {2,1}, {1,1} };
      static std::pair<int,int> quad[]     = { {0,2}, {1,2}, {3,2}, {2,2}, {2,1}, {1,1}, {3,1}, {0,1} };
      static std::pair<int,int> tetra[]    = { {0,3}, {1,3}, {2,3}, {3,3}, {0,2}, {2,2}, {1,2}, {3,2}, {4,2}, {5,2} };
      static std::pair<int,int> hexa[]     = { {0,3}, {1,3}, {3,3}, {2,3}, {4,3}, {5,3}, {7,3}, {6,3},
                                               {6,2}, {5,2}, {7,2}, {4,2}, {10,2}, {9,2}, {11,2}, {8,2}, {0,2}, {1,2}, {3,2}, {2,2} };
      if (type.isTriangle())       return triangle[vtkPoint];
      if (type.isQuadrilateral())  return quad[vtkPoint];
      if (type.isTetrahedron())    return tetra[vtkPoint];
      if (type.isHexahedron())     return hexa[vtkPoint];
      return std::make_pair(-1,-1); // never get here
    }
  
    // ----------------------------------------------------------
    
    namespace 
    {
      constexpr char const base64chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    }
    
    /**
     * \brief A class for encoding binary data in base64.
     * 
     * On destruction of the writer (or by calling flush()), the terminal padding characters are written.
     * 
     * For base64 encoding see https://en.wikipedia.org/wiki/Base64.
     */
    class Base64Writer
    {
    public:
      Base64Writer(std::ostream& out_): out(out_) {}
      
      /**
       * \brief Encode binary data.
       */
      template <class T>
      Base64Writer& operator << (T const& t)
      {
        // copy data to character buffer
        unsigned char const* p = reinterpret_cast<unsigned char const*>(&t);
        std::copy(p,p+sizeof(T),std::back_inserter(buffer));
        encode();
        return *this;
      }
      
      /**
       * \brief Writes terminal padding characters.
       * 
       * After calling this method, a new sequence of data can be encoded.
       */
      void flush()
      {
        if (!buffer.empty())
        {
          assert(buffer.size() < 3);
          // buffer may contain one or two bytes. Fill up the remaining bits with zero and the remaining characters with '='.
          out << base64chars[(buffer[0] & 0xfc) >> 2];
          if (buffer.size()==1)
            out << base64chars[(buffer[0] & 0x03) << 4] << "==";
          else // buffer.size()==2
            out << base64chars[((buffer[0] & 0x03) << 4) + ((buffer[1] & 0xf0) >> 4)]
                << base64chars[(buffer[1] & 0x0f) << 2] << '=';
          buffer.clear();
        }
      }
      
      ~Base64Writer()
      {
        flush();
      }
      
    private:
      std::ostream& out;
      std::vector<unsigned char> buffer;
      
      void encode()
      {
        // encode the blocks three bytes of the buffer into four characters written to the stream. The encoded bytes are removed
        // from the buffer.
        auto bufp = begin(buffer);
        for (int i=0; i<buffer.size()/3; ++i, bufp+=3)
          out << base64chars[(bufp[0] & 0xfc) >> 2]
              << base64chars[((bufp[0] & 0x03) << 4) + ((bufp[1] & 0xf0) >> 4)]
              << base64chars[((bufp[1] & 0x0f) << 2) + ((bufp[2] & 0xc0) >> 6)]
              << base64chars[bufp[2] & 0x3f];
        buffer.erase(begin(buffer),bufp);
      }
      
    };
    
    //----------------------------------------------------------------------------------------------------
    
    // a streaming writer for data array tags, uses ASCII inline format
    template<class T>
    VTKDataArrayWriter<T>::VTKDataArrayWriter(std::ostream& s_, IoOptions::OutputType outputType_, std::string const& name, int ncomps_,
                                              size_t entries, int indentCount_, int precision_)
    : s(s_), outputType(outputType_), ncomps(ncomps_), counter(0), numPerLine(12), indentCount(indentCount_), precision(precision_)
    {
      if (ncomps==2) // 2D vectorial data shall be written as 3D with z component 0
        ncomps = 3;
      
      indent(s,indentCount); s << "<DataArray type=\"" << boost::fusion::at_key<T>(vtkTypeNames) << "\" Name=\"" << name << "\" ";
      s << "NumberOfComponents=\"" << ncomps << "\" ";
      s << "format=\"" << (outputType==IoOptions::ascii? "ascii": "binary") << "\">" << std::endl; 
      
      if (outputType==IoOptions::ascii)
        indent(s,indentCount);
      else 
      {
        base64 = std::make_unique<Base64Writer>(s);
        // VTK binary base64 is structured as |header|data| where header is a UInt32, base64-encoded separately, and 
        // gives the number of data bytes before base64 encoding. data follows immediately after the header. No white
        // space is allowed within header or data, as VTK reads the data apparently without parsing as a plain block 
        // of data. The header shall start immediately after a newline following the closing '>' of the DataArray tag.
        // see http://public.kitware.com/pipermail/paraview/2005-April/001391.html
        uint32_t bytes = sizeof(T)*ncomps*entries;
        
        *base64 << bytes;
        base64->flush();
      }
    }
    
    template <class T>
    void VTKDataArrayWriter<T>::write(T data)
    {
      if (outputType==IoOptions::ascii)
      {
        s << std::setprecision(precision) << static_cast<PrintType<T>>(data) << " ";
        ++counter;
        if (counter%numPerLine==0) {s << std::endl; indent(s,indentCount);}
      }
      else
        *base64 << data;
    }
    
    //! write one data element to output stream
    template <class T>
    template <int m>
    void VTKDataArrayWriter<T>::write (Dune::FieldVector<T,m> const& data)
    {
      for (int i=0; i<ncomps; ++i)
        write(i<m? data[i]: 0);
    }
    
    //! finish output; writes end tag
    template<class T>
    VTKDataArrayWriter<T>::~VTKDataArrayWriter ()
    {
      base64.reset(); // writes remaining binary data, if any.
      
      if (outputType==IoOptions::ascii && counter%numPerLine!=0) 
      {
        s << std::endl; 
        indent(s,indentCount);
      }
      s << "</DataArray>" << std::endl;
    }
    
    // explicit instantiaion of class
    template class VTKDataArrayWriter<double>;
    template class VTKDataArrayWriter<float>;
    template class VTKDataArrayWriter<int>;
    template class VTKDataArrayWriter<unsigned char>;
    
    template void VTKDataArrayWriter<double>::template write<1>(Dune::FieldVector<double,1> const& data);
    template void VTKDataArrayWriter<double>::template write<2>(Dune::FieldVector<double,2> const& data);
    template void VTKDataArrayWriter<double>::template write<3>(Dune::FieldVector<double,3> const& data);
    template void VTKDataArrayWriter<float>::template write<1>(Dune::FieldVector<float,1> const& data);
    template void VTKDataArrayWriter<float>::template write<2>(Dune::FieldVector<float,2> const& data);
    template void VTKDataArrayWriter<float>::template write<3>(Dune::FieldVector<float,3> const& data);
    
    // ----------------------------------------------------------
    
  }
}



#ifdef UNITTEST

#include <iostream>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include <fem/functionspace.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk2.hh>
#include <utilities/gridGeneration.hh>

using namespace Kaskade;


int main(void)
{
  using Grid = Dune::UGGrid<2>;
  Dune::FieldVector<double,2> x0(0), dx(1);
  GridManager<Grid> gridManager(createRectangle<Grid>(x0,dx,1.0));
  
  typedef Grid::LeafGridView LeafView;
  LeafView leafView = gridManager.grid().leafGridView();

  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView>>;
  H1Space h1Space2(gridManager,leafView,2);
  
  using L2Space = FEFunctionSpace<DiscontinuousLagrangeMapper<double,LeafView>>;
  L2Space l2Space2(gridManager,leafView,2);
  L2Space l2Space0(gridManager,leafView,0);

  // Create a scalar discontinuos piecewise quadratic function
  L2Space::Element_t<1> f(l2Space2);
  for (int i=0; i<f.size(); ++i)
    f[i] = i;
  
  H1Space::Element_t<1> g(h1Space2); 
  g = f;
  
  L2Space::Element_t<1> f0(l2Space0);
  f0 = f;
  
  
  writeVTK(f,"testVTKconforming2-l2",IoOptions().setOrder(2).setDataMode(IoOptions::conforming),"f");
  writeVTK(f,"testVTKnonconforming2-l2",IoOptions().setOrder(2).setDataMode(IoOptions::nonconforming),"f");
  writeVTK(g,"testVTKconforming2-h1",IoOptions().setOrder(2).setDataMode(IoOptions::conforming),"g");
  writeVTK(g,"testVTKnonconforming2-h1",IoOptions().setOrder(2).setDataMode(IoOptions::nonconforming),"g");
  writeVTK(g,"testVTKconforming1-l2",IoOptions().setOrder(1).setDataMode(IoOptions::conforming),"g");
  writeVTK(f,"testVTKnonconforming1-l2",IoOptions().setOrder(1).setDataMode(IoOptions::nonconforming),"f");
  writeVTK(g,"testVTKconforming1-h1",IoOptions().setOrder(1).setDataMode(IoOptions::conforming),"g");
  writeVTK(g,"testVTKnonconforming1-h1",IoOptions().setOrder(1).setDataMode(IoOptions::nonconforming),"g");
  writeVTK(f0,"testVTKconforming1-l0",IoOptions().setOrder(1).setDataMode(IoOptions::conforming),"f0");

  
  return 0;
}

#endif


