/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2009-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef IOBASE_HH
#define IOBASE_HH

#include <string>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/zip_view.hpp>
#include <boost/fusion/include/vector.hpp>

/**
 * @file
 * @brief  Output of mesh and solution for visualization software
 * @author Martin Weiser
 *
 * This file contains common basis functionality needed to write mesh
 * data and solution data for visualization software,e.g., in
 *
 * - writeVTK: write data to VTK XML file, e.g. for visualization in Paraview software
 * - TimeSeriesVTKWriter: write a (time-) sequence of variable sets into VTK XML files
 * - writeAMIRAFile: write data to file in AmiraMesh format, e.g. for visualization in AMIRA software
 *
 */

//---------------------------------------------------------------------
//---------------------------------------------------------------------

namespace Kaskade
{
  struct UnaryTrue {
    template <class T>
    bool operator()(T const&) const { return true; }
  };


  /**
   * \cond internals
   */
  namespace IoDetail {


    template <class GridView,
              class Writer,
              class Filter,
              class RAIter>
    struct AddDataTo {};


  } // End of namespace IoDetail
  /**
   * \endcond
   */

  /**
   * \ingroup IO
   * \brief options for VTK/AMIRA output 
   *
   * IoOptions are used to specify details for output,e.g., ascii or binary format. 
   * For ascii format, output precision (i.e. number of significant digits to be written) can be specified; 
   * if not set explicitly, six digits will be used.
   * \code
   * writeVTK(variableSet,"temperaturefile",IoOptions().setOutputType(IoOptions::ascii).setPrecision(10)); 
   * \endcode
   */
  struct IoOptions {
    /**
     * \brief Determines text or binary output format.
     * Currently this is only used by VTK output.
     */
    enum OutputType {
      /** Output data arrays in ascii. */
      ascii,
      /** Output data arrays in base64 encoded binary. */
      binary
    };
    
    /** 
     * \brief The data mode determines the continuity structure of the output.
     */
    enum DataMode {
      /**
       * Cells share a node at a common vertex. This means that point data (written for 
       * linear and quadratic FE) is necessarily continuous.
       */ 
      conforming, 
      /**
       * Joint vertices are duplicated, i.e. each cell has an independent set of nodes (in effect, 
       * the grid is just a bunch of disconnected cells). Thus, also point data can be discontinuous. 
       * The drawback is that the cell connectivity is lost and the file size increased.
       */
      nonconforming, 
      /**
       * Tells the file writer to choose the most appropriate version automatically.
       */
      inferred
    };
    
    enum Info {
      none, summary, detail
    };
    
    /** 
     *\brief Defaults: outputType=ascii, dataMode=conforming, info=none, 6 decimal digits, order=1 
     */
    IoOptions()
    : outputType(ascii), dataMode(conforming), info(none), precision(6), order(1)
    {};
    
    /**
     * \brief Sets the polynomial order of output.
     * 
     * VTK file formats currently support order 0 (discontinuous piecewise constant), 1 (linear continuous or discontinuous), 
     * and 2 (quadratic continuous or discontinuous). Other orders are not supported.
     */
    IoOptions& setOrder(int ord) { order = ord; return *this; }
    
    /**
     * \brief Sets the output type to either ascii or binary.
     * 
     * In VTK files, both define XML file types. The coefficient data is either written as ASCII floating point
     * values or as base64 encoded binary data.
     */
    IoOptions& setOutputType(OutputType out) { outputType = out; return *this; }

    /**
     * \brief Sets the number of decimal digits to be written for numbers in ASCII output mode.
     */
    IoOptions& setPrecision(int pre) { precision = pre; return *this; }
    
    /**
     * \brief Sets the data mode (conforming or nonconforming).
     */
    IoOptions& setDataMode(DataMode mode) { dataMode = mode; return *this; }
    
    OutputType outputType;
    DataMode dataMode;
    Info info;
    int precision;
    int order;
  };

  extern IoOptions ioOptions_default;
  //---------------------------------------------------------------------
  
  /**
   * \ingroup IO
   * \brief creates a zero-padded string representation of the given number
   * 
   * Use this to construct sequences of filenames.
   * 
   * \param n the number to represent as string
   * \param places the total number of digits
   */
  std::string paddedString(int n, int places=3);
  
  //---------------------------------------------------------------------

  /**
   * \ingroup IO
   * Writes a subset of variables contained in the variable set \a vars
   * with names taken from the \a description to an Amira or VTK XML
   * file with the name specified by \a filename. The subset is
   * specified by the STL predicate \a filter working on the variable
   * names.
   *	@param gridView GridView, i.e. LeafGridView
   *	@param writer writer for mesh data and solution data, e.g., LeafAmiraMeshWriter, VTKWriter
   *  @param vars data of variables
   *  @param filename name of output file
   *	@param filter filter for selecting variables to be written
   *  @param options defines the format for the output,e.g., ascii/binary
   */
  template <class GridView,
            class Writer,
            class VariableSet,
            class Filter,
            class Options>
  void writePartialFile(GridView const& gridView,
                        Writer& writer,
                        VariableSet const& vars,
                        std::string const& filename,
                        Filter const& filter,
                        Options const& options)
  {
    using namespace boost::fusion;

    // Bug in boost::fusion? The following line does not work!
    // for_each(zip(description.variables,vars.data),AddDataTo<VTKWriter>(writer));

    typedef vector<typename VariableSet::Descriptions::Variables const&,
                   typename VariableSet::Functions const&> ZipVector;
    zip_view<ZipVector> zipView(ZipVector(typename VariableSet::Descriptions::Variables(),vars.data));
    for_each(zipView,IoDetail::AddDataTo<GridView,Writer,Filter,std::string const*>(gridView,writer,filter,&vars.descriptions.names[0]));

    writer.write(filename.c_str(),options);
  }

} // namespace Kaskade

#endif
