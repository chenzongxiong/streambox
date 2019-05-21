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
#ifndef CHECK_DERIVATIVE_HH
#define CHECK_DERIVATIVE_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <boost/mpl/size.hpp>
#include <boost/lexical_cast.hpp>
#include "algorithm/dune_bridge.hh"
#include "tools/linalg/scalarproducts.hh"

/// Some helper's for a correct indentation of XML-code
namespace XMLHelper{

  static int xml_level;

  struct XMLBlock{
    XMLBlock(std::ofstream &f, const char* s, const char* d="") : file(f), str(s), data(d), offset("")
    {
      xml_level += 2;
      getOffset();
      writeBegin();
    }
    XMLBlock(std::ofstream &f, std::string &s, std::string d="") : file(f), str(s), data(d), offset("")
    { 
      xml_level += 2;
      getOffset();
      writeBegin();
    }
    ~XMLBlock()
    {
      file << offset.c_str() << "</" << str.c_str() << ">" << std::endl;
      xml_level -= 2;
    }
    
  private:
    void getOffset(){
      for(int i=0; i<xml_level; ++i) offset.append(" ");
    }

    void writeBegin(){
      file << offset.c_str() << "<" << str.c_str() << data.c_str() << ">" << std::endl;
    }

    std::ofstream& file;
    std::string const str, data;
    std::string offset;
  };
  
}

/// Class that checks the derivatives of a functional at a linearization point.
/**
 * The analytic derivatives are compared with a finite difference approximation. Note that this
 * just gives you a hint (but a good one if you don't use stupid parameters for the tolerance and
 * the increment) about the correctness of your analytically calculated derivatives.
 * Keep in mind that it is a heuristic.
 */
template <class Functional, class EvaluationPoint, class SparseInt=int>
                                                                   class DerivativeChecker
                                                                   {
                                                                   public:
  typedef Bridge::KaskadeLinearization<Functional, EvaluationPoint> Linearization;
  typedef typename Functional::Scalar Scalar;
  enum{ idOfFirstBlock = 0,
        idOfLastBlock  = boost::mpl::size<typename Functional::AnsatzVars::Variables>::type::value };

  /// Constructor
  /**
   * \param functional instance of the Functional
   * \param x_ linearization point
   * \param t error tolerance (point-wise, default=1.0e-6)
   * \param incr increment for finite differences (default=1.0e-9)
   */
  DerivativeChecker(Functional &functional, EvaluationPoint &x_, Scalar const t = 1.0e-6, Scalar const incr = 1.0e-9) :
      x(x_), tolerance(t), increment(incr), linearization(functional, x)
  {
    checkFirstDerivative();
    checkSecondDerivative();
  }


  /// Check first derivative.
  /**
   * This function is called by the constructor. You only need to call this function if you have
   * changed the tolerance.
   */
  void checkFirstDerivative(){

    std::cout << "checking first derivative" << std::endl;
    int const numRows = linearization.rows(idOfFirstBlock, idOfLastBlock);
    std::cout << "number of rows: " << numRows << std::endl;
    std::cout << "assembling...";;

    // analytic solution
    std::vector<Scalar> analytic_solution;
    // storage for the evaluation of the derivative at x and the at x+increment*e_i
    Scalar distorted_value;

    linearization.flush();
    // get reference solution
    linearization.getRHSBlocks(analytic_solution, idOfFirstBlock, idOfLastBlock);
    // get undistorted rhs
    Scalar const reference_value = linearization.getValue();

    EvaluationPoint copyOf_x = linearization.getX();
    std::vector<Scalar> xAsVector(numRows);
    copyOf_x.write(xAsVector.begin());

    // reserve storage for finite difference solution
    firstDerivativeError.resize(xAsVector.size(),0);

    // get finite difference gradient
    for(size_t i=0; i<xAsVector.size(); ++i)
    {
      xAsVector[i] += increment;
      copyOf_x.read(xAsVector.begin());
      linearization.setX(copyOf_x);
      linearization.flush();
      distorted_value = linearization.getValue();
      xAsVector[i] -= increment;

      // store finite difference
      firstDerivativeError[i] = (distorted_value - reference_value)/increment;
      std::cout << ".";
    }
    std::cout << "done." << std::endl;
    for(size_t i=0; i<firstDerivativeError.size(); ++i)
      firstDerivativeError[i] -= analytic_solution[i];

    Scalar infinityError = InfinityNorm(firstDerivativeError);
    firstDerivativeOk = infinityError < tolerance;
    if(firstDerivativeOk)
      std::cout << "\n first derivative seems to be trustworthy!\n" << std::endl;
    else
      std::cout << "\n first derivative does not seem to be trustworthy!\n" << std::endl;
  }


  /// Check second derivative.
  /**
   * This function is called by the constructor. You only need to call this function if you have
   * changed the tolerance.
   */
  void checkSecondDerivative(){

    std::cout << "checking second derivative:" << std::endl;
    // Define linearization at x
    int const numRows = linearization.rows(idOfFirstBlock, idOfLastBlock);
    std::cout << "number of rows: " << numRows << std::endl;
    std::cout << "assembling...";


    // sparse storage for the analytic solution
    MatrixAsTriplet<Scalar,SparseInt> analytic_solution;
    // storage for the evaluation of the derivative at x and the at x+increment*e_i
    std::vector<Scalar> reference_value, distorted_value;
    linearization.flush();
    // get reference solution
    linearization.getMatrixBlocks(analytic_solution, idOfFirstBlock, idOfLastBlock, idOfFirstBlock, idOfLastBlock);
    // get undistorted rhs
    linearization.getRHSBlocks(reference_value, idOfFirstBlock, idOfLastBlock);

    // get evaluation point
    EvaluationPoint copyOf_x = linearization.getX();
    std::vector<Scalar> xAsVector(numRows);
    copyOf_x.write(xAsVector.begin());

    // get finite difference gradient matrix
    for(size_t i=0; i<xAsVector.size(); ++i)
    {
      xAsVector[i] += increment;
      copyOf_x.read(xAsVector.begin());
      linearization.setX(copyOf_x);
      linearization.flush();
      linearization.getRHSBlocks(distorted_value, idOfFirstBlock, idOfLastBlock);
      xAsVector[i] -= increment;

      // store finite difference
      for(size_t k=0; k<numRows; ++k)
      {
        secondDerivativeError.addEntry(k, i, (distorted_value[k] - reference_value[k])/increment);
      }
      std::cout << ".";
    }
    std::cout << "done." << std::endl;
    analytic_solution *= -1;
    secondDerivativeError += analytic_solution;


    Scalar infinityError = InfinityNorm(secondDerivativeError);
    secondDerivativeOk = tolerance > infinityError;
    if(secondDerivativeOk)
      std::cout << "\n second derivative seems to be trustworthy!\n" << std::endl;
    else
      std::cout << "\n second derivative does not seem to be trustworthy!\n" << std::endl;
  }

  /// get block
  MatrixAsTriplet<Scalar> getBlocksD2Error(int firstBlockId, int lastBlockId) const
  {
    assert(firstBlockId < lastBlockId);
    int offset = 0;
    if(firstBlockId > 0) offset = linearization.rows(0,firstBlockId);
    int const lastBlockIndex = offset + linearization.rows(firstBlockId, lastBlockId);

    MatrixAsTriplet<Scalar> blocks;
    for(int i=0; i<secondDerivativeError.size(); ++i)
    {
      if(secondDerivativeError.ridx[i] >= offset && secondDerivativeError.ridx[i] < lastBlockIndex &&
         secondDerivativeError.cidx[i] >= offset && secondDerivativeError.cidx[i] < lastBlockIndex)
        blocks.addEntry(secondDerivativeError.ridx[i], secondDerivativeError.cidx[i], secondDerivativeError.data[i]);
    }
    return blocks;
  }

  /// Only errors below tolerance will be displayed.
  void setTolerance(Scalar const tol){ tolerance = tol; }

  /// Set increment for finite differences
  void setIncrement(Scalar const incr){ increment = incr; }

  /// print error in first derivative (only values below tolerance)
  void printD1Error() const {
    std::cout << "components in D1 over tolerance:\n";
    for(size_t i=0; i<firstDerivativeError.size(); ++i)
      if(fabs(firstDerivativeError[i]) > tolerance)
        std::cout << i << ": " << firstDerivativeError[i] << std::endl;
  }

  /// print error in first derivative
  void printFullD1Error() const {
    std::cout << "components in D1:\n" << firstDerivativeError << std::endl;
  }

  /// print error in second derivative (only values below tolerance)
  void printD2Error() const {
    std::cout << "components in D2 over tolerance:\n";
    secondDerivativeError.print(std::cout,tolerance);
  }

  /// print error in second derivative
  void printFullD2Error() const {
    std::cout << "components in D2:\n";
    secondDerivativeError.print(std::cout);
  }

  /// write error in first derivative to .vtu file
  void d1ErrorToVTK()  {

    using std::endl;
    using std::string;

    size_t const numberOfCells = firstDerivativeError.size(),
                 numberOfPoints = (1+numberOfCells)*2;
    size_t const cell_type = 8,
    offset = 4;

    XMLHelper::xml_level = 0;
    typedef typename XMLHelper::XMLBlock XMLBlock;

    std::ofstream file("d1error.vtu");
    file << "<?xml version=\"1.0\"?>" << endl;
    { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
      { XMLBlock xgrid(file, "UnstructuredGrid");
        { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" NumberOfCells=\"" +
                                                      boost::lexical_cast<string>(numberOfCells) + "\"";
          XMLBlock xpiece(file, "Piece", str.c_str());
          { XMLBlock xcd(file, "CellData");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << fabs(firstDerivativeError[cell]) << " ";
                if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
              }
            }
          }
          { XMLBlock xpoints(file, "Points");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells+1; ++cell){
                if(cell%2==0) file << "    ";
                file << cell << " 0 0 " << cell << " 1 0 ";
                if(cell%2==1 || cell==numberOfCells) file << endl;
              }
            }
          }
          { XMLBlock xcells(file, "Cells");
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%3==0) file << "    ";
                file << 2*cell << " " << 2*cell+1 << " " << 2*cell+2 << " " << 2*cell+3 << " " << endl;
                if(cell%3==11 || cell==(numberOfCells-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << cell*offset << " ";
                if(cell%12==1 || cell==(numberOfCells-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << cell_type << " ";
                if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
              }
            }
          }
        }
      }
    }
  }

  /// write error in second derivative to .vtu file
  void d2ErrorToVTK()  {

    using std::endl;
    using std::string;

    /// store index pairs
    std::vector<std::pair<SparseInt,SparseInt> > indices(secondDerivativeError.data.size());
    for(size_t i=0; i<indices.size(); ++i)
      indices[i] = std::make_pair(secondDerivativeError.ridx[i], secondDerivativeError.cidx[i]);

    SparseInt const numberOfRows = secondDerivativeError.nrows()+1,
    numberOfColumns = secondDerivativeError.ncols()+1,
    numberOfPoints = numberOfRows*numberOfColumns,
    numberOfCells = (numberOfRows-1)*(numberOfColumns-1);
    size_t const cell_type = 8, // = vtk_pixel
    offset = 4;

    XMLHelper::xml_level = 0;
    typedef typename XMLHelper::XMLBlock XMLBlock;

    std::ofstream file("d2error.vtu");
    file << "<?xml version=\"1.0\"?>" << endl;
    { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
      { XMLBlock xgrid(file, "UnstructuredGrid");
        { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" NumberOfCells=\"" +
                                             boost::lexical_cast<string>(numberOfCells) + "\"";
          XMLBlock xpiece(file, "Piece", str.c_str());
          { XMLBlock xcd(file, "CellData");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << fabs(firstDerivativeError[cell]) << " ";
                if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
              }
            }
          }
          { XMLBlock xp(file, "Points");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells+1; ++cell){
                if(cell%2==0) file << "    ";
                file << cell << " 0 0 " << cell << " 1 0 ";
                if(cell%2==1 || cell==numberOfCells) file << endl;
              }
            }
          }
          { XMLBlock xcell(file, "Cells");
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%3==0) file << "    ";
                file << 2*cell << " " << 2*cell+1 << " " << 2*cell+2 << " " << 2*cell+3 << " " << endl;
                if(cell%3==11 || cell==(numberOfCells-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << cell*offset << " ";
                if(cell%12==1 || cell==(numberOfCells-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << "    ";
                file << cell_type << " ";
                if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
              }
            }
          }
        }
      }
    }
  }

private:
  EvaluationPoint & x;
  Scalar tolerance, increment;
  Linearization linearization;
  std::vector<Scalar> firstDerivativeError;
  MatrixAsTriplet<Scalar,SparseInt> secondDerivativeError;
  bool firstDerivativeOk, secondDerivativeOk;
};
#endif /* CHECK_DERIVATIVE_HH */
