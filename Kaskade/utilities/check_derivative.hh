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
#include <string>
#include <utility>
#include <vector>

#include <boost/mpl/size.hpp>
#include <boost/lexical_cast.hpp>

#include "algorithm/dune_bridge.hh"
#include "utilities/linalg/scalarproducts.hh"

namespace Kaskade{
/// Some helper's for a correct indentation of XML-code
namespace XMLHelper{

  static int const xml_increment=2;
  static int xml_level;

  std::string getXMLIndent()
  {
    std::string result;
    for(int i=0; i<xml_level; ++i) result.append(" ");
    return result;
  }

  struct XMLBlock{
    XMLBlock(std::ofstream &f, const char* s, const char* d="") : file(f), str(s), data(d), offset("")
    {
      xml_level += xml_increment;
      getOffset();
      writeBegin();
    }
    XMLBlock(std::ofstream &f, std::string &s, std::string d="") : file(f), str(s), data(d), offset("")
    { 
      xml_level += xml_increment;
      getOffset();
      writeBegin();
    }
    ~XMLBlock()
    {
      file << offset.c_str() << "</" << str.c_str() << ">" << std::endl;
      xml_level -= xml_increment;
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

namespace DerivativeCheck
{
  namespace Private
  {
    template <class SparseInt>
    void writeConnectivity(std::ofstream &file, SparseInt row, SparseInt col, SparseInt numRows)
    {
      int const upperLeftVertexId = row+col*numRows;
      file << upperLeftVertexId << " " << upperLeftVertexId+1 << " " << upperLeftVertexId+numRows << " " << upperLeftVertexId+numRows+1 << " ";
    }
  }

  template <class Assembler, class Functional,int firstBlock=0, int lastBlock = Functional::AnsatzVars::noOfVariables>
  typename Functional::Scalar d1(Assembler& assembler, Functional const& f, typename Functional::AnsatzVars::VariableSet& x, typename Functional::Scalar tolerance = 1e-6, typename Functional::Scalar increment = 1e-12, bool toFile=false, std::string const& filename = std::string("d1error"))
  {
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<firstBlock,lastBlock>::type CoefficientVector;

    std::cout << "checking first derivative" << std::endl;
    std::cout << "assembling...";

    // storage for the evaluation of the derivative at x and the at x+increment*e_i
    Scalar distorted;

    // get reference solution
    assembler.assemble(linearization(f,x), Assembler::RHS | Assembler::VALUE );

    CoefficientVector analytic(assembler.template rhs<firstBlock,lastBlock>());
    std::vector<Scalar> analyticD1(assembler.template nrows<firstBlock,lastBlock>());
    analytic.write(analyticD1.begin());

    // get undistorted rhs
    Scalar const reference_value = assembler.functional();

    std::vector<Scalar> xAsVector(assembler.nrows());
    x.write(xAsVector.begin());

    // reserve storage for finite difference solution
    std::vector<Scalar> firstDerivativeError(analyticD1.size(),0.);

    size_t offset = (0!=firstBlock) ? assembler.template nrows<0,firstBlock>() : 0;
    size_t end = (lastBlock!=Functional::AnsatzVars::noOfVariables) ? (offset + assembler.template nrows<firstBlock,lastBlock>()) : xAsVector.size();

    // get finite difference gradient
    for(size_t i=offset; i<end; ++i)
    {

      xAsVector[i] += increment;
      x.read(xAsVector.begin());

      assembler.assemble(linearization(f,x),Assembler::VALUE);
      distorted = assembler.functional();
      xAsVector[i] -= increment;

      // store finite difference
      firstDerivativeError[i-offset] = (distorted - reference_value)/increment;
      if(i%100==0) std::cout << ".";
    }
    // set x back to initial state
    x.read(xAsVector.begin());
    std::cout << "done." << std::endl;
    for(size_t i=0; i<firstDerivativeError.size(); ++i)
        firstDerivativeError[i] -= analyticD1[i];

    Scalar infinityError = LinAlg::InfinityNorm()(firstDerivativeError);
    bool errorBelowTolerance = infinityError < tolerance;
    if(errorBelowTolerance)
      std::cout << "\n first derivative seems to be trustworthy!\nInfinityError: " << infinityError << std::endl << std::endl;
    else
      std::cout << "\n first derivative does not seem to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;

    if(toFile) vectorToVTK(firstDerivativeError);

    return infinityError;
  }

  template <class Assembler, class Functional, int firstRow=0, int lastRow=Functional::TestVars::noOfVariables, int firstCol=0, int lastCol=Functional::AnsatzVars::noOfVariables>
  typename Functional::Scalar d2(Assembler& assembler, Functional const& f, typename Functional::AnsatzVars::VariableSet const& x, typename Functional::Scalar increment = 1e-12, typename Functional::Scalar tolerance = 1e-6, bool toFile = false, std::string const& savefilename = std::string("d2error"))
  {
    typedef typename Functional::Scalar Scalar;
    typedef typename Functional::TestVars::template CoefficientVectorRepresentation<firstRow,lastRow>::type  RowVector;
    typedef typename Functional::AnsatzVars::template CoefficientVectorRepresentation<firstCol,lastCol>::type ColumnVector;

    std::cout << "checking second derivative:" << std::endl;
    std::cout << "assembling...";

    // sparse storage for the analytic solution
    typedef MatrixAsTriplet<Scalar> SparseMatrix;
    SparseMatrix secondDerivativeError;

    // get reference solution
    assembler.assemble(linearization(f,x), Assembler::RHS, Assembler::MATRIX);
    SparseMatrix analytic_solution = assembler.template get<SparseMatrix,firstRow,lastRow,firstCol,lastCol>(false);
    RowVector reference_value(assembler. template rhs<firstRow,lastRow>());
    RowVector copyOfRefVal(reference_value);

    // get evaluation point
    std::vector<Scalar> xAsVector(assembler.nrows());
    x.write(xAsVector.begin());

    size_t offset = (firstRow != 0) ? assembler.template nrows<0,firstRow>() : 0;
    size_t end = (lastRow!=Functional::TestVars::noOfVariables) ? (offset + assembler.template nrows<firstRow,lastRow>()) : xAsVector.size();

    // get finite difference gradient matrix
    for(size_t i=offset; i<end; ++i)
    {
      copyOfRefVal = reference_value;
      xAsVector[i] += increment;
      x.read(xAsVector.begin());

      assembler.assemble(linearization(f,x), Assembler::RHS);
      RowVector distorted_value(assembler.template rhs<firstRow,lastRow>());

      xAsVector[i] -= increment;

      std::vector<Scalar> column(end-offset);
      copyOfRefVal -= distorted_value;
      copyOfRefVal.write(column.begin());

      secondDerivativeError.addColumn(column,i);
      // store finite difference
      //      for(size_t k=0; k<numRows; ++k)
      //      {
      //        secondDerivativeError.addEntry(k, i, (distorted_value[k] - reference_value[k])/increment);
      //      }
      if(i%100==0) std::cout << ".";
    }
    // set x back to initial state
    x.read(xAsVector.begin());
    std::cout << "done." << std::endl;
    analytic_solution *= -1;
    secondDerivativeError += analytic_solution;


    Scalar infinityError = LinAlg::InfinityNorm()(secondDerivativeError);
    if( (tolerance > infinityError) )
      std::cout << "\nSecond derivative seems to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;
    else
      std::cout << "\nSecond derivative does not seem to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;

    return infinityError;
  }


  template <class Vector>
  void vectorToVTK(Vector const& vec, std::string savefilename = std::string("d1error"))
  {
    using std::endl;
    using std::string;

    size_t const numberOfRows = vec.size(),
                 numberOfPoints = (1+numberOfRows)*2;
    size_t const cell_type = 8,
                 offset = 4;

    XMLHelper::xml_level = 0;
    typedef typename XMLHelper::XMLBlock XMLBlock;

    savefilename.append(".vtu");
    std::ofstream file(savefilename.c_str());
    file << "<?xml version=\"1.0\"?>" << endl;
    { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
      { XMLBlock xgrid(file, "UnstructuredGrid");
        { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" numberOfRows=\"" +
          boost::lexical_cast<string>(numberOfRows) + "\"";
          XMLBlock xpiece(file, "Piece", str.c_str());
          { XMLBlock xcd(file, "CellData");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                // invert row order for visualization tools
                file << fabs(vec[numberOfRows-1-cell]) << " ";
                if(cell%12==11 || cell==(numberOfRows-1) ) file << endl;
              }
            }
          }
          { XMLBlock xpoints(file, "Points");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows+1; ++cell){
                if(cell%2==0) file << "    ";
                file << cell << " 0 0 " << cell << " 1 0 ";
                if(cell%2==1 || cell==numberOfRows) file << endl;
              }
            }
          }
          { XMLBlock xcells(file, "Cells");
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%3==0) file << "    ";
                file << 2*cell << " " << 2*cell+1 << " " << 2*cell+2 << " " << 2*cell+3 << " " << endl;
                if(cell%3==11 || cell==(numberOfRows-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                file << cell*offset << " ";
                if(cell%12==1 || cell==(numberOfRows-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                file << cell_type << " ";
                if(cell%12==11 || cell==(numberOfRows-1) ) file << endl;
              }
            }
          }
        }
      }
    }
  }

  template <class Scalar, class SparseInt>
  void matrixToVTK(MatrixAsTriplet<Scalar,SparseInt> const& matrix, std::string savefilename = std::string("d2error"))
    {
      using std::cout;
      using std::endl;
      using std::string;

      cout << "writing file";

      /// store index pairs
      std::vector<std::pair<SparseInt,SparseInt> > indices(matrix.data.size());
      for(size_t i=0; i<indices.size(); ++i)
        indices[i] = std::make_pair(matrix.ridx[i], matrix.cidx[i]);

      SparseInt const numberOfRows = matrix.nrows()+1,
          numberOfColumns = matrix.ncols()+1,
          numberOfPoints = numberOfRows*numberOfColumns,
          numberOfCells = (numberOfRows-1)*(numberOfColumns-1);
      size_t const cell_type = 8, // = vtk_pixel
          offset = 4;

      XMLHelper::xml_level = 0;
      string indent;
      typedef typename XMLHelper::XMLBlock XMLBlock;

      savefilename.append(".vtu");
      std::ofstream file(savefilename.c_str());
      file << "<?xml version=\"1.0\"?>" << endl;
      { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
        { XMLBlock xgrid(file, "UnstructuredGrid");
          { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" NumberOfCells=\"" +
            boost::lexical_cast<string>(numberOfCells) + "\"";
            XMLBlock xpiece(file, "Piece", str.c_str());
            { XMLBlock xcd(file, "CellData");
              { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
                int cell = 0;
                indent = XMLHelper::getXMLIndent();
                for(SparseInt col=0; col<numberOfColumns-1; ++col){
                  for(SparseInt row=0; row<numberOfRows-1; ++row){
                    if(cell%12==0) file << indent.c_str();
                    typename std::vector<std::pair<SparseInt,SparseInt> >::const_iterator iter = std::find(indices.begin(), indices.end(), std::make_pair(numberOfRows-2-row,col)),
                        iterStart = indices.begin();
                    if(iter!=indices.end())
                      file << fabs(matrix.data[std::distance(iterStart, iter)]) << " ";
                    else
                      file << 0.0 << " ";
                    if(cell%12==11 || cell==numberOfCells-1) file << endl;
                    ++cell;
                  }
                }
              }
            }
            cout << ".";
            { XMLBlock xp(file, "Points");
              { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
                int cell = 0;
                for(SparseInt col=0; col<numberOfColumns; ++col){
                  for(SparseInt row=0; row<numberOfRows; ++row){
                    if(cell%4==0) file << indent.c_str();
                    file << row << " " << col << " 0 ";
                    if(cell%4==3 || cell==(numberOfCells-1)) file << endl;
                    ++cell;
                  }
                }
              }
            }
            cout << ".";
            { XMLBlock xcell(file, "Cells");
              { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
                int cell = 0;
                for(SparseInt col=0; col<numberOfColumns-1; ++col){
                  for(SparseInt row=0; row<numberOfRows-1; ++row){
                    if(cell%3==0) file << indent.c_str();
                    Private::writeConnectivity(file, row, col, numberOfRows);
                    if(cell%3==2 || cell==(numberOfCells-1)) file << endl;
                    ++cell;
                  }
                }
              }
              cout << ".";
              { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
                for(SparseInt cell=0; cell<numberOfCells; ++cell){
                  if(cell%12==0) file << indent.c_str();
                  file << (cell+1)*offset << " ";
                  if(cell%12==11 || cell==(numberOfCells-1)) file << endl;
                }
              }
              cout << ".";
              { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
                for(SparseInt cell=0; cell<numberOfCells; ++cell){
                  if(cell%12==0) file << indent.c_str();
                  file << cell_type << " ";
                  if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
                }
              }
            }
          }
        }
      }
      cout << "done" << endl;
    }

} // end of namespace DerivativeCheck

/// Class that checks the derivatives of a functional at a linearization point.
/**
 * The analytic derivatives are compared with a finite difference approximation. Note that this
 * just gives you a hint (but a good one if you don't use stupid parameters for the tolerance and
 * the increment) about the correctness of your analytically calculated derivatives.
 * Keep in mind that it is a heuristic.
 *
 * \param checkD1 enable/disable check of first derivative (i.e. disable it if d0() has not been implemented
 *                correctly in boundary or domain cache (optional, default=true(enabled))
 * \param SparseInt integer type used for sparse matrix indices
 */
template <class Functional, bool checkD1=true, class SparseInt=int>
class DerivativeChecker
{
  typedef VariationalFunctionalAssembler<LinearizationAt<Functional> > Assembler;
  typedef typename Functional::AnsatzVars VariableSet;
  typedef typename VariableSet::VariableSet EvaluationPoint;
  typedef typename VariableSet::template CoefficientVectorRepresentation<> CoefficientVectorRepresentation;
  typedef typename CoefficientVectorRepresentation::type CoefficientVector;

public:
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
  DerivativeChecker(Assembler& assembler_, Functional const& f_, EvaluationPoint const& x_, VariableSet const& variableSet_, Scalar const tol = 1.0e-6, Scalar const incr = 1.0e-9) :
    assembler(assembler_), f(f_), x(x_), variableSet(variableSet_), tolerance(tol), increment(incr)
  {
    if(checkD1) checkFirstDerivative();
    checkSecondDerivative();
  }


  /// Check first derivative.
  /**
   * This function is called by the constructor. You only need to call this function if you have
   * changed the tolerance.
   */
  void checkFirstDerivative()
  {
    std::cout << "checking first derivative" << std::endl;
    std::cout << "assembling...";

    // storage for the evaluation of the derivative at x and the at x+increment*e_i
    Scalar distorted;

    // get reference solution
    assembler.assemble(linearization(f,x), Assembler::RHS | Assembler::VALUE );

    CoefficientVector analytic(assembler.rhs());
    std::vector<Scalar> analyticD1(assembler.nrows());
    analytic.write(analyticD1.begin());

    // get undistorted rhs
    Scalar const reference_value = assembler.functional();

    std::vector<Scalar> xAsVector(analyticD1.size());
    x.write(xAsVector.begin());

    // reserve storage for finite difference solution
    firstDerivativeError.resize(xAsVector.size(),0);

    EvaluationPoint copyOfX(x);
    // get finite difference gradient
    for(size_t i=0; i<xAsVector.size(); ++i)
    {
      xAsVector[i] += increment;
      copyOfX.read(xAsVector.begin());

      assembler.assemble(linearization(f,copyOfX),Assembler::VALUE);
      Scalar distorted_value = assembler.functional();
      xAsVector[i] -= increment;

      // store finite difference
      firstDerivativeError[i] = (distorted_value - reference_value)/increment;
      if(i%100==0) std::cout << ".";
    }
    std::cout << "done." << std::endl;
    for(size_t i=0; i<firstDerivativeError.size(); ++i)
      firstDerivativeError[i] -= analyticD1[i];

    Scalar infinityError = LinAlg::InfinityNorm()(firstDerivativeError);
    firstDerivativeOk_ = infinityError < tolerance;
    if(firstDerivativeOk_)
      std::cout << "\n first derivative seems to be trustworthy!\nInfinityError: " << infinityError << std::endl << std::endl;
    else
      std::cout << "\n first derivative does not seem to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;
  }


  /// Check second derivative.
  /**
   * This function is called by the constructor. You only need to call this function if you have
   * changed the tolerance.
   */
  void checkSecondDerivative()
  {
    std::cout << "checking second derivative:" << std::endl;
    // Define linearization at x
    std::cout << "assembling...";

    // sparse storage for the analytic solution
    typedef MatrixAsTriplet<Scalar,SparseInt> SparseMatrix;

    // get reference solution
    assembler.assemble(linearization(f,x), Assembler::RHS, Assembler::MATRIX);
    SparseMatrix analytic_solution = assembler.template get<SparseMatrix>(false);
    CoefficientVector reference_value(assembler.rhs());
    CoefficientVector copyOfRefVal(reference_value);

    // get evaluation point
    std::vector<Scalar> xAsVector(assembler.nrows());
    x.write(xAsVector.begin());

    EvaluationPoint copyOfX(x);
    // get finite difference gradient matrix
    for(size_t i=0; i<xAsVector.size(); ++i)
    {
      copyOfRefVal = reference_value;
      xAsVector[i] += increment;
      copyOfX.read(xAsVector.begin());

      assembler.assemble(linearization(f,copyOfX), Assembler::RHS);
      CoefficientVector distorted_value(assembler.rhs());

      xAsVector[i] -= increment;

      std::vector<Scalar> column(xAsVector.size());
      copyOfRefVal -= distorted_value;
      copyOfRefVal.write(column.begin());

      secondDerivativeError.addColumn(column,i);
      // store finite difference
//      for(size_t k=0; k<numRows; ++k)
//      {
//        secondDerivativeError.addEntry(k, i, (distorted_value[k] - reference_value[k])/increment);
//      }
      if(i%100==0) std::cout << ".";
    }
    std::cout << "done." << std::endl;
    analytic_solution *= -1;
    secondDerivativeError += analytic_solution;


    Scalar infinityError = LinAlg::InfinityNorm()(secondDerivativeError);
    if( (secondDerivativeOk_ = tolerance > infinityError) )
      std::cout << "\nSecond derivative seems to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;
    else
      std::cout << "\nSecond derivative does not seem to be trustworthy!\nInfinity error: " << infinityError << std::endl << std::endl;
  }

  bool firstDerivativeOk() const { return firstDerivativeOk_; }

  bool secondDerivativeOk() const { return secondDerivativeOk_; }

//  /// get block of d1Error
//  std::vector<Scalar> getBlocksD1Error(int const firstBlockId, int const lastBlockId) const
//  {
//    assert(firstBlockId < lastBlockId);
//    std::pair<int,int> ids = getRowIds(firstBlockId, lastBlockId);
//    std::vector<Scalar> result;
//    for(int row = ids.first; row<ids.second; ++row)
//      result.push_back(firstDerivativeError[row]);
//    return result;
//  }
//
//  /// get block of d2Error
//  MatrixAsTriplet<Scalar> getBlocksD2Error(int const firstBlockId, int const lastBlockId) const
//  {
//    assert(firstBlockId < lastBlockId);
//    std::pair<int,int> ids = getRowIds(firstBlockId, lastBlockId);
//
//    MatrixAsTriplet<Scalar> blocks;
//    for(int i=0; i<secondDerivativeError.size(); ++i)
//    {
//      if(secondDerivativeError.ridx[i] >= ids.first && secondDerivativeError.ridx[i] < ids.second &&
//          secondDerivativeError.cidx[i] >= ids.first && secondDerivativeError.cidx[i] < ids.second)
//        blocks.addEntry(secondDerivativeError.ridx[i], secondDerivativeError.cidx[i], secondDerivativeError.data[i]);
//    }
//    return blocks;
//  }

  /// Only errors below tolerance will be displayed.
  void setTolerance(Scalar const tol){ tolerance = tol; }

  /// Set increment for finite differences
  void setIncrement(Scalar const incr){ increment = incr; }

  /// print error in first derivative (only values below tolerance)
  void printD1Error() const
  {
    std::cout << "components in D1 over tolerance:\n";
    for(size_t i=0; i<firstDerivativeError.size(); ++i)
      if(fabs(firstDerivativeError[i]) > tolerance)
        std::cout << i << ": " << firstDerivativeError[i] << std::endl;
  }

  /// print error in first derivative
  void printFullD1Error() const
  {
    std::cout << "components in D1:\n" << firstDerivativeError << std::endl;
  }

  /// print error in second derivative (only values below tolerance)
  void printD2Error() const
  {
    std::cout << "components in D2 over tolerance:\n";
    secondDerivativeError.print(std::cout,tolerance);
  }

  /// print error in second derivative
  void printFullD2Error() const
  {
    std::cout << "components in D2:\n";
    secondDerivativeError.print(std::cout);
  }

  /// write d2 to .vtu file
  void d2ToVTK(){
    assembler.assemble(linearization(f,x),Assembler::MATRIX);
    // sparse storage for the analytic solution
    MatrixAsTriplet<Scalar,SparseInt> analytic_solution = assembler.template get<MatrixAsTriplet<Scalar,SparseInt>>(false);
    // get reference solution
    matrixToVTK(analytic_solution);
  }

  /// write error in first derivative to .vtu file
  void d1ErrorToVTK(int const steps) const { vectorToVTK(firstDerivativeError, steps); }

//  /// write block of error in first derivative to .vtu file
//  void d1BlockErrorToVTK(int const firstBlockId, int const lastBlockId) const { vectorToVTK(getBlocksD1Error(firstBlockId,lastBlockId)); }

  /// write Vector vector to .vtu file
  /**
   * Vector must provide the function size() and  operator[](size_type n)
   */
  template <class Vector>
  void vectorToVTK(Vector &vector, int const steps=0) const {

    using std::endl;
    using std::string;

    size_t const numberOfRows = vector.size(),
        numberOfPoints = (1+numberOfRows)*2;
    size_t const cell_type = 8,
        offset = 4;

    XMLHelper::xml_level = 0;
    typedef typename XMLHelper::XMLBlock XMLBlock;

    string filename("d1error");
    filename.append(boost::lexical_cast<string>(steps)).append(".vtu");

    std::ofstream file(filename.c_str());
    file << "<?xml version=\"1.0\"?>" << endl;
    { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
      { XMLBlock xgrid(file, "UnstructuredGrid");
        { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" numberOfRows=\"" +
          boost::lexical_cast<string>(numberOfRows) + "\"";
          XMLBlock xpiece(file, "Piece", str.c_str());
          { XMLBlock xcd(file, "CellData");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                // invert row order for visualization tools
                file << fabs(vector[numberOfRows-1-cell]) << " ";
                if(cell%12==11 || cell==(numberOfRows-1) ) file << endl;
              }
            }
          }
          { XMLBlock xpoints(file, "Points");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows+1; ++cell){
                if(cell%2==0) file << "    ";
                file << cell << " 0 0 " << cell << " 1 0 ";
                if(cell%2==1 || cell==numberOfRows) file << endl;
              }
            }
          }
          { XMLBlock xcells(file, "Cells");
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%3==0) file << "    ";
                file << 2*cell << " " << 2*cell+1 << " " << 2*cell+2 << " " << 2*cell+3 << " " << endl;
                if(cell%3==11 || cell==(numberOfRows-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                file << cell*offset << " ";
                if(cell%12==1 || cell==(numberOfRows-1)) file << endl;
              }
            }
            { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfRows; ++cell){
                if(cell%12==0) file << "    ";
                file << cell_type << " ";
                if(cell%12==11 || cell==(numberOfRows-1) ) file << endl;
              }
            }
          }
        }
      }
    }
  }

  /// write error in second derivative to .vtu-file
  void d2ErrorToVTK(int const steps) const { matrixToVTK(secondDerivativeError,steps); }

//  /// write block of second derivative to .vtu-file
//  void d2BlockErrorToVTK(int firstBlockId, int lastBlockId) const { matrixToVTK(getBlocksD2Error(firstBlockId,lastBlockId)); }

  /// write, possibly sparse, matrix in triplet format to .vtu-file
  void matrixToVTK(MatrixAsTriplet<Scalar,SparseInt> const& matrix, int const steps=0) const
  {
    using std::cout;
    using std::endl;
    using std::string;

    cout << "writing file";

    /// store index pairs
    std::vector<std::pair<SparseInt,SparseInt> > indices(matrix.data.size());
    for(size_t i=0; i<indices.size(); ++i)
      indices[i] = std::make_pair(matrix.ridx[i], matrix.cidx[i]);

    SparseInt const numberOfRows = matrix.nrows()+1,
        numberOfColumns = matrix.ncols()+1,
        numberOfPoints = numberOfRows*numberOfColumns,
        numberOfCells = (numberOfRows-1)*(numberOfColumns-1);
    size_t const cell_type = 8, // = vtk_pixel
        offset = 4;

    XMLHelper::xml_level = 0;
    string indent;
    typedef typename XMLHelper::XMLBlock XMLBlock;

    string filename("d2error");
    filename.append(boost::lexical_cast<string>(steps)).append(".vtu");

    std::ofstream file(filename.c_str());
    file << "<?xml version=\"1.0\"?>" << endl;
    { XMLBlock xvtk(file, "VTKFile", " type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"");
      { XMLBlock xgrid(file, "UnstructuredGrid");
        { string str = string(" NumberOfPoints=\"") + boost::lexical_cast<string>(numberOfPoints) + " \" NumberOfCells=\"" +
          boost::lexical_cast<string>(numberOfCells) + "\"";
          XMLBlock xpiece(file, "Piece", str.c_str());
          { XMLBlock xcd(file, "CellData");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"d2error\" NumberOfComponents=\"1\" format=\"ascii\"");
              int cell = 0;
              indent = XMLHelper::getXMLIndent();
              for(SparseInt col=0; col<numberOfColumns-1; ++col){
                for(SparseInt row=0; row<numberOfRows-1; ++row){
                  if(cell%12==0) file << indent.c_str();
                  typename std::vector<std::pair<SparseInt,SparseInt> >::const_iterator iter = std::find(indices.begin(), indices.end(), std::make_pair(numberOfRows-2-row,col)),
                      iterStart = indices.begin();
                  if(iter!=indices.end())
                    file << fabs(matrix.data[std::distance(iterStart, iter)]) << " ";
                  if(cell%12==11 || cell==numberOfCells-1) file << endl;
                  ++cell;
                }
              }
            }
          }
          cout << ".";
          { XMLBlock xp(file, "Points");
            { XMLBlock xda(file, "DataArray", " type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\"");
              int cell = 0;
              for(size_t col=0; col<numberOfColumns; ++col){
                for(size_t row=0; row<numberOfRows; ++row){
                  if(cell%4==0) file << indent.c_str();
                  file << row << " " << col << " 0 ";
                  if(cell%4==3 || cell==(numberOfCells-1)) file << endl;
                  ++cell;
                }
              }
            }
          }
          cout << ".";
          { XMLBlock xcell(file, "Cells");
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
              int cell = 0;
              for(size_t col=0; col<numberOfColumns-1; ++col){
                for(size_t row=0; row<numberOfRows-1; ++row){
                  if(cell%3==0) file << indent.c_str();
                  writeConnectivity(file, row, col, numberOfRows);
                  if(cell%3==2 || cell==(numberOfCells-1)) file << endl;
                  ++cell;
                }
              }
            }
            cout << ".";
            { XMLBlock xda(file, "DataArray", " type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << indent.c_str();
                file << (cell+1)*offset << " ";
                if(cell%12==11 || cell==(numberOfCells-1)) file << endl;
              }
            }
            cout << ".";
            { XMLBlock xda(file, "DataArray", " type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
              for(size_t cell=0; cell<numberOfCells; ++cell){
                if(cell%12==0) file << indent.c_str();
                file << cell_type << " ";
                if(cell%12==11 || cell==(numberOfCells-1) ) file << endl;
              }
            }
          }
        }
      }
    }
    cout << "done" << endl;
  }

private:
//  /// get ids of first and last row in [firstBlockId,lastBlockId]
//  std::pair<int,int> getRowIds(int const firstBlockId, int const lastBlockId) const
//  {
//    int const firstIdx = (firstBlockId > 0) ? linearization.rows(0,firstBlockId) : 0;
//    int const lastIdx = firstIdx + linearization.rows(firstBlockId, lastBlockId);
//    return std::make_pair(firstIdx, lastIdx);
//  }

  void writeConnectivity(std::ofstream &file, SparseInt row, SparseInt col, SparseInt numRows) const
  {
    int const upperLeftVertexId = row+col*numRows;
    file << upperLeftVertexId << " " << upperLeftVertexId+1 << " " << upperLeftVertexId+numRows << " " << upperLeftVertexId+numRows+1 << " ";
  }

  Assembler& assembler;
  Functional const& f;
  EvaluationPoint x;
  VariableSet const& variableSet;
  Scalar tolerance, increment;
  std::vector<Scalar> firstDerivativeError;
  MatrixAsTriplet<Scalar,SparseInt> secondDerivativeError;
  bool firstDerivativeOk_, secondDerivativeOk_;
};
} // end of namespace Kaskade
#endif /* CHECK_DERIVATIVE_HH */
