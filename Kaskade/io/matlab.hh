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

#ifndef MATLAB_HH
#define MATLAB_HH

#include <fstream>
#include <string>

#include "linalg/triplet.hh"

namespace Kaskade
{
  /**
   * \ingroup IO
   *
   * Writes the assembled matrix and right hand side of a
   * VariationalFunctionalAssembler to the given file. The contents is a
   * Matlab function that can be executed. The file extension .m is
   * appended automatically.
   *
   * \param[in] assembler a VariationalFunctionalAssembler containing the matrix and the right hand side
   * \param[in] basename the file name to be written to. A ".m" suffix is appended automatically.
   * \param[in] precision the number of valid decimal digits written for each scalar
   */
  template <class Assembler>
  void writeToMatlab(Assembler const& assembler, std::string const& basename, int precision=10)
  {
    std::string fname = basename + ".m";
    std::ofstream f(fname.c_str());
    f.precision(precision);

    f << "function [A,b] = " << basename << '\n'
        << " b = [\n";

    assembler.toSequence(0,Assembler::TestVariableSet::noOfVariables,
        std::ostream_iterator<typename Assembler::Scalar>(f,"\n"));
    f << "\n];\ndata = [\n";

    typedef MatrixAsTriplet<typename Assembler::Scalar> Matrix;

    Matrix A = assembler.template get<Matrix>(false);

    for (size_t i=0; i<A.ridx.size(); ++i)
      f << A.ridx[i]+1 << ' ' << A.cidx[i]+1 << ' ' << A.data[i] << '\n';
    f << "];\nA = spconvert(data);\n";
  }
}


#endif
