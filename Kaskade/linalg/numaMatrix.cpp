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
#include <iostream>

#include "dune/common/fvector.hh"

#include "linalg/numaMatrix.hh"

using namespace Kaskade;

#ifdef UNITTEST 
#include <boost/timer/timer.hpp>
 
int main(void)
{ 
  using namespace std;
  
  std::cout << "8x8 matrix a_ij = i+j\n";
  NumaDenseMatrix<Dune::FieldVector<double,1>> mat(8,8);
  for (auto ri=mat.begin(); ri!=mat.end(); ++ri)
    for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
      *ci = ri.index()+ci.index();
    
  std::cout << mat << "\nFrobenius norm: " << mat.frobenius_norm() << "\n";
  
  std::cout << "Timing test of 100 additions of two 60000x100 matrices\n";
  NumaDenseMatrix<Dune::FieldMatrix<double,1,1>> m2(60000,100);
  for (auto ri=m2.begin(); ri!=m2.end(); ++ri)
    for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
      *ci = ri.index()+ci.index();
  boost::timer::cpu_timer timer;
  for (int i=0; i<100; ++i)
  {
    m2 += m2; // these two are a no-op
    m2 /= 2;
  }
  std::cout << timer.format() << "\n";
  // subtract original data - should be zero.
  for (auto ri=m2.begin(); ri!=m2.end(); ++ri) 
    for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
      *ci -= ri.index()+ci.index();
  std::cout << "Frobenius norm (should be 0): " << m2.frobenius_norm() << "\n";
  
  
  NumaVector<double> vec(12);
  for (auto i=begin(vec); i!=end(vec); ++i)
    *i = i.index();
  std::cout << "vec 12: " << vec << "\ndot product (should be 506): " << vec.dot(vec) << "\n2-norm (should be 22.4944): " << vec.two_norm() << "\n";
  
  vec.axpy(1.0,vec);
  std::cout << "twice vec: " << vec << "\n";
  
  return 0;
}

#endif
