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

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/timer.hpp>

#include "rangecoder.hh"

template <class UInt>
void encode(std::vector<double> const& data,
	    Alphabet<UInt> const& alphabet,
	    std::ostream& out) {
  UInt const N=alphabet.size();
  RangeEncoder<UInt> encoder(out);

  for (size_t i=0; i<data.size(); ++i) {
    UInt s = static_cast<UInt>(floor(N*data[i]));
    s = std::min(s,N-1);
    encodeSymbol(encoder,alphabet,s);
  }
  // std::cerr << "compression rate: "
  //           << (double)data.size()*sizeof(double)/encoder.size()
  //           << '\n';
  // std::cerr.flush();
}

template <class UInt>
void basic(std::vector<double> const& data,
	    std::ostream& out) {
  out.write((char const*)(&data[0]),data.size()*sizeof(double));
}

 
template <class UInt>
double test(size_t n, char const* fname, bool renc=true) {
  std::vector<double> data(n);
  for (size_t i=0; i<n; ++i) {
    data[i] = (random()%1000000)/1000000.0;
    data[i] = std::pow(data[i],3.0);
  }

  UInt N = 1000;
  std::vector<UInt> count(N);
  for (size_t i=0; i<N; ++i) 
    count[i] = (UInt) floor(20*std::pow((i+.5)/N,-0.6667));

  
  Alphabet<UInt> alphabet(count.begin(),count.end());
  

  std::ofstream out(fname,std::ios::binary);  
  boost::timer timer;
  if (renc)
    encode(data,alphabet,out);
  else
    basic<UInt>(data,out);
  return timer.elapsed();
}

  


int main(int argc, char const* argv[]) {
  size_t n=100000;
  for (int i=0; i<45; ++i) {
    std::cout << n << ' ' << test<unsigned long>(n,"/dev/null")
	      << ' ' << test<unsigned long>(n,"/tmp/nix2")
	      << ' ' << test<unsigned long>(n,"/tmp/nix2",false)
	      << '\n';
    std::cout.flush();
    
    n = 6*n/5;
  }
  
    
  return 0;
}

