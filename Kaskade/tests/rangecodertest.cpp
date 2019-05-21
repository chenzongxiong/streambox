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

#include <cassert>
#include <cctype>
#include <iostream>

#include <boost/timer.hpp>

#include "rangecoder.hh"

int main(int argc, char const* argv[]) {
  typedef unsigned long UInt;
  
  std::vector<UInt> count(256,1);
  Alphabet<UInt> alphabet(count.begin(),count.end());

  assert(argc>1);
  
  
  if (argv[1][0]=='e') {
    RangeEncoder<UInt> encoder(std::cout);
    int i = 0;
    while (std::cin) {
      unsigned char s = std::cin.get();
      encodeSymbol(encoder,alphabet,s);
      
      ++count[s];
      ++i;
      if (i>2*alphabet.size()) {
	alphabet.update(count.begin(),count.end());
	i=0;
      }
    }
  } else {
    RangeDecoder<UInt> decoder(std::cin);
    try {
      int i = 0;
      while (std::cin) {
	UInt s = decodeSymbol(decoder,alphabet);
	std::cout.put(static_cast<unsigned char>(s));

	count[s]++;
	++i;
	if (i>2*alphabet.size()) {
	  alphabet.update(count.begin(),count.end());
	  i=0;
	}
      }
    } catch (...) {
    }
    
  }

  return 0;
}

