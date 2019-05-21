/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef VECTOR_TOOLS_HH_
#define VECTOR_TOOLS_HH_

#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>

namespace Kaskade
{
  template <class Arg, class PrintArg,
             class enable = typename std::enable_if<std::is_arithmetic<Arg>::value,void>::type>
  void printVector(std::vector<Arg> const& vec)
  {
    std::for_each(vec.begin(), vec.end(), [](Arg const& arg){ std::cout << arg << std::endl; });
  }

  template <class Arg, class PrintArg>
  void printVector(std::vector<Arg> const& vec, PrintArg printArg)
  {
    std::for_each(vec.begin(),vec.end(),printArg);
  }

  void printVector(std::vector<Dune::FieldVector<double,3> > const& vec)
  {
    printVector(vec, [](Dune::FieldVector<double,3> const& e) { std::cout << e << std::endl; });
  }

  void printVector(std::vector<std::vector<unsigned int> > const& vec)
  {
    printVector(vec, [](std::vector<unsigned int> const& t)
    {
      for(size_t i=0; i<t.size(); ++i) std::cout << t[i] << " ";
      std::cout << std::endl;
    });
  }

  /// Insert v into vec, if not yet existent in vec.
  /**
   * \param v element to be inserted
   * \param vec
   * \return if v was already an element of vec return its index, else return the index of the inserted v (i.e. vec.size()-1)
   */
  template <class Scalar, int dim>
  size_t addToVector(Dune::FieldVector<Scalar,dim> const& v, std::vector<Dune::FieldVector<Scalar,dim>>& vec)
  {
    auto tmp = std::find_if(vec.begin(), vec.end(), [&v](Dune::FieldVector<Scalar,dim> const& arg) -> bool
    {
      return ( (v-arg).two_norm() < 1e-6 );
    });
    if(tmp==vec.end())
    {
      vec.push_back(v);
      return vec.size()-1;
    }
    else return std::distance(vec.begin(),tmp);
  }
}

#endif
