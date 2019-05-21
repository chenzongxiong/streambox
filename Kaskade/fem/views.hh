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

#ifndef UTILITIES_VIEWS_HH
#define UTILITIES_VIEWS_HH

#include <iterator>

namespace Kaskade
{
/**
 * \brief DEPRECATED. Use boost::iterator_range instead.
 */
template <class It>
class RangeView 
{
public:
  typedef It iterator;
  typedef It const_iterator;
  typedef typename std::iterator_traits<It>::value_type value_type;



  RangeView(It first_, It last_): first(first_), last(last_) {}

  const iterator begin() const { return first; }
  const iterator end()   const { return last; }

  bool empty() const { return first==last; }
  size_t size() const { return std::distance(first,last); }
  
  

  value_type const& operator[](int i) const 
  {
    iterator it = first;
    std::advance(it,i);
    return *it;
  }

private:
  iterator first, last;
};

  
/**
 * \brief Convenience function for constructing range views on the fly.
 */
template <class It>
RangeView<It> rangeView(It first, It last) { return RangeView<It>(first,last); }
} // end of namespace Kaskade

#endif
