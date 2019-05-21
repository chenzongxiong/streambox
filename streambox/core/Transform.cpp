/*
 * Transform.cpp
 *
 *  Created on: Aug 7, 2016
 *      Author: xzl
 *
 *  Purdue University, 2016
 */


#include "Transforms.h"
#include "WinKeyReducer/WinKeyReducer.h"
#include "WinSum/WinSumBase.h"

#if 0 /* moved to WindowKeyedReducer.cpp */
// a default, trivial reducer implementation.
// the the sum of all values (belonging to the same key)
template<>
std::pair<long, long> WindowKeyedReducer<std::pair<long, long>>::do_reduce
  (long const & key, ValueContainer<long> const & vcontainer) {

  VV("vcontainer has %lu elements...", vcontainer.vals.size());
  auto bcontainer_ptr = vcontainer.vals[0];
  VV("bcontainer0 %ld %ld ...", (*bcontainer_ptr)[0], (*bcontainer_ptr)[1]);

  auto & end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}

template<>
std::pair<string, long> WindowKeyedReducer<std::pair<string, long>>::do_reduce
  (string const & key, ValueContainer<long> const & vcontainer) {

  VV("vcontainer has %u elements...", vcontainer.vals.size());
  auto bcontainer_ptr = vcontainer.vals[0];
  VV("bcontainer0 %ld %ld ...", (*bcontainer_ptr)[0], (*bcontainer_ptr)[1]);

  auto & end = vcontainer.cend(); /* avoid calling it in each iteration */

  long sum = 0;
  for (auto it = vcontainer.cbegin(); it != end; ++it)
    sum += *it;

  return make_pair(key, sum);
}
#endif
