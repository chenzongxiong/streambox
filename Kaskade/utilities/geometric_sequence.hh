/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GEOMETRIC_SEQUENCE_HH_
#define GEOMETRIC_SEQUENCE_HH_

#include <cmath>
#include <numeric>
#include <iterator>
#include <utility>

#include "utilities/power.hh"

/**
 * \brief Estimates parameters of geometric sequence.
 * 
 * Given a (finite) sequence \f$ a_i, i=0,\dots , n-1 \f$, this function estimates the parameters \f$ c,q\f$ of a
 * geometric sequence \f$ c q^i\f$, such that the least squares error \f$ \sum_{i=0}^{n-1} (ln(cq^i) - \ln a_i)^2 \f$ is
 * minimal. \f$ n \f$ has to be at least 2. The values \f$ a_i \f$ have to be positive.
 * 
 * \tparam Iterator an input iterator type with a scalar floating point value type
 */
template <class Iterator>
std::pair<typename std::iterator_traits<Iterator>::value_type,typename std::iterator_traits<Iterator>::value_type>
estimateGeometricSequence(Iterator first, Iterator last) {
  
  typedef typename std::iterator_traits<Iterator>::value_type Scalar;
  
  // Form the normal equations A^T A x = A^T b. As A = [ ones(n,1) (0:n-1)^T ], we have
  //          /n   n(n-1)/2              \          / ln(a_0)+...+ln(a_{n-1})               \   .
  // A^T A =  |                          |  and b = |                                       | 
  //          \ n(n-1)/2  n(n-1)(2n-1)/6 /          \ ln(a_1)+2ln(a_2)+...+(n-1)ln(a_{n-1}) /   such that
  //
  //                   2   / 2n-1    -3    \   .
  // (A^T A)^{-1} = ------ |               |
  //                n(n+1) \  -3   6/(n-1) /   .
  
  // compute right hand side
  Scalar b[2] = {0.0, 0.0};
  int n = 0;
  for ( ; first!=last; ++first) {
    Scalar lna = std::log(*first);
    b[0] += lna;
    b[1] += n*lna;
    ++n;
  }
  assert(n>=2);
  
  // solve system by multiplication with (A^T A)^{-1}
  return std::make_pair(std::exp(2*((2*n-1)*b[0]-3*b[1])/n/(n+1)),
                        std::exp( 2*(-3*b[0]+6*b[1]/(n-1))/n/(n+1)));
}

/**
 * \brief Estimates errors of truncating geometrically convergent sequences.
 * 
 * Assume there is a geometrically convergent series \f$ (x_i)_{i\in N} \f$ with data available for \f$ i=0,\dots , n \f$. The 
 * differences satisfy \f$ a_i = x_{i+1} - x_i \approx c q^i \f$. Then we estimate the truncation error
 * \f$ x_k - x_\infty \approx \sum_{i=k}^{n-1} a_i + \frac{cq^n}{1-q}. \f$
 * This value is computed here, given the \f$ n \f$ difference values \f$ a_i \f$.
 * 
 * \tparam Iterator an input iterator with scalar floating point value type.
 * \param k defines for which iterate \f$ x_k \f$ the error is estimated. \f$ 0\le k \le n \f$.
 */
template <class Iterator>
typename std::iterator_traits<Iterator>::value_type estimateGeometricTruncationError(Iterator first, Iterator last, int k) {
  typedef typename std::iterator_traits<Iterator>::value_type Scalar;

  int const n = std::distance(first,last);
  assert(0<=k && k<=n);
  
  std::pair<Scalar,Scalar> cq = estimateGeometricSequence(first,last);
  if (cq.second<1) {
    std::advance(first,k);
    return std::accumulate(first,last,0.0) + cq.first*power(cq.second,n)/(1-cq.second);
  } else
    return std::numeric_limits<Scalar>::max(); // TODO: throw an exception?
}
  

#endif /* GEOMETRIC_SEQUENCE_HH_ */