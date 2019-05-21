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

#if !defined(FIXFUSION_H)
#define FIXFUSION_H

#include <boost/mpl/range_c.hpp>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

template <class Seq1, class Seq2>
boost::fusion::zip_view<boost::fusion::vector<Seq1&,Seq2&> > zip2(Seq1& s1, Seq2& s2)
{
  return  boost::fusion::zip_view< boost::fusion::vector<Seq1&,Seq2&> >(boost::fusion::vector<Seq1&,Seq2&>(s1,s2));
}

template <class Seq1, class Seq2, class Func>
void for_each2(Seq1 const& s1, Seq2& s2, Func const& func)
{
  using namespace boost::fusion;
  zip_view<vector<Seq1 const&,Seq2&> > zipped(vector<Seq1 const&,Seq2&>(s1,s2));
  for_each(zipped,func);
}

/**
 * \brief A half-open integer-range based for-loop.
 */
template <int first, int last, class Functor>
void for_range(Functor& f) {
  boost::fusion::for_each(typename boost::mpl::range_c<int,first,last>::type(),f);
}

/**
 * \brief Constructs an iterator range from two static boost::fusion iterators
 */
template <class First, class Last>
boost::fusion::iterator_range<First,Last> make_iterator_range(First first, Last last)
{
  return boost::fusion::iterator_range<First,Last>(first,last);
}

namespace boost {
  namespace fusion {
    namespace result_of {
      
      
      /**
       * \brief A meta-function defining a boost::fusion iterator range from given sequence indices
       */
      template <int First, int Last, class Sequence>
      struct make_range
      {
        typedef iterator_range<typename advance_c<typename begin<Sequence>::type,First>::type,
        typename advance_c<typename begin<Sequence>::type,Last>::type> type;
      };
      
      
    } // End of namespace result_of
  } // End of namespace fusion
} // End of namespace Kaskade

/**
 * \brief Constructs a boost::fusion view that picks a certain range from a forward sequence
 */
template <int First, int Last, class Sequence>
typename boost::fusion::result_of::make_range<First,Last,Sequence>::type 
make_range(Sequence& s)
{
  using namespace boost::fusion;

  return make_iterator_range(advance_c<First>(begin(s)),advance_c<Last>(begin(s)));
}


#endif
