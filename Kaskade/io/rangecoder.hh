/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef RANGECODER_HH
#define RANGECODER_HH

/**
 * \file
 * \brief Range coder for fast entropy coding.
 * \author Martin Weiser
 */

#include <algorithm>
#include <cassert>
#include <exception>
#include <ios>
#include <limits>
#include <numeric>
#include <vector>


/**
 * \ingroup entropycoding
 * \brief Base class for entropy coding with range encoder and decoder.
 *
 * The coding works on sequences \f$ (s_i)_i \f$ of symbols.  The
 * symbols \f$ s_i \f$ are elements of an alphabet \f$ S_i \f$, which
 * may be a different one for any position \f$ i \f$.
 *
 * Each symbol \f$ s \f$ is a pair of values \f$
 * (s_{\text{low}},s_{\text{high}}) \f$ with \f$ s_{\text{low}} <
 * s_{\text{high}} \f$. Symbols represent half-open integer
 * ranges. Symbols from the same alphabet are non-overlapping, i.e.,
 * \f$ s,t\in S \Rightarrow s_{\text{high}} \le t_{\text{low}} \vee
 * t_{\text{high}} \le s_{\text{low}} \f$.
 *
 * Each alphabet \f$ S \f$ is a setof nonoverlapping symbols such that
 * \f$ s_{\text{high}} \le \text{maxRange} \f$ for all \f$ s\in S \f$.
 *
 * Note that the range coder is not self-terminating. That means that
 * spurious trailing symbols can be extracted. For a correct
 * termination of the symbol sequence on decoding, you either need to
 * transmit the length of the symbol sequence beforehand or define a
 * special EOF symbol.
 *
 * The range coder algorithm is derived from the implementation by
 * Dmitry Subbotin. See also https://en.wikipedia.org/wiki/Range_encoding.
 *
 * \tparam UInt an unsigned integral type with more than 16 bits (usually 32 or 64 bit
 *         usigned integer).
 */
template <class UInt>
class RangeCoder
{
protected:
  // The number of binary digits in our unsigned integral type UInt.
  static int  const digits = std::numeric_limits<UInt>::digits;

  // Mask for shifting of ranges during normalization. All numbers
  // smaller than top have 8 zero leading bits.
  //
  // Smaller value of top, e.g. 1 << digits-16, would allow to shift
  // the range by a larger amount and output short's instead of char's
  // for improved performance. A dawback of this is that bottom needs
  // to be decreased as well, which leads to a small maxRange and
  // hence a low resolution of the alphabet's probability table.
  static UInt const top    = static_cast<UInt>(1) << (digits-8);

  // The smallest value of range that is possible. The class invariant
  // is range>=bottom. In contrast to top, the shift of bottom does
  // not need to be a multiple of 8. Smaller values of bottom lead to
  // a lower resolution of the alphabet's probability table (and hence
  // a worse approximation of symbol probabilities, in particular for
  // very different probabilities), however, it also saves space due
  // to fewer carry compensations.
  static UInt const bottom = static_cast<UInt>(1) << (digits-16);
  
public:
  /**
   * \brief Maximal total range of alphabets.
   */
  static const UInt maxRange = bottom;

  /**
   * \brief Number of processed encoded bytes.
   *
   * When encoding, this is the number of bytes already written to the
   * output stream. When decoding, this is the number of bytes already
   * read from the input stream.
   */
  size_t size() const {
    return count;
  }

protected:
  RangeCoder()
    : low(0), range(std::numeric_limits<UInt>::max()), count(0)
  {
    assert(std::numeric_limits<UInt>::is_specialized);
    assert(!std::numeric_limits<UInt>::is_signed);
    assert(digits > 16);
  }
  
  UInt low, range;
  size_t count;
};

//---------------------------------------------------------------------

/**
 * \ingroup entropycoding
 * \brief Entropy coding with range encoder.
 *
 * A range encoder. It encodes a sequence of symbols \f$ (s_i)_i \f$.
 * See base class RangeCoder for a definition of symbols and
 * alphabets.
 *
 * Note that the encoding is completed only after destruction of the
 * encoder, since the destructor writes the remaining characters to
 * the output stream.
 *
 * \tparam UInt an unsigned integral type (usually 32 or
 *         64 bit usigned integer).
 */
template <class UInt>
class RangeEncoder:public RangeCoder<UInt>
{
  static UInt const bottom = RangeCoder<UInt>::bottom;
  static UInt const top = RangeCoder<UInt>::top;
  static UInt const digits = RangeCoder<UInt>::digits;
  
public:
  static UInt const maxRange = RangeCoder<UInt>::maxRange;

  RangeEncoder(std::ostream& out_)
    : out(out_), fill(0)
  {}
  
  /**
   * \brief Destructor writes the remaining data to the output stream.
   */
  ~RangeEncoder() {
    for(int i=0; i<digits/8; i++) 
      put();
    out.write((char const*)buf,fill);
  }

  /**
   * \brief Encodes one symbol.
   * 
   * Encodes the symbol that covers the half-open range
   * [s.first,s.second[ of the totalRange. totalRange must not exceed
   * maxRange.  Characters may be written to the output stream. If the
   * write operation is not successful, this method throws an
   * exception.
   *
   * See also the convenience function \ref encodeSymbol().
   */
  void push(std::pair<UInt,UInt> s, UInt totalRange) 
  {
    // Check parameters.
    assert(s.first<s.second);
    assert(totalRange <= maxRange);
    
    UInt& range = this->range;
    UInt& low = this->low;

    // Check class invariant
    assert(range >= bottom);
    assert(range-1 <= std::numeric_limits<UInt>::max()-low);
    
    
    range /= totalRange;
    low += s.first*range;
    range *= s.second-s.first;

    // Perform range normalization
    while ( (low^(low+range)) < top
	    || (range<bottom && ((range= -low & (bottom-1)),true)) ) 
      put();       // Output of highest byte and shifting of range
    
    // Check class invariant
    assert(range >= bottom);
    assert(range-1 <= std::numeric_limits<UInt>::max()-low);
  }

private:
  std::ostream& out;

  // Putting one char at a time into the std::ostream is quite
  // inefficient. We store the bytes in a small buffer and feed the
  // buffer into the ostream once it is full. A buffer size of 64
  // seems to be a good compromise between performance improvement and
  // memory consumption.
  static int const bufsize = 64;
  unsigned char buf[bufsize];
  int fill;  

  void put() {
    UInt& range = this->range;
    UInt& low = this->low;

    if (fill==bufsize) {
      out.write((char const*)buf,bufsize);
      if (!out)
	throw std::ios_base::failure("write error");
      fill = 0;
      this->count += bufsize;
    }

    // Store eight most significant bits (top) in stream buffer.
    buf[fill] = low >> (digits-8);
    ++fill;
    
    // shift window
    range <<= 8;
    low <<= 8;
  }
};

//---------------------------------------------------------------------

/**
 * \ingroup entropycoding
 * \brief Entropy coding with range decoder.
 *
 * A range decoder. It decodes to a sequence of symbols \f$ (s_i)_i
 * \f$.  See base class RangeCoder for a definition of symbols and
 * alphabets.
 *
 * \tparam UInt an unsigned integral type (usually 32 or
 *         64 bit usigned integer).
 */
template <class UInt>
class RangeDecoder:public RangeCoder<UInt>
{
  static UInt const bottom = RangeCoder<UInt>::bottom;
  static UInt const top = RangeCoder<UInt>::top;
  static UInt const digits = RangeCoder<UInt>::digits;

  UInt code;
  std::istream& in;

public:
  static UInt const maxRange = RangeCoder<UInt>::maxRange;

  /**
   * Constructor. Note that reading the input stream starts already on
   * construction. Keep this in mind if you decide to transmit the
   * length of the symbol sequence before the encoded characters.
   */
  RangeDecoder(std::istream& input_)
    : code(0), in(input_)
  {
    // Throw ios_base::failure if end of stream is reached.
    in.exceptions(std::ios_base::eofbit | in.exceptions());
    
    // Fill range
    for(int i=0; i<digits/8; i++) {
      code = (code << 8) | in.get();
      if (!in) throw std::ios_base::failure("read error");
      ++this->count;
    }
  }

  /**
   * Returns a value \f$ c \f$ that is contained in the half-open
   * range of the current symbol \f$ s \f$: \f$ s_{\text{low}} \le c <
   * s_{\text{high}} \f$.
   * 
   * See also the convenience function decodeSymbol().
   */
  UInt front(UInt totalRange) const {
    UInt const& range = this->range;
    UInt const& low = this->low;

    return (code-low)/(range/totalRange);
  }

  /**
   * Removes the current symbol from the sequence. The given symbol
   * must mach the current symbol, i.e., s.first <= front(totalRange) <
   * s.second.
   *
   * Characters may be read from the input stream. If this operation
   * fails (e.g. because of end of stream is reached), an exception is
   * thrown.
   *
   * If the same number of symbols are pop'ed as have previously been
   * push'ed, the same number of encoded characters is read from the
   * input stream as has previously been written to the output
   * stream. Note that this does NOT imply that only as many symbols
   * can be pop'ed as have been push'ed before (it only means that at
   * least as many symbols can be decoded as have been
   * encoded). Spurious trailing symbols may be encountered before EOF
   * is reached.
   */
  void pop(std::pair<UInt,UInt> s, UInt totalRange) {
    UInt& range = this->range;
    UInt& low = this->low;

    assert(s.first<=front(totalRange)
	   && front(totalRange)<s.second);
    
    range /= totalRange;
    low += s.first*range;
    range *= s.second-s.first;

    // Perform range normalization
    while ( (low^(low+range))<top
	    || (range<bottom && ((range= -low & (bottom-1)),true)) ) {
      code = code<<8 | in.get();
      range <<= 8;
      low <<= 8;
      ++this->count;
    }
  }
};


//---------------------------------------------------------------------
//---------------------------------------------------------------------

/**
 * \ingroup entropycoding
 * \brief A simple alphabet of symbols with frequencies to be used with the range coder.
 *
 * The symbols are consecutive integers 0..(size()-1). The half-open
 * ranges associated to the symbols are derived from given frequencies
 * of the symbols.
 */
template <class UInt>
class Alphabet {
public:
  /**
   * Constructs an alphabet from a sequence of frequencies of
   * symbols. Note that for each symbol the frequency must be larger
   * than 0.
   */
  template <class InIter>
  Alphabet(InIter first, InIter last)
    : cumFreq(std::distance(first,last)+1) 
  {
    cumFreq[0] = 0;
    std::partial_sum(first,last,cumFreq.begin()+1,std::plus<UInt>());
  }

  /**
   * \brief Modifies the symbols' half-open ranges to represent the symbols'
   * frequencies given in the range [first,last[. 
   * 
   * Note that for each symbol the frequency must be larger than 0. The length of the
   * given sequence must match the size of the alphabet.
   */
  template <class InIter>
  void update(InIter first, InIter last) 
  {
    assert(std::distance(first,last)==size()); // check that all symbols are covered
    assert(std::find(first,last,0)==last);     // check that all frequencies are greater than zero
    std::partial_sum(first,last,cumFreq.begin()+1,std::plus<UInt>()); // compute cumulative distribution
  }

  /**
   * \brief Returns the total range, i.e. the maximum upper range bound of any symbol in the alphabet.
   */
  UInt totalRange() const 
  {
    return cumFreq.back();
  }

  /**
   * \brief Returns the number of symbols in the alphabet.
   */
  UInt size() const 
  {
    return cumFreq.size()-1;
  }

  /**
   * \brief Returns the symbol that contains the given value.
   */
  UInt symbol(UInt value) const 
  {
    typename std::vector<UInt>::const_iterator i
      = std::upper_bound(cumFreq.begin(),cumFreq.end(), value);
    assert(i==cumFreq.end() || *i>value);
    assert(i!=cumFreq.begin());
    return i-cumFreq.begin()-1;
  }

  /**
   *  \brief Returns the symbol's half-open cumulative frequency range.
   */
  std::pair<UInt,UInt> range(UInt symbol) const 
  {
    assert(symbol<size());
    return std::make_pair(cumFreq[symbol],cumFreq[symbol+1]);
  }
  

private:
  // cumulated frequencies for the symbols
  std::vector<UInt> cumFreq;
};

//---------------------------------------------------------------------
//---------------------------------------------------------------------

/**
 * \ingroup entropycoding
 * \brief A convenience function that encodes a symbol in a range encoder.
 * \relates RangeEncoder
 */
template <class UInt, class Symbol>
void encodeSymbol(RangeEncoder<UInt>& encoder, Alphabet<UInt> const& alphabet,
		  Symbol const& s) {
  encoder.push(alphabet.range(s),alphabet.totalRange());
}


/**
 * \ingroup entropycoding
 * \brief A convenience function that retrieves and pops the current symbol from a range decoder.
 * \relates RangeDecoder
 */
template <class UInt>
UInt decodeSymbol(RangeDecoder<UInt>& decoder, Alphabet<UInt> const& alphabet) {
  UInt s = alphabet.symbol(decoder.front(alphabet.totalRange()));
  decoder.pop(alphabet.range(s),alphabet.totalRange());
  return s;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------



#endif
