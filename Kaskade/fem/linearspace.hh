/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LINEARSPACE_HH
#define LINEARSPACE_HH

#include <cmath>
#include <type_traits>

#include <boost/version.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>

#include "linalg/crsutil.hh"

namespace Kaskade
{
  /// \internal 
  // forward declaration
  template <class> class VariableSet;

  namespace LinearSpace_Detail
  {

    template <class Space0, class Space2, int id> struct Copy;

    template <class Space, int id, class VSDescriptions>
    struct Copy<VariableSet<VSDescriptions>,Space,id>
    {
      static void apply(const VariableSet<VSDescriptions>& from, Space& to)
      {
        boost::fusion::at_c<id>(to.data) = boost::fusion::at_c<id>(from.data).coefficients();
        Copy<VariableSet<VSDescriptions>,Space,id-1>::apply(from,to);
      }
    };

    template <class Space, class VSDescriptions>
    struct Copy<VariableSet<VSDescriptions>,Space,0>
    {
      static void apply(const VariableSet<VSDescriptions>& from, Space& to)
      {
        boost::fusion::at_c<0>(to.data) = boost::fusion::at_c<0>(from.data).coefficients();
      }
    };
  }
  /// \endinternal 

  /**
   * \ingroup linalgbasic 
   * \brief A product space of linear spaces. 
   * 
   * Access to the components of the product space is provided by the the boost::fusion
   * sequence data:
   * \code
   * LinearProductSpace<double,boost::fusion::vector<Dune::FieldVector<double,2>,Dune::FieldVector<double,3>>> x;
   * Dune::FieldVector<double,3>& comp1 = boost::fusion::at_c<1>(x.data);
   * \endcode
   * 
   * Arithmetic operations are delegated componentwise to the heterogeneous subspaces.
   * 
   *
   * \tparam Scalar the scalar field type of the linear space
   * \tparam Seq a boost::fusion sequence type defining the variables of the cartesian product
   */
  template <class Scalar_, class Seq>
  class LinearProductSpace
  {
  private:
    typedef LinearProductSpace<Scalar_,Seq> Self;

    // forward declarations (definitions at end of file)
    struct Add;
    struct Sub;
    struct Assign;
    struct ScalarMult;
    struct Axpy;
    struct ScalarProduct;
    struct ScalarProductHack;
    template <typename> struct ReadBlock;
    template <typename> struct WriteBlock;

  public:
    /// scalar type
    typedef Scalar_ Scalar;
    typedef Scalar field_type;
    
    /// boost::fusion::vector of element vectors.
    typedef Seq Sequence;
    typedef Seq Functions;


    // use with care: default constructed subvectors may not be properly initialized (e.g. size?)
    // LinearProductSpace() {}
    LinearProductSpace() = delete;

    /// Copy constructor
    LinearProductSpace(LinearProductSpace<Scalar,Sequence> const& y): data(y.data) {}

    /// Copy from VariableSet
    template <class VSDescriptions>
    explicit LinearProductSpace(const VariableSet<VSDescriptions>& y)
    {
      *this = y;
    }

    /// Constructor with explicit data
    template <class S>
    explicit LinearProductSpace(S const& init): data(init)
    {}

    template <typename... Args>
    explicit LinearProductSpace(boost::fusion::transform_view<Args...> const& data_) : data(data_)
    {
      applyUnaryOp(Assign(0));
    }

    /**
     *\brief Number of scalar degrees of freedom.
     * This computes the total number of degrees of freedom, i.e. the sum of the dimensions of 
     * the components.
     */
    size_t dim() const
    {
      return boost::fusion::accumulate(data,0,[](auto dim, auto const& v) { return dim+v.dim(); });
    }


    /// Assignment from the same type
    Self& operator=(Self const& y) { if (this!=&y) data = y.data; return *this; }

    /// Assignment from a different (hopefully compatible) type
    template <class VSDescriptions>
    Self& operator=(VariableSet<VSDescriptions> const& y)
    {
      using namespace boost::fusion;
      static_assert(result_of::size<Sequence>::type::value == VSDescriptions::noOfVariables,"Numbers of variables do not match.");
      LinearSpace_Detail::Copy<VariableSet<VSDescriptions>,LinearProductSpace,
                               result_of::size<typename VSDescriptions::RepresentationData>::type::value-1>::apply(y,*this);
      return *this;
    } // need not check for self assignment here

    /// Assignment from a different (hopefully compatible) type
    template <class OtherSeq>
    Self& operator=(LinearProductSpace<Scalar,OtherSeq> const& y) 
    { 
      using namespace boost::fusion;
      static_assert(result_of::size<Sequence>::type::value == result_of::size<OtherSeq>::type::value, "Numbers of variables do not match.");
      data = y.data;
      return *this; 
    } // need not check for self assignment here

    /// Assignment from a (hopefully compatible) FunctionSpaceElement
    template <class FunctionSpace, int m>
    Self& operator=(FunctionSpaceElement<FunctionSpace,m> const& fse) 
    { 
      data = fse.coefficients(); 
      return *this; 
    }

    /// Assignment to constant scalar
    Self& operator=(Scalar a) 
    { 
      return applyUnaryOp(Assign(a));
    }
    
    /// Assignment from a (hopefully compatible) boost::fusion iterator range 
    template <class First, class Last>
    Self& operator=(boost::fusion::iterator_range<First,Last> range)
    {
      static_assert(boost::fusion::result_of::size<Sequence>::type::value == boost::fusion::result_of::distance<First,Last>::type::value,
                    "Numbers of variables do not match.");
      data = range;
      return *this;
    }
    
    //
    
    /// Scaling
    Self& operator*=(Scalar a) { return applyUnaryOp(ScalarMult(a)); }

    /// In place addition
    template <class SequenceY>
    Self& operator+=(LinearProductSpace<Scalar,SequenceY> const& y) {
      return applyBinaryOp(Add(),y); }

    /// In place subtraction
    Self& operator-=(Self const& y) { return applyBinaryOp(Sub(),y); }

    /// this <- this + a*y
    Self& axpy(Scalar a, Self const& y) { return applyBinaryOp(Axpy(a),y); }

    /// Scalar product
    field_type operator*(Self const& y) const {
      // For some unknown reason, this does not work... (will work in more recent fusion version).
      // return boost::fusion::accumulate(boost::fusion::zip(data,y.data),0,ScalarProduct());
      // We therefore use the following ugly hack workaround
      using namespace boost::fusion;
      zip_view<vector<Sequence const&,Sequence const&> > zipped(vector<Sequence const&,Sequence const&>(data,y.data));
      ScalarProductHack sp;
      for_each(zipped,sp);
      return sp.s;
    }
    
    /// Scalar product
    field_type dot(Self const& y) const { return *this * y; } 

    /**
     * \brief DEPRECATED use vectorFromSequence instead
     * 
     * Reads the coefficients sequentially from an input iterator. 
     * The InIterator's value type must be convertible to field_type.
     */
    template <class InIterator>
    void read(InIterator i) 
    { 
      vectorFromSequence(*this,i);
    }

    /**
     * \brief Reads the coefficients of the subrange [rbegin,rend[ sequentially from an input iterator. 
     * The InIterator's value type must be convertible to field_type.
     */
    template <class InIterator>
    void read(int rbegin, int rend, InIterator i) { boost::fusion::for_each(data,ReadBlock<InIterator>(rbegin,rend,i)); }

    /**
     * \brief DEPRECATED use vectorToSequence instead
     * 
     * Writes the coefficients sequentially to an output iterator (flattening).
     * The field_type must be convertible to the OutIterator's value type.
     */
    template <class OutIterator>
    void write(OutIterator i) const 
    { 
      vectorToSequence(*this,i);
    }

    /**
     * \brief Writes the coefficients of the subrange [rbegin,rend[ sequentially to an output iterator. 
     * The InIterator's value type must be convertible to field_type.
     */
    template <class DataOutIter>
    void write(int rbegin, int rend, DataOutIter i) const
    {
      boost::fusion::for_each(data,WriteBlock<DataOutIter>(rbegin,rend,i));
    }

    /// Euclidean Norm
    double two_norm() const { return std::sqrt((*this)*(*this)); }
    
    /// Data
    Sequence data;

  private:

    // applies a binary operator to the elements of a pair
    template <class BinaryOp>
    struct PairOp
    {
      PairOp(BinaryOp const& op_): op(op_) {}

      template <class Pair> void operator()(Pair p) const {
        op(boost::fusion::at_c<0>(p),boost::fusion::at_c<1>(p));
      }

    private:
      BinaryOp const& op;
    };

    template <class UnaryOp>
    Self& applyUnaryOp(UnaryOp const& op) {
      boost::fusion::for_each(data,op);
      return *this;
    }

    template <class BinaryOp, class SpaceY>
    Self& applyBinaryOp(BinaryOp const& op, SpaceY const& y) {
      using namespace boost::fusion;
      typedef typename SpaceY::Sequence SequenceY;
      zip_view<vector<Sequence&,SequenceY const&> > zipped(vector<Sequence&,SequenceY const&>(data,y.data));
      for_each(zipped,PairOp<BinaryOp>(op));
      return *this;
    }

    struct Add { template <class T, class S> void operator()(T& a, S const& b) const { a += b; } };
    struct Sub { template <class T> void operator()(T& a, T const& b) const { a -= b; } };
    struct Assign {
      Assign(Scalar s_): s(s_) {}
      template <class T> void operator()(T& a) const { a = s; }
    private: Scalar s;
    };


    struct ScalarMult {
      ScalarMult(Scalar s_): s(s_) {}
      template <class Element> void operator()(Element& e) const { e *= s; }
    private: Scalar s;
    };

    struct Axpy
    {
      Axpy(Scalar s_): s(s_) {}
      template <class T> void operator()(T& a, T const& b) const {

        // Funny runtime error could only be avoided by this funny code.
        // Probably a compiler bug.
        // assert(a.size()==b.size());
        std::cout << "";


        a.axpy(s,b);
      }

    private: Scalar s;
    };

    struct ScalarProduct
    {
      template <class T> struct result {};

      template <class Pair, class T> struct result<ScalarProduct(Pair,T)> { typedef field_type type; };
      template <class Pair>
      field_type operator()(Pair const& pair, field_type res) const {
        using namespace boost::fusion;

        return res ;//+ at_c<0>(pair)*at_c<1>(pair);
      }
    };

    // Ugly hack -- see above
    struct ScalarProductHack
    {
      ScalarProductHack(): s(0) {}
      template <class Pair> void operator()(Pair const& p) const { s += boost::fusion::at_c<0>(p) * boost::fusion::at_c<1>(p); }
      mutable field_type s;
    };

    template <class DataOutIter>
    struct WriteBlock
    {


      WriteBlock(int& rbegin_, int& rend_, DataOutIter& out_): rbegin(rbegin_), rend(rend_), out(out_) {}
      template <class VectorBlock> void operator()(VectorBlock const& v) const
      {
        if (rbegin<=0 && rend>0)
          vectorToSequence(v,out);
        --rbegin;
        --rend;
      }
    private:
      int& rbegin;
      int& rend;
      DataOutIter& out;

    };

    template <class DataInIter>
    struct ReadBlock
    {
      ReadBlock(int& rbegin_, int& rend_, DataInIter& in_): rbegin(rbegin_), rend(rend_), in(in_) {}
      
      template <class VectorBlock> void operator()(VectorBlock& v) const
      {
        if (rbegin<=0 && rend>0)
          vectorFromSequence(v,in);
        --rbegin;
        --rend;
      }
    private:
      int& rbegin;
      int& rend;
      DataInIter& in;

    };

  };
  
  /**
   * \ingroup linalgbasic
   * \brief Provides access to the m-th component of a product space.
   * 
   * This simplifies the access to individual components of a product space, e.g., a variable set of finite
   * element functions. Instead of 
   * \code
   * boost::fusion::at_c<m>(x.data)
   * \endcode
   * one can write
   * \code
   * component<m>(x)
   * \endcode
   * 
   * \tparam m the index of the component (nonnegative)
   * \tparam Scalar a scalar type
   * \tparam Sequence a boost::fusion sequence type
   * 
   * \return a reference to the m-th component
   * 
   * \relates LinearProductSpace
   */
  template <int m, class Scalar, class Sequence>
  typename boost::fusion::result_of::at_c<Sequence const,m>::type component(LinearProductSpace<Scalar,Sequence> const& x) { return boost::fusion::at_c<m>(x.data); }

  template <int m, class Scalar, class Sequence>
  typename boost::fusion::result_of::at_c<Sequence,m>::type component(LinearProductSpace<Scalar,Sequence>& x) { return boost::fusion::at_c<m>(x.data); }
  

  /**
   * \ingroup linalgbasic
   * \brief writes the coefficients of a vector to a flat scalar sequence
   * \related LinearProductSpace<Scalar,Seq>
   */
  template <class Scalar, class Seq, class OutIter>
  OutIter vectorToSequence(LinearProductSpace<Scalar,Seq> const& v, OutIter i)
  {
    boost::fusion::for_each(v.data,[&](auto& x) { i = vectorToSequence(x,i); }); 
    return i;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief reads the coefficients of a vector from a flat scalar sequence
   * \related LinearProductSpace<Scalar,Seq>
   */
  template <class Scalar, class Seq, class InIter>
  InIter vectorFromSequence(LinearProductSpace<Scalar,Seq>& v, InIter i)
  {
    boost::fusion::for_each(v.data,[&](auto& x) { i = vectorFromSequence(x,i); }); 
    return i;
  }
  
} // end of namespace Kaskade


#endif
