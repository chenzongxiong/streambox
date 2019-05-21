/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef NUMAMATRIX_HH
#define NUMAMATRIX_HH

#include <algorithm>
#include <cassert>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>

#include "linalg/dynamicMatrix.hh"
#include "utilities/threading.hh"

namespace Kaskade {
  
  // forward declaration
  template <class Entry> class NumaDenseMatrix;
  
  /**
   * \cond internals
   */
  namespace NumaDenseMatrixDetail
  {
    
    template <class Entry>
    class ColumnIterator: public boost::iterator_facade<ColumnIterator<Entry>,Entry,boost::random_access_traversal_tag>
    {
    public:
      ColumnIterator() = default;
      
      ColumnIterator(Entry* base_, size_t stride_, size_t idx_=0)
      : base(base_), stride(stride_), idx(idx_) {}
      
      size_t index() const { return idx; }
      
    private:
      friend class boost::iterator_core_access;
      
      bool equal(ColumnIterator<Entry> const& p) const
      {
        return idx==p.idx;
      }
      
      void increment() { ++idx; }
      void decrement() { --idx; }
      void advance(typename ColumnIterator<Entry>::difference_type n) { idx += n; }
      template <class E>
      typename ColumnIterator<Entry>::difference_type distance_to(ColumnIterator<E> const& i) const { return i.idx - idx; }
      
      
      Entry& dereference() const { return base[idx*stride]; }
      
      Entry* base;
      size_t stride;
      size_t idx;
    };
    
    // ----------------------------------------------------------------------------------------------
    
    template <class Entry>
    class Row 
    {
    public:
      typedef ColumnIterator<Entry> iterator;
      typedef ColumnIterator<Entry const> const_iterator;
      
      Row() = default;
      
      Row(Entry* first_, size_t stride_, std::size_t cols_)
      : first(first_), stride(stride_), cols(cols_) {}
      
      iterator begin() { return iterator(first,stride,0); }
      iterator end()   { return iterator(first,stride,cols); }
      const_iterator begin() const { return const_iterator(first,stride,0); }
      const_iterator end()   const { return const_iterator(first,stride,cols); }
      
      Entry&       operator[](size_t c)       { return first[stride*c]; }
      Entry const& operator[](size_t c) const { return first[stride*c]; }
      
    private:
      Entry* first;
      size_t stride;
      size_t cols;
    };
    
    // ----------------------------------------------------------------------------------------------

    template <class Entry>
    class RowIterator: public boost::iterator_facade<RowIterator<Entry>,
                                                     typename std::conditional<std::is_const<Entry>::value,
                                                                               Row<typename std::remove_const<Entry>::type> const,
                                                                               Row<Entry>>::type,
                                                     boost::random_access_traversal_tag>
    {
      typedef typename std::remove_const<Entry>::type NonConstEntry;
      
      typedef typename std::conditional<std::is_const<Entry>::value,
                                        NumaDenseMatrix<NonConstEntry> const,
                                        NumaDenseMatrix<Entry>>::type Matrix;
      
      typedef boost::iterator_facade<RowIterator<Entry>,
                                    typename std::conditional<std::is_const<Entry>::value,
                                                              Row<NonConstEntry> const,
                                                              Row<Entry>>::type,
                                    boost::random_access_traversal_tag> Base;
      
      
    public:
      RowIterator() = default;
      
      RowIterator(Matrix& mat_, size_t idx_)
      : mat(&mat_) 
      {
        update(idx_);
      }
      
      size_t index() const { return idx; }
      
    private:
      friend class boost::iterator_core_access;
      using DifferenceType = Base::difference_type;
      
      bool equal(RowIterator<Entry> const& other) const { return idx==other.idx;}
      void increment() { update(idx+1); }
      void decrement() { update(idx-1); }
      void advance(typename RowIterator<Entry>::difference_type n) { update(idx+n); }
      template <class E> DifferenceType distance_to(RowIterator<E> const& i) const { return static_cast<DifferenceType>(i.idx) - idx; }
      typename Base::reference dereference() const { return row; } 
      
      // moves row as to represent row with index newidx.
      void update(size_t newidx)
      {
        idx = newidx;
        if (idx<mat->N()) // to not trip over the edge when arriving at end()
        {
          size_t chunk = uniformWeightRange(idx,mat->chunks.size(),mat->N());
          size_t rowStart = uniformWeightRangeStart(chunk,mat->chunks.size(),mat->N());
          size_t stride = uniformWeightRangeStart(chunk+1,mat->chunks.size(),mat->N()) - rowStart; // that's "LDA"
          row = typename Base::value_type(const_cast<NonConstEntry*>(&mat->chunks[chunk][idx-rowStart]),stride,mat->M());
        }
      }
      
      Matrix* mat;
      size_t idx; // row index
      mutable typename Base::value_type row; // dereference is const method, but can "point to" mutable row value type..
    };

    // ----------------------------------------------------------------------------------------------

    template <class T, int n>
    typename Dune::FieldTraits<T>::field_type frobenius_product(Dune::FieldVector<T,n> const& x, Dune::FieldVector<T,n> const& y) 
    { 
      typename Dune::FieldTraits<T>::field_type result = 0;
      for (int i=0; i<n; ++i)
        result += frobenius_product(x[i],y[i]);
      return result;
    }
    
    template <class T, int n, int m>
    typename Dune::FieldTraits<T>::field_type frobenius_product(Dune::FieldMatrix<T,n,m> const& x, Dune::FieldMatrix<T,n,m> const& y) 
    { 
      typename Dune::FieldTraits<T>::field_type result = 0;
      for (int i=0; i<n; ++i)
        for (int j=0; j<m; ++j)
          result += frobenius_product(x[i][j],y[i][j]);
      return result;
    }
    
    double frobenius_product(double x, double y) { return x*y; }
    float  frobenius_product(float x, float y)  { return x*y; }
    
    
    
    template <class T, int n>
    typename Dune::FieldTraits<typename Dune::FieldTraits<T>::field_type>::real_type frobenius_norm2(Dune::FieldVector<T,n> const& x) { return x.two_norm2(); }
    
    template <class T, int n, int m>
    typename Dune::FieldTraits<typename Dune::FieldTraits<T>::field_type>::real_type frobenius_norm2(Dune::FieldMatrix<T,n,m> const& x) { return x.frobenius_norm2(); }
    
    double frobenius_norm2(double x) { return x*x; }
    float  frobenius_norm2(float x)  { return x*x; }
    
    
    template <class T, int n>
    typename Dune::FieldTraits<typename Dune::FieldTraits<T>::field_type>::real_type one_norm(Dune::FieldVector<T,n> const& x) { return x.one_norm(); }
    
    template <class T, int n, int m>
    typename Dune::FieldTraits<typename Dune::FieldTraits<T>::field_type>::real_type one_norm(Dune::FieldMatrix<T,n,m> const& x) { return x.one_norm(); }
    
    double one_norm(double x) { return std::abs(x); }
    float  one_norm(float x)  { return std::abs(x); }
    
    template <class T, int n>
    void axpy(Dune::FieldVector<T,n>& y, typename Dune::FieldTraits<T>::field_type const& a, Dune::FieldVector<T,n> const& x) { y.axpy(a,x); }
    
    template <class T, int n, int m>
    void axpy(Dune::FieldMatrix<T,n,m>& y, typename Dune::FieldTraits<T>::field_type const& a, Dune::FieldMatrix<T,n,m> const& x) { y += a*x; }
    
    void axpy(double& y, double a, double x) { y += a*x; }
    void axpy(float& y,  float a,  float x) { y += a*x; }
    
    // ----------------------------------------------------------------------------------------------
    
    // forward declaration
    template <class Entry> class VectorIterator;

    /**
     * \ingroup linalgbasic
     * \brief A base class for dense matrix and vector classes.
     * 
     * This distributes the entries across several NUMA nodes by storing a contiguous block of rows (chunk) on each node.
     * In each chunk, the entries are stored column-major, i.e. as a LAPACK dense matrix.
     */
    template <class Entry>
    class NumaDenseBase
    {
    protected:
      typedef std::vector<Entry,NumaAllocator<Entry>> Chunk;  // the actual data storage type
      typedef NumaDenseBase<Entry> Self;                      
      
    public:
      typedef size_t      size_type;
      typedef Entry       value_type;
      
      /**
       * \brief the type of the scalar field.
       */
      typedef typename Dune::FieldTraits<Entry>::field_type     field_type;
      typedef typename Dune::FieldTraits<field_type>::real_type real_type;
      
      /**
       * \brief the type of matrix entries.
       */
      typedef Entry block_type;
      
      /**
       * \name Construction and assignment
       * @{
       */
      
      /**
       * \brief Creates a 0 x 0 matrix.
       */
      NumaDenseBase (): NumaDenseBase(0,0)
      {
      } 
      
      /**
       * \brief Copy constructor.
       */
      NumaDenseBase (NumaDenseBase const& other)
      : rows(other.rows), cols(other.cols), chunks(other.chunks.size())
      {
        parallelForNodes([&](int i, int n) // perform data copy on the respective
        {                                  // NUMA node
          chunks[i] = other.chunks[i];
        });
      }
      
      /**
       * \brief Move constructor.
       */
      NumaDenseBase(NumaDenseBase&& a) = default;
      
      /**
       * \brief Creates an r x c matrix with given entries.
       */
      NumaDenseBase (size_type r, size_type c, value_type v = value_type()) :
      rows(r), cols(c)
      {
        size_type n = NumaThreadPool::instance().nodes();
        chunks.reserve(n);
        for (size_type i=0; i<n; ++i)
          chunks.emplace_back((uniformWeightRangeStart(i+1,n,r)-uniformWeightRangeStart(i,n,r))*c,v,NumaAllocator<Entry>(i));
      }
      
      /**
       * \brief Copy assignment.
       */
      NumaDenseBase& operator=(NumaDenseBase const& other)
      {
        assert(chunks.size() == other.chunks.size()); // this can only be violated if the number of NUMA nodes changes...
        rows = other.rows;
        cols = other.cols;
        parallelForNodes([&](int i, int n) // perform data copy on the respective
        {                                  // NUMA node
          chunks[i] = other.chunks[i];
        });
      }
      
      /**
       * \brief Assignment from scalar. 
       */
      Self& operator=(field_type const& a)
      {
        for_each([&](Entry& e) { e = a; }); 
      }
      
      /// @}

    protected:
      /**
       * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
       * 
       * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::DynamicMatrix interface.
       */
      void resize (size_type r, size_type c)
      {
        // Check whether we need to change anything. Note that due to row blocking,
        // the check rows*cols==n*m is NOT sufficient (the distribution of rows on 
        // NUMA nodes may change nevertheless).
        if (rows==r && cols==c)
          return;
        
        rows = r;
        cols = c;
        for (int i=0; i<chunks.size(); ++i) // TODO: parallelize
          chunks[i].resize((uniformWeightRangeStart(i+1,chunks.size(),r)-uniformWeightRangeStart(i,chunks.size(),r))*c);
      }

    public:
      /**
       * \brief Return the number of rows. 
       */
      size_type N() const { return rows; }
      
    protected:
      /**
       * \brief Return the number of rows. 
       */
      size_type M() const { return cols; }
      
      
      /**
       * \brief scalar product
       */
      field_type dot(Self const& other) const
      {
        assert(rows==other.rows && cols==other.cols);
        std::vector<field_type> dots(chunks.size());
        parallelForNodes([&](int i, int n) 
        {
          assert(n==chunks.size());
          field_type sum = 0;
          for (size_t k=0; k<chunks[i].size(); ++k)
            sum += NumaDenseMatrixDetail::frobenius_product(chunks[i][k],other.chunks[i][k]);
          dots[i] = sum;
        });
        return std::accumulate(dots.begin(),dots.end(),0);
      }
      
      /**
       * \brief Performs a (parallel) accumulation procedure.
       * 
       * Note that this groups the operations in a different way than the std::accumulate, and hence 
       * gives only the same results if the operation F is associative.
       */
      template <class Result, class F>
      Result accumulate(F f, Result const& init) const
      {
        std::vector<Result> rs(chunks.size());
        parallelForNodes([&](int i, int n) 
        {
          assert(n==this->chunks.size());
          Result r = init;
          for (auto const& e: chunks[i])
            r = f(r,e);
          rs[i] = r;
        });
        return std::accumulate(rs.begin()+1,rs.end(),rs[0],f);
      }
      
      /**
       * \brief Applies the given functor to all matrix entries.
       */
      template <class F>
      void for_each(F const& f) 
      {
        parallelForNodes([&](int i, int n) 
        {
          assert(n==this->chunks.size());
          std::for_each(begin(chunks[i]),end(chunks[i]),f);
        });
      }
      
      /**
       * \brief Applies the given functor to pairs of corresponding matrix elements.
       */
      template <class F>
      void for_each2(F const& f, Self const& other) 
      {
        assert(rows==other.rows && cols==other.cols);
        parallelForNodes([&](int i, int n) 
        {
          assert(n==this->chunks.size());
          for (size_t j=0; j<chunks[i].size(); ++j)
            f(chunks[i][j],other.chunks[i][j]);
        });
      }
      
    private:
      friend class NumaDenseMatrixDetail::RowIterator<Entry>;
      friend class NumaDenseMatrixDetail::RowIterator<Entry const>;
      friend class NumaDenseMatrixDetail::VectorIterator<Entry>;
      friend class NumaDenseMatrixDetail::VectorIterator<Entry const>;
      
      size_type rows, cols;      // the matrix dimensions
      std::vector<Chunk> chunks; // for each NUMA node a separate array to store the entries
    };
    
  }
  /**
   * \endcond 
   */
  
  //-------------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic
   * \brief A dense matrix class tailored towards NUMA architectures.
   * 
   * The class satisfies the Dune ISTL recursive matrix interface similar to Dune::Matrix. The data storage, however,
   * is organized in such a way as to benefit from NUMA architectures: the rows of the matrix are distributed evenly 
   * over the nodes. This implies that a good load balancing can only be achieved for matrices with (significantly) 
   * more rows than NUMA nodes. Avoid flat, long matrices.
   * 
   * There is no guarantee about placement or ordering of matrix elements, in particular none about the 
   * contiguity of element storage.
   * 
   * \tparam Entry the entry type (usually a Dune::FieldMatrix<n,m,double>)
   * 
   * \TODO: Currently, the implementation is mainly intended to gather a bag of coefficient vectors in one place,
   * rather than to provide a fully-fledged matrix class.
   */
  template <class Entry>
  class NumaDenseMatrix: public NumaDenseMatrixDetail::NumaDenseBase<Entry>
  {
    typedef NumaDenseMatrix<Entry> Self;
    typedef NumaDenseMatrixDetail::NumaDenseBase<Entry> Base;
    
  public:
    typedef size_t      size_type;
    typedef Entry       value_type;
    
    /**
     * \brief the type of the scalar field.
     */
    typedef typename Dune::FieldTraits<Entry>::field_type         field_type;
    typedef typename Dune::FieldTraits<field_type>::real_type     real_type;
    
    typedef NumaDenseMatrixDetail::RowIterator<Entry>       RowIterator;
    typedef NumaDenseMatrixDetail::RowIterator<Entry const> ConstRowIterator;

    typedef typename RowIterator::value_type      row_type;
    typedef typename ConstRowIterator::value_type const_row_type;
    
    /**
     * \name Construction, assignment, and resize.
     * @{
     */
    
    /**
     * \brief Creates a 0 x 0 matrix.
     */
    NumaDenseMatrix () = default;

    /**
     * \brief Copy constructor.
     */
    NumaDenseMatrix (NumaDenseMatrix const& a) = default;

    /**
     * \brief Move constructor.
     */
    NumaDenseMatrix (NumaDenseMatrix&& a) = default;

    /**
     * \brief Creates an r x c matrix with given entries.
     */
    NumaDenseMatrix (size_type r, size_type c, value_type v = value_type()) 
    : Base(r,c,v)
    {}
    
    /**
     * \brief Copy assignment.
     */
    NumaDenseMatrix& operator=(NumaDenseMatrix const& a) = default;
    
    /**
     * \brief Assignment from scalar. 
     */
    Self& operator=(field_type const& a)
    {
      Base::operator=(a);
      return *this;
    }
     
    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     * 
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::DynamicMatrix interface.
     */
    void resize (size_type r, size_type c)
    {
      Base::resize(r,c);
    }

    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     * 
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::Matrix interface.
     */
    void setSize (size_type r, size_type c)
    {
      resize(r,c);
    }
    
    /// @}
    
    /**
     * \name Entry access
     * @{
     */
    
    /**
     * \brief Get iterator to first row. 
     */
    RowIterator begin() { return RowIterator(*this,0); }
    
    /**
     * \brief Get iterator to one beyond last row. 
     */
    RowIterator end() { return RowIterator(*this,this->N()); }
    
    /**
     * \brief Get const iterator to first row. 
     */
    ConstRowIterator begin() const { return ConstRowIterator(*this,0); }
    
    /**
     * \brief Get const iterator to one beyond last row. 
     */
    ConstRowIterator end() const { return ConstRowIterator(*this,this->N()); }
    
    /**
     * \brief The subscript operator. 
     * 
     * Rows do not exist as persistent data structures in the matrix. Hence the subscript operator
     * returns a row object (in fact a view onto data held by the matrix) by value.
     */
    row_type  operator[] (size_type row) { return *RowIterator(*this,row); }
    
    /**
     * \brief The const subscript operator. 
     * 
     * Rows do not exist as persistent data structures in the matrix. Hence the subscript operator
     * returns a row object (in fact a view onto data held by the matrix) by value.
     */
    const_row_type const operator[] (size_type row) const { return *ConstRowIterator(*this,row); }
    
    /// @}
    
    /**
     * \brief Return the number of columns. 
     */
    size_type M() const { return Base::M(); }
    
    /**
     * \name Linear algebra operations.
     * @{
     */
    
    /**
     * \brief Multiplication with a scalar. 
     */
    Self& operator*= (field_type const& a)
    {
      this->for_each([&](Entry& e) { e *= a; });
      return *this;
    }
    
    /**
     * \brief Division by a scalar.
     */
    Self& operator/= (field_type const& scalar)
    {
      return (*this) *= (1/scalar);
    }
    
    /**
     * \brief Matrix addition.
     * Both matrices have to have the same sizes.
     */
    Self& operator+= (Self const& b)
    {
      this->for_each2([](Entry& y, Entry const& x) { y += x; },b);
      return *this;
    }
   
    /**
     * \brief Matrix subtraction.
     */
    Self& operator-= (Self const& b)
    {
      this->for_each2([](Entry& y, Entry const& x) { y -= x; },b);
      return *this;
    }
    
    /// @}
    
    /**
     * \brief frobenius norm: sqrt(sum over squared frobenius norms of entries) 
     */
    real_type frobenius_norm () const { return std::sqrt(frobenius_norm2()); }
    
    /**
     * \brief frobenius norm squared: sum over squared frobenius norms of entries
     */
    real_type frobenius_norm2 () const
    {
      return this->accumulate([](real_type s, Entry const& e) { return s+NumaDenseMatrixDetail::frobenius_norm2(e); }, real_type(0)); 
    }
  };
  
  /**
   * \relates NumaDenseMatrix
   */
  template <class Entry>
  std::ostream& operator<<(std::ostream& out, NumaDenseMatrix<Entry> const& mat)
  {
    for (auto ri=mat.begin(); ri!=mat.end(); ++ri)
    {
      std::copy(ri->begin(),ri->end(),std::ostream_iterator<Entry>(out," "));
      out << '\n';
    }
    return out;
  }
  
  // --------------------------------------------------------------------------
  
  // forward declaration
  template <class Entry> class NumaVector;
  
  /**
   * \cond internals
   */
  namespace NumaDenseMatrixDetail
  {
    
    template <class Entry>
    class VectorIterator: public boost::iterator_facade<VectorIterator<Entry>,Entry,boost::random_access_traversal_tag>
    {
      typedef typename std::remove_const<Entry>::type NonConstEntry;
      
      typedef typename std::conditional<std::is_const<Entry>::value,
                                        NumaVector<NonConstEntry> const,
                                        NumaVector<Entry>>::type Vector;
    public:
      VectorIterator() = default;
      
      VectorIterator(Vector& vec_, typename Vector::size_type idx_)
      : vec(&vec_), idx(idx_) {}
      
      size_t index() const { return idx; }
      
    private:
      friend class boost::iterator_core_access;
      
      bool equal(VectorIterator<Entry> const& p) const
      {
        return idx==p.idx;
      }
      
      void increment() { ++idx; }
      void decrement() { --idx; }
      void advance(typename VectorIterator<Entry>::difference_type n) { idx += n; }
      template <class E>
      typename VectorIterator<Entry>::difference_type distance_to(VectorIterator<E> const& i) const { return i.idx - idx; }
      Entry& dereference() const 
      { 
        size_t chunk = uniformWeightRange(idx,vec->chunks.size(),vec->N());
        size_t chunkStart = uniformWeightRangeStart(chunk,vec->chunks.size(),vec->N());
        return vec->chunks[chunk][idx-chunkStart];
      }
      
      Vector* vec;
      typename Vector::size_type idx;
    };
  }
  /**
   * \endcond
   */
  
  /**
   * \ingroup linalgbasic
   * \brief A vector class tailored towards NUMA architectures.
   * This vector distributes its entries in blocks of approximately equal size to the different nodes of NUMA machines
   * in order to exploit the larger accumulated memory bandwidth for fast vector operations. This is most beneficial
   * for large vectors exceeding the cache of the CPUs. Modern processors (2015-01-01) have up to 6MB 3rd level cache,
   * so this is intended for vectors of more than 50,000 doubles, say.
   * 
   * \tparam Entry the type of vector entries (usually Dune::FieldVector<double,n>)
   */
  template <class Entry>
  class NumaVector: public NumaDenseMatrixDetail::NumaDenseBase<Entry>
  {
    typedef NumaVector<Entry> Self;
    typedef NumaDenseMatrixDetail::NumaDenseBase<Entry> Base;
    
   public:
    typedef typename Dune::FieldTraits<Entry>::field_type         field_type;
    typedef typename Dune::FieldTraits<field_type>::real_type     real_type;
    typedef typename Base::size_type                              size_type;
    typedef typename Base::value_type                             value_type;
    typedef NumaDenseMatrixDetail::VectorIterator<Entry>          iterator;
    typedef NumaDenseMatrixDetail::VectorIterator<Entry const>    const_iterator;
    

    /**
     * \name Construction, assignment, and shape changes.
     * @{
     */

    /**
     * \brief Creates a 0 vector.
     */
    NumaVector(): Base(0,1) 
    {} 

    /**
     * \brief Copy constructor.
     */
    NumaVector(NumaVector const& a) = default;

    /**
     * \brief Move constructor.
     */
    NumaVector (NumaVector&& a) = default;

    /**
     * \brief Creates an r vector with given entries.
     */
    NumaVector (size_type r, value_type v = value_type()) :
      Base(r,1,v)
    {}
    
    /**
     * \brief Copy assignment.
     */
    NumaVector& operator=(NumaVector const& a) = default;
    
    /**
     * \brief Resizes the vector to r entries, leaving the entries in an undefined state.
     * 
     * On resize, all iterators are invalidated. This method adheres to the Dune::BlockVector interface.
     */
    void resize (size_type r) { Base::resize(r,1); }
    
    /// @}
    
    /**
     * \name Entry access
     * @{
     */
    /**
     * \brief iterator to first entry. 
     */
    iterator begin() { return iterator(*this,0); }
    
    /**
     * \brief iterator to one beyond last entry. 
     */
    iterator end() { return iterator(*this,this->N()); }
    
    /**
     * \brief const iterator to first row. 
     */
    const_iterator begin() const { return const_iterator(*this,0); }
    
    /**
     * \brief const iterator to one beyond last entry. 
     */
    const_iterator end() const { return const_iterator(*this,this->N()); }
    
    /**
     * \brief The subscript operator. 
     */
    Entry&  operator[] (size_type row) { return *iterator(*this,row); }
    
    /**
     * \brief The const subscript operator. 
     */
    Entry const& operator[] (size_type row) const { return *const_iterator(*this,row); }
    
    /// @}
     
    /**
     * \name Linear algebra opterations.
     * The computations are performed in parallel on the NUMA nodes.
     * @{
     */
    
    /**
     * \brief \f$ y \leftarrow a y \f$
     */
    Self& operator*= (field_type const& a)
    {
      this->for_each([&](Entry& e) { e *= a; });
      return *this;
    }
    
    /**
     * \brief \f$ y \leftarrow \frac{1}{a} y \f$
     */
    Self& operator/= (field_type const& a)
    {
      return (*this) *= (1/a);
    }
    
    /**
     * \brief \f$ y \leftarrow y+x \f$
     * 
     * Both vectors have to have the same size.
     */
    Self& operator+= (Self const& x)
    {
      this->for_each2([](Entry& yi, Entry const& xi) { yi += xi; },x);
      return *this;
    }
   
    /**
     * \brief \f$ y \leftarrow y - x \f$
     * 
     * Both vectors have to have the same size.
     */
    Self& operator-= (Self const& x)
    {
      this->for_each2([](Entry& yi, Entry const& xi) { yi -= xi; },x);
      return *this;
    }

    /**
     * \brief \f$ y \leftarrow ax + y \f$
     */
    Self& axpy(field_type const& a, Self const& x)
    {
      this->for_each2([&](Entry& yi, Entry const& xi) { NumaDenseMatrixDetail::axpy(yi,a,xi); },x);
      return *this;
    }
    
    /// @}
    
    /**
     * \name Norms and scalar product.
     * 
     * The computations are performed in parallel on the NUMA nodes.
     * @{
     */
    
    /**
     * \brief euclidean norm: sqrt(sum over squared euclidean norms of entries) 
     */
    real_type two_norm () const { return std::sqrt(two_norm2()); }
    
    /**
     * \brief euclidean norm squared: sum over squared euclidean norms of entries
     */
    real_type two_norm2 () const 
    { 
      return this->accumulate([](real_type s, Entry const& e) { return s+NumaDenseMatrixDetail::frobenius_norm2(e); }, real_type(0)); 
    }
    
    /**
     * \brief one-norm: sum over absolute values of entries
     */
    real_type one_norm() const 
    { 
      return this->accumulate([](real_type s, Entry const& e) { return s+NumaDenseMatrixDetail::one_norm(e); }, real_type(0)); 
    }
    
    /**
     * \brief 
     */
    field_type dot(Self const& other) const { return Base::dot(other); }
    
    /// @}
  };
  
  /**
   * \relates NumaVector
   */
  template <class Entry>
  std::ostream& operator<<(std::ostream& out, NumaVector<Entry> const& vec)
  {
    using namespace std;
    copy(begin(vec),end(vec),ostream_iterator<Entry>(out," "));
    return out;
  }
}

#endif
