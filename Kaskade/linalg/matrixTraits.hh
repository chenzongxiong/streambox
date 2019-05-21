/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LINALG_MATRIX_TRAITS
#define LINALG_MATRIX_TRAITS

// forward declarations
namespace Dune
{
  template <class Scalar, int n, int m>                          class FieldMatrix;
  template <class Scalar, int n>                                 class FieldVector;
  template <class Entry, class Allocator>                        class BlockVector;
}

/**
 * \file Provides forward declaration of NumaBCRSMatrix and a couple of traits for corresponding type computations.
 */

namespace Kaskade 
{
  // forward declaration
  template <class Entry, class Index> class NumaBCRSMatrix;

  // --------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic 
   * \brief Defines domain and range types for matrix classes.
   */
  template <class Matrix>
  struct MatrixTraits
  {
    /**
     * \brief The natural domain type.
     * 
     * This is the type of vectors that can be multiplied from the right to matrices, e.g. by the mv method or 
     * by operator*, depending on the matrix type.
     * 
     * There may be further supported vector types, but this one is the "natural" one.
     */
    using NaturalDomain = void;
    
    /**
     * \brief The natural range type.
     * 
     * This is the type of vectors that are the result of a matrix-vector multiplication with the natural
     * domain type. 
     */
    using NaturalRange = void;
  };
  
  template <class Scalar, int n, int m>
  struct MatrixTraits<Dune::FieldMatrix<Scalar,n,m>>
  {
    using NaturalDomain = Dune::FieldVector<Scalar,m>;
    using NaturalRange = Dune::FieldVector<Scalar,n>;
  };
  
  template <class Entry, class Index>
  struct MatrixTraits<NumaBCRSMatrix<Entry,Index>>
  {
    using NaturalDomain = Dune::BlockVector<typename MatrixTraits<Entry>::NaturalDomain>;
    using NaturalRange  = Dune::BlockVector<typename MatrixTraits<Entry>::NaturalRange>;
  };
}

#endif