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

#ifndef CONVERTER_HH
#define CONVERTER_HH

#include <cassert>

#include "dune/common/fmatrix.hh"

#include "fem/fixdune.hh"
#include "fem/variables.hh"

namespace Kaskade
{
  /**
   * \brief A Converter for scalar shape functions that do not change their
   * value under transformation from reference to actual element geometry. 
   * 
   * If grid and world dimension are the same, the restriction of ansatz functions \f$ \varphi \f$ on the actual element are defined 
   * in terms of the shape functions \f$ \phi \f$ on the reference elements as
   * \f[ \varphi(x) = \phi(B^{-1}(x-b)), \quad \xi = B^{-1}(x-b). \f]
   * Therefore the derivatives and gradients are transformed as
   * \f[ \varphi'(x) = \phi'(\xi)B^{-1}, \quad \nabla \varphi(x) = B^{-T}\nabla\phi(\xi). \f]
   * Hessians (derivatives of gradients) are transformed as
   * \f[ H\varphi(x) = B^{-T}H_\phi(\xi)B^{-1}. \f]
   * This is used as a converter for scalar Lagrange and hierarchical spaces.
   * 
   * If the world dimension is larger than the grid dimension, the gradient is defined as
   * the gradient of the function extending as a constant in normal direction. This is 
   * just 
   * \f[ \varphi'(x) = \phi'(\xi)B^{+}, \quad \nabla \varphi(x) = B^{+T}\nabla\phi(\xi), \f]
   * where \f$ B \in {\bf R}^{w\times d}\f$ .
   * \todo: The Hessian depends on the curvature of the grid - this is not taken into account yet.
   *
   * Template parameters:
   * \tparam Cell the type of codim 0 entity on which the shape functions are defined.
   */
  template <class Cell>
  class ScalarConverter
  {
    static int const dim = Cell::dimension;
    static int const dimw = Cell::Geometry::coorddimension;

  public:
    ScalarConverter(): cell_(0) {}

    ScalarConverter(Cell const& cell): cell_(&cell) {}

    /**
     * \brief Indicates that following accesses refer to the given cell.
     */
    void moveTo(Cell const& cell) { cell_ = &cell; }

    void setLocalPosition(Dune::FieldVector<typename Cell::Geometry::ctype,dim> const& xi)
    {
      assert(cell_);
      inv = cell_->geometry().jacobianInverseTransposed(xi);
    }

    /// Applies the transformation \f$ \psi(x) \f$ to shape function value.
    template <class Scalar>
    Dune::FieldMatrix<Scalar,1,1> global(Dune::FieldMatrix<Scalar,1,1> const& sf) const { return sf; }

    /**
     * \brief Applies the transformation \f$ \psi \f$ to shape function value, gradient, and Hessian.
     */
    template <class Scalar>
    VariationalArg<Scalar,dimw,1> global(VariationalArg<Scalar,dim,1> const& sf, int deriv=1) const
    {
      VariationalArg<Scalar,dimw,1> phi;
      phi.value = sf.value;
      if (deriv >= 1) phi.derivative[0] = inv * sf.derivative[0];
      if (deriv >= 2) phi.hessian[0] = inv * sf.hessian[0] * transpose(inv);
      return phi;
    }
    
    template <class Scalar>
    Dune::FieldMatrix<Scalar,1,1> local(Dune::FieldMatrix<Scalar,1,1> const& glob) const { return glob; }


  private:
    Cell const*                                      cell_;
    Dune::FieldMatrix<typename Cell::Geometry::ctype,dimw,dim> inv;
  };
} /* end of namespace Kaskade */


#endif
