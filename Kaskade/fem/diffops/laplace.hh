/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef LAPLACE_OP_HH
#define LAPLACE_OP_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "fem/fixdune.hh"
#include "utilities/linalg/scalarproducts.hh"

namespace Kaskade {

  /**
   * \ingroup diffops
   * \brief Convenience class for handling diffusion terms in elliptic/parabolic equations.
   * 
   * This class defines the variational formulation of the diffusion operator \f$ -\mathrm{div}(\kappa \nabla \cdot)\f$
   * for scalar-valued functions, i.e. the variational functional \f$ u \mapsto \nabla \frac{1}{2}u^T \kappa \nabla u \f$.
   */
  template <class Scalar, int dim, bool isotropic=true>
  class ScalarLaplace {
  public:
    /**
     * \brief Constructor
     * The diffusion constant \f$ kappa \f$ is initialized to 1, the linearization point \f$u_0\f$ to \f$ \nabla u_0 = 0 \f$.
     */
    ScalarLaplace(): kappa(1), du0(0)
    {
    }
    
    /**
     * \brief Specifies a scalar diffusion constant.
     * \param kappa the diffusion constant (has to be positive and real)
     */
    Laplace& setDiffusionTensor(Scalar kappa_)
    {
      kappa = kappa_;
      assert(kappa > 0); // TODO: does this work for complex types? Maybe real(kappa) > 0 && imag(kappa) == 0?
      return *this;
    }
    
    /**
     * \brief Specifies the gradient of the linearization point of the Laplace variational functional.
     */
    Laplace& setLinearizationPoint(Dune::FieldVector<Scalar,dim> const& du0_)
    {
      du0 = du0_;
      return *this;
    }
    
    /**
     * \brief Returns the value of the Laplace variational functional.
     * The value returned is \f[ \frac{\kappa}{2} |\nabla u_0|^2, \f] where the gradient \f$ \nabla u_0 \f$ of the
     * linearization point \f$ u_0 \f$ is set via setLinearizationPoint.
     */
    Scalar d0() const
    {
      return 0.5*kappa*(du0 * du0);
    }
    
    Dune::FieldVector<Scalar,1> d1(VariationalArg<Scalar,dim> const& arg) const
    {
      return kappa*(arg.derivative*du0);
    }
    
    Dune::FieldMatrix<Scalar,1,1> d2(VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return kappa*(arg1.derivative[0]*arg2.derivative[0]);
    }
    
  private:
    Scalar kappa;                        // diffusion constant
    Dune::FieldVector<Scalar,dim> du0;   // gradient of linearization point
  };

  /**
   * \ingroup diffops
   * \brief Convenience class for handling diffusion terms in elliptic/parabolic equations.
   *
   * This class defines the variational formulation of the diffusion operator \f$ -\mathrm{div}(\kappa \nabla \cdot)\f$
   * for scalar-valued functions, i.e. the variational functional \f$ u \mapsto \nabla \frac{1}{2}u^T \kappa \nabla u \f$.
   */
  template <class Scalar, int dim, int components, bool isotropic=true>
  class Laplace {
  public:
    /**
     * \brief Constructor
     * The diffusion constant \f$ kappa \f$ is initialized to 1, the linearization point \f$u_0\f$ to \f$ \nabla u_0 = 0 \f$.
     */
    Laplace(): kappa(1), du0(0)
    {
    }

    void setDiffusionTensor(Scalar kappa_)
    {
      // TODO: shall we enforce kappa > 0? And kappa real?
      kappa = kappa_;
    }

    void setLinearizationPoint(Dune::FieldMatrix<Scalar,components,dim> const& du0_)
    {
      du0 = du0_;
    }

    Scalar d0() const
    {
      return 0.5*kappa*sp(du0, du0);
    }

    Scalar d1(VariationalArg<Scalar,dim,components> const& arg) const
    {
      return kappa*sp(arg.derivative, du0);
    }

    Scalar d2(VariationalArg<Scalar,dim,components> const &arg1, VariationalArg<Scalar,dim,components> const &arg2) const
    {
      return kappa*sp(arg1.derivative, arg2.derivative);
    }

  private:
    // diffusion constant
    Scalar kappa;
    Dune::FieldMatrix<Scalar,components,dim> du0;
    LinAlg::EuclideanScalarProduct sp;
  };

  // partial specialization for anisotropic case
  template <class Scalar, int dim>
  class ScalarLaplace<Scalar,dim,false>
  {
  public:
    /**
     * \brief Constructor
     * The diffusion constant \f$ kappa \f$ is initialized to 1, the linearization point \f$u_0\f$ to \f$ \nabla u_0 = 0 \f$.
     */
    ScalarLaplace(): kappa(unitMatrix<Scalar,dim,dim>()), du0(0)
    {
    }
    
    void setDiffusionTensor(Dune::FieldMatrix<Scalar,dim,dim> kappa_)
    {
      kappa = kappa_;
      // TODO: shall we enforce kappa spd?
      assert(kappa==transpose(kappa));
    }
    
    void setLinearizationPoint(Dune::FieldVector<Scalar,dim> const& du0_)
    {
      du0 = du0_;
    }
    
    Scalar d0() const
    {
      return 0.5*(du0 * (kappa*du0));
    }
    
    Dune::FieldVector<Scalar,1> d1(VariationalArg<Scalar,dim> const& arg) const
    {
      return arg.derivative*(kappa*du0);
    }
    
    Dune::FieldMatrix<Scalar,1,1> d2(VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return arg1.derivative*(kappa*arg2.derivative[0]);
    }
    
  private:
    Dune::FieldMatrix<Scalar,dim,dim> kappa;
    Dune::FieldVector<Scalar,dim> du0;
  };

}

#endif
