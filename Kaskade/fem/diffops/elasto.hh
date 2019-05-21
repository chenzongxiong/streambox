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

#ifndef ELASTO_HH
#define ELASTO_HH

#include <cassert>

#include "dune/common/fmatrix.hh"
#include "dune/grid/config.h"


#include "fem/fixdune.hh"
#include "fem/integration.hh"

#include "linalg/determinant.hh"

#include "utilities/power.hh"

namespace Kaskade 
{
  namespace Elastomechanics
  {
    /**
     * \ingroup diffopsElasto
     * \defgroup stationaryElasticity Stationary elasticity
     * \brief Building blocks for stationary elasticity problems (linear elasticity and hyperelasticity)
     */
    
    /**
     * \ingroup diffopsElasto
     * \defgroup viscoPlasticity Viscoplasticity
     * \brief Building blocks for viscoplastic problems
     */
    
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * \ingroup stationaryElasticity
     * \brief Material parameters for isotropic linearly elastic materials.
     */
    class ElasticModulus
    {
    public:
      /// Default constructor uses \f$ \lambda=\mu=1 \mathrm{Pa} \f$ (pretty useless "unit square" material).
      ElasticModulus(): lambda(1), mu(1) {}
      
      /// Constructor expects first and second Lame parameters \f$ \lambda, \mu \f$.
      ElasticModulus(double lambda_, double mu_): lambda(lambda_), mu(mu_) {}
      
      /// Define material properties in terms of Young's modulus \f$ E \f$ and Poisson's ratio \f$ \nu \f$.
      ElasticModulus& setYoungPoisson(double E, double nu);
      
      /// Returns the shear modulus \f$ G \f$, also known as Lame's second parameter \f$ \mu \f$.
      double shear() const { return mu; }
      
      /// Returns the bulk modulus \f$ K \f$.
      double bulk() const { return lambda + 2*mu/3; }
      
      /// Returns the Young's modulus \f$ E \f$.
      double young() const { return mu*(3*lambda+2*mu)/(lambda+mu); }
      
      /// Returns the first Lame parameter \f$ \lambda \f$.
      double lame() const { return lambda; }
      
      /// Returns Poisson's ratio \f$ \nu \f$.
      double poisson() const { return lambda/2/(lambda+mu); }
      
      /**
       * \brief Returns the material parameters of the given materials.
       * 
       * Throws a LookupException if the material is not available in the data base.
       */
      static ElasticModulus const& material(std::string const& name);
      
      /**
       * \brief A map from known material names to elastic moduli.
       * 
       * Use this to inquire the list of known materials
       */
      static std::map<std::string,ElasticModulus> const& materials();
      
    private:
      double lambda, mu;
    };
    
    /**
     * \ingroup stationaryElasticity
     * \brief Mass density of several materials in kg/m^3.
     * 
     * Throws a LookupException if the material is not available in the data base.
     */
    double massDensity(std::string const& name);
    
    /**
     * \ingroup stationaryElasticity
     * \brief A map from known material names to their mass densities.
     * 
     * Use this to inquire the list of known materials
     */
    std::map<std::string,double> const& massDensities();
    
    /**
     * \ingroup viscoPlasticity
     * \brief Yield strengths of several materials in Pa.
     * 
     * Throws a LookupException if the material is not available in the data base.
     */
    double yieldStrength(std::string const& name);
    
    
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * \cond internals
     */
    namespace ElastoDetails
    {
      constexpr int dim(int n) { return (sqrti(1+8*n)-1)/2; }
    }
    /**
     * \endcond
     */
    
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * \ingroup diffopsElasto
     * \brief Returns the Voigt vector representation of a symmetric matrix (with or without doubled off-diagonal entries).
     * 
     * The result is 
     * \f[ \begin{bmatrix} a_{11} & a_{12} \\ * & a_{22} \end{bmatrix} \rightarrow \begin{bmatrix} a_{11} \\ a_{22} \\ s a_{12} \end{bmatrix}, \qquad 
     * \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ * & a_{22} & a_{23} \\ * & * & a_{33} \end{bmatrix} \rightarrow 
     * \begin{bmatrix} a_{11} \\ a_{22} \\ a_{33} \\ sa_{23} \\ s a_{13} \\ sa_{12} \end{bmatrix} \f]
     * 
     * Note that traditional Voigt notation treats the doubling of off-diagonal entries differently for stress and strain tensors
     * (\f$ s=2\f$ for strains but \f$ s=1 \f$ for stresses, such that \f$\sigma:\epsilon = \sigma^T \epsilon\f$). 
     * 
     * The inverse operation is provided as \ref unpack.
     * 
     * \param e the symmetric matrix
     * \param s the off-diagonal factor (defaults to 2)
     */
    template <class Scalar, int n>
    Dune::FieldVector<Scalar,n*(n+1)/2> pack(Dune::FieldMatrix<Scalar,n,n> const& e, Scalar s=2.0);
    
    /**
     * \ingroup diffopsElasto
     * \brief Returns the symmetric matrix representation of a Voigt vector with or without doubled off-diagonal entries.
     * 
     * This is the inverse operation of \ref pack.
     * 
     * \param e the symmetric matrix
     * \param s the off-diagonal factor (defaults to 2)
     */
    template <class Scalar, int n>
    Dune::FieldMatrix<Scalar,ElastoDetails::dim(n),ElastoDetails::dim(n)> unpack(Dune::FieldVector<Scalar,n> const& c, Scalar s=2.0);
    
    /**
     * \ingroup diffopsElasto
     * \brief The deviatoric part of a tensor.
     * 
     * This computes the deviatoric part \f$ \mathrm{dev}(s) = s - \frac{I_1}{d} I \f$ of the tensor \f$ s \f$, where
     * \f$ I_1 = \mathrm{tr}(s) \f$ is the first invariant.
     */
    template <int d, class Scalar>
    Dune::FieldMatrix<Scalar,d,d> deviatoricPart(Dune::FieldMatrix<Scalar,d,d> const& s)
    {
      return s - (trace(s)/d) * unitMatrix<Scalar,d>();
    }
    
    // --------------------------------------------------------------------------------------------
    
    /**
     * \ingroup stationaryElasticity
     * \brief A class for computing determinants and their derivatives.
     * 
     * This class supports the evaluation of \f$ |I+A|-1 \f$ for \f$ A \in \mathbb{K}^{d\times d} \f$
     * in a numerically stable way also for small matrices close to zero.
     */
    template <int dim, class Scalar>
    class DetIpm1
    {
    public:
      using Matrix = Dune::FieldMatrix<Scalar,dim,dim>;
    
      DetIpm1(Matrix const& A_): A(A_) {}
      
      Scalar d0() const;
      Scalar d1(Matrix const& dA) const;
      Scalar d2(Matrix const& dA1, Matrix const& dA2) const;
    private:
      Matrix A;
    };
    
    /**
     * \ingroup stationaryElasticity
     * \brief Numerically stable evaluation of \f$ (x+1)^p -1  \f$
     * 
     * The direct evaluation is numerically instable around zero due to cancellation of 
     * leading digits. This class performs a Taylor expansion around zero and achieves 
     * a stable evaluation.
     */
    class Pstable
    {
    public:
      /**
       * \brief Constructor
       * \param p exponent
       * \param x argument, x > -1 has to hold
       */
      Pstable(double p_, double x_);
      double d0() const { return f; }
      double d1(double dx) const { return df*dx; }
      double d2(double dx1, double dx2) const { return ddf*dx1*dx2; }
    private:
      double f, df, ddf;
    };
    
    // --------------------------------------------------------------------------------------------
    
    /**
     * \ingroup stationaryElasticity
     * \brief A class for shifted invariants of a tensor.
     * 
     * For a tensor \f$ A\in \mathbb{K}^{d\times d} \f$, the shifted invariants \f$ i_k \f$ are just the usual invariants
     * of \f$ I+A \f$ shifted by an additive constant. For 2D, they are
     * \f[ i_1 = I_1(I+A)-2 = \mathrm{tr}(I+A)-2 = \mathrm{tr}(A), \quad i_2 = I_2(I+A)-1 = \mathrm{det}(I+A)-1 \f]
     * and for 3D 
     * \f[ i_1 = I_1(I+A)-3 = \mathrm{tr}(I+A)-3 = \mathrm{tr}(A), \quad 
     *     i_2 = I_2(I+A)-3 = 2\mathrm{tr}(A) + I_2(A) = 2\mathrm{tr}(A) + \frac{1}{2}((\mathrm{tr}(A))^2-\mathrm{tr}(A^2)), \quad 
     *     i_3 = I_3(I+A) = \mathrm{det}(I+A)-1. \f]
     * 
     * Implementing material laws in terms of shifted invariants of \f$ 2E \f$ instead of in the usual invariants of \f$ C \f$ 
     * has the advantage of numerical stability in the vicinity of the reference configuration.
     * 
     * \tparam dim spatial dimension \f$ d \f$
     * \tparam Scalar the type of the underlying field \f$ K \f$
     */
    template <int dim, class Scalar=double>
    class ShiftedInvariants
    {
    public:
      using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
      using Invariants = Dune::FieldVector<Scalar,dim>;
      
      ShiftedInvariants(Tensor const& A_): A(A_), det(A) {}
      Invariants d0() const;
      Invariants d1(Tensor const& dA) const;
      Invariants d2(Tensor const& dA1, Tensor const& dA2) const;
      
    public:
      Tensor A;
      DetIpm1<dim,Scalar> det;
    };
  
    // --------------------------------------------------------------------------------------------

    /**
     * \ingroup stationaryElasticity
     * \brief Linearized (right) Green-Lagrange strain tensor, the workhorse of linear elastomechanics.
     *
     * This implements the linearized Green-Lagrange strain tensor
     * \f[ \epsilon(u_x) = \frac{1}{2} (u_x + u_x^T). \f]
     * This is an approximation of the (full) Green-Lagrange strain tensor linearized around the undeformed reference
     * configuration, and hence is valid only for small strains. In particular it violates rotational invariance.
     * 
     * \tparam Scalar the scalar field type, usually double
     * \tparam dim the spatial dimension, usually 2 or 3
     * \tparam byValue if true (default), the displacement derivative \f$ u_x \f$ is stored by value, otherwise by reference
     */
    template <class Scalar, int dim, bool byValue=true>
    class LinearizedGreenLagrangeTensor
    {
    public:
      /**
       * \brief The type of displacement derivatives.
       */
      typedef Dune::FieldMatrix<Scalar,dim,dim> Tensor;
      
      /**
       * \brief Constructor.
       * \param du The spatial derivative of the displacement.
       */
      explicit LinearizedGreenLagrangeTensor(Tensor const& du_): du(du_) {}
      
      /**
       * \brief Default constructor.
       * This initializes the displacement derivative to 0, i.e. to the reference configuration.
       */
      LinearizedGreenLagrangeTensor(): du(0) {}
      
      
      void setDisplacementDerivative(Tensor const& du_) [[gnu::deprecated("create a new LinearizedGreenLagrangeTensor if needed")]]
      {
        du = du_;
      }
      
      /**
       * \brief The linearized Green-Lagrange strain tensor \f$ \epsilon = \frac{1}{2} (u_x + u_x^T) \f$.
       */
      Tensor d0() const { return 0.5*(du+transpose(du)); }
      
      /**
       * \brief The first derivative of the linearized Green-Lagrange strain tensor in direction dv.
       *
       * This is \f$ \epsilon_{u_x}(u_x) v_x = \epsilon(v_x) \f$, as \f$ \epsilon \f$ is a linear operator.
       */
      Tensor d1(Tensor const& dv) const { return 0.5*(dv+transpose(dv)); }
      
      /**
       * \brief The second derivative of the linearized Green-Lagrange strain tensor in direction dv,dw.
       *
       * This is \f$ \epsilon_{u_x,u_x}(u_x) [v_x,w_x] =0 \f$, as \f$ \epsilon \f$ is a linear operator.
       */
      Tensor d2(Tensor const&, Tensor const&) const { return Tensor(0.0); }
      
    private:
      std::conditional_t<byValue,Tensor,Tensor const&> du;
    };
    
    /**
     * \ingroup stationaryElasticity
     * \brief Full (right) Green-Lagrange strain tensor, the workhorse of hyperelasticity.
     *
     * The right Green Lagrange strain tensor for a displacement \f$ u \f$ is defined as
     * \f[ E = \frac{1}{2}(u_x + u_x^T + u_x^T u_x) = \frac{1}{2}(F^TF - I) = \frac{1}{2}(C-I). \f]
     * 
     * \tparam Scalar the scalar field type, usually double
     * \tparam dim the spatial dimension, usually 2 or 3
     * \tparam byValue if true (default), the displacement derivative \f$ u_x \f$ is stored by value, otherwise by reference
     */
    template <class Scalar, int dim, bool byValue=true>
    class GreenLagrangeTensor
    {
    public:
      typedef Dune::FieldMatrix<Scalar,dim,dim> Tensor;
      
      /**
       * \brief Constructor.
       * \param du The spatial derivative of the displacement.
       */
      GreenLagrangeTensor(Tensor const& du_): du(du_) {}
      
      /**
       * \brief Default constructor.
       * This initializes the displacement derivative to 0, i.e. to the reference configuration.
       */
      GreenLagrangeTensor(): du(0) {}
      
      void setDisplacementDerivative(Tensor const& du_) [[gnu::deprecated("create a new LinearizedGreenLagrangeTensor if needed")]]
      {
        du = du_;
      }
      
      /**
       * \brief The linearized Green-Lagrange strain tensor \f$ E = \frac{1}{2} (u_x + u_x^T + u_x^T u_x) \f$.
       */
      Tensor d0() const { return 0.5*(du+transpose(du)+transpose(du)*du); }
      
      /**
       * \brief The first derivative of the linearized Green-Lagrange strain tensor in direction dv.
       *
       * This is \f$ \E_{u_x}(u_x) v_x = \frac{1}{2}(v_x + v_x^T + v_x^T u_x + u_x^T v_x). \f$
       */
      Tensor d1(Tensor const& dv) const { return 0.5*(dv+transpose(dv)+transpose(dv)*du + transpose(du)*dv); }
      
      /**
       * \brief The second derivative of the linearized Green-Lagrange strain tensor in direction dv,dw.
       *
       * This is \f$ E_{u_x,u_x}(u_x) [v_x,w_x] = \frac{1}{2}(v_x^T w_x + w_x^T v_x). \f$
       */
      Tensor d2(Tensor const& dv, Tensor const& dw) const { return 0.5*(transpose(dv)*dw + transpose(dw)*dv); }
      
    private:
      std::conditional_t<byValue,Tensor,Tensor const&> du;
    };
    
    
    /**
     * \ingroup stationaryElasticity
     * \brief Isochoric part of the full (right) Green-Lagrange strain tensor, used in compressible hyperelasticity.
     *
     * The isochoric right Green Lagrange strain tensor for a strain tensor \f$ E \f$ is defined as
     * \f[ \bar E = \frac{1}{2}(\bar C - I), \quad \bar C = (z+1) C, \quad z = |C|^{-1/d}-1. \f]
     * From that we infer the representation
     * \f[ \bar E = (z+1) E + \frac{z}{2} I, \f]
     * which is preferable as it is numerically stable in the vicinity of the reference configuration. \f$ z \f$
     * is computed in a stable way as \f[ z = \exp\left((-\frac{1}{d}\log( |I+2E|-1 + 1)\right)-1 , \f]
     * where the C++ library functions std::expm1 and std::log1p as well as \ref detIpm1 can be used for stable
     * evaluation.
     * 
     * 
     * \tparam dim the spatial dimension, usually 2 or 3
     * \tparam Scalar the scalar field type, defaults to double
     */
    template <int dim, class Scalar=double>
    class IsochoricGreenLagrangeTensor
    {
    public:
      typedef Dune::FieldMatrix<Scalar,dim,dim> Tensor;
      
      /**
       * \brief Constructor.
       * \param E The Green-Lagrange strain tensor.
       */
      IsochoricGreenLagrangeTensor(Tensor const& E_): E(E_), det(2*E), z(-1.0/dim,det.d0())
      {
      }
      
      /**
       * \brief Default constructor.
       * This initializes the strain tensor to 0, i.e. to the reference configuration.
       */
      IsochoricGreenLagrangeTensor(): E(0), z(0.0,0.0) {}
      
      /**
       * \brief The isochoric linearized Green-Lagrange strain tensor \f$ \bar E \f$.
       */
      Tensor d0() const { 
        return E + z.d0()*(E + 0.5*unitMatrix<Scalar,dim>()); 
      }
      
      /**
       * \brief The first derivative of the linearized Green-Lagrange strain tensor in direction de.
       *
       */
      Tensor d1(Tensor const& dE) const 
      { 
        return dE + z.d1(det.d1(2*dE))*(E + 0.5*unitMatrix<Scalar,dim>()) + z.d0()*dE;
      }
      
      /**
       * \brief The second derivative of the linearized Green-Lagrange strain tensor in direction dv,dw.
       *
       */
      Tensor d2(Tensor const& dE1, Tensor const& dE2) const 
      { 
        return z.d2(det.d1(2*dE1,2*dE2))*(E + 0.5*unitMatrix<Scalar,dim>()) + z.d1(det.d1(2*dE1))*dE2 + z.d1(det.d1(2*dE2))*dE1; 
      }
      
      /**
       * \brief The shifted determinant of the original Green Lagrange strain tensor.
       */
      DetIpm1<dim,Scalar> const& determinant() const
      {
        return det;
      }
      
    private:
      Tensor E;
      DetIpm1<dim,Scalar> det;
      Pstable z;
    };
    
    // --------------------------------------------------------------------------------------------
    
    
    /**
     * \ingroup stationaryElasticity
     * \brief A function view that provides on the fly computed strain tensors of displacemnts.
     * 
     * \tparam Displacement A finite element function type with vectorial value type
     * \tparam StrainTensor A strain tensor class, usually GreenLagrangeTensor or LinearizedGreenLagrangeTensor.
     */
    template <class Displacement, class StrainTensor>
    class StrainView
    {
    public:
      using Space = typename Displacement::Space;
      using ValueType = typename StrainTensor::Tensor;
      
      
      /**
       * \brief Constructor.
       * 
       * \param u the displacement (has to exist during the lifetime of the StrainView).
       */
      StrainView(Displacement const& u_)
      : u(u_) {}
      
      Space const& space() const
      {
        return u.space();
      }
      
      ValueType value(typename Space::Evaluator const& evaluator) const
      {
        return StrainTensor(u.derivative(evaluator)).d0();
      }
      
    private:
      Displacement const& u;
    };
    
    /**
     * \ingroup stationaryElasticity
     * \brief A convenience function for creating function views for linearized Green-Lagrange strain tensors.
     */
    template <class Displacement>
    auto makeLinearizedGreenLagrangeStrainView(Displacement const& u)
    {
      return StrainView<Displacement,LinearizedGreenLagrangeTensor<typename Displacement::Scalar,Displacement::Space::GridView::dimension>>(u);
    }
    
    /**
     * \ingroup stationaryElasticity
     * \brief A convenience function for creating function views for full Green-Lagrange strain tensors.
     */
    template <class Displacement>
    auto makeGreenLagrangeStrainView(Displacement const& u)
    {
      return StrainView<Displacement,GreenLagrangeTensor<typename Displacement::Scalar,Displacement::Space::GridView::dimension>>(u);
    }
    
    // ---------------------------------------------------------------------------------------------------------
    
    
    
    // ---------------------------------------------------------------------------------------------------------
    
    
    /// \cond internals
    namespace Elasto_Details {
      
      template <class Real, int n> 
      Real maxOrientationPreservingStepsize(Dune::FieldMatrix<Real,n,n> const& A, Dune::FieldMatrix<Real,n,n> const& dA, Real eps);
      
      template <class Function>
      struct Limiter {
        static int const dim = Function::components;
        
        Limiter(Function const& y_, Function const& dy_, double epsilon_)
        : y(y_), dy(dy_), epsilon(epsilon_), step(1), I(unitMatrix<typename Function::Scalar,dim>()) { mindet = std::numeric_limits<double>::max(); }
        
        template <class Cell, class QuadraturePoint, class Evaluator>
        void operator()(Cell const&, QuadraturePoint const&, Evaluator const& eval) {
          auto yx = y.derivative(eval);
          auto dyx = dy.derivative(eval);
          auto tmp = determinant(I+yx);
          
          if (tmp <= 0)
            std::cerr << "#### AIEEE: given determinant nonpositive: " << tmp << " !\n";
          
          mindet = std::min(mindet,tmp);
          while (step>0 && determinant(I+yx+step*dyx)<epsilon*tmp)
            step = 0.7*step - 1e-4;
          
          for (int i=1; i<10; ++i)
            if (determinant(I+yx+i*step/10*dyx)<epsilon*tmp)
              std::cerr << "##### Aieee: at step " << step << " * " << i << "/10 we have " << determinant(I+yx+i*step/10*dyx) << " < " << epsilon*tmp << "\n";
        }
        
        template <class Evaluator>
        int order(Evaluator const& eval) const { return y.order(eval); }
        
        Function const& y, dy;
        double epsilon, step, mindet;
        Dune::FieldMatrix<typename Function::Scalar,dim,dim> I;
      };
      
    }
    /// \endcond
    
    /** \ingroup stationaryElasticity
     * \brief Computes the maximal orientation preserving stepsize.
     * 
     * Given two FE functions \f$ y, \delta y \f$ with \f$ \det (I+y_x) > 0 \f$, computes the maximum step size \f$ \alpha \le 1 \f$ such
     * that \f$ \det(I+y_x+\alpha \delta y_x) \ge \epsilon \det (I+y_x) \f$.
     * 
     * \param epsilon (<1)
     */
    template <class Function>
    double orientationPreservingStepsize(Function const& y, Function const& dy, double epsilon)
    {
      
      Elasto_Details::Limiter<Function> g(y,dy,epsilon);
      forEachQP(g,y.space());
                     std::cout << "minimal determinant: " << g.mindet << "\n";
      return std::max(g.step,0.0);
    }
    
  // ---------------------------------------------------------------------------------------------------------
  }
  
  // TODO: remove this backwards compatibility hack as soon as all users switched over
  using namespace Elastomechanics;
  
}


  

#endif
