/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FEM_DIFFOPS_ELASTOVARIATIONALFUNCTIONALS_HH
#define FEM_DIFFOPS_ELASTOVARIATIONALFUNCTIONALS_HH

#include "fem/diffops/materialLaws.hh"

namespace Kaskade 
{
  namespace Elastomechanics
  {
    
    /**
     * \ingroup stationaryElasticity
     * \brief General base class for variational functionals defined in terms of hyperelastic stored energies.
     * 
     * This defines the interface for stored energies working on displacement derivatives \f$ W = W(u_x) \f$. 
     * The interface is intended to be used for defining variational functionals.
     * 
     * The class has as a member an object of the hyperelastic energy defining the stored energy \f$ W(E) \f$ 
     * depending on the Green-Lagrange strain tensor.
     * 
     * \tparam HyperelasticEnergy the hyperelastic energy class 
     * \tparam StrainTensor the strain tensor class, usually one of GreenLagrangeTensor or LinearizedGreenLagrangeTensor
     * 
     * Example for linear elastomechanics of steel:
     * \code
     * HyperelasticVariationalFunctional<MaterialLaws::StVenantKirchhoff<dim>,LinearizedGreenLagrangeTensor<double,dim>> energy(ElasticModulus::material("steel"));
     * energy.setLinearizationPoint(du);
     * \endcode
     */
    template <class HyperelasticEnergy, class StrainTensor>
    class HyperelasticVariationalFunctional
    {
    public:
      using Scalar = typename HyperelasticEnergy::Scalar;
      static int const dim = HyperelasticEnergy::dim;
      using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
      
      /**
       * \brief Constructor.
       * 
       * The constructor arguments are forwarded to the hyperelastic stored energy constructor.
       * The strain tensor is default initialized.
       */
      template <typename... Args>
      HyperelasticVariationalFunctional(Args&&... args): energy(std::forward<Args>(args)...)
      {
      }
      
      /**
       * \brief Constructor.
       */
      HyperelasticVariationalFunctional(HyperelasticEnergy const& energy_, StrainTensor const& strain_)
      : du(0), energy(energy_), strain(strain_)
      { 
      }
      
      /**
       * \brief Defines the displacement derivative around which to linearize.
       * 
       * This method shall be called when the evaluation of a stored energy and its derivatives
       * is intended for a different displacement derivative, in particular in the "evaluateAt" 
       * method of a domain cache.
       * 
       * \param du0 the current displacement derivative \f$ u_x \f$
       */
      void setLinearizationPoint(Dune::FieldMatrix<Scalar,dim,dim> const& du0)
      {
        du = du0;
        strain = StrainTensor(du);
        energy.setLinearizationPoint(strain.d0());
      }
      
      /**
       * \brief Computes the hyperelastic stored energy density.
       */
      Scalar d0() const
      {
        return energy.d0();
      }
      
      /**
       * \brief Computes the directional derivative of the hyperelastic stored energy density.
       */
      Dune::FieldVector<Scalar,dim> d1(VariationalArg<Scalar,dim> const& arg) const
      {
        // TODO: I have doubts that this is the most efficient implementation...
        Dune::FieldVector<Scalar, dim> ret;
        
        // step trough all spatial dimensions and test with vectorial test function 
        // with values only in that direction.
        for (int i=0; i<dim; ++i) {
          Dune::FieldMatrix<Scalar,dim,dim> dv(0);
          dv[i] = arg.derivative[0];
          
          Dune::FieldMatrix<Scalar,dim,dim> dE = strain.d1(dv);
          ret[i] = energy.d1(dE);
        }
                
        return ret;
      }
      
      /**
       * \brief Computes the second directional derivative of the hyperelastic stored energy density.
       * 
       * For scalar functions \f$ \phi, \psi \f$, the values and derivatives at a certain point are provided in arg1 and arg 2,
       * this computes the \f$ d\times d \f$ matrix \f$ A \f$ with \f$ A_{ij} = \epsilon((u_i)_x)^T \mathcal{C} \epsilon((v_j)_x) \f$,
       * where \f$ u_i = \phi e_i, v_j = \psi e_j \f$ and \f$ e_k \f$ is the \f$ k \f$-th unit vector.
       * 
       * This is exactly what is to be returned from the d2 method of a a variational functional's domain cache implementing the Lame-Navier
       * equations.
       */
      Tensor d2(VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
      {
        // TODO: I have doubts that this is the most efficient implementation...
        Dune::FieldMatrix<Scalar, dim, dim> ret;
        
        // In a double loop, step trough all spatial dimensions and test with vectorial test functions 
        // with values only in that direction.
        // TODO: due to symmetry, can we skip half the computation?
        for (int i=0; i<dim; ++i) {
          Dune::FieldMatrix<Scalar,dim,dim> dv1(0);
          dv1[i] = arg1.derivative[0];
          Dune::FieldMatrix<Scalar,dim,dim> dE1 = strain.d1(dv1);
          
          for (int j=0; j<dim; ++j) {
            Dune::FieldMatrix<Scalar,dim,dim> dv2(0);
            dv2[j] = arg2.derivative[0];
            Dune::FieldMatrix<Scalar,dim,dim> dE2 = strain.d1(dv2);
            Dune::FieldMatrix<Scalar,dim,dim> ddE = strain.d2(dv1,dv2);
            
            // according to product rule...
            ret[i][j] = energy.d2(dE1,dE2) + energy.d1(ddE);
          }
        }
        return ret;
      }
      
      /**
       * \brief Provides access to the stored energy density.
       * 
       * The stored energy function has the linearization (and evaluation) point given by
       * the latest call to the setLinearizationPoint() method.
       */
      HyperelasticEnergy& getHyperelasticEnergy() 
      {
        return energy;
      }
      
      /**
       * \brief Returns the 2nd Piola-Kirchhoff stress tensor \f$ S = \frac{\partial W}{\partial E} \f$ corresponding to the current displacement derivative.
       */
      Tensor stress() const { return energy.stress(); }
      
      /**
       * \brief Returns the Cauchy stress tensor \f$ \sigma = \frac{1}{J} F S F^T \f$ corresponding to the current displacement derivative.
       */
      Tensor cauchyStress() const 
      { 
        Tensor F = du; 
        for (int i=0; i<dim; ++i) F[i][i] += 1; // F = I + ux
        
        return F*stress()*transpose(F) / F.determinant(); 
      }
      
    protected:
      Tensor             du;
      HyperelasticEnergy energy;
      StrainTensor       strain;
    };
    
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * \ingroup stationaryElasticity
     * \brief Convenience class for handling linear elastomechanics.
     * 
     * This class defines the variational formulation of the Lame-Navier operator \f$ u \mapsto -2\mu\Delta u - \lambda \nabla \mbox{div} u\f$
     * for vector-valued displacements \f$ u \f$, i.e. the variational functional 
     * \f$ u \mapsto \frac{\lambda}{2} \mbox{tr}(E)^2 + \mu(E:E) \f$ with \f$ E = \frac{1}{2}(u_x + u_x^T) \f$.
     */
    template <int dim_, class Scalar_=double>
    class LameNavier: public HyperelasticVariationalFunctional<MaterialLaws::StVenantKirchhoff<dim_,Scalar_>,LinearizedGreenLagrangeTensor<Scalar_,dim_>>
    {
      using Energy = MaterialLaws::StVenantKirchhoff<dim_,Scalar_>;
      using Base = HyperelasticVariationalFunctional<Energy,LinearizedGreenLagrangeTensor<Scalar_,dim_>>;
      
      
    public:
      using Scalar = Scalar_;
      
    private:
      Scalar lambda, mu;
      
    public:
      static int const dim = dim_;
      
      /**
       * \brief Constructor.
       * Both Lame constants are initialized to 1, the linearization point to \f$ (u_0)_x = 0 \f$.
       */
      LameNavier(): LameNavier(ElasticModulus(1,1)) {}
      
      LameNavier(ElasticModulus const& moduli): Base(moduli), lambda(moduli.lame()), mu(moduli.shear())
      {
      }
      
      /**
       * \brief Define material properties.
       * 
       * This invalidates previously set linearization points.
       */
      void setElasticModuli(ElasticModulus const& modulus)
      {
        lambda = modulus.lame();
        mu = modulus.shear();
        this->energy = Energy(modulus);
      }
      
      /**
       * \brief Retrieve current material properties.
       */
      ElasticModulus getElasticModuli()
      {
        return ElasticModulus(lambda,mu);
      }
      
      /**
       * \brief Returns the matrix \f$ C \f$ mapping the strain tensor to the stress tensor: \f$ \sigma = C \epsilon \f$ 
       * 
       * In this contex, the symmetric matrices \f$ \sigma \f$ and \f$ \epsilon \f$ are interpreted as vectors in \f$ \mathbb{R}^{d(d+1)/2} \f$:
       * \f[ \sigma = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{12} \\ \sigma_{13} \\ \sigma_{23} \end{bmatrix} \quad\mbox{or}\quad
       *      \sigma = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\  \sigma_{12}  \end{bmatrix} \f]
       */
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> strainToStressMatrix() const
      {
        // The implementation follows Braess, Finite Elemente, Chapter VI.3
        Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> C(0.0);
        
        if (dim==3)
        {
          C[0][0] = lambda+2*mu; C[0][1] = lambda; C[0][2] = lambda;
          C[1][0] = lambda; C[1][1] = lambda+2*mu; C[1][2] = lambda;
          C[2][0] = lambda; C[2][1] = lambda; C[2][2] = lambda+2*mu;
          for (int i=3; i<6; ++i)
            C[i][i] = 2*mu;
        }
        else if (dim==2)       
        {
          C[0][0] = lambda+2*mu; C[0][1] = lambda; 
          C[1][0] = lambda; C[1][1] = lambda+2*mu; 
          C[2][2] = 2*mu;
        }
        return C;
      }
      
      /**
       * \brief Computes the symmetric matrix times vector product, when the matrix is compactly stored in a dim*(dim+1)/2-vector.
       * 
       * This is useful, e.g., for computing normal stresses \f$ \sigma n \f$.
       */
      static Dune::FieldVector<Scalar,dim> mv(Dune::FieldVector<Scalar,dim*(dim+1)/2> const& A, Dune::FieldVector<Scalar,dim> const& x)
      {
        Dune::FieldVector<Scalar,dim> result; 
        
        if (dim==3)
        {
          result[0] = A[0]*x[0] + A[3]*x[1] + A[4]*x[2];
          result[1] = A[3]*x[0] + A[1]*x[1] + A[5]*x[2];
          result[2] = A[4]*x[0] + A[5]*x[1] + A[2]*x[2];
        }
        else if (dim==2)
        {
          result[0] = A[0]*x[0] + A[2]*x[1];
          result[1] = A[2]*x[0] + A[1]*x[1];
        }
        else if (dim==1)
        {
          abort();
        }
        
        return result;
      }
      
      /**
       * \brief returns the matrix \f$ E \f$ mapping the displacement derivative to the strain tensor: \f$ \epsilon = E u_x \f$
       * 
       * In this context, the symmetric matrix \f$ \epsilon \f$ is interpreted as vector in \f$ \mathbb{R}^{d(d+1)/2} \f$:
       * \f[ \epsilon = \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ \epsilon_{12} \\ \epsilon_{13} \\ \epsilon_{23} \end{bmatrix} \quad\mbox{or}\quad
       *      \epsilon = \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\  \epsilon_{12}  \end{bmatrix} \f]
       * The displacement derivative \f$ u_x \f$ is interpreted as vector in \f$ \mathbb{R}^{d^2} \f$ (column-major format):
       * \f[ \frac{\partial u_i}{\partial x_j} = (u_x)_{i+dj} \f]
       */
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*dim> displacementDerivativeToStrainMatrix() const
      {
        Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*dim> E(0.0);
        
        if (dim==3)
        {
          E[0][0] = 2; // diagonal elements of epsilon
          E[1][4] = 2;
          E[2][8] = 2;
          E[3][1] = 1; E[3][3] = 1; // off-diagonal elements of epsilon
          E[4][2] = 1; E[4][6] = 1;
          E[5][5] = 1; E[5][7] = 1;
        }
        else if (dim==2)
        {
          E[0][0] = 2;
          E[1][3] = 2;
          E[2][1] = 1; E[2][2] = 1;
        } else if (dim==1)
          E[0][0] = 2;
        
        E *= 0.5;
        return E;
      }
    };
    
    // ---------------------------------------------------------------------------------------------------------
    
    /**
     * \ingroup stationaryElasticity
     * \brief Convenience class for handling orthotropic materials in linear elastomechanics.
     * 
     * This class defines the variational formulation of the linear elasticity for orthotropic materials
     * for vector-valued displacements \f$ u \f$, i.e. the variational functional 
     * \f$ u \mapsto \frac{1}{2} \epsilon : \mathcal{C}\epsilon \f$ with \f$ \epsilon = \frac{1}{2}(u_x + u_x^T) \f$.
     * In an appropriately rotated coordinate system, the stiffness matrix, representing the tensor \f$\mathcal{C}\f$ by interpreting 
     * \f$ \epsilon, \sigma \f$ as vectors, has the form 
     * \f[ C = \begin{bmatrix} * & * & \\ * & * & \\ & & * \end{bmatrix}, \quad C = \begin{bmatrix} * & * & * & & & \\ * & * & * & & & \\ * & * & * & & & \\
     *     & & & * & & \\ & & & & * & \\ & & & & & * \end{bmatrix} \f]
     * depending on the spatial dimension (see <a href="https://en.wikipedia.org/wiki/Orthotropic_material">Wikipedia</a>).
     */
    template <class Scalar, int dim>
    class OrthotropicLameNavier {
    public:
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat  The material parameters.
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat is a matrix containing the material parameters in the following way:
       * \f[ \begin{bmatrix} E_1 & G_{12}  \\ \nu_{12} & E_2 \end{bmatrix}, \quad 
       *     \begin{bmatrix} E_1 & G_{12} & G_{13} \\ \nu_{12} & E_2 & G_{23} \\ \nu_{13} & \nu_{23} & E_3 \end{bmatrix}, \f]
       * where \f$ E_i \f$ is the elastic modulus in material coordinate direction \f$ i \f$, \f$ \nu_{ij} \f$ is the Poissons ratio
       * for contraction in \f$ j \f$ direction due to strain in \f$ i \f$ direction, and \f$ G_{ij} \f$ is the shear modulus in
       * \f$ ij \f$ plane.
       */
      OrthotropicLameNavier(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, Dune::FieldMatrix<Scalar,dim,dim> const& mat_);
      
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat1   material parameters
       * \param mat2   material parameters
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat1 is a \f$ d \times d \f$ matrix containing the top left part of the stiffness matrix \f$ C \f$, and \a mat2
       * is a \f$ d(d-1)/2 \f$ vector containing the diagonal of the bottom right part.
       * \f[ C = \begin{bmatrix} \mathrm{mat1} &  \\  & \mathrm{diag}(\mathrm{mat2}) \end{bmatrix}, \f]
       */
      OrthotropicLameNavier(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, 
                            Dune::FieldMatrix<Scalar,dim,dim> const& mat1, Dune::FieldVector<Scalar,dim*(dim-1)/2> const& mat2);
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with given material parameters.
       */
      OrthotropicLameNavier(ElasticModulus const& p);
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with Lame parameters set to 1.
       */
      OrthotropicLameNavier()
      : OrthotropicLameNavier(ElasticModulus()) {}
      
      /**
       * \brief Returns the stiffness tensor.
       *
       * The stiffness tensor \f$ C \f$ is represented in matrix form for Voigt notation (i.e. the mapping from Voigt strain vector to Voigt stress vector). 
       * In 3D this reads
       * \f[ \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix} = C
       *     \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} \end{bmatrix} \f]
       */
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& stiffnessMatrix() const { return stiffness; }
      
      /**
       * \brief Defines the displacement derivative around which to linearize.
       */
      void setLinearizationPoint(Dune::FieldMatrix<Scalar,dim,dim> const& du0_)
      {
        du0 = du0_;
      }
      
      /**
       * \brief Computes the elastic energy \f$ \frac{1}{2} \epsilon(u_x) : \mathcal{C}\epsilon(u_x) \f$.
       */
      Scalar d0() const;
      
      /**
       * \brief Computes the elastic energy \f$ \frac{1}{2} \epsilon(u_x) : \mathcal{C}\epsilon(u_x) \f$.
       * 
       * \a arg is assumed to be a scalar variational argument.
       */
      Dune::FieldVector<Scalar,dim> d1(VariationalArg<Scalar,dim> const& arg) const;
      
      /**
       * \brief Computes the local Hessian of the variational functional.
       * For scalar functions \f$ \phi, \psi \f$, the values and derivatives at a certain point are provided in arg1 and arg 2,
       * this computes the \f$ d\times d \f$ matrix \f$ A \f$ with \f$ A_{ij} = \epsilon((u_i)_x)^T \mathcal{C} \epsilon((v_j)_x) \f$,
       * where \f$ u_i = \phi e_i, v_j = \psi e_j \f$ and \f$ e_k \f$ is the \f$ k \f$-th unit vector.
       * 
       * This is exactly what is to be returned from the d2 method of a a variational functional's domain cache implementing the Lame-Navier
       * equations.
       */
      Dune::FieldMatrix<Scalar,dim,dim> d2(VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const;
      
      
    private:
      Dune::FieldMatrix<Scalar,dim,dim> orth; // the orthonormal matrix mapping material to global coordinates
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> stiffness;
      Dune::FieldMatrix<Scalar,dim,dim> du0; // the linearization point of the displacement derivative
    };
    
    // ---------------------------------------------------------------------------------------------------------
    
    
  }
}

#endif