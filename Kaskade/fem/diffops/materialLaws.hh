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

#ifndef FEM_DIFFOPS_MATERIALLAWS_HH
#define FEM_DIFFOPS_MATERIALLAWS_HH

#include <cassert>

#include "dune/common/fmatrix.hh"

#include "fem/diffops/elasto.hh"

namespace Kaskade
{
  namespace Elastomechanics
  {
    /**
     * \ingroup diffopsElasto
     * \brief A namespace containing various material laws.
     * 
     * Hyperelastic material laws define a stored energy density \f$ W(E) \f$ in terms of the
     * Green-Lagrange strain tensor \f$ E \f$.
     */
    namespace MaterialLaws 
    {
      /**
       * \ingroup stationaryElasticity
       * \brief Base class for hyperelastic material laws, providing default implementations of the stress and the tangent stiffness tensor \f$ C \f$.
       * 
       * The stiffness tensor computation is injected by CRTP:
       * \code
       * class MyMaterial: public MaterialLawBase<MyMaterial> {...};
       * \endcode
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law implementation (derived class)
       */
      template <int dim_, class Scalar_, class MaterialLaw>
      class MaterialLawBase
      {
      public:
        using Scalar = Scalar_;
        static int const dim = dim_;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>; 
        using VoigtTensor = Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2>;
        
        /**
         * \brief Returns the second Piola-Kirchhoff stress tensor \f$ S = W'(E) \f$ corresponding to the current strain \f$ E \f$.
         */
        Tensor stress() const
        {
          Tensor sigma;
          
          for (int i=0; i<dim; ++i)
            for (int j=0; j<=i; ++j)
            {
              Tensor E(0); 
              E[i][j] = 1;
              sigma[i][j] = static_cast<MaterialLaw const*>(this)->d1(E);
            }
            
          for (int i=0; i<dim-1; ++i)
            for (int j=i+1; j<dim; ++j)
              sigma[i][j] = sigma[j][i];
            
          return sigma;
        }
        
        /**
         * \brief Returns the tangent stiffness tensor \f$ C(E) = W''(E)\f$ mapping strain tensor variations to stress tensor (2nd Piola-Kirchhoff) variations.
         * 
         * In this context, the symmetric matrices \f$ \sigma \f$ and \f$ \epsilon \f$ are interpreted as \f$d(d+1)/2\f$-vectors in Voigt notation as
         * returned by pack. Note that strain tensor off-diagonal entries appear with factor two in the vector.
         * 
         * Now, a matrix \f$ C\in\mathbb{R}^{(d+1)d/2\times(d+1)d/2} \f$ is returned such that for infinitesimal changes \f$\delta\epsilon\f$ the 
         * corresponding change in \f$ \sigma(\epsilon+\delta\epsilon) = \sigma(\epsilon)+\delta\sigma \f$ is given as \f$ \delta\sigma = C \delta\epsilon \f$. 
         * Note that due to the symmetry of the stored energy Hessian \f$ W''(E) \f$, the Voigt representation of the stiffness tensor is a
         * symmetric matrix.
         * 
         * Note that is generic implementation usually exhibits a suboptimal performance. Specialized implementations can provide a significantly
         * higher performance. For linear material laws, caching of the (constant) result might be an option.
         */
        VoigtTensor strainToStressMatrix() const
        {
          int const n = dim*(dim+1)/2;
          Dune::FieldMatrix<Scalar,n,n> C;
          
          // For the Voigt notation it holds that x^T C y = W''[unpack(x),unpack(y)]. 
          // Due to symmetry of C, we only compute the lower half directly...
          for (int i=0; i<n; ++i)
          {
            Dune::FieldVector<Scalar,n> x;
            x[i] = 1;
            auto X = unpack(x);
            
            for (int j=0; j<=i; ++j)
            {
              Dune::FieldVector<Scalar,n> y;
              y[j] = 1;
              auto Y = unpack(y);
              
              C[i][j] = static_cast<MaterialLaw const*>(this)->d2(X,Y);
            }
          }
          
          // ... and copy it to the upper half.
          for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
              C[i][j] = C[j][i];
            
            return C;
        }
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup stationaryElasticity
       * \brief Adaptor for hyperelastic material laws, providing an easy way to formulate incompressible material laws in terms of 
       * the invariants of the Cauchy-Green strain tensor \f$ C \f$.
       * 
       * This class essentially implements the composition \f$ W(E) = \hat W(I_1,I_2,I_3) \f$ with \f$ I_i = I_i(I+2E) \f$, computing
       * derivatives by the chain rule. Note that most hyperelastic material laws given in terms of the invariants are actually designed
       * for incompressible material, i.e. \f$ I_d = 1 \f$, and should be used directly only in incompressible computations. For compressible
       * hyperelasticity consider using \ref CompressibleInvariantsMaterialLaw.
       * 
       * Model of HyperelasticMaterialLaw.
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law class, a model of InvariantsMaterialConcept
       */
      template <class Material>
      class InvariantsMaterialLaw: public MaterialLawBase<Material::dim,typename Material::Scalar,InvariantsMaterialLaw<Material>>
      {
      public:
        using Scalar = typename Material::Scalar;
        static int const dim = Material::dim;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Constructor.
         * \param E the current Green-Lagrange strain tensor (possibly linearized).
         */
        template <class ... Args>
        InvariantsMaterialLaw(Tensor const& E, Args... args): invariants(2*E), material(args...) {}
        
        /**
         * \name Hyperelastic material law interface
         * @{
         */
        
        /**
         * \brief Defines a new evaluation/linearization point.
         * \param e the strain tensor
         */
        void setLinearizationPoint(Tensor const& e)
        {
          invariants = ShiftedInvariants<dim,Scalar>(2*e);
          material.setLinearizationPoint(invariants.d0());
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(E) \f$.
         */
        Scalar d0() const
        {
          return material.d0();
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(E)E_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          return 2*material.d1(invariants.d1(e1));
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(E)E_1 E_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          return 4*(material.d2(invariants.d1(e1),invariants.d1(e2)) + material.d1(invariants.d2(e1,e2)));
        }
        
        /**
         * @}
         */
      private:
        ShiftedInvariants<dim,Scalar> invariants;
        Material                      material;
      };
      
      /**
       * \ingroup stationaryElasticity
       * \brief Adaptor for hyperelastic material laws, providing an easy way to formulate compressible material laws in terms of 
       * the invariants of the isochoric part \f$ \bar C = C/(det C)^(1/d) \f$ of the Cauchy-Green strain tensor and a penalization
       * of the deviatoric part \f$ I_d \f$.
       * 
       * This class essentially implements the law \f$ W(E) = \bar W(\bar E) + D(I_d)\f$, where \f$ \bar E = \frac{E}{z} + \frac{1-z}{2z} I \f$
       * is the isochoric part of the strain tensor and \f$ z = I_d^{1/d} \f$, computing derivatives by the chain rule.
       * 
       * Model of HyperelasticMaterialLaw.
       *
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law class, a model of InvariantsMaterialConcept
       * \tparam DeviatoricPenalty
       */
      template <class Material, class DeviatoricPenalty>
      class CompressibleInvariantsMaterialLaw
      : public MaterialLawBase<Material::dim,typename Material::Scalar,CompressibleInvariantsMaterialLaw<Material,DeviatoricPenalty>>
      {
      public:
        using Scalar = typename Material::Scalar;
        static int const dim = Material::dim;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        template <class ... Args>
        CompressibleInvariantsMaterialLaw(Tensor const& E, Args... args): bE(E), material(bE.d0(),args...), deviator(bE.determinant().d0()) {}
        
        /**
         * \brief Defines a new evaluation/linearization point.
         * \param e the strain tensor
         */
        void setLinearizationPoint(Tensor const& e)
        {
          bE = IsochoricGreenLagrangeTensor<dim,Scalar>(e);
          deviator = DeviatoricPenalty(bE.determinant().d0());
          material.setLinearizationPoint(bE.d0());
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(\epsilon) \f$.
         */
        Scalar d0() const
        {
          return material.d0() + deviator.d0();
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(\epsilon)\epsilon_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          return material.d1(bE.d1(e1)) + deviator.d1(bE.determinant().d1(e1));
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          return material.d2(bE.d1(e1),bE.d1(e2)) + material.d1(bE.d2(e1,e2))
                 + deviator.d2(bE.determinant().d1(e1),bE.determinant().d1(e2)) + deviator.d1(bE.determinant().d2(e1,e2));
        }
        
      private:
        IsochoricGreenLagrangeTensor<dim,Scalar> bE;
        InvariantsMaterialLaw<Material>          material;
        DeviatoricPenalty                        deviator;
      };
      
      // ---------------------------------------------------------------------------------------------------------
      // ---------------------------------------------------------------------------------------------------------
      
      
      /**
       * \ingroup stationaryElasticity
       * \brief The St. Venant-Kirchhoff material, foundation of linear elastomechanics.
       * 
       * The stored energy density is defined as
       * \f[ W(E) = \frac{\lambda}{2} (\mathop\mathrm{tr} E)^2 + \mu E:E. \f]
       * 
       * The St. Venant-Kirchhoff material, characterized by any two of the five elastic moduli, is most
       * useful for linear elasticity. For finite strain elasticity, its use is not recommended, as it
       * is not polyconvex and does not preserve orientation of the deformation.
       * 
       * This is a model of the \ref HyperelasticMaterialLaw concept.
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the type of real numbers to use
       */
      template <int dim, class Scalar=double>
      class StVenantKirchhoff: public MaterialLawBase<dim,Scalar,StVenantKirchhoff<dim,Scalar>>
      {
      public:
        
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Pretty useless default constructor.
         */
        StVenantKirchhoff(): lambda(1), mu(1), e(0), tre(0)
        {}
        
        /**
         * \brief Constructor.
         * \param moduli the material parameters
         */
        StVenantKirchhoff(ElasticModulus const& moduli): lambda(moduli.lame()), mu(moduli.shear()), e(0), tre(0)
        {
        }
        
        /**
         * \brief Constructor.
         * \param moduli the material parameters
         * \param e the initial Green-Lagrange strain tensor
         * 
         * This is equivalent to:
         * \code
         * StVenantKirchhoff material(moduli);
         * setLinearizationPoint(e);
         * \endcode
         */
        StVenantKirchhoff(ElasticModulus const& moduli, Tensor const& e_): lambda(moduli.lame()), mu(moduli.shear()), e(e_), tre(trace(e))
        {
        }
        
        /**
         * \name Hyperelastic material law interface
         * @{
         */
        
        /**
         * \brief Defines the linearization/evaluation point for subsequent calls to d0, d1, d2.
         */
        void setLinearizationPoint(Tensor const& e_)
        {
          e = e_;
          tre = trace(e);
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(\epsilon) \f$.
         */
        Scalar d0() const
        {
          return lambda*tre*tre/2 + mu*contraction(e,e);
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(\epsilon)\epsilon_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          return lambda*tre*trace(e1) + 2*mu*contraction(e,e1);
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          return lambda*trace(e1)*trace(e2) + 2*mu*contraction(e1,e2);
        }
        
        /**
         * @}
         */
        
        private:
          Scalar lambda, mu; // Lame parameters
          Tensor e;          // the strain tensor
          Scalar tre;        // trace of strain tensor
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup stationaryElasticity
       * \brief Mooney-Rivlin material law formulated in terms of the (shifted) invariants \f$ i_1, i_2, i_3 \f$ of the doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * The Mooney-Rivlin hyperelastic material energy for incompressible materials is defined in terms of the Cauchy-Green strain tensor \f$ C = I+2E \f$ as
       * \f[ W = C_1(I_1-3) + C_2(I_2-3) = C_1 i_1 + C_2 i_2. \f]
       * The material parameters \f$ 2(C_1+C_2) \f$ correspond to the shear modulus \f$ \mu \f$ of linear elasticity. 
       * 
       * The implementation is numerically stable in the vicinity of the reference configuration.
       * 
       * Use the \ref IncompressibleInvariantsMaterialLaw adaptor for creating a hyperelastic material law for incompressible situations, i.e. \f$ I_3 = 1 \f$, 
       * or the \ref CompressibleInvariantsMaterialLaw adaptor in case of compressible materials.
       * 
       * Model of InvariantsMaterialConcept.
       *
       * \see Kaskade::InvariantsMaterialLaw
       */
      template <int dimension>
      class MooneyRivlin
      {
      public:
        /**
         * \brief The scalar field type (usually double).
         */
        using Scalar = double;
        
        /**
         * \brief The spatial dimension (2 or 3).
         */
        static int const dim = dimension;
        
        /**
         * \brief A 3-dimensional vector type.
         */
        using Invariants = Dune::FieldVector<Scalar,3>;
        
        /**
         * \brief Constructor.
         * \param c1 material parameter (factor for I1)
         * \param c2 material parameter (factor for I2)
         */
        MooneyRivlin(double c1_, double c2_): c1(c1_), c2(c2_)
        {}
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Invariants const& i)
        {
          i1 = i[0];
          i2 = i[1];
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(I) \f$.
         */
        Scalar d0() const
        {
          return c1*i1 + c2*i2;
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(I)d_1I \f$.
         */
        Scalar d1(Invariants const& di1) const
        {
          return c1*di1[0] + c2*di1[1];
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(I)\d_1I d_2I \f$.
         */
        Scalar d2(Invariants const& di1, Invariants const& di2) const
        {
          return 0;
        }
        
      private:
        double i1, i2, c1, c2;
      };
      
      
      /**
       * \ingroup stationaryElasticity
       * \brief Neo-Hookean material law formulated in terms of the shifted invariants \f$ i_1, i_2, i_3 \f$ of the doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * The Neo-Hookean hyperelastic material energy for incompressible materials is defined in terms of the Cauchy-Green strain tensor \f$ C \f$ as
       * \f[ W = \frac{G}{2} (I_1-d) = \frac{G}{2} i_1. \f]
       * The material parameter \f$ G \f$ corresponds to the shear modulus of linear elasticity. The material law is a special case of the more general
       * MooneyRivlin material law.
       * 
       * The implementation is numerically stable in the vicinity of the reference configuration.
       * 
       * Use the \ref InvariantsMaterialLaw adaptor for creating a hyperelastic material law, 
       * or the \ref CompressibleInvariantsMaterialLaw adaptor in case of compressible materials.
       * 
       * Model of InvariantsMaterialConcept.
       * 
       * \see \ref InvariantsMaterialLaw
       */
      template <int dimension>
      class NeoHookean: public MooneyRivlin<dimension>
      {
      public:
        using typename MooneyRivlin<dimension>::Invariants;

        /**
         * \brief Constructor.
         * \param moduli the elastic moduli of the material (only shear modulus is relevant)
         */
        NeoHookean(ElasticModulus const& moduli): MooneyRivlin<dimension>(moduli.shear()/2,0)
        {}
      };
      
      /**
       * \ingroup stationaryElasticity
       * \brief Blatz-Ko material law for rubber foams in terms of the shifted invariants \f$ i_1, i_2, i_3 \f$ of the doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * In contrast to many other material laws given in terms of the invariants, the Blatz-Ko material model for foams is explicitly formulated for compressible
       * materials. The stored energy density is \f[ W = \frac{\mu}{2} \left(i_1 + \frac{2}{a} ( (1+i_3)^{-a/2}-1 ) \right), \f] with \f$ \mu \f$ the shear modulus
       * and \f$ a = 2\nu / (1-2\nu) \f$ related to Poisson's ratio.
       * 
       * Use the \ref InvariantsMaterialLaw adaptor for creating a hyperelastic material law.
       * 
       * Model of InvariantsMaterialConcept.
       * 
       * \see \ref InvariantsMaterialLaw
       */
      template <int dimension>
      class BlatzKo
      {
      public:
        /**
         * \brief The scalar field type (usually double).
         */
        using Scalar = double;
        
        /**
         * \brief The spatial dimension (2 or 3).
         */
        static int const dim = dimension;
        
        /**
         * \brief A 3-dimensional vector type.
         */
        using Invariants = Dune::FieldVector<Scalar,3>;
        
        /**
         * \brief Constructor.
         * \param moduli The elastic moduli.
         */
        BlatzKo(ElasticModulus const& moduli): mu(moduli.shear()), a(2*moduli.poisson()/(1-2*moduli.poisson()))
        {}
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Invariants const& i)
        {
          i1 = i[0];
          i3 = i[dimension-1];
          pstable = Pstable(-a/2,i3);
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(I) \f$.
         */
        Scalar d0() const
        {
          return mu*(i1 + 2*pstable.d0()/a)/2;
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(I)d_1I \f$.
         */
        Scalar d1(Invariants const& di1) const
        {
          return mu*(di1[0] + 2*pstable.d1(di1[dimension-1])/a)/2;
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(I)\d_1I d_2I \f$.
         */
        Scalar d2(Invariants const& di1, Invariants const& di2) const
        {
          return mu*pstable.d2(di1[dimension-1],di2[dimension-1])/a;
        }
        
      private:
        double i1, i3, mu, a;
        Pstable pstable;
      };
      
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup viscoPlasticity
       * \brief An adaptor for using hyperelastic stored energies for viscoplasticity.
       * 
       * The stored energy density for viscoplasticity is defined as \f$ W^{vp}(E) = W(E-E^{vp}). \f$ 
       * The viscoplastic strain \f$ E^{vp} \f$ acts as an offset, such that the material is stress-free
       * in a deformed configuration. In time-dependent viscoplastic models, the viscoplastic strain 
       * usually evolves according to some law depending on the current stress.
       * 
       * Model of HyperelasticMaterialLaw.
       * 
       * \tparam HyperelasticEnergy the material law of elastic response
       */
      template <class HyperelasticEnergy>
      class ViscoPlasticEnergy: public HyperelasticEnergy
      {
      public:
        using Scalar = typename HyperelasticEnergy::Scalar;
        static int const dim = HyperelasticEnergy::dim;
        
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Default constructor.
         */
        ViscoPlasticEnergy() = default;
        
        /**
         * \brief Constructor.
         * \param evp The viscoplastic stain.
         * The remaining parameters are forwarded to the hyperelastic energy density constructor.
         */
        template<typename... Args>
        ViscoPlasticEnergy(Tensor const& evp_, Args&&... args): HyperelasticEnergy(std::forward<Args>(args)...), evp(evp_)
        {}
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Tensor const& e)
        {
          HyperelasticEnergy::setLinearizationPoint(e-evp);
        }
                
      private:
        Tensor             evp;
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup viscoPlasticity
       * \brief Computes the von Mises equivalent stress \f$ \sigma_v \f$.
       * 
       * The equivalent von Mises stress is defined as 
       * \f[ 2\sigma_v^2 = (\sigma_{11}-\sigma_{22})^2 + (\sigma_{22}-\sigma_{33})^2 + (\sigma_{33}-\sigma_{11})^2 + 6 (\sigma_{23}^2+\sigma_{31}^2 + \sigma_{12}^2) \f]
       * in terms of the Cauchy stress \f$ \sigma \f$.
       * 
       * \tparam dim the spatial dimension. For d=2, plane stress is assumed.
       * \tparam Scalar the scalar field type of the stress tensor, usually double
       * 
       * Assuming there are finite element functions lambda and mu giving the material parameters, and u the displacement, a scalar finite element function sv of
       * von Mises equivalent stress values can be obtained as follows:
       * \code
       * interpolateGlobally<PlainAverage>(sv,makeFunctionView(u.space(), [&] (auto const& evaluator)
       * {
       *   ElasticModulus modulus(lambda.value(evaluator.cell(),evaluator.xloc()),
       *                          mu.value(evaluator.cell(),evaluator.xloc()));
       *   HyperelasticVariationalFunctional<StVenantKirchhoff<dim>,GreenLagrangeTensor<double,dim>>  material(modulus);
       *   material.setLinearizationPoint(u.derivative(evaluator));
       *   return Dune::FieldVector<double,1>(vonMisesStress(material.cauchyStress()));
       * }));
       * \endcode
       */
      template <int dim, class Scalar>
      Scalar vonMisesStress(Dune::FieldMatrix<Scalar,dim,dim> const& stress);

      // ---------------------------------------------------------------------------------------------------------

      /**
       * \ingroup viscoPlasticity
       * \brief A simple Duvaut-Lions flow rule with \f$ J_2 \f$ (von Mises) yield surface.
       * 
       * The Duvaut-Lions model of viscoplasticity defines the flow rate, i.e. the time derivative of the viscoplastic strain,
       * as \f[ \dot \epsilon^{\rm vp} = \tau^{-1} C^{-1} (\sigma - P\sigma), \f]
       * where \f$ P \f$ is the closest point projector onto the admissible stress state. In \f$ J_2 \f$ viscoplasticity, the 
       * admissible stresses are characterized by \f$ \|\sigma\| \le \sqrt{2/3}\,\sigma_Y \f$, and the boundary of that set 
       * is known as von Mises yield surface.
       * 
       * \param cauchyStress the true stress tensor
       * \param C the strain to stress matrix in Voigt notation
       * \param tau the relaxation time (>0)
       * \param sigmaY the yield strength (>=0)
       * \return the flow rate \f$ \dot \epsilon^{\rm vp} \f$ in Voigt notation
       */
      template <int dim, class Scalar=double>
      Dune::FieldVector<Scalar,dim*(dim+1)/2> 
      duvautLionsJ2Flow(Dune::FieldMatrix<Scalar,dim,dim> const& cauchyStress, Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& C, double tau, double sigmaY);
      
      // ---------------------------------------------------------------------------------------------------------
      
    }
  }
}

#endif