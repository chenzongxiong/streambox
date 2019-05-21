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

#include "elastoVariationalFunctionals.hh"
  
namespace Kaskade 
{
  namespace Elastomechanics
  {
    template <class Scalar, int dim>
    OrthotropicLameNavier<Scalar,dim>::OrthotropicLameNavier(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, Dune::FieldMatrix<Scalar,dim,dim> const& mat)
    : orth(orth_), du0(0)
    {
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> compliance(0);
      if (dim==3)
      {
        // see http://en.wikipedia.org/wiki/Orthotropic_material
        compliance[0][0] = 1 / mat[0][0];             // 1/E_11
        compliance[0][1] = - mat[1][0] / mat[0][0];   // - nu_12 / E_11 = - nu_21 / E_22
        compliance[0][2] = - mat[2][0] / mat[0][0];   // - nu_13 / E_11 = - nu_31 / E_33
        compliance[1][0] = compliance[0][1];          // symmetry
        compliance[1][1] = 1 / mat[1][1];
        compliance[1][2] = - mat[2][1] / mat[1][1];   // - nu_23 / E_22 = - nu_32 / E_33
        compliance[2][0] = compliance[0][2];          // symmetry
        compliance[2][1] = compliance[1][2];          // symmetry
        compliance[2][2] = 1 / mat[2][2];
        
        compliance[3][3] = 1 / mat[1][2];             // 1 / G_23
        compliance[4][4] = 1 / mat[0][2];             // 1 / G_13
        compliance[5][5] = 1 / mat[0][1];             // 1 / G_12
      } 
      else if (dim==2)
      {
        compliance[0][0] = 1 / mat[0][0];             // 1/E_11;
        compliance[0][1] = - mat[1][0] / mat[0][0];   // - nu_12 / E_11 = - nu_21 / E_22
        compliance[1][0] = compliance[0][1];          // symmetry
        compliance[1][1] = 1 / mat[1][1];             // 1/E_22
        
        compliance[2][2] = 1/mat[0][1];               // 1 / G_12
      }
      
      // stiffness tensor is inverse of compliance tensor
      stiffness = compliance;
      stiffness.invert();
    }
    // explicit instantiation
    template OrthotropicLameNavier<double,2>::OrthotropicLameNavier(Dune::FieldMatrix<double,2,2> const& orth_, Dune::FieldMatrix<double,2,2> const& mat_);
    template OrthotropicLameNavier<double,3>::OrthotropicLameNavier(Dune::FieldMatrix<double,3,3> const& orth_, Dune::FieldMatrix<double,3,3> const& mat_);
    
    template <class Scalar, int dim>
    OrthotropicLameNavier<Scalar,dim>::OrthotropicLameNavier(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, 
                                                             Dune::FieldMatrix<Scalar,dim,dim> const& mat1, 
                                                             Dune::FieldVector<Scalar,dim*(dim-1)/2> const& mat2)
    : orth(orth_), stiffness(0), du0(0)
    {
      for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
          stiffness[i][j] = mat1[i][j];
        for (int i=0; i<dim*(dim-1)/2; ++i)
          stiffness[dim+i][dim+i] = mat2[i];
    }
    // explicit instantiation
    template OrthotropicLameNavier<double,2>::OrthotropicLameNavier(Dune::FieldMatrix<double,2,2> const& orth_, Dune::FieldMatrix<double,2,2> const& mat1, Dune::FieldVector<double,1> const& mat2);
    template OrthotropicLameNavier<double,3>::OrthotropicLameNavier(Dune::FieldMatrix<double,3,3> const& orth_, Dune::FieldMatrix<double,3,3> const& mat1, Dune::FieldVector<double,3> const& mat2);
    
    
    template <class Scalar, int dim>
    OrthotropicLameNavier<Scalar,dim>::OrthotropicLameNavier(ElasticModulus const& p)
    : orth(unitMatrix<Scalar,dim>()), stiffness(0), du0(0) 
    {
      Scalar lambda = p.lame(), mu = p.shear();
      if (dim==2)
      {
        stiffness[0][0] = stiffness[1][1] = 2*mu + lambda;
        stiffness[1][0] = stiffness[0][1] = lambda;
        stiffness[2][2] = mu;
      }
      if (dim==3)
        abort();
    }
    // explicit instantiation
    template OrthotropicLameNavier<double,2>::OrthotropicLameNavier(ElasticModulus const& p);
    template OrthotropicLameNavier<double,3>::OrthotropicLameNavier(ElasticModulus const& p);
    
    
    
    template <class Scalar, int dim>
    Scalar OrthotropicLameNavier<Scalar,dim>::d0() const
    {
      auto eps = pack( LinearizedGreenLagrangeTensor<Scalar,dim>(du0*orth).d0() );
      auto sigma = stiffness * eps;
      return 0.5 * (sigma*eps); // remember that eps contains doubled off-diagonal entries. Hence the scalar product 
    }                           // of the Voigt vectors is just the tensor contraction.
    // explicit instantiation
    template double OrthotropicLameNavier<double,2>::d0() const;
    template double OrthotropicLameNavier<double,3>::d0() const;
    
    template <class Scalar, int dim>
    Dune::FieldVector<Scalar,dim> OrthotropicLameNavier<Scalar,dim>::d1(VariationalArg<Scalar,dim> const& arg) const
    {
      // TODO: I have doubts that this is the most efficient implementation...
      Dune::FieldVector<Scalar, dim> ret;
      LinearizedGreenLagrangeTensor<Scalar,dim> E(du0*orth);
      auto epsU = pack(E.d0());
      auto sigma = stiffness*epsU;
      
      for (int i=0; i<dim; ++i) {
        Dune::FieldMatrix<Scalar,dim,dim> dv(0);
        dv[i] = arg.derivative[0];
        auto epsV = pack( E.d1(dv*orth) );
        
        ret[i] = epsV*sigma; // remember that epsV contains doubled off-diagonal entries. Hence the scalar product 
      }                      // of the Voigt vectors is just the tensor contraction.
      return ret;
    }
    // explicit instantiation
    template Dune::FieldVector<double,2> OrthotropicLameNavier<double,2>::d1(VariationalArg<double,2> const &arg) const;
    template Dune::FieldVector<double,3> OrthotropicLameNavier<double,3>::d1(VariationalArg<double,3> const &arg) const;
    
    template <class Scalar, int dim>
    Dune::FieldMatrix<Scalar,dim,dim> OrthotropicLameNavier<Scalar,dim>::d2(VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      // TODO: I have doubts that this is the most efficient implementation...
      Dune::FieldMatrix<Scalar, dim, dim> ret;
      LinearizedGreenLagrangeTensor<Scalar,dim> E(du0*orth);
      
      for (int i=0; i<dim; ++i) {
        Dune::FieldMatrix<Scalar,dim,dim> dv(0);
        dv[i] = arg1.derivative[0];
        auto epsV = pack( E.d1(dv*orth) );
        auto sigmaV = stiffness*epsV;
        
        for (int j=0; j<dim; ++j) {
          Dune::FieldMatrix<Scalar,dim,dim> dw(0);
          dw[j] = arg2.derivative[0];
          auto epsW = pack( E.d1(dw*orth) );
          
          ret[i][j] = sigmaV*epsW; // remember that epsW contains doubled off-diagonal entries. Hence the scalar product 
        }                          // of the Voigt vectors is just the tensor contraction.
      }
      return ret;
    }
    // explicit instantiation
    template Dune::FieldMatrix<double,2,2> OrthotropicLameNavier<double,2>::d2(VariationalArg<double,2> const &arg1, VariationalArg<double,2> const &arg2) const;
    template Dune::FieldMatrix<double,3,3> OrthotropicLameNavier<double,3>::d2(VariationalArg<double,3> const &arg1, VariationalArg<double,3> const &arg2) const;
    
    // ---------------------------------------------------------------------------------------
    
  }
}



#ifdef UNITTEST

#include "fem/shapefunctioncache.hh"
#include <random>

using namespace Kaskade;
using namespace Kaskade::Elastomechanics;

template <int dim>
void checkDerivatives()
{
  using Matrix = Dune::FieldMatrix<double,dim,dim>;
  std::random_device rd;
  std::uniform_real_distribution<> rand(-0.1,0.1);
  
  Matrix A;
  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      A[i][j] = rand(rd);
    
  std::cout << "Checking shifted invariants derivatives at A=\n" << A << "\n";
    
    
  // Create shifted invariants
  ShiftedInvariants<dim> inv(A);
  std::cout << "value is: " << inv.d0() << "\n";
  
  // Check first derivative
  std::cout << "First derivative\n";
  for(int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
    {
      Matrix B(0);
      B[i][j] = 1e-7;
      ShiftedInvariants<dim> inv2(A+B);
      std::cout << "(" << i << "," << j << "): " << inv.d1(B) << " [should be " << inv2.d0()-inv.d0() << "]\n";
    }
     
  // Check second derivative
  std::cout << "Second derivative\n";
  for(int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      for (int k=0; k<dim; ++k)
        for (int l=0; l<dim; ++l)
        {
          Matrix B(0);
          B[i][j] = 1e-7;
          ShiftedInvariants<dim> inv2(A+B);
          Matrix C(0);
          C[k][l] = 1;
          std::cout << "(" << i << "," << j << ")(" << k << "," << l << "): " << inv.d2(B,C) << " [should be " << inv2.d1(C)-inv.d1(C) << "]\n";
        }
}


int main(void)
{
  
  checkDerivatives<2>();
  checkDerivatives<3>();
  
  int const dim = 2;
  
  ElasticModulus moduli;
  LameNavier<dim> lameNavier; lameNavier.setElasticModuli(moduli);
  OrthotropicLameNavier<double,dim> oLameNaiver(moduli);
  
  std::cout << "stiffness:\n" << oLameNaiver.stiffnessMatrix() << "\n";
  
  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
    {
      Dune::FieldMatrix<double,dim,dim> ux(0);
      ux[i][j] = 1;
      std::cout << "ux = \n" << ux << "\n";
      lameNavier.setLinearizationPoint(ux);
      oLameNaiver.setLinearizationPoint(ux);
      std::cout << "d0 = " << lameNavier.d0() << " <-> " << oLameNaiver.d0() << "\n";
      for (int k=0; k<dim; ++k)
      {
        Dune::FieldVector<double,1> v(0);
        Dune::FieldMatrix<double,1,dim> dv(0); dv[0][k] = 1;
        VariationalArg<double,dim> arg1(v,dv);
        std::cout << "d1 = " << lameNavier.d1(arg1) << " <-> " << oLameNaiver.d1(arg1) << "\n";
      }
      std::cout << "---------------------------------------------------------\n";
    }
}

#endif