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

#include <map>

#include "fem/diffops/elasto.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"

namespace Kaskade 
{
  namespace Elastomechanics
  {
    // Writes a symmetric matrix as vector in Voigt notation with doubled off-diagonal entries. 
    // First the diagonal entries, then the doubled off-diagonal entries.
    template <class Scalar, int n>
    Dune::FieldVector<Scalar,n*(n+1)/2> pack(Dune::FieldMatrix<Scalar,n,n> const& e, Scalar s)
    {
      Dune::FieldVector<Scalar,n*(n+1)/2> c;
      
      for (int i=0; i<n; ++i)
        c[i] = e[i][i];
      
      if (n==3)
      {
        c[3] = s*e[1][2];
        c[4] = s*e[0][2];
        c[5] = s*e[0][1];
      }
      else if (n==2)
        c[2] = s*e[0][1];
      
      return c;
    }
    
    // explicit instantiaion of function
    template Dune::FieldVector<double,3> pack(Dune::FieldMatrix<double,2,2> const& e, double s);
    template Dune::FieldVector<double,6> pack(Dune::FieldMatrix<double,3,3> const& e, double s);
    
    template <class Scalar, int n>
    Dune::FieldMatrix<Scalar,ElastoDetails::dim(n),ElastoDetails::dim(n)> unpack(Dune::FieldVector<Scalar,n> const& v, Scalar s)
    {
      int const d = ElastoDetails::dim(n);
      Dune::FieldMatrix<Scalar,d,d> mat;
      
      for(int i=0; i<d; ++i) mat[i][i] = v[i];
      
      if(d==2) 
        mat[1][0] = mat[0][1] = v[2]/s;
      else if(d==3)
      {
        mat[1][0] = mat[0][1] = v[5]/s;
        mat[2][0] = mat[0][2] = v[4]/s;
        mat[1][2] = mat[2][1] = v[3]/s;
      }
      
      return mat;
    }
    
    template Dune::FieldMatrix<double,2,2> unpack(Dune::FieldVector<double,3> const& v, double s);
    template Dune::FieldMatrix<double,3,3> unpack(Dune::FieldVector<double,6> const& v, double s);
    
    
    
    // ---------------------------------------------------------------------------------------
    
    namespace
    {
      
      // general data from Wikipedia
      std::map<std::string,ElasticModulus> materialData{
        { "A36 steel",     ElasticModulus( 86.0e9, 79.3e9) }, // American structural steel
        { "aluminum",      ElasticModulus( 61.2e9, 25.5e9) },
        { "copper",        ElasticModulus( 45.0e9, 47.0e9) },
        { "glass",         ElasticModulus(460.0e6, 26.2e9) },
        { "magnesium",     ElasticModulus( 31.2e9, 17.0e9) },
        { "polyethylene",  ElasticModulus(151.4e6,117.0e6) },
        { "rubber",        ElasticModulus( 15.0e6,300.0e3) },
        { "steel",         ElasticModulus(107.1e9, 79.3e9) },
        { "titanium",      ElasticModulus( 81.9e9, 41.4e9) },
        { "bone",          ElasticModulus( 9.81e9, 6.54e9) },
        { "quartz",        ElasticModulus( 1.56307e10, 3.03419e+10) },
        { "polycarbonate", ElasticModulus( 3.57143e+09, 8.92857e+08) },
        { "polyamid12",    ElasticModulus( 8.65385e+08, 5.76923e+08) }
      };
      
      // general data from Wikipedia
      std::map<std::string,double> massDensityMap{
        { "A36 steel",     7800.0 }, // American structural steel
        { "aluminum",      2700.0 },
        { "concrete",      2100.0 },
        { "copper",        8920.0 },
        { "cork",           500.0 },
        { "glass",         2550.0 },
        { "magnesium",     1738.0 },
        { "polyethylene",   940.0 },
        { "rubber",         940.0 },
        { "silicon",       2330.0 },
        { "steel",         7870.0 },
        { "titanium",      4500.0 },
        { "bone",          2000.0 },
        { "quartz",        2204.0 },
        { "polycarbonate", 1200.0 },
        { "polyamid12",    1010.0 }
      };
      
      // general data from Wikipedia
      std::map<std::string,double> yieldStrengthMap {
        { "A36 steel",   250.0e6 } // American structural steel
      };
      
    }
    
    // ---------------------------------------------------------------------------------------
    
    std::map<std::string,double> const& massDensities() { return massDensityMap; }
    
    double massDensity(std::string const& name)
    {
      auto it = massDensityMap.find(name);
      if (it == massDensityMap.end())
        throw LookupException("Material " + name + " not found in mass density materials data base.",__FILE__,__LINE__);
      return it->second;
    }
    
    double yieldStrength(std::string const& name)
    {
      auto it = yieldStrengthMap.find(name);
      if (it == yieldStrengthMap.end())
        throw LookupException("Material " + name + " not found in yield strength materials data base.",__FILE__,__LINE__);
      return it->second;
    }
    
    // ---------------------------------------------------------------------------------------
    
    ElasticModulus& ElasticModulus::setYoungPoisson(double E, double nu)
    {
      assert(-1<nu && nu<0.5);
      assert(E>0);
      lambda = E*nu/(1+nu)/(1-2*nu);
      mu = E/2/(1+nu);
      return *this;
    }
    
    
    ElasticModulus const& ElasticModulus::material(std::string const& name)
    {
      auto it = materialData.find(name);
      if (it == materialData.end())
        throw LookupException("Material " + name + " not found in elastic materials data base.",__FILE__,__LINE__);
      return it->second;
    }
    
    std::map<std::string,ElasticModulus> const& ElasticModulus::materials() 
    {
      return materialData;
    }
    
    // ---------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------

    
    template <int d, class Scalar>
    Scalar DetIpm1<d,Scalar>::d0() const
    {
      if (d==2)
        return trace(A) + A.determinant();
      if (d==3)
        return A[2][2] + (A[1][1] + A[0][0]*(1+A[1][1]))*(1+A[2][2])
        + A[0][1]*A[1][2]*A[2][0]    // rule of Sarrus - subtracted 1 in the first (expanded) product
        + A[0][2]*A[1][0]*A[2][1]
        - A[2][0]*(1+A[1][1])*A[0][2]
        - A[2][1]*A[1][2]*(1+A[0][0])
        - (1+A[2][2])*A[1][0]*A[0][1];
    }
    
    template <int d, class Scalar>
    Scalar DetIpm1<d,Scalar>::d1(Dune::FieldMatrix<Scalar,d,d> const& B) const
    {
      if (d==2)
        return trace(B) + A[0][0]*B[1][1] + B[0][0]*A[1][1] - A[1][0]*B[0][1] - B[1][0]*A[0][1];
      if (d==3)
        return B[0][0]*(1+A[1][1])*(1+A[2][2]) + (1+A[0][0])*B[1][1]*(1+A[2][2]) + (1+A[0][0])*(1+A[1][1])*B[2][2]
        + B[0][1]*A[1][2]*A[2][0] + A[0][1]*B[1][2]*A[2][0]  + A[0][1]*A[1][2]*B[2][0]         
        + B[0][2]*A[1][0]*A[2][1] + A[0][2]*B[1][0]*A[2][1] + A[0][2]*A[1][0]*B[2][1]
        - B[2][0]*(1+A[1][1])*A[0][2] - A[2][0]*B[1][1]*A[0][2] - A[2][0]*(1+A[1][1])*B[0][2]
        - B[2][1]*A[1][2]*(1+A[0][0]) - A[2][1]*B[1][2]*(1+A[0][0]) - A[2][1]*A[1][2]*B[0][0]
        - B[2][2]*A[1][0]*A[0][1]- (1+A[2][2])*B[1][0]*A[0][1] - (1+A[2][2])*A[1][0]*B[0][1];
    }
    
    template <int d, class Scalar>
    Scalar DetIpm1<d,Scalar>::d2(Dune::FieldMatrix<Scalar,d,d> const& B, Dune::FieldMatrix<Scalar,d,d> const& C) const
    {
      if (d==2)
        return C[0][0]*B[1][1] + B[0][0]*C[1][1] - C[1][0]*B[0][1] - B[1][0]*C[0][1];
      if (d==3)
        return B[0][0]*C[1][1]*(1+A[2][2]) + B[0][0]*(1+A[1][1])*C[2][2] + C[0][0]*B[1][1]*(1+A[2][2]) + (1+A[0][0])*B[1][1]*C[2][2] + C[0][0]*(1+A[1][1])*B[2][2] + (1+A[0][0])*C[1][1]*B[2][2]
        + B[0][1]*(C[1][2]*A[2][0]+A[1][2]*C[2][0]) + B[1][2]*(C[0][1]*A[2][0]+A[0][1]*C[2][0])  + (C[0][1]*A[1][2]+A[0][1]*C[1][2])*B[2][0]         
        + B[0][2]*(C[1][0]*A[2][1]+A[1][0]*C[2][1]) + B[1][0]*(C[0][2]*A[2][1]+A[0][2]*C[2][1]) + (C[0][2]*A[1][0]+A[0][2]*C[1][0])*B[2][1]
        - B[2][0]*(C[1][1]*A[0][2]+(1+A[1][1])*C[0][2]) - B[1][1]*(C[2][0]*A[0][2]+A[2][0]*C[0][2]) - (C[2][0]*(1+A[1][1])+A[2][0]*C[1][1])*B[0][2]
        - B[2][1]*(C[1][2]*(1+A[0][0])+A[1][2]*C[0][0]) - B[1][2]*(C[2][1]*(1+A[0][0])+A[2][1]*C[0][0]) - (C[2][1]*A[1][2]+A[2][1]*C[1][2])*B[0][0]
        - B[2][2]*(C[1][0]*A[0][1]+A[1][0]*C[0][1]) - B[1][0]*(C[2][2]*A[0][1]+(1+A[2][2])*C[0][1]) - (C[2][2]*A[1][0]+(1+A[2][2])*C[1][0])*B[0][1];
    }
    
    // explicit instantiations
    template class DetIpm1<2,double>;
    template class DetIpm1<3,double>;
    
    // ---------------------------------------------------------------------------------------
    
    
    Pstable::Pstable(double p, double x)
    { 
      assert(x>-1); 
      if (std::abs(x) < 5e-6) // according to gnuplot graphs, this should be a reasonable threshold for double precision
      {
        f = p*x + p*(p-1)*x*x/2; // Taylor expansion around 0
        df = p + p*(p-1)*x;
        ddf = p*(p-1);
      }
      else
      {
        f = std::pow(x+1,p)-1;
        df = p*std::pow(x+1,p-1);
        ddf = p*(p-1)*std::pow(x+1,p-2);
      }
    }
    
    // ---------------------------------------------------------------------------------------
    
    template <int dim, class Scalar>
    typename ShiftedInvariants<dim,Scalar>::Invariants ShiftedInvariants<dim,Scalar>::d0() const
    {
      Invariants i;
      i[0] = trace(A);
      if (dim>1)
        i[dim-1] = det.d0();
      if (dim==3)
        i[1] = 2*i[0] + A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2] - A[0][1]*A[1][0] - A[0][2]*A[2][0] - A[1][2]*A[2][1];
      return i;
    }
    
    template <int dim, class Scalar>
    typename ShiftedInvariants<dim,Scalar>::Invariants ShiftedInvariants<dim,Scalar>::d1(Tensor const& dA) const
    {
      Invariants i;
      i[0] = trace(dA);
      if (dim>1)
        i[dim-1] = det.d1(dA);
      if (dim==3)
        i[1] = 2*i[0] + dA[0][0]*A[1][1] + A[0][0]*dA[1][1] + dA[0][0]*A[2][2] + A[0][0]*dA[2][2] + dA[1][1]*A[2][2]  + A[1][1]*dA[2][2] 
              - dA[0][1]*A[1][0] - A[0][1]*dA[1][0] - dA[0][2]*A[2][0] - A[0][2]*dA[2][0] - dA[1][2]*A[2][1] - A[1][2]*dA[2][1];
      return i;
    }
    
    template <int dim, class Scalar>
    typename ShiftedInvariants<dim,Scalar>::Invariants ShiftedInvariants<dim,Scalar>::d2(Tensor const& dA1, Tensor const& dA2) const
    {
      Invariants i;
      i[0] = 0;
      if (dim>1)
        i[dim-1] = det.d2(dA1,dA2);
      if (dim==3)
        i[1] = dA1[0][0]*dA2[1][1] + dA2[0][0]*dA1[1][1] + dA1[0][0]*dA2[2][2] + dA2[0][0]*dA1[2][2] + dA1[1][1]*dA2[2][2]  + dA2[1][1]*dA1[2][2] 
              - dA1[0][1]*dA2[1][0] - dA2[0][1]*dA1[1][0] - dA1[0][2]*dA2[2][0] - dA2[0][2]*dA1[2][0] - dA1[1][2]*dA2[2][1] - dA2[1][2]*dA1[2][1];
      return i;
    }
    
    // explicit instantiations
    template class ShiftedInvariants<2,double>;
    template class ShiftedInvariants<3,double>;
  }
}



