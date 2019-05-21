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

#include <cmath>

#include "fem/diffops/materialLaws.hh"
#include "utilities/power.hh"

namespace Kaskade
{
  namespace Elastomechanics {
    
    namespace MaterialLaws
    {
      
      template <int dim, class Scalar>
      Scalar vonMisesStress(Dune::FieldMatrix<Scalar,dim,dim> const& s)
      {
        if (dim==2)
          return std::sqrt(0.5 * ( square(s[0][0]-s[1][1]) + square(s[0][0]) + square(s[1][1]) + 6*square(s[0][1]) ));
        if (dim==3)
          return std::sqrt(0.5 * ( square(s[0][0]-s[1][1]) + square(s[0][0]-s[2][2]) + square(s[1][1]-s[2][2]) 
                                   + 6*(square(s[0][1])+square(s[0][2])+square(s[1][2])) ));
        // never get here
        return -1;
      }
        
      // explicit instantiation
      template double vonMisesStress(Dune::FieldMatrix<double,2,2> const& stress);
      template double vonMisesStress(Dune::FieldMatrix<double,3,3> const& stress);
            
      // ---------------------------------------------------------------------------------------------------------------------------------
      
      template <int dim, class Scalar>
      Dune::FieldVector<Scalar,dim*(dim+1)/2> 
      duvautLionsJ2Flow(Dune::FieldMatrix<Scalar,dim,dim> const& cauchyStress, Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& C, double tau, double sigmaY)
      {
        using ReturnType = Dune::FieldVector<Scalar,dim*(dim+1)/2>;
        
        assert(tau>0);
        assert(sigmaY>=0);
        
        auto s = deviatoricPart(cauchyStress); 
        auto cSigma = s.frobenius_norm2();
        
        if (cSigma <= 2*sigmaY*sigmaY/3) // admissible stress state
          return ReturnType(0.0);
        
        // compute (s-barS)/tau with barS = Ps = sqrt(2/3)sigmaY/|s|
        s *= (1-std::sqrt(2.0/3/cSigma)*sigmaY) / tau;
        
        // compute deps = C^{-1} (s-barS)/tau
        ReturnType deps;      
        C.solve(deps,pack(s,1.0));
        
        return deps;
      }
      
      // explicit instantiation
      template Dune::FieldVector<double,3> 
      duvautLionsJ2Flow(Dune::FieldMatrix<double,2,2> const& cauchyStress, Dune::FieldMatrix<double,3,3> const& C, double tau, double sigmaY);
      template Dune::FieldVector<double,6> 
      duvautLionsJ2Flow(Dune::FieldMatrix<double,3,3> const& cauchyStress, Dune::FieldMatrix<double,6,6> const& C, double tau, double sigmaY);
    }
  }
}
