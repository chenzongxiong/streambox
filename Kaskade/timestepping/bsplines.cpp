#include <iostream>

#include "bsplines.hh"

namespace Kaskade {
  
  template <class Scalar>
  void extendedNodeSet(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar>& tau, int p)
  {
    int const n = t.size()-1;
    tau.resize(n+p+2);
    
    // Compute the extended node set such as to have a reasonable, possibly symmetric set of basis functions.
    // To this extent, we have to add p+1 nodes symmetrically at the left and right interval boundaries.
    // For even degree p, the basis functions are centered around interval *midpoints*, hence we then form 
    // the dual grid.
    double const spacer = 1.0/n;
    if (p%2 == 0) // even degree
    {
      int const k = p/2+1;
      for (int i=0; i<k; ++i)
        tau[i] = t[0] - (k-i)*spacer;
      for (int i=0; i<n; ++i)
        tau[i+k] = (t[i] + t[i+1])/2; // dual grid
      for (int i=0; i<k; ++i)
        tau[n+k+i] = t[n]+(i+1)*spacer;
    } 
    else // odd degree
    {
      int const k = (p+1)/2;
      for (int i=0; i<k; ++i)
        tau[i] = t[0] - (k-i)*spacer;
      for (int i=0; i<=n; ++i)
        tau[i+k] = t[i];
      for (int i=1; i<=k; ++i)
        tau[n+k+i] = t[n]+i*spacer;
    }
  }
  
  // Evaluates the B-spline basis functions belonging to the time grid t at points s
  // TODO: The implementation based on B-spline basis functions has the drawback that the corresponding Lagrangian
  //       basis (as evaluated in bsplineEvaluationMatrix) has global support. Look for Lagrangian spline bases with 
  //       local support.
  template <class Scalar>
  void bsplineBasisEvaluationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, 
				    int p, Dune::DynamicMatrix<Scalar>& e)
  {
    int const n = t.size()-1;
    int const m = s.size()-1;
    
    assert(p>=0);
    
    // Compute the extended node set such as to have a reasonable, possibly symmetric set of basis functions.
    Dune::DynamicVector<Scalar> tau(n+p+2); // contains the extended nodes
    extendedNodeSet(t,tau,p);
    Dune::DynamicVector<Scalar> b(n+p+2); // contains the basis functions values 
    
    e.resize(m+1,n+1);
    
    for (int i=0; i<=m; ++i) 
    {
      Scalar const x = s[i];
      // initialize the basis function's values for order 0
      for (int j=0; j<=n+p+1; ++j)
	b[j] = (tau[j]<=x && x<tau[j+1]) ? 1: 0; 
      
      // increase order one by one, using the de Boor/Cox/Mansfield recursion formula
      for (int q=1; q<=p; ++q)
	for (int j=0; j<=n+p-q; ++j)
	  b[j] = (x-tau[j])/(tau[j+q]-tau[j]) * b[j] + (tau[j+q+1]-x)/(tau[j+q+1]-tau[j+1]) * b[j+1];
	
      // copy result into matrix
      for (int j=0; j<=n; ++j)
	e[i][j] = b[j];
    }
  }
  
  
  // Evaluates the B-spline basis functions' derivatives belonging to the time grid t at points s
  template <class Scalar>
  void bsplineBasisDifferentiationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, int p,
					 Dune::DynamicMatrix<Scalar>& e)
  {
  }
  
  
  template <class Scalar>
  void bsplineDifferentiationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, int p, 
				    Dune::DynamicMatrix<Scalar>& d)
  {
    // Let At be the basis evaluation matrix at t and As the basis differentiation matrix at s.
    // Then As (At)^(-1) gives the desired result, transfering the values at ti to the basis coefficients
    // and then to the derivatives at si
    Dune::DynamicMatrix<Scalar> at, as;
    bsplineBasisEvaluationMatrix(t,t,p,at);
    bsplineBasisDifferentiationMatrix(t,s,p,as);
    at.invert();
    as.rightmultiply(at);
    d = as;
  }
  
  template <class Scalar>
  void bsplineEvaluationMatrix(Dune::DynamicVector<Scalar> const& t, Dune::DynamicVector<Scalar> const& s, int p,
			       Dune::DynamicMatrix<Scalar>& e)
  {
    // Let At be the basis evaluation matrix at t and As the basis evaluation matrix at s.
    // Then As (At)^(-1) gives the desired result, transfering the values at ti to the basis coefficients
    // and then to the values at si
    Dune::DynamicMatrix<Scalar> at, as;
    bsplineBasisEvaluationMatrix(t,t,p,at);
    bsplineBasisEvaluationMatrix(t,s,p,as);
    at.invert();
    as.rightmultiply(at);
    e = as;
  }
  
  // explicit instantiation
  template void bsplineDifferentiationMatrix<double>(Dune::DynamicVector<double> const& t, Dune::DynamicVector<double> const& s,
						     int p, Dune::DynamicMatrix<double>& d);
  template void bsplineEvaluationMatrix<double>(Dune::DynamicVector<double> const& t, Dune::DynamicVector<double> const& s, 
						int p, Dune::DynamicMatrix<double>& e);
};


#ifdef UNITTEST

#include <fstream>
#include <iostream>

int main(void)
{
  Dune::DynamicVector<double> t(11); // 11 equidistant interpolation nodes
  for (int i=0; i<t.size(); ++i)
    t[i] = static_cast<double>(i) / (t.size()-1);
  
  Dune::DynamicVector<double> s(101); // 101 equidistant interpolation nodes
  for (int i=0; i<s.size(); ++i)
    s[i] = static_cast<double>(i) / (s.size()-1);
  
  Dune::DynamicMatrix<double> e;
  
  Kaskade::bsplineBasisEvaluationMatrix(t,s,0,e);
  std::ofstream out0("bsplines0.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out0 << s[j] << ' ' << e[j][i] << '\n';
    out0 << "\n\n";
  }
  out0.close();
    
  Kaskade::bsplineEvaluationMatrix(t,s,1,e);
  std::ofstream out1("bsplines1.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out1 << s[j] << ' ' << e[j][i] << '\n';
    out1 << "\n\n";
  }
  out1.close();
    
  Kaskade::bsplineEvaluationMatrix(t,s,2,e);
  std::ofstream out2("bsplines2.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out2 << s[j] << ' ' << e[j][i] << '\n';
    out2 << "\n\n";
  }
  out2.close();
    
  Kaskade::bsplineEvaluationMatrix(t,s,3,e);
  std::ofstream out3("bsplines3.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out3 << s[j] << ' ' << e[j][i] << '\n';
    out3 << "\n\n";
  }
    

  Kaskade::bsplineEvaluationMatrix(t,s,4,e);
  std::ofstream out4("bsplines4.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out4 << s[j] << ' ' << e[j][i] << '\n';
    out4 << "\n\n";
  }
    

  Kaskade::bsplineEvaluationMatrix(t,s,5,e);
  std::ofstream out5("bsplines5.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out5 << s[j] << ' ' << e[j][i] << '\n';
    out5 << "\n\n";
  }
    
  Kaskade::bsplineEvaluationMatrix(t,s,6,e);
  std::ofstream out6("bsplines6.gnu");
  for (int i=0; i<t.size(); ++i) 
  {
    for (int j=0; j<s.size(); ++j)
      out6 << s[j] << ' ' << e[j][i] << '\n';
    out6 << "\n\n";
  }
    

  return 0;
}

#endif
