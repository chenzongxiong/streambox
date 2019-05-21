#include <cmath>
#include <iostream>

#include "dune/geometry/type.hh"
#include "dune/geometry/quadraturerules.hh"

#include "utilities/power.hh"

double fac(int n) 
{
  double f = 1;
  for (int i=2; i<=n; ++i)
    f *= i;
  return f;
}



int main(void) 
{
  using namespace Kaskade;
  // Tetrahedral integration
  Dune::GeometryType gt(Dune::GeometryType::simplex,3);

  typedef Dune::QuadratureRule<double,3> Rule;

  std::cout.precision(16);
  
  for (int p=0; p<15; ++p) {
          
    Rule const& rule = Dune::QuadratureRules<double,3>::rule(gt,p);
    size_t nQuadPos = rule.size();
    std::cout << "\n Integrating order " << p << " with " << nQuadPos << " points\n";
    

    // Integrate over all monomials.
    for (int a=0; a<=p; ++a) 
      for (int b=0; b<=p-a; ++b)
        for (int c=0; c<=p-a-b; ++c) {
          // This is the exact value:
          double exact = fac(a)*fac(b)*fac(c)/fac(a+b+c+3);
          
          double num = 0;
          for (size_t g=0; g<nQuadPos; ++g) {
            // pos of integration point
            Dune::FieldVector<double,3> const& quadPos = rule[g].position();
            double w = rule[g].weight();
            
            num += w * power(quadPos[0],a) * power(quadPos[1],b) * power(quadPos[2],c);
          }
          
          std::cout << "integ(x^" << a << "*y^" << b << "*z^" << c << ") -> " << num << " (should be: " << exact
                    << ", relative error " << std::abs(num-exact)/exact << ")\n";
        }
  }
  
          

  return 0;
}
