// Access on acml include file in dependence of the platform HTC-Linux or Mac OS.
// ZIBHTC is to set in the Makefile.Local (copy of Makefile-htc.Local) on HTC platform.
#ifdef ZIBHTC
#include "acml.h"
#endif
#if defined(Darwin) || defined(Cygwin) || defined(Atlas)
#include "linalg/acml_to_stdlapack.hh"
#endif




#include <iostream>

#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

#include "linalg/simpleLAPmatrix.hh"
#include "timestepping/sdc.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade {
  
  void eulerIntegrationMatrix(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix& Shat)
  {
    SDCTimeGrid::RealVector const& p = grid.points();
    int n = p.size()-1;
    
    Shat.resize(n,n+1,0);
    for (int i=0; i<n; ++i)
      Shat[i][i+1] = p[i+1]-p[i];
  }
  
  void luIntegrationMatrix(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix& Shat)
  {
    SDCTimeGrid::RealMatrix const&  S = grid.integrationMatrix();
    int n = S.N(); // number of rows
    
    // Compute transpose of S
    SLAPMatrix<double> A(n,n);
    for (int j=0; j<n; ++j)
      for (int i=0; i<n; ++i)
        A(i,j) = S[j][i+1];
        
    // Compute LU decomposition. Do this here by hand, for we need a factorization *without*
    // pivoting in order to retain the forward sweep structure of SDC (i.e., Shat is lower triangular),
    // and LAPACK does not provide an LU decomposition without pivoting. Note that performance is a non-issue here.
    // TODO: Check the impact of S = P L D D^{-1} U factorization with Shat = P L D on SDC.
    //       Probably a P != I corresponds to a kind of non-forward "sweep" that might actually 
    //       be better (due to pivoting). The reason is that pivoting minimizes the entries in L,
    //       which are the ones that remain in the iteration matrix I - Shat^{-1} * S. But it is unclear
    //       whether pivoting would result in P != I for the integration matrices.
    for (int k=0; k<n-1; ++k)
    {
      assert(A(k,k) != 0);
      for (int i=k+1; i<n; ++i)
      {
        // subtract row k from row i
        double a = A(i,k)/A(k,k); 
        for (int j=k+1; j<n; ++j)
          A(i,j) -= a*A(k,j);
        // do not store L entry as it is not needed: A(i,k) = a; 
      }
    }
    
    // extract the transpose of U factor as Shat and copy the very first column (corresponding to u(t0) which is fixed).
    Shat.resize(n,n+1,0);
    for (int i=0; i<n; ++i)
    {
      Shat[i][0] = S[i][0];
      for (int j=0; j<=i; ++j)
        Shat[i][j+1] = A(j,i);
    }
  }

  namespace {
    // Computes the n+1 Lobatto nodes on [-1,1]
    void computeLobattoNodes(int n, Dune::DynamicVector<double>& t)
    {
      t.resize(n+1);
      
      // The interior Lobatto points are eigenvalues of the symmetric tridiagonal matrix
      // with A_i,i+1 = (i+1)*(i+3)/(2*i+3)/(2*i+5) and every other uppper diagonal
      // entries zero. See, e.g., J. Shen, T. Tang and L. Wang, Spectral Methods:
      // Algorithms, Analysis and Applications, Springer Series in Compuational
      // Mathematics, 41, Springer, 2011.
      //
      // Of course, this works only for n at least 3.
     // if (n<=2)
      //{
	t[0] = -1;
	t[n] = 1;
	
	if (n==2)
	  t[1] = 0;
	//return t;
      //}
      
      // Build up matrix. Symmetric tridiagonal matrices are stored in two arrays, one 
      // for diagonal and one for subdiagonal.
      std::vector<double> diag(n-1,0.0), subdiag(n-2);
      for (int i=0; i<n-2; ++i)
        subdiag[i] = sqrt((i+1)*(i+3)/(2.*i+3)/(2.*i+5));
  
      // compute eigenvalues
      int info;
      dstev('N',n-1,&diag[0],&subdiag[0],nullptr,1,&info);
      assert(info>=0);
      if (info>0)
        throw LinearAlgebraException("Lapack routine dstev failed to converge",__FILE__,__LINE__);
      
      
      // write result
      std::copy(diag.begin(),diag.end(),&t[1]);
    }
    
    // Computes the n Radau nodes on ]-1,1], and stores them together with first point -1 in t
    void computeRadauNodes(int n, Dune::DynamicVector<double>& t)
    {
      t.resize(n+1);
      
      // The interior Radau points are the zeros of the Jacobi polynomials P^(1,0) and are obtained as 
      // eigenvalues of the symmetric tridiagonal matrix with A_00 = -1/3, A_jj = -1/(2j+1)/(2j+3), A_j(j-1) = 2j(j+1)/(2j+1)/sqrt(2j(2j+2))
      // See Deuflhard/Bornemann 6.3.2 and Wikipedia on Jacobi polynomials.
      
      std::vector<double> diag(n-1,-1.0/3), subdiag(std::max(0,n-2)); // take care of the special case n=1
      for (int i=1; i<n-1; ++i)
      {
        diag[i]      = -1.0/(2.0*i+1)/(2.0*i+3);
        subdiag[i-1] = 2.*i*(i+1)/(2.0*i+1)/std::sqrt(2*i*(2.0*i+2));
      }
      
      // compute eigenvalues
      int info;
      dstev('N',n-1,&diag[0],&subdiag[0],nullptr,1,&info);
      assert(info>=0);
      if (info>0)
        throw LinearAlgebraException("Lapack routine dstev failed to converge",__FILE__,__LINE__);
      
      // write result
      t[0] = -1;
      std::copy(diag.begin(),diag.end(),&t[1]);
      t[n] = 1;
    }
    
    // Evaluates the integral over the product of linear factors with zeros given in t
    // on the interval [a,b]
    double integrateLinearFactors(double a, double b, Dune::DynamicVector<double> const& t)
    {
      // This is done by representing the polynomial in a monomial basis w.r.t. x = t-a
      // and integrating this analytically.
      // let p(x) = \sum_{k=0}^n \alpha_k x^k. Then (t-t_i)*p(x) = (x+a-t_i)*p(x)
      // = \sum_{k=0}^{n} (a-t_i)*\alpha_k*x^k + \sum_{k=1}^{n+1} \alpha_{k-1} x^k.
      // Or, more explicitly: (a-t_i)*\alpha_0 + \sum_{k=1}^n (a-t_i*\alpha_k + \alpha_{k-1}) x^k
      // + \alpha_{n}x^{n+1}.
     int const n = t.size();
     Dune::DynamicVector<double> alpha(n+1);
     alpha[0] = 1.0;
     for (int i=0; i<n; ++i)
     {
	//double const c = c-t[i];
	double const c = a-t[i];
	alpha[i+1] = alpha[i];
	for (int k=i; k>=1; k--)
	  alpha[k] = c*alpha[k] + alpha[k-1];
	alpha[0] = alpha[0]*c;
      }
      
      // Now integrate the polynomial p over [0,b-a]
      double p = 1.0;
      double integral = 0;
      for (int i=0; i<=n; ++i) 
      {
	p *= b-a;
	integral += p/(i+1)*alpha[i];
      }
      
      // TODO: Performance improvement by not treating all Lagrange polynomials one by one
      //       but simultaneously.
      
      return integral;
      
    }
    
    // computes the integrals of the Lagrangian polynomials for the nodes s_j
    // over the intervals [t_i,t_{i+1}]
    void computeLagrangeIntegrationMatrix(Dune::DynamicVector<double> const& t,
                                          Dune::DynamicVector<double> const& s,
					  Dune::DynamicMatrix<double>& integ)
    {
      int const nt = t.size()-1; // number of intervals over which to integrate
      int const ns = s.size();   // number of Lagrangian polynomials which to integrate
      integ.resize(nt,ns);
      
      Dune::DynamicVector<double> lf(ns-1); // contains zeros of Lagrangian polynomials (linear factors)
      
      // step through the Lagrangian basis functions (through the interpolation nodes)
      for (int j=0; j<ns; ++j)
      {
	// The Lagrange polynomials consists of linear factors at the nodes
	// *except* s_j. The normalization constant is c.
	double c = 1.0;
	for (int k=0; k<ns-1; ++k)
	{
	  lf[k] = k<j? s[k]: s[k+1];
	  c *= s[j] - lf[k];
	}
	
	// Now integrate
	// i for time steps, intergration interval
	for (int i=0; i<nt; ++i)
	  integ[i][j] = integrateLinearFactors(t[i],t[i+1],lf)/c;
      }
    
      
    }
    
    // Evaluates p'(t) with p(x) = (x-lf[0])*...*(x-lf[n]).
    double differentiateLinearFactors(double t, Dune::DynamicVector<double> const& lf)
    {
      double d = 0;
      // straightforward product rule -- somewhat inefficient, but works
      for (int i=0; i<lf.size(); ++i)
      {
        double s = 1;
        for (int j=0; j<lf.size(); ++j)
          if (j!=i)
            s *= (t-lf[j]);
        d += s;
      }
      return d;
    }
    
    // evaluates the derivatives of the Lagrangian polynomials for the nodes s_j at points t_i
    void computeLagrangeDifferentiationMatrix(Dune::DynamicVector<double> const& t,
                                              Dune::DynamicVector<double> const& s,
                                              Dune::DynamicMatrix<double>& diff)
    {
      int const nt = t.size();   // number of evaluation points
      int const ns = s.size();   // number of Lagrangian polynomials which to integrate
      diff.resize(nt,ns);
      
      Dune::DynamicVector<double> lf(ns-1); // contains zeros of Lagrangian polynomials (linear factors)
      
      // step through the Lagrangian basis functions (through the interpolation nodes)
      for (int j=0; j<ns; ++j)
      {
        // The Lagrange polynomials consists of linear factors at the nodes
        // *except* s_j. The normalization constant is c.
        double c = 1.0;
        for (int k=0; k<ns-1; ++k)
        {
          lf[k] = k<j? s[k]: s[k+1];
          c *= s[j] - lf[k];
        }
        
        // Now differentiate
        // i for time steps, intergration interval
        for (int i=0; i<nt; ++i)
          diff[i][j] = differentiateLinearFactors(t[i],lf)/c;
      }
    }
    
    // Evaluates all Lagrangian interpolation polynomials to nodes t (n+1 of them) at the points in x (m+1 of them).
    // p is resized to (m+1) x (n+1).
    void computeLagrangeEvaluationMatrix(Dune::DynamicVector<double> const& t,
					 Dune::DynamicVector<double> const& x,
					 Dune::DynamicMatrix<double>& p)
    {
      int const n = t.size()-1;
      int const m = x.size()-1;
      
      p.resize(m+1,n+1);
      
      // step through all evaluation points
     for (int i=0; i<=m; ++i)
     { 
       // step through all interpolation nodes
       for (int j=0; j<=n; ++j){
         p[i][j] =1 ;
         for (int k=0; k<n; ++k)
         {
           //double tk = k<j? k: k+1;
           double tk = k<j? t[k]: t[k+1];
           
           p[i][j] *= (x[i]-tk)/(t[j]-tk);
         }
       }
     }
     // TODO: Use barycentric interpolation formula for performance (see Higham or Trefethen)
    }
  }
  
  
  LobattoTimeGrid::LobattoTimeGrid(int n, double a, double b)
  {
    assert(n>=1);
    computeLobattoNodes(n,pts);
    
    // rescale from [-1,1] to [a,b]
    for (int i=0; i<pts.size(); ++i)
      pts[i] = a + (pts[i]+1)*(b-a)/2;
    computeLagrangeIntegrationMatrix(pts,pts,integ);
    
    // interpolation nodes for differentiation are all, including 
    // the initial point 0.
    computeLagrangeDifferentiationMatrix(pts,pts,diff);
  }
  
  void LobattoTimeGrid::refine(Dune::DynamicMatrix<double>& p)
  {
    Dune::DynamicVector<double> newPts(pts.size()+1);
    computeLobattoNodes(pts.size(),newPts);
    // rescale from [-1,1] to [a,b]
    for (int i=0; i<newPts.size(); ++i)
      newPts[i] = pts[0] + (newPts[i]+1)*(pts[pts.size()-1]-pts[0])/2;
    computeLagrangeEvaluationMatrix(pts,newPts,p);
    pts = newPts;
    computeLagrangeIntegrationMatrix(pts,pts,integ);
  }
  
  SDCTimeGrid::RealMatrix LobattoTimeGrid::interpolate(SDCTimeGrid::RealVector const& x) const 
  {
    SDCTimeGrid::RealMatrix p;
    computeLagrangeEvaluationMatrix(points(),x,p);
    return p;
  }
  
  
  
  
  
  

  RadauTimeGrid::RadauTimeGrid(int n, double a, double b)
  {
    assert(n>=1);
    computeRadauNodes(n,pts);
    
    // rescale from [-1,1] to [a,b]
    for (int i=0; i<pts.size(); ++i)
      pts[i] = a + (pts[i]+1)*(b-a)/2;
    
    computeMatrices();
  }
  
  void RadauTimeGrid::refine(Dune::DynamicMatrix<double>& p)
  {
    Dune::DynamicVector<double> newPts(pts.size()+1);
    computeRadauNodes(pts.size(),newPts);
    // rescale from [-1,1] to [a,b]
    for (int i=0; i<newPts.size(); ++i)
      newPts[i] = pts[0] + (newPts[i]+1)*(pts[pts.size()-1]-pts[0])/2;
    computeLagrangeEvaluationMatrix(pts,newPts,p);
    pts = newPts;
    
    computeMatrices();
  }

  SDCTimeGrid::RealMatrix RadauTimeGrid::interpolate(SDCTimeGrid::RealVector const& x) const 
  {
    SDCTimeGrid::RealMatrix p;
    computeLagrangeEvaluationMatrix(points(),x,p);
    return p;
  }
  
  void RadauTimeGrid::computeMatrices()
  {
    int n = pts.size()-1;
    
    // interpolation nodes for integration are the true collocation points,
    // i.e. without initial point 0!
    RealVector nodes(n);
    for (int i=0; i<n; ++i)
      nodes[i] = pts[i+1];
    Dune::DynamicMatrix<double> rInteg;
    computeLagrangeIntegrationMatrix(pts,nodes,rInteg);
    // extend the integration matrix by leading zero column
    integ.resize(rInteg.N(),rInteg.M()+1);
    for (int i=0; i<integ.N(); ++i)
    {
      integ[i][0] = 0;
      for (int j=0; j<rInteg.M(); ++j)
        integ[i][j+1] = rInteg[i][j];
    }    
    
    // interpolation nodes for differentiation are all, including 
    // the initial point 0.
    computeLagrangeDifferentiationMatrix(pts,pts,diff);
  }

  namespace
  {
    double radauFlatN2D[2][2][2] = { { { 2.9291377, 0}, 
                                    {-1.4964389, 1.3083054} },
                                  { { 2.9213062, 0}, 
                                    {-1.4391869, 1.4727835} } };
    double radauFlatN2S[2][2][2] = { { {1.2602552, 0}, 
                                    {0.50174357, 0.57526967} },
                                  { {1.2843820, 0}, 
                                    {0.4954965, 0.51313156} } };

    double radauFlatN3D[2][3][3] = { { { 6.4261711,  0,          0}, 
                                       {-2.216832,   0.65451625, 0},
                                       {-0.2433249, -4.9218363,  2.1510773} },
                                     { {6.8756085, 0, 0}, 
                                       {-1.7188632, 2.034571, 0},
                                       {-0.60311261, -2.7236651, 2.7852899} } };
    double radauFlatN3S[2][3][3] = { { {1.3433028, 0, 0}, 
                                    {0.35672933, 0.47919602, 0},
                                    {-0.11823505, 0.28189846, 0.70369694} },
                                  { {2.2422795, 0, 0}, 
                                    {-0.2889824, 0.41983286, 0},
                                    {-0.26880299, 0.53050208, 0.44346113} } };

    double radauFlatN4D[2][4][4] = { { {11.208892,     0,           0,         0}, 
                                       {-3.0833747,    3.0367841,   0,         0},
                                       { 0.15793988,  -1.3403872,   1.5893103, 0},
                                       {-0.057982995, -0.65127407, -9.1627462, 3.4076761} },
                                     { {11.433291,       0,           0,         0 }, 
                                       {-3.517097,       2.0931311,   0,         0},
                                       {-0.00017653478, -1.8014578,   2.6269088, 0},
                                       {-0.39581998,    -0.19426875, -4.5863355, 4.6702588} } };
    double radauFlatN4S[2][4][4] = { { { 1.378026,     0, 0, 0 }, 
                                       { 0.37291834,   0.90022911, 0,        0},
                                       { 0.052159577,  1.4801239,  1.33701,  0},
                                       {-0.23067997,  -0.65234572, 3.175459, 0.8859888} },
                                     { { 3.0086247,  0,          0,          0}, 
                                       {-0.86590462, 0.46578692, 0,          0},
                                       { 1.4302848,  0.82388991, 0.40783232, 0},
                                       { 0.72369989, 0.08540081, 0.54316245, 0.45714202} } };
//     double radauFlatN5D[][][] = { { {}, 
//                                     {},
//                                     {},
//                                     {},
//                                     {} },
//                                   { {}, 
//                                     {},
//                                     {},
//                                     {},
//                                     {} } };
//     double radauFlatN5S[][][] = { { {}, 
//                                     {},
//                                     {},
//                                     {},
//                                     {} },
//                                   { {}, 
//                                     {},
//                                     {},
//                                     {},
//                                     {} } };
                                    
    template <int n>                                
    std::pair<SDCTimeGrid::RealMatrix,SDCTimeGrid::RealMatrix>
    getData(double const d[n][n], double const s[n][n])
    {
      SDCTimeGrid::RealMatrix dm(n,n), sm(n,n);
      for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
        {
          dm[i][j] = d[i][j];
          sm[i][j] = s[i][j];
        }
      return std::make_pair(dm,sm);
    }
  };
  
  std::pair<SDCTimeGrid::RealMatrix,SDCTimeGrid::RealMatrix> 
  optimizedIntegrationMatrix(SDCTimeGrid const& grid, int k, IntegrationMatrixOptimizationWeight w)
  {
    int const n = grid.points().size();
    
    if (dynamic_cast<RadauTimeGrid const*>(&grid))
    {
      switch (n)
      {
        case 2: return getData<2>(radauFlatN2D[k%2],radauFlatN2S[k%2]);
        case 3: return getData<3>(radauFlatN3D[k%2],radauFlatN3S[k%2]);
        case 4: return getData<4>(radauFlatN4D[k%2],radauFlatN4S[k%2]);
        default: break;
      }
    }
    
    // nothing found - alarm
    throw LookupException("No precomputed optimized differentiation/integration matrices found.",__FILE__,__LINE__);
  }

  std::pair<SDCTimeGrid::RealMatrix,SDCTimeGrid::RealMatrix> 
  optimizedIntegrationMatrix(SDCTimeGrid const& grid, int k, SDCTimeGrid::RealMatrix const& Df, SDCTimeGrid::RealMatrix const& Sf, 
                             IntegrationMatrixOptimizationWeight w)
  {
    try 
    {
      return optimizedIntegrationMatrix(grid,k,w);
    }
    catch (LookupException const& ex)
    {
      return std::make_pair(Df,Sf);
    }
  }
}

//#define UNITTEST
#ifdef UNITTEST

#include <algorithm>
#include <iostream>
#include <vector>

int main(void) 
{
  using namespace Kaskade;
  
  int const n = 4;
  
  std::cerr.precision(16);
 
  std::cerr << "Testing Lobatto time grid\n";
  {
    LobattoTimeGrid grid(n,0,1);
    
    std::cerr << "time points are: " << grid.points() << '\n';
    
    std::cerr << "checking integration\n";
    // check that the integral over 1 is 1
    Dune::DynamicMatrix<double> const& A = grid.integrationMatrix();
    double one = 0;
    for (int i=0; i<n; ++i)
      for (int j=0; j<=n; ++j)
        one += A[i][j];
    std::cerr << "1.0 = " << one << " (approximately)\n";
    
    std::cerr << "checking prolongation\n";
    // check that refinement and prolongation work
    Dune::DynamicMatrix<double>  P(n+2,n+1);
    grid.refine(P);
    for (int i=0; i<=n+1; ++i)
    {
      double one = 0;
      for (int j=0; j<=n; ++j)
        one += P[i][j];
      std::cerr << "1.0 = " << one << " (approximately)\n";
    }
  }
  
  std::cerr << "\nTesting Radau time grid\n";
  {
    RadauTimeGrid grid(n,0,1);
    
    std::cerr << "time points are: " << grid.points() << '\n';
    
    std::cerr << "checking integration\n";
    // check that the integral over 1 is 1
    Dune::DynamicMatrix<double> const& A = grid.integrationMatrix();
    std::cerr << "integration matrix is \n" << A ;
    double one = 0;
    for (int i=0; i<n; ++i)
      for (int j=0; j<=n; ++j)
        one += A[i][j];
    std::cerr << "1.0 = " << one << " (approximately)\n";
    
    std::cerr << "checking prolongation\n";
    // check that refinement and prolongation work
    Dune::DynamicMatrix<double>  P(n+2,n+1);
    grid.refine(P);
    for (int i=0; i<=n+1; ++i)
    {
      double one = 0;
      for (int j=0; j<=n; ++j)
        one += P[i][j];
      std::cerr << "1.0 = " << one << " (approximately)\n";
    }
  }
  
  
  std::cout << "\nTesting SDC contraction\n";
  {
    // use Radau grid on [0,1]
    RadauTimeGrid grid(4,0,1);
    
    // LU based DIRK sweep
    SDCTimeGrid::RealMatrix Shat;
    luIntegrationMatrix(grid, Shat);
//     eulerIntegrationMatrix(grid, Shat);
    
    std::cerr << "grid: " << grid.points() << "\n";
    std::cerr << "Integration matrix:\n" << Shat;
    
    typedef Dune::FieldMatrix<double,1,1> Matrix;
    typedef Dune::FieldVector<double,1> Vector;
    int n = grid.points().N();
    std::vector<Vector> ui (n,1.0), rUi(n), Mdu(n), du(n);
    std::vector<Matrix> rDu(n,0.0); // no reaction term - zero derivative
    
    // "solver" is just division by scalar value
    auto solver = [](Matrix A, Vector& x, Vector b) { x[0] = b[0]/A[0][0]; };
    
    auto printout = [](std::vector<Vector> const& xs) { for (auto const& x: xs) {std::cout.width(7); std::cout << x[0] << "  "; } };
    auto printoutScalar = [](std::vector<double> const& xs) { for (auto const& x: xs) {std::cout.width(7); std::cout << x << "  "; } };

    // several steps SDC on the Dahlquist test equation y' = lambda y, starting at y0=1
    for (double lambda=-1e-2; lambda>-1e8; lambda*=1.1)
//     double lambda = -1e4;
    {
      std::vector<double> duNorm;
      
      std::cout.precision(3);
      int k;
      for (k=0; k<30; ++k)
      {
        // compute right hand side
        for (int i=0; i<ui.size(); ++i)        rUi[i] = lambda*ui[i];
        for (int i=1; i<ui.size(); ++i)        Mdu[i-1] = ui[i-1]-ui[i];
        duNorm.push_back( sdcIterationStep2(grid,Shat,solver, Matrix(1.0), Matrix(lambda),rUi,rDu,Mdu,du) );
        for (int i=0; i<ui.size(); ++i)        ui[i] += du[i];
        printout(du); std::cout << "\n";
        if (duNorm.back() < 1e-11) {
          std::cout << "l: " << lambda << "    "; printoutScalar(duNorm); std::cout << "\n";
          break;
        }
      }
//       printout(ui);      std::cout << "\n";
      std::cout << "lambda: " << lambda << "    contraction estimate: " << std::pow(duNorm[duNorm.size()-1]/duNorm[0],1.0/(duNorm.size()-1)) << " " << k << "\n";
    }
  }
  
  return 0;
}


#endif

