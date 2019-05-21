#ifndef ALIEV_HH
#define ALIEV_HH

#include <cmath>
#include <tr1/tuple>

#include "fem/fetransfer.hh"
#include "fem/fixdune.hh"
#include "fem/variables.hh"

template <bool c>
struct if_ 
{
  template <class A, class B>
  static A const& value(A const& a, B const& b) { return a; }
};

template <>
struct if_<false> 
{
  template <class A, class B>
  static B const& value(A const& a, B const& b) { return b; }
};

/**
 * A class that defines the weak formulation of the Aliev-Panfilov equation
 * \f[ \begin{bmatrix}\dot u\\ \dot v\end{bmatrix}
 *     = \begin{bmatrix} \Delta u + f(u,v) \\ g(u,v)  \end{bmatrix} \f] with Neumann boundary conditions.
 * The reaction terms are given by \f$ f(u,v)=\epsilon^{-1} u(1-u)(u-\frac{v+b}{a}) \f$ and \f$ g(u,v) = u-v \f$.
 */
template <class RType, class VarSet>
class AlievPanfilovEquation
{
  typedef AlievPanfilovEquation<RType,VarSet> Self;
  
public:
  typedef RType  Scalar;
  typedef VarSet OriginVars;
  typedef VarSet AnsatzVars;
  typedef VarSet TestVars;
  static ProblemType const type = WeakFormulation;

  typedef typename AnsatzVars::Grid::template Codim<0>::Entity Cell;
  typedef Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> Position;

private:
  static int const uIdx = 0;
  static int const vIdx = 1;
  static int const uSIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,uIdx>::type::spaceIndex;
  static int const vSIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,vIdx>::type::spaceIndex;

public:
  class DomainCache 
  {
  public:
    DomainCache(Self const& F_,
                typename AnsatzVars::VariableSet const& vars_,int flags=7):
      F(F_), vars(vars_)
    {}

    void moveTo(Cell const& entity) { e = &entity; }

    template <class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;

      u = at_c<uIdx>(vars.data).value(at_c<uSIdx>(evaluators));
      // Manipulation:
      if (u < 0.0) u = 0.0;
      
      du = at_c<uIdx>(vars.data).gradient(at_c<uSIdx>(evaluators))[0];
      assert(finite(u));
      assert(finite(du[0]) && finite(du[1]));

      v = at_c<vIdx>(vars.data).value(at_c<vSIdx>(evaluators));

      theta = F.theta(u);
      dtheta = F.dtheta(u);

      assert(finite(theta));
      assert(finite(dtheta));
      
      ex = F.excitation(e->geometry().global(x));
    }

    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      if (row==uIdx) return -F.kappa*dtheta*(du*arg.gradient[0]) + f(theta,v)*arg.value + ex*arg.value;
      if (row==vIdx) return g(theta,v)*arg.value;
      assert("wrong index\n"==0);
      return 0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      if (row==uIdx && col==uIdx) return -F.kappa*dtheta*(arg1.gradient[0]*arg2.gradient[0])
                                         +  F.reactionImplicitFactor * fu(theta,v)*dtheta*arg1.value*arg2.value;
      if (row==uIdx && col==vIdx) return fv(theta,v)*arg1.value*arg2.value;
      if (row==vIdx && col==uIdx) return gu(theta,v)*dtheta*arg1.value*arg2.value;
      if (row==vIdx && col==vIdx) return gv(theta,v)*arg1.value*arg2.value;

      assert("wrong index\n"==0);
      return 0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    b2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      if (row==uIdx && col==uIdx) return dtheta*arg1.value[0]*arg2.value[0];
      else if (row==col)          return arg1.value[0]*arg2.value[0];
      else                        return 0;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;

    Cell const* e;

    Scalar f(Scalar u, Scalar v)  const { return -(F.ga*u*(u-F.a)*(u-1.0) + u*v); }
    Scalar fu(Scalar u, Scalar v) const { Scalar  Fu = -F.ga*( (2*u-F.a)*(u-1.0) + u*(u-F.a) ) - v; 
    						  if (Fu > 0.0) Fu=0.0; return Fu;}
    Scalar fv(Scalar u, Scalar v) const { return -u; }
    
    Scalar g(Scalar u, Scalar v)  const { return 0.25*(F.eps1+F.mu1*v/(u+F.mu2))*(-v-F.gs*u*(u-F.a-1.0)); }
    Scalar gu(Scalar u, Scalar v) const { return 0.25* 
    						  		( - F.mu1*v/((u+F.mu2)*(u+F.mu2))*(-v-F.gs*u*(u-F.a-1)) 
    						  		  - (F.eps1+F.mu1*v/(u+F.mu2))*(F.gs*(2*u-F.a-1)) ); }
    Scalar gv(Scalar u, Scalar v) const { Scalar Gv =  0.25*( F.mu1/(u+F.mu2)*(-v-F.gs*u*(u-F.a-1)) - (F.eps1+F.mu1*v/(u+F.mu2)) ); 
    						  if (Gv > 0.0) {
    						  	 //std::cout << "Gv = " << Gv << std::endl;
    						  	 Gv = 0.0;
    						  }
    						  return Gv;}

    Scalar theta, dtheta;
    Scalar u, v, ex;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> du;
  };

  class BoundaryCache 
  {
  public:
    typedef typename AnsatzVars::Grid::template Codim<0>::Entity::LeafIntersectionIterator FaceIterator;

    BoundaryCache(Self const&,
                  typename AnsatzVars::VariableSet const&,int flags=7)
    {}

    void moveTo(FaceIterator const& entity) {}

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators) { }

    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m> d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return 0; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return 0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    b2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return 0;
    }
  };

  class Scaling 
  {
  public:
    Scaling(Self const& self_): self(&self_) {}

    template <class Vals>
    void scale(Cell const&, Position const&, Vals& vals) const 
    {
      double u = boost::fusion::at_c<uIdx+2>(vals);
      double dTheta = self->dtheta(u);
      boost::fusion::at_c<0>(vals) *= dTheta;

      double theta = self->theta(u);
      boost::fusion::at_c<2>(vals) = theta;
      
    }

  private:
    Self const* self;
  };

  Scalar zetaInitial(Scalar x) const {
    assert(finite(x));
    if (x>1) return zeta(1)+(x-1)*dzetaInitial(1);
    if (x<0) return zeta(0)+x*dzetaInitial(0);
    
    int i; Scalar l, r; idx(x,zetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<0>(zetaData[i])+r*std::tr1::get<0>(zetaData[i+1]);

    // Define zeta in terms of theta. zeta is only used for transforming the initial value.
    Scalar dresult = -(theta(result)-x)/dtheta(result);
    do {
      result += dresult;
      dresult = -(theta(result)-x)/dtheta(result);
    } while(dresult>=1e-10);
    
    
    assert(finite(result));
    return result;
  }

  Scalar dzetaInitial(Scalar x) const {
    assert(finite(x));
    x = std::max(0.0,std::min(x,1.0));
    
    int i; Scalar l, r; idx(x,zetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<1>(zetaData[i])+r*std::tr1::get<1>(zetaData[i+1]);

    // Define zeta' in terms of theta and theta'
    
    assert(finite(result));
    return result;
  }

  Scalar thetaInitial(Scalar x) const {
    
    assert(finite(x));
    if (x>1) return theta(1)+(x-1)*dtheta(1);
    if (x<0) return theta(0)+x*dtheta(0);
    
    int i; Scalar l, r; idx(x,thetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<0>(thetaData[i])+r*std::tr1::get<0>(thetaData[i+1]);
    assert(finite(result));
    return result;
  }
  

  Scalar zeta(Scalar x) const {
    
    //if (time() > 0.015) return x;			// switch on if no scaling
    return x; 								// if zeta(v) is identity
    
    assert(finite(x));
    if (x>1) return zeta(1)+(x-1)*dzeta(1);
    if (x<0) return zeta(0)+x*dzeta(0);
    
    int i; Scalar l, r; idx(x,zetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<0>(zetaData[i])+r*std::tr1::get<0>(zetaData[i+1]);

    // Define zeta in terms of theta. zeta is only used for transforming the initial value.
    Scalar dresult = -(theta(result)-x)/dtheta(result);
    do {
      result += dresult;
      dresult = -(theta(result)-x)/dtheta(result);
    } while(dresult>=1e-10);
    
    
    assert(finite(result));
    return result;
  }
  
  Scalar dzeta(Scalar x) const {
    
    //if (time() > 0.015) return 1.0;			// switch on if no scaling
    return 1.0;									// if zeta(v) is identity

    assert(finite(x));
    x = std::max(0.0,std::min(x,1.0));
    
    int i; Scalar l, r; idx(x,zetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<1>(zetaData[i])+r*std::tr1::get<1>(zetaData[i+1]);

    // Define zeta' in terms of theta and theta'
    
    assert(finite(result));
    return result;
  }
  
  
  Scalar theta(Scalar x) const {
    
    //if (time() > 0.015) return x;		// switch on if no scaling
    return x;							// if theta(v) is identity
    
    assert(finite(x));
    if (x>1) return theta(1)+(x-1)*dtheta(1);
    if (x<0) return theta(0)+x*dtheta(0);
    
    int i; Scalar l, r; idx(x,thetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<0>(thetaData[i])+r*std::tr1::get<0>(thetaData[i+1]);
    assert(finite(result));
    return result;
  }
  
  Scalar dtheta(Scalar x) const {

    //if (time() > 0.015) return 1.0;			// switch on if no scaling
    return 1.0;									// if theta(v) is identity

    assert(finite(x));
    x = std::max(0.0,std::min(x,1.0));
    
    int i; Scalar l, r; idx(x,thetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<1>(thetaData[i])+r*std::tr1::get<1>(thetaData[i+1]);
    // Scalar result = (std::tr1::get<0>(thetaData[i+1])-std::tr1::get<0>(thetaData[i]))*(thetaData.size()-1);
    assert(finite(result));
    return result;
  }
  
  Scalar ddtheta(Scalar x) const {
    assert(finite(x));
    if (x<0 || x>1) return 0;
    
    int i; Scalar l, r; idx(x,thetaData.size(),i,l,r);
    Scalar result = l*std::tr1::get<2>(thetaData[i])+r*std::tr1::get<2>(thetaData[i+1]);
    assert(finite(result));
    return result;
  }
  

  void idx(Scalar x, int n, int& id, Scalar& l, Scalar& r) const 
  {
    assert(0<=x && x<=1);
    id = static_cast<int>(std::floor((n-1)*x));
    if (id<0 || id>n-1) {
      std::cerr << "Assertion: x=" << x << " id=" << id << "\n";
      abort();
    }
    
    if (id==n-1) --id;

    r = (n-1)*x - id;
    l = 1-r;
  }


  Scalar dthetaMax() const { return maxDTheta; }
  
  
public:

  void newTime(Scalar dt) { t += dt; }
  Scalar time() const { return t; }
  void time(Scalar tnew)  { t=tnew; }
  Scaling scaling() const { return Scaling(*this); }
  void temporalEvaluationRange(double t0, double t1) { tau = t1-t0; }
  
  

  template <int row>
  struct D1 
  {
    static bool const present   = true;
    static bool const constant  = false;
    static int const derivatives = 1;
  };

  template <int row, int col>
  struct D2 
  {
    static bool const present   = ((row==col) || (row==vIdx && col==uIdx));
    								//((row==col) || (row==vIdx && col==uIdx));			
    								// row==col;		
    								// true;   
    static bool const symmetric = row==col;
    static bool const constant  = false;
    static bool const lumped = false;
    static int const derivatives = 1;
  };

  template <int row, int col>
  struct B2 
  {
    static bool const present   = row==col;
    static bool const symmetric = true;
    static bool const constant  = true;		// Ponosca: row==vIdx;
    static bool const lumped = false;
  };

  AlievPanfilovEquation(Scalar xi_, Scalar kappa_, Scalar a_, Scalar gs_, Scalar ga_, Scalar mu1_, Scalar mu2_, Scalar eps1_,
  				  Scalar v_rest_, Scalar v_peak_,
                  std::string const& zetaFile, std::string const& thetaFile):
    			  t(0), xi(xi_), kappa(kappa_), a(a_), gs(gs_), ga(ga_), 
    			  mu1(mu1_), mu2(mu2_), eps1(eps1_), v_rest(v_rest_), v_peak(v_peak_),
    			  maxDTheta(1), reactionImplicitFactor(1) {
    if (zetaFile.size()>0) {
      reactionImplicitFactor = 0.7;
      std::ifstream in(zetaFile.c_str());
      std::string line;

      do {
        std::getline(in,line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;
        std::istringstream in_line(line);
        Scalar x, z, dz;
        in_line >> x >> z >> dz;
        zetaData.push_back(std::tr1::make_tuple(z,dz));
      } while(!in.eof());
    } else {
      zetaData.push_back(std::tr1::make_tuple(0,1));
      zetaData.push_back(std::tr1::make_tuple(1,1));
    }
    
    if (thetaFile.size()>0) {
      std::ifstream in(thetaFile.c_str());
      std::string line;

      do {
        std::getline(in,line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;
        std::istringstream in_line(line);
        Scalar x, t;
        in_line >> x >> t;
        thetaData.push_back(std::tr1::make_tuple(t,0.0,0.0));
      } while(!in.eof());

      double tau = 1.0/(thetaData.size()-1);

      using namespace std::tr1;
      
      for (int i=1; i<thetaData.size()-1; ++i) {
        get<1>(thetaData[i]) = (get<0>(thetaData[i+1])-get<0>(thetaData[i-1]))/(2*tau);
        maxDTheta = std::max(maxDTheta,get<1>(thetaData[i]));
      }
      get<1>(thetaData[0]) = get<1>(thetaData[1]);
      get<1>(thetaData.back()) = get<1>(thetaData[thetaData.size()-2]);
      
      for (int i=1; i<thetaData.size()-1; ++i) {
        get<2>(thetaData[i]) = (get<1>(thetaData[i+1])-get<1>(thetaData[i-1]))/(2*tau);
      }
      get<2>(thetaData[0]) = 0;
      get<2>(thetaData.back()) = 0;
        
    } else {
      thetaData.push_back(std::tr1::make_tuple(0,1,0));
      thetaData.push_back(std::tr1::make_tuple(1,1,0));
    }
  }

  void check() const 
  {
    std::ofstream check("check.gnu");
    for (int i=0; i<=12000; ++i) {
      double h = 1e-5;
      double x = -.1 + 1.2*i/12000.0;
      check << x << ' ' << zeta(x) << ' ' << theta(x) << ' ' << dtheta(x) << ' ' << (theta(x+h)-theta(x-h))/(2*h) << ' ' 
            << ddtheta(x) << ' ' << (dtheta(x+h)-dtheta(x-h))/(2*h) << '\n';
    }
  }
  

  /**
   * Given an initial value, this is transfered to a properly scaled
   * finite element iterate.
   */
  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& u0, typename AnsatzVars::VariableSet& u) const 
  {
    interpolateGloballyWeak<Volume>(boost::fusion::at_c<row>(u.data),
                                      ScaledFunction<WeakFunctionView>(row==uIdx,u0,*this));
  }
  

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 0;
    else
      return 2*shapeFunctionOrder+1;
  }

private:
  Scalar t, tau;  
  Scalar xi, kappa, a, gs, ga, mu1, mu2, eps1, v_rest, v_peak;
  std::vector<std::tr1::tuple<Scalar,Scalar> >    zetaData;
  std::vector<std::tr1::tuple<Scalar,Scalar,Scalar> > thetaData;
  Scalar maxDTheta;
  double reactionImplicitFactor;

  Scalar excitation(Position const& x) const 
  {
  	Scalar Iapp = 0.0;
  	
	  // slab
	  
//     if(t<0.015) 
//       	if ( (x[0]>=1.4) && (x[0]<=1.6) && (x[1]>=1.4) && (x[1]<=1.6) && (x[2]>=1.05) && (x[2]<=1.1) ) {
//       		Iapp=10000.0*t;
//     		if (Iapp > 100.0) Iapp = 100.0;
//       	}
      

	// ventricles
  	if (t <= 1.0) 
  	{
  		double rX = x[0] - (-0.9); 
  		double rY = x[1] - 1.7; 
  		double rZ = x[2] - 1.1;
  		double r = sqrt(rX*rX + rY*rY + rZ*rZ);
   	 	if ( r < 1.0 )	{
   	 		Iapp = 1.0;
    	}
  	}
   
//   	if (t <= 0.014) 
//   	{
//   		double rX = x[0] - (-0.9); 
//   		double rY = x[1] - 1.7; 
//   		double rZ = x[2] - 1.1;
//   		double r = sqrt(rX*rX + rY*rY + rZ*rZ);
//    	 	if ( r < 1.0 )	{
//    	 		Iapp=10000.0*t;
//     		if (Iapp > 100.0) Iapp = 100.0;
//     	}
//   	}
   
   	if ( (t > 225.0) && (t < 226.0) )
   	{
   		double rX = x[0] - (0.0); 
  		double rY = x[1] - 0.5; 
   		double rZ = x[2] - 1.3;
   		double r = sqrt(rX*rX + rY*rY + rZ*rZ);
    		if ( r < 1.0 )  Iapp = 1.0;
   	}
      
    return Iapp;
  }

  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    typedef typename WeakFunctionView::Scalar Scalar;
    static int const components = WeakFunctionView::components;
    typedef Dune::FieldVector<Scalar,components> ValueType;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& u0_,
                   AlievPanfilovEquation<RType,AnsatzVars> const& f_): doScaling(doScaling_), u0(u0_), f(f_) {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::ctype,Cell::dimension> const& localCoordinate) const 
    {
      ValueType u = u0.value(cell,localCoordinate);
      if (doScaling)
        for (int i=0; i<components; ++i) {
        //std::cout << "u, zeta(u): " << u[i] << ", " << f.zeta(u[i]) << std::endl;
          u[i] = f.zeta(u[i]);
        }
      return u;
    }

  private:
    bool doScaling;
    WeakFunctionView const& u0;
    AlievPanfilovEquation<RType,AnsatzVars> const& f;
  };
  
    
};
#endif
