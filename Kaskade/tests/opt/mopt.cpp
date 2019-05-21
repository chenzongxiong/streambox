#include "dune/grid/config.h"
#include "dune/common/fmatrix.hh"
#include "dune/istl/bvector.hh"
#include "linalg/simpleLAPmatrix.hh"	
#include "linalg/umfpack_solve.hh"
#include "linalg/discrete_solver.hh"
#include "algorithm/discrete_bridge.hh"
#include "linalg/tcg.hh"
#include "algorithm/optimization.hh"
#include "testfunctions.hh"
#include "composite_step_linalg.hh"

//typedef UNCMINVector dvec;
//typedef SLAPMatrix<double,1> dmat;

using namespace Kaskade;

template<class VectorImpl>
class IPChartProjection : public AbstractChart
{
public:
  IPChartProjection(int dimFree_,int dimPrimal_, double fractionToTheBoundary_=0.01) : dimFree(dimFree_), dimPrimal(dimPrimal_), fractionToTheBoundary(fractionToTheBoundary_) {}
  
  void addPerturbation(AbstractVector& newIterate, AbstractVector const& perturbation, AbstractLinearization const& lin,
		        std::vector<std::shared_ptr<AbstractVector> > basis = std::vector<std::shared_ptr<AbstractVector> >()) const
  {
    std::vector<double> p,it,ni;
    dynamic_cast<Bridge::Vector<VectorImpl>const& >(perturbation).write(p);   
    dynamic_cast<Bridge::Vector<VectorImpl>const& >(lin.getOrigin()).write(it);   

    ni.resize(it.size());
    
    for(int i=0; i<dimFree; ++i)
    {
      ni[i] = it[i]+p[i];
    }
    for(int i=dimFree; i<dimPrimal;++i)
    {
      ni[i] = std::max(it[i]+p[i],fractionToTheBoundary*it[i]);
    }
    
    for(int i=dimPrimal; i<it.size(); ++i)
      ni[i]=it[i];
    
    dynamic_cast<Bridge::Vector<VectorImpl>& >(newIterate).read(ni);   
  }
  private:
  int dimFree,dimPrimal;
  double fractionToTheBoundary;
};

template<class VectorImpl>
class IPChartDamping : public AbstractChart
{
public:
  IPChartDamping(int dimFree_,int dimPrimal_, double fractionToTheBoundary_=0.01) : dimFree(dimFree_), dimPrimal(dimPrimal_), fractionToTheBoundary(fractionToTheBoundary_) {}
  
  void addPerturbation(AbstractVector& newIterate, AbstractVector const& perturbation, AbstractLinearization const& lin,
		        std::vector<std::shared_ptr<AbstractVector> > basis = std::vector<std::shared_ptr<AbstractVector> >()) const
  {
    std::vector<double> p,it,ni;
    dynamic_cast<Bridge::Vector<VectorImpl>const& >(perturbation).write(p);   
    dynamic_cast<Bridge::Vector<VectorImpl>const& >(lin.getOrigin()).write(it);   

    ni.resize(it.size());
    
    double alpha(1.0);
    
    for(int i=dimFree; i<dimPrimal; ++i)
    {
      if(it[i]+p[i] < fractionToTheBoundary*it[i])
        alpha = std::min((fractionToTheBoundary-1.0)*it[i]/p[i],alpha);
    }
    
    for(int i=0; i<dimPrimal; ++i)
    {
      ni[i] = it[i]+alpha*p[i];
    }
    
    for(int i=dimPrimal; i<it.size(); ++i)
      ni[i]=it[i];
    
    dynamic_cast<Bridge::Vector<VectorImpl>& >(newIterate).read(ni);   
  }
  private:
  int dimFree,dimPrimal;
  double fractionToTheBoundary;
};


template<typename Scalar, typename Vector>
class EuclideanScalarProduct
{
public:

  EuclideanScalarProduct() {}

  double operator()(Vector const& v, Vector const& w)
  {
    Scalar sum(0.0);
    for(int i=0; i<v.size(); ++i)
       if(v.getRole(i)=="primal") sum +=v[i][0]*w[i][0];
    return sum;
  }

  template<class Matrix, class Functional>
  void metricTensor(Matrix& m, Vector const& v, Functional& f) const
  {
    for(int i=0; i<f.dimPrimal(); ++i)
      m[i][i] = 1.0;
  }

template<class Linearization>
  void setOrigin(Linearization const&l)
  {
  }

private:
};


template<typename Scalar, typename PB, class Vector>
class IPScaledScalarProduct
{
public:
  IPScaledScalarProduct(double exponent_=1.0) : exponent(exponent_), idx0(0), idx1(0) {}

  double operator()(Vector const& v, Vector const& w)
  {
    Scalar sum(0.0);
    //std::cout << "---------------  Scalar Product  ------------" << std::endl;
    for(int i=0; i<dimx; ++i)
    {
      sum += s[i]*v[i][0]*w[i][0];
     // std::cout << s[i] << ".." << v[i][0] << "." << w[i][0] << std::endl;
    }
    return sum;
  }

  
  template<class Matrix, class Functional>
  void metricTensor(Matrix& m, Vector const& v, Functional& f) const
  {
    double muscaling=1.0;
    if(exponent != 0.0) muscaling = std::pow(f.Mu(),exponent);
    for(int i=0; i<f.dimPrimal()-f.nIneq(); ++i)
      m[i][i] = 1.0;
    for(int i=f.dimPrimal()-f.nIneq(); i<f.dimPrimal(); ++i)
      m[i][i] = muscaling/v[i]/v[i];
  }

  template<class Linearization>
  void setOrigin(Linearization const&l)
  {
    typename Linearization::Implementation::Implementation const& linI(l.getLinImpl().getLinImpl());
    s.resize(linI.dimPrimal(),0.0);
    dimx=linI.dimPrimal();
    idx0=linI.dimPrimal()-linI.nIneq();
    idx1=linI.dimPrimal();
    nIneq=linI.nIneq();
    nEq=linI.nEq();
    for(int i=0; i<idx0;++i)
      s[i]=1.0;
    
    Vector const& x(l.getOriginImpl());

    mu = linI.Mu();
    
    double muscaling=1.0;
    if(exponent != 0.0) muscaling = std::pow(mu,exponent);
    
    for(int i=idx0;i<idx1;++i)
    {
      s[i]= muscaling/x[i]/x[i];
    }
  }

private:
  double exponent;
  std::vector<Scalar> s;
  int idx0,idx1,nIneq,nEq,dimx;
  double mu;
};


struct pullbackExp
{
  static double logdd(double x)
  {
    return 0.0;
  }

  static double slack(double x, double dx)
  {
    return x*std::exp(dx/x);
  }
  static double slackdd(double x)
  {
    return 1.0/x;
  }
};

struct pullbackPoly
{
  static double logdd(double x)
  {
    return 0.0;
  }

  static double slack(double x, double dx)
  {
    double arg=dx/x;
    return std::max(x*(1+arg+0.5*arg*arg+0.5*arg*arg*arg),1e-16*x);
  }
  static double slackdd(double x)
  {
    return 1.0/x;
  }
};


struct pullbackE2
{
  static double logdd(double x)
  {
    return 1/(x*x);
  }

  static double slack(double x, double dx)
  {
    return x*exp(dx/x-0.5*dx*dx/x/x); 
  }
  static double slackdd(double x)
  {
    return 0.0;
  }
};

struct pullbackTrivial
{
  static double logdd(double x)
  {
    return 1.0/(x*x);
  }

  static double slack(double x, double dx)
  {
    return x+dx; 
  }
  static double slackdd(double x)
  {
    return 0.0;
  }
};

template<class Vector>
class VectorWithRoles : public Vector
{
public:
   
   std::string getRole(int component) const
   {
      if(component < dimX) return "primal";
      else return "dual";
   }
   
   template<class Objective, class Constraints>
   VectorWithRoles(Objective const& obj, Constraints const& con) : Vector(obj.size()+2*con.nIneq()+con.nEq()), dimX(obj.size()+con.nIneq()) {}

   VectorWithRoles(int size, int dimX_) : Vector(size), dimX(dimX_) {}
   
   VectorWithRoles() : Vector() {}
   
   VectorWithRoles<Vector> operator=(double d) { Vector::operator=(d); return *this; }
private:
  int dimX;
   
};


/**
 * 
 * 
 * template <class Vct>
 * class ConceptConstraints 
 * {
 * public:
 *   typedef double Scalar; 
 *   typedef Vct Vector;
 *   typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
 *   
 *   int nIneq() const { return "Number of Inequality Constraints" }
 *   int nEq() const{ return "Number of Equality Constraints"; } 
 *
 *   If you have inequality constraints and equality constraints:
 *   define the inequality constraints first, then the equality constraints
 * 
 *   double d0(Vector const& x, int n) const
 *   {
 *     return "Value of the n-th constraint function at the position x";
 *   }
 *   void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
 *   {
 *     rhs = "1st Derivative of the n-th constraint function at the position x";
 *   }
 * 
 *   void d2(Matrix& m, Vector const& x, int n) const
 *   {
 *     m = "2nd Derivative of the n-th constraint function at the position x";
 *   }
 * };
 * 
 **/

/**
 * template <class Vct>
 * class ConceptObjective
 * {
 * public:
 *   typedef double Scalar; 
 *   typedef Vct Vector;
 *   typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
 * 
 *   bool inDomain(Vector const& x) const 
 *   { 
 *      return "x is in the domain of definition" 
 *   }
 * 
 *   int size() const{ return "Number of variables; } 
 *   double d0(Vector const& x) const
 *   {
 *     return "Function value at the position x"
 *   }
 *   void d1(std::vector<Scalar>& rhs, Vector const& x) const
 *   {
 *     return "1st derivative at the position x"
 *   }
 * 
 *   void d2(Matrix& m, Vector const& x) const
 *   {
 *     return "2nd derivative at the position x"
 *   }
 * };
 **/


/** In this class objective and constraints are put together. Inequality constraints are handled by slack variables and a barrier term:
*
*   The problem min f(x)  s.t. g(x) <= 0  c(x) = 0 is transformed into
* 
*   min f(x) s.t. g(x)+s=0,  c(x)=0  s >= 0  (s : slack variables) and then into
* 
*   min f(x)+mu b(s) s.t. g(x)+s=0  c(x)=0  (b: barrier term, mu: barrier parameter)
* 
*   Template parameters:
* 
*     Objective : class that defines an objective functional, should provide functionality as ConceptObjective (see above)
* 
*     Constraints : class that defines (in)equality constraints, should provide functionality as ConceptConstraints (see above)
* 
*     ScalarProduct: class that defines the scalar product to be used
* 
*     PullBack: ....
**/
template<class Objective, class Constraints, class ScalarProduct, class PullBack>
class ConstrainedProblem
{
public:
/// Scalar type (usually double) 
  typedef typename Objective::Scalar Scalar;
/// Vector type for iterates (for residuals we use std::vector) 
  typedef typename Objective::Vector Vector;
/// Matrix type  (for second derivatives)
  typedef typename Objective::Matrix Matrix;

/// number of inequalities involved -> for each inequality we get a slack variable  and a dual variable
  int nIneq() const {return con.nIneq();}
  
/// number of equalities involved -> for each inequality we get a dual variable
  int nEq() const {return con.nEq();}
/// total number of constraints  
  int nConstraints() const{ return con.nEq()+con.nIneq(); } 

/// total size  
  int size() const{ return dimPrimal()+dimDual(); }  
/// number of primal variables  
  int dimPrimal() const { return obj.size()+con.nIneq(); } 
/// number of free (non-slack) primal variables  
  int dimFree() const { return obj.size(); }
/// number of dual variables  
  int dimDual() const { return nConstraints(); }
 
/// domain of definition 
  bool inDomain(Vector const& x) const 
  { 
    bool inDom = obj.inDomain(x); 
    for(int i=dimFree(); i<dimPrimal(); ++i)
      inDom = inDom && (x[i]>0.0);
    return inDom;
  };

/// barrier parameter  
  double Mu() const { return mu; }

  ConstrainedProblem(const Objective& obj_,const  Constraints& con_,const ScalarProduct& scp_, double mu_, bool tangstep_) 
    : con(con_), obj(obj_), scp(scp_),  x(obj,con), tangstep(tangstep_), mu(mu_)
  {
    mb.setSize(size(),size()); 
    x=0.0;
  }

/// point of evaluation  
  void setOrigin(Vector const& x_) { 
    x = x_; 
  }

/// function value f(x)
  double d0() const
  {
    double sum(obj.d0(x));
    for(int i=dimFree(); i<dimPrimal();++i)
    {
      sum -= mu*std::log(x[i]);
    }
    return sum;
  }

/** residual of the KKT-system:
 *  f'(x)+mu b'(s) + p_I (g'(x)+I) + p_E c'(x)
 *  g(x) + s
 *  c(x)
 * 
 * here  p_I (g'(x)+I) + p_E c'(x) is provided by addpcs
 * 
 */
  void d1(std::vector<Scalar>& rhs) const
  {
      std::vector<Scalar> rhs0(dimFree());
// f'(x)      
      obj.d1(rhs0,x);
      for(int i=0; i<dimFree(); ++i) rhs[i]=rhs0[i];
// + mu b'      (barrier term for inequality constraints
      for(int i=dimFree();i<dimPrimal();++i) rhs[i]=-mu/x[i];
// + p c'(x)      
      for(int i=0; i<nConstraints();++i) addpcs(rhs,i,x[i+dimPrimal()]);
// c(x)      
      for(int i=0;i<nConstraints();++i) rhs[i+dimPrimal()]=con.d0(x,i);
// + s          (slack variables)
      for(int i=0;i<nIneq();++i) rhs[i+dimPrimal()]+=x[i+dimFree()];
  }

/** derivative of KKT system
 * 
 *   f''(x)+mu b''(s)+p_I g''(x)+p_E c''(x)   (g'(x)+I)*     c'(x)*
 *   g'(x)+I                                  0              0  
 *   c'(x)                                    0              0
 */
  void d2(Matrix& m) const
  {
    
// Hessian of Equality constrained problem
// and linearized equality constraints
    for(int i=0; i<m.N(); ++i)
      for(int j=0; j<m.M();++j)
        m[i][j]=0.0;

    for(int i=0; i<mb.N(); ++i)
      for(int j=0; j<mb.M();++j)
        mb[i][j]=0.0;

    if(tangstep)
    {
      obj.d2(mb,x);
      std::vector<Scalar> c(dimFree());

// f''(x)      
       for(int i=0; i<dimFree();++i)
        for(int j=0; j<dimFree();++j)
          m[i][j]=mb[i][j];
// ... + mu b''(s)
      for(int i=dimFree(); i<dimPrimal();++i)
        m[i][i]+= mu*PullBack::logdd(x[i]);
// ...+ p_I g''(x)+p_E c''(x)    
      for(int i=0; i<nConstraints();++i)
        addpcss(m,i,x[i+dimPrimal()]);
    }
    else
    {
      scp.metricTensor(m,x,*this);
    }
  
// Linearization of slack equations   
    std::vector<Scalar> c(dimFree());
    for(int i=0; i<nConstraints(); ++i)
    {
      int ip=i+dimPrimal();
      con.d1(c,x,i);
      for(int j=0; j<dimFree();++j)
      {
        m[ip][j]=c[j];
        m[j][ip]=c[j];
      }
// Slack variables
      if(i<nIneq())
      {
        m[ip][dimFree()+i]= 1.0;
        m[dimFree()+i][ip]= 1.0;
      }
    }
   }

private:  
  
// add n-th component of pc'(x)  
  void addpcs(std::vector<double>& m, int n, double p) const
  {
    std::vector<double> tmp(x.size(),0.0);
    con.d1(tmp,x,n);
    if(n<nIneq())
      tmp[dimFree()+n]+= 1.0;
    for(int i=0; i<m.size();++i)
    m[i]+=tmp[i]*p;
  }

  
// add n-th component of pc''(x)  
  void addpcss(Matrix& m, int n, double p) const
  {
    mb*=0.0;
    con.d2(mb,x,n);
    if(n<nIneq())
      mb[dimFree()+n][dimFree()+n]+= PullBack::slackdd(x[dimFree()+n]);
    mb*=p;
    m+=mb;
  }

  const Constraints& con;
  const Objective& obj;
  const ScalarProduct& scp;
  mutable Matrix mb;
  mutable Vector x;
  bool tangstep;
public:
  double mu;
};


template<class Functional>
void checkDerivatives(Functional& f,typename Functional::Vector& x) 
{
  f.setOrigin(x);
  double f0=f.d0();
  double h=1e-8;
  std::vector<double> df(f.size()),Df(f.size());
  typename Functional::Matrix ddf(f.size(),f.size()),DDf(f.size(),f.size());
  std::cout << "Function:" << std::endl;
  std::cout << f0 << std::endl;
  std::cout << "Gradient:" << std::endl;
  f.d1(Df);
  for(int i=0; i<f.size();++i)
  {
    std::cout << Df[i] << " ";
  } 
  std::cout << std::endl;
  std::cout << "Hessian:" << std::endl;
  f.d2(DDf);
  for(int i=0; i<f.size();++i)
  {
    for(int j=0; j<f.size();++j)
      std::cout << DDf[j][i] << " ";
    std::cout << std::endl;
  } 
  std::cout << "Checking Gradient:" << std::endl;
  for(int i=0; i<x.size();++i)
  {
    x[i]+=h;
    f.setOrigin(x);
    df[i]=(f.d0()-f0)/h;
    x[i]-=h;
    std::cout << df[i]-Df[i] << " ";
  } 
  std::cout << std::endl;
  std::cout << "Checking Hessian:" << std::endl;
  for(int i=0; i<x.size();++i)
  {  
    x[i]+=h;
    f.setOrigin(x);
    f.d1(df);
  for(int j=0; j<f.size();++j)
  {
    ddf[j][i]=(df[j]-Df[j])/h;
    std::cout << ddf[j][i]-DDf[j][i] << " ";
  }
  std::cout << std::endl;
  x[i]-=h;
  } 
}

int main()
{	
  typedef pullbackTrivial pullback; 
  
  typedef VectorWithRoles<Dune::BlockVector<Dune::FieldVector<double,1> > > DomainVector;

  std::vector<double> coeffs(1,1.0);
  
  typedef LinearFunction<DomainVector> Functional;
  Functional obj(coeffs);
  
  typedef WBConstraints<DomainVector> Constraints;
  Constraints con(1,1);

  double mu0 = 1.0;

  typedef IPScaledScalarProduct<double,pullback,DomainVector> MyScalarProduct;
  MyScalarProduct es;


  typedef ConstrainedProblem<Functional,Constraints,MyScalarProduct,pullback> MyConstrainedProblem;
  Bridge::Functional<MyConstrainedProblem,DomainVector,DomainVector> normalStep(obj,con,es,mu0,false);
  Bridge::Functional<MyConstrainedProblem,DomainVector,DomainVector> tangentialStep(obj,con,es,mu0,true);
  
  DomainVector x_0(obj,con);

  x_0[0]=-1.0;
  x_0[1]=2.0;
  x_0[2]=2.0;
  x_0[3]=2.0;
  x_0[4]=2.0;

  Bridge::Vector<DomainVector> x(x_0);
  
  x.print();

  Bridge::LocalScalarProduct<MyScalarProduct,DomainVector,MyConstrainedProblem> normL(Bridge::makeUP(new MyScalarProduct));
  Bridge::LocalScalarProduct<MyScalarProduct,DomainVector,MyConstrainedProblem> normC(Bridge::makeUP(new MyScalarProduct));
  typedef UMFFactorization<double> MySolver;

  typedef Bridge::DirectInnerSolver<MySolver,DomainVector> NormalSolver;
  typedef Bridge::SimpleTangentialSolver<MySolver,DomainVector> TangentialSolver;

  NormalSolver normalSolver;
  TangentialSolver tangentialSolver;

  OptimizationParameters p(1e-3,50);
  p.sigma=0.4;
  p.Thetad=0.9;
  p.gamma=1.0;
  p.alpha=1.0;
  p.fracBdr=0.1;
  
  IPChartDamping<DomainVector> chart(obj.size(),obj.size()+con.nIneq(),p.fracBdr);

  int verbosity = 2;
  
  Optimization compositeStepSolver(normalSolver,tangentialSolver, normL, normC, chart, p, verbosity);

  compositeStepSolver.reportOnIteration(2);
  {
    p.desiredAccuracy=1e-12;//*std::pow(addrpf->Mu(),0.5);
    compositeStepSolver.solve(&normalStep, &tangentialStep, x);
  }
  return 0;
}
