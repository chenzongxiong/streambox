#ifndef TESTFUNCTIONS_HH
#define TESTFUNCTIONS_HH

// Constraints


template <class Vct>
class CircleConstraints 
{
public:
  CircleConstraints(double r_=1.0,double x0_=0.0,double x1_=0.0) : r(r_),x0(x0_),x1(x1_) {};
 
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  int nEq() const{ return 1; } 
  int nIneq() const { return 0; }

  double d0(Vector const& x, int n) const
  {
    return std::pow(x[0]-x0,2)+std::pow(x[1]-x1,2)-r*r;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    rhs[0]=2*(x[0]-x0);
    rhs[1]=2*(x[1]-x1);
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    m[0][0]=2;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=2;
  }
private:
  double r,x0,x1;
};

template <class Vct>
class LinearConstraints
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  

  int nEq() const{ return 1; } 
  int nIneq() const { return 0; }
  double d0(Vector const& x, int n) const
  {
    return x[1];
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    rhs[0]=0;
    rhs[1]=1;
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    m[0][0]=0;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=0;
  }
};

template <class Vct>
class CubicConstraints
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  int nEq() const{ return 1; } 
  int nIneq() const { return 0; }
  double d0(Vector const& x, int n) const
  {
    return x[0]+x[1]+std::pow(x[1],3)/6.0-1.0;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    rhs[0]=1;
    rhs[1]=1+std::pow(x[1],2)/2.0;
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    m[0][0]=0;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=x[1];
  }
};


template <class Vct>
class NoConstraints 
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  int nEq() const{ return 0; } 
  int nIneq() const { return 0; }
  
  double d0(Vector const& x, int n) const
  {
    return 0;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    rhs[0]=0;
    rhs[1]=0;
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    m[0][0]=0;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=0;
  }
};

//--------------------------------------------------------------------------------------------------------
// Functionals


template <class Vct>
class Himmelblau
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  int dimx() const { return 2;} 
  bool inDomain(Vector const& x) const { return true; };

  int nEq() const{ return 0; } 
  int size() const{ return dimx(); } 
  double d0(Vector const& x) const
  {
    return std::pow(x[0]*x[0]+x[1]-11,2)+std::pow(x[0]+x[1]*x[1]-7,2);
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=(x[0]*x[0]+x[1]-11)*4*x[0]+2*(x[0]+x[1]*x[1]-7);
    rhs[1]=2*(x[0]*x[0]+x[1]-11)+(x[0]+x[1]*x[1]-7)*4*x[1];
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=12*x[0]*x[0]+4*x[1]+2-11*4;
    m[0][1]=4*x[0]+4*x[1];
    m[1][0]=4*x[0]+4*x[1];
    m[1][1]=12*x[1]*x[1]+4*x[0]+2 - 7*4;
  }
};

template <class Vct>
class Rosenbrock
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  bool inDomain(Vector const& x) const { return true; };

  int size() const{ return 2; } 
  double d0(Vector const& x) const
  {
    return 100*std::pow(x[1] - x[0]*x[0],2) + std::pow(1 - x[0],2);
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=400*std::pow(x[0],3) - 400*x[0]*x[1]-2 + 2*x[0];
    rhs[1]=200*x[1] - 200*x[0]*x[0];
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=1200*x[0]*x[0] - 400*x[1] + 2;
    m[0][1]=-400*x[0];
    m[1][0]=-400*x[0];
    m[1][1]=200;
  }
};

template <class Vct>
class QuadraticOptFunction
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  QuadraticOptFunction(double alpha_,double beta_) : alpha(alpha_), beta(beta_) {}

  int size() const{ return 2; } 
  double d0(Vector const& x) const
  {
    return alpha*x[0]*x[0]+beta*x[1]*x[1]+0.1*x[1]*x[1]*x[1];
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=2*alpha*x[0];
    rhs[1]=2*beta*x[1]+3*0.1*x[1]*x[1];
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=2*alpha;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=2*beta+6*0.1*x[1];
  }
private:
  double alpha,beta;
};

template <class Vct>
class MarathosFunction
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;

  bool inDomain(Vector const& x) const { return true; };
  int size() const{ return 2; } 
  double d0(Vector const& x) const
  {
    return 2*((x[0])*(x[0])+x[1]*x[1]-1)-x[0];
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=4*(x[0])-1;
    rhs[1]=4*x[1];
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=4;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=4;
  }
private:
  double alpha,beta;
};

template <class Vct>
class BCCSSFunction
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  bool inDomain(Vector const& x) const { return true; };

  int size() const{ return 2; } 
  double d0(Vector const& x) const
  {
    return 0.5*x[0]*x[0]+2*x[1];
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=x[0];
    rhs[1]=2;
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=1;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=0;
  }
private:
  double alpha,beta;
};

template<class Functional, class ScalarProduct>
class MetricFunctional
{
public:
  typedef typename Functional::Scalar Scalar;
  
  typedef typename Functional::Vector Vector;
  typedef typename Functional::Matrix Matrix;
  bool inDomain(Vector const& x) const { return f.inDomain(x); };

  int size() const{ return f.size(); } 

  MetricFunctional(Functional &f_, ScalarProduct &s_) : f(f_), s(s_) {}

  double d0(Vector const& x) const
  {
    return f.d0(x);
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    f.d1(rhs,x);
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m *= 0.0;
    s.metricTensor(m,x,f);
  }
private:
  Functional& f;
  ScalarProduct& s;
};

template <class Vct>
class LinearFunction
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  LinearFunction(std::vector<double> const& coeffs_) : coeffs(coeffs_) {}

  bool inDomain(Vector const& x) const { return true; };

  int size() const { return coeffs.size(); } 

  double d0(Vector const& x) const
  {
    double sum(0.0);
    for(int i=0; i<coeffs.size();++i)
      sum+=coeffs[i]*x[i];
    return sum;
  }

  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    for(int i=0; i<coeffs.size();++i)
    {
      rhs[i]=coeffs[i];
    }
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m*=0.0;
  }
private:
  std::vector<double> const& coeffs;
};


template <class Vct>
class ConstFunction
{
public:
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;

  bool inDomain(Vector const& x) const { return true; };

  int size() const { return 1;} 
  double d0(Vector const& x) const
  {
    return 0;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x) const
  {
    rhs[0]=0;
    rhs[1]=0;
  }

  void d2(Matrix& m, Vector const& x) const
  {
    m[0][0]=0;
    m[0][1]=0;
    m[1][0]=0;
    m[1][1]=0;
  }
};

template <class Vct>
class WBConstraints
{
public:
  double s,t;

  WBConstraints(double a_, double b_) : s(1.0),t(1.0),a(a_), b(b_) {};
 
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  int nEq() const { return 0; } 
  int nIneq() const { return 2; }

  double d0(Vector const& x, int n) const
  {
      if(n==0)
        return s*(-a-x[0]*x[0]);
      if(n==1)
        return t*(b-x[0]);
      return 0;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    if(n==0)
      rhs[0]=-s*2*x[0];
    if(n==1)
      rhs[0]=-t*1;
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    if(n==0)
      m[0][0]=-2*s;
    if(n==1)
      m[0][0]=0;
  }
private:
  double a,b;
};


template <class Vct>
class BoxConstraints
{
public:
  double s,t;

  BoxConstraints(double a_, double b_) : s(1.0),t(1.0),a(a_), b(b_) {};
 
  typedef double Scalar; 
  typedef Vct Vector;
  typedef Dune::Matrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
  
  int nEq() const{ return 0; } 
  int nIneq() const { return 2; }

  double d0(Vector const& x, int n) const
  {
      if(n==0)
        return s*(-a+x[0]);
      if(n==1)
        return t*(b-x[0]);
      return 0;
  }
  void d1(std::vector<Scalar>& rhs, Vector const& x, int n) const
  {
    if(n==0)
      rhs[0]= s*1;
    if(n==1)
      rhs[0]=-t*1;
  }

  void d2(Matrix& m, Vector const& x, int n) const
  {
    if(n==0)
      m[0][0]=0;
    if(n==1)
      m[0][0]=0;
  }
private:
  double a,b;
};


#endif
