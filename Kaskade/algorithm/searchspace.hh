#ifndef SEARCHSPACE_HH
#define SEARCHSPACE_HH

#include <memory>
#include "abstract_interface.hh"
#include "opt_interface.hh"
#include "algorithm_base.hh"
#include "istl/matrix.hh"
#include "linalg/simpleLAPmatrix.hh"
#include <vector>
#include "dune/common/fmatrix.hh"


namespace Kaskade
{

/**
 * @file 
 * @brief Some model functions for various purposes
 * @author Anton Schiela
 *
 */

class SearchSpaceCreator
{
public:
  SearchSpaceCreator(AbstractTangentialSolver* tangentialSolver_,AbstractNormalSolver* normalSolver_): 
    tangentialSolver(tangentialSolver_),normalSolver(normalSolver_), normalactive(false) {}; 


  virtual ~SearchSpaceCreator() {};
  int getMaxDimension() const
  {
    int dim=0;
    if(tangentialSolver) dim+=tangentialSolver->nSolutionVectors();
    if(normalSolver) dim+=1;
    return dim;
  };

  AbstractNormalSolver* getNormalSolver() { return normalSolver; }
  AbstractTangentialSolver* getTangentialSolver() { return tangentialSolver; }

  int getDimension() const
  {
    return currentdimension;
  };

  bool normalFactorizationPresent() {return normalactive;}

  bool hasEqualityConstraints() { return (normalSolver!=0); }
  bool hasNontrivialFunctional() { return (tangentialSolver!=0); }

  virtual void computeBasisVectors(std::vector<AbstractFunctionSpaceElement *>& basisVectors, 
                                   AbstractFunctionSpaceElement& iterate, 
                                   AbstractLinearization& normalLinearization,
                                   AbstractFunctional& tangentialFunctional,
                                   std::unique_ptr<AbstractLinearization>& tangentialLinearization,
                                   AbstractNorm const& norm,
                                   double ThetaAim,
                                   double omegaC,
                                   double omegaL,
                                   double omegaH,
                                   int report,
                                   double& nu0,
                                   double& normNormal);
  void getLinearCombination(std::vector<AbstractFunctionSpaceElement *>& basisVectors, std::vector<double>const & coefficients, AbstractFunctionSpaceElement& result) const;

private:
  
  AbstractTangentialSolver* tangentialSolver; 
  AbstractNormalSolver* normalSolver; 
  bool normalactive;
  int currentdimension;
};



/// For optimization with cubic upper bounds method
class CUBThetaModelFunction
{
public:

  typedef UNCMINVector dvec;
  typedef SLAPMatrix<double,1> dmat;

  virtual ~CUBThetaModelFunction(){};

  CUBThetaModelFunction(AbstractLinearization const& linT, AbstractLinearization const& linN, std::vector<AbstractFunctionSpaceElement* >& basis, AbstractScalarProduct& scL, AbstractScalarProduct& scC, SearchSpaceCreator& sp);

  CUBThetaModelFunction(AbstractLinearization const& linT, 
                        AbstractLinearization const& linN,  
                        std::vector<AbstractFunctionSpaceElement *>& basis, 
                        AbstractScalarProduct& scL, 
                        AbstractScalarProduct& scC,
                        AbstractTangentialSolver const* tS,          
                        SearchSpaceCreator& sp, 
                        bool useHessianAsNorm);


template<class Paras>
void reparametrize(double omegaL_,double omegaC_, double gamma_, Paras p)
  {
    omegaL=omegaL_;
    omegaC=omegaC_;
    gamma=gamma_;
    aimContraction = p.ThetaAim/p.safetyFactor;
    normalContraction = p.ThetaNormal/p.safetyFactor;
    alpha = p.alpha;
  }

  // minimize functional
  void getMinimalCoefficients(std::vector<double>& coeff);
  
// Interface to UNCMIN
virtual double evalModel(std::vector<double> & iterate)
  {
    double d= eval(iterate,0.0,0.0,0.0);
//    d -= 0.5*(computeblf(normblf,iterate,iterate)-computeblf(norm2blf,iterate,iterate));
    return d;
  }

  virtual double evalRegModel(std::vector<double> & iterate, double omegaL)
  {
    double d= eval(iterate,omegaL,0.0,0.0);
//    d -= 0.5*(computeblf(normblf,iterate,iterate)-computeblf(norm2blf,iterate,iterate));
    return d;
  }


virtual double eval(std::vector<double> & iterate, double oL, double oC, double g, double soscale=1.0,double foscale=1.0,double alpha=1.0)
  {
    if(fixednu && iterate.size()<sz) iterate.push_back(lopt);
    double f(0.0);
    if(tangentialstep)
    {
      // f'dx
      f += foscale*dot(gradientcoeff, iterate);
      // f''dx^2 blf = bilinear form
      f += soscale*0.5*computeblf(hessianblf,iterate,iterate);
      double np3=std::pow(std::fabs(computeblf(normMblf,iterate,iterate)),1.5);
      f += oL/6.0*np3;
    }
    f *= alpha;

    // barrier term instead of trust region
    if(normalstep && g != 0.0) 
      f -= g*std::log(compute_barrier_arg(iterate,oC));
    return f;
  }
	
  // Function to minimize
  double f_to_minimize(dvec &p)
  {
    std::vector<double> iterate(p.data);
    double f=eval(iterate,omegaL,omegaC,gamma,1.0,1.0,alpha);
    return f;
  }
	
  // Gradient of function to minimize
virtual void gradient(dvec &p, dvec &g) 
  {
    std::vector<double> iterate(p.data);
    if(fixednu && iterate.size()<sz) iterate.push_back(1.0);
    std::vector<double> v(sz),w(sz),w2(sz);

    MatMultiply(hessianblf,iterate,v);

    MatMultiply(normMblf,iterate,w);
    MatMultiply(normCblf,iterate,w2);

    double f=std::pow(std::fabs(computeblf(normMblf,iterate,iterate)),0.5);
    double normCiter=std::pow(std::fabs(computeblf(normCblf,iterate,iterate)),0.5);

    if(tangentialstep) 
      for(int i=0; i<g.data.size();++i)
        g.data[i]=alpha*(gradientcoeff[i]+v[i]+0.5*omegaL*f*w[i]);
    else
      for(int i=0; i<g.data.size();++i)
        g.data[i]=0.0;

    if(normalstep)
    {
      double b=-gamma/compute_barrier_arg(iterate,omegaC);
      for(int i=0; i<g.data.size();++i)
        g.data[i]-=b*omegaC/2.0*w2[i]/normCiter;
    }
  }
	
  // Hessian not used in this example
  void hessian(dvec /* &x */, dmat /* &h */) {}
	
  // Indicates analytic gradient is used
virtual int HasAnalyticGradient() {return 1;}
	
  // Indicates analytic hessian is not used
  int HasAnalyticHessian() {return 0;}
	
  // Any real vector will contain valid parameter values
virtual int ValidParameters(dvec &x) {
    std::vector<double> iterate(x.data);
    if(fixednu && iterate.size()<sz) iterate.push_back(lopt);
    return (!normalstep) || compute_barrier_arg(iterate,omegaC)>0;
  }
	
  // Dimension of problem 
  int dim() {return sz-(fixednu==true && normalstep);}


  double normC(std::vector<double>& coeff)
  {
    return std::sqrt(std::fabs(computeblf(normCblf,coeff,coeff)));
  }

  double normL(std::vector<double>& coeff)
  {
    return std::sqrt(std::fabs(computeblf(normLblf,coeff,coeff)));
  }

  double normH(std::vector<double>& coeff)
  {
    return std::sqrt(std::fabs(computeblf(normHblf,coeff,coeff)));
  }

protected:

  void MatMultiply(SLAPMatrix<double>& mat, std::vector<double>& in, std::vector<double>& out)
  {
    if(in.size() < mat.cols() || out.size() < mat.rows())
    {
      std::vector<double> in2(mat.cols(),0.0), out2(mat.rows(),0.0);
      for(int i=0; i<in.size();++i)
        in2[i]=in[i];
      MatrixMultiply(mat,in2,out2);
      for(int i=0; i<out.size();++i)
        out[i]=out2[i];
    }
    else
      MatrixMultiply(mat,in,out);
  }

  double computeblf(SLAPMatrix<double>& blf, std::vector<double>& left, std::vector<double>& right)
  {
    std::vector<double> Ar(left.size());
    MatMultiply(blf,right, Ar);
    return dot(Ar,left);
  }

  // Function to minimize
  double compute_barrier_arg(std::vector<double> &iterate, double oC)
  {
    double normIter=std::pow(std::fabs(computeblf(normCblf,iterate,iterate)),0.5);
    double barg(aimContraction);
    barg -= oC/2 * normIter;
    return barg;
  }

  double dot(std::vector<double>const& v1,std::vector<double>const& v2)
  {
    double sum(0.0);
    for(int i=0; i<v1.size(); ++i)
      sum+=v1[i]*v2[i];
    return sum;
  } 
  void debugInfo();

  void plausibility();

  int sz, nidx,tidx0,tidx1;
  SLAPMatrix<double> hessianblf, normHblf, ablf, normLblf, normCblf, normMblf;
  std::vector<double> gradientcoeff;
  double omegaL, omegaC;
  double gamma;
  double aimContraction, normalContraction, alpha; 
  double scalef,lopt;
  bool fixednu;
  bool const normalstep;
  bool const tangentialstep;
};

}  // namespace Kaskade
#endif


