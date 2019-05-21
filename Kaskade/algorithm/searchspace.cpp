#include "./opt_aux/src/include/Uncmin.h"
#include <dune/grid/config.h>
#include "linalg/simpleLAPmatrix.hh"	
#include "searchspace.hh"    
#include <cstdio>	// for printf
//#include "composite_step.hh"


namespace Kaskade
{

typedef UNCMINVector dvec;
typedef SLAPMatrix<double,1> dmat;


void SearchSpaceCreator::computeBasisVectors(std::vector<AbstractFunctionSpaceElement *>& basisVectors, 
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
                                   double& normNormal)
  {    
    assert(basisVectors.size()>=getMaxDimension());
    currentdimension=0;
    if(normalSolver)
    {
      normalactive = true;

// Important: here the iterate is modified, namely the adjoint state, as p influences the second derivative of the Lagrangian

      // basis of overall problems: basisVectors

      normalSolver->solve(*(basisVectors[currentdimension]),iterate,normalLinearization);
      if(report>=2) { basisVectors[currentdimension]->print("   dn:");}
      normNormal=norm(*basisVectors[currentdimension]);
      // nu0: predictor normal damping
      nu0 = std::min(1.0,2.0*ThetaAim/(omegaC*normNormal));

      std::cout << nu0 << std::endl;

      currentdimension+=1;
    }
    else
    {
      std::cout << "Unconstrained Problem" << std::endl;
      nu0=1.0;
      normNormal = 0.0;
    }

    if(tangentialSolver){
      if(!tangentialLinearization.get())
      {
        tangentialLinearization=tangentialFunctional.getLinearization(iterate);
      } 


      currentdimension+=tangentialSolver->solve(basisVectors,*tangentialLinearization,normalLinearization,currentdimension/*,ThetaAim,omegaC,omegaL,omegaH*/,nu0);
        
      if(report >= 2)
      for(int i=normalactive; i<currentdimension; ++i)
      {
        basisVectors[i]->print("dt:");
      }
    }

    // basisVectors dn, dt1, dt2 -> dt1, dt2, ..., dn
    if(normalSolver)
    {
      std::vector<AbstractFunctionSpaceElement *> bV(basisVectors);
      for(int i=0; i<currentdimension-1;++i) basisVectors[i]=bV[i+1];
      basisVectors[currentdimension-1]=bV[0]; // normal step to the end

    }
  }

void SearchSpaceCreator::getLinearCombination(std::vector<AbstractFunctionSpaceElement *>& basisVectors, std::vector<double>const & coefficients, AbstractFunctionSpaceElement& result) const
  {
    assert(getDimension() > 0);
    assert(basisVectors.size() >= getDimension());
    assert(coefficients.size() >= getDimension());
    result *= 0.0;
    for(int i=0; i<getDimension();++i)
    {
      for(int j=0; j<result.nComponents();++j)
        // do not change dual variables
        if(result.getRole(j)!="dual")
        {
          result.axpy(coefficients[i],*(basisVectors[i]),j);      
        }
    }
  }

// in+f''*x
void ddxpyprimal(AbstractLinearization const& lin, AbstractFunctionSpaceElement& result, AbstractFunctionSpaceElement const& in)
{
  for(int i=0; i<result.nComponents(); ++i)
    for(int j=0; j<in.nComponents(); ++j)
    {
      // do not change
      if(result.getRole(i) != "dual" && in.getRole(j) != "dual")
        lin.ddxpy(result,in,i,i+1,j,j+1);
    }
}

CUBThetaModelFunction::CUBThetaModelFunction(AbstractLinearization const& linT, AbstractLinearization const& linN, std::vector<AbstractFunctionSpaceElement *>& basis, 
                                   AbstractScalarProduct& scL, AbstractScalarProduct& scC, 
                                   SearchSpaceCreator& sp):
  sz(sp.getDimension()),
  nidx(sz-sp.hasEqualityConstraints()),
  tidx0(0),
  tidx1(nidx),
  hessianblf(sz, sz),
  normHblf(sz,sz),
  ablf(sz,sz),
  normLblf(sz, sz),
  normCblf(sz, sz),
  normMblf(sz, sz),
  gradientcoeff(sz),
  fixednu(false),
  normalstep(sp.hasEqualityConstraints()),
  tangentialstep(sp.hasNontrivialFunctional())
{
  std::unique_ptr<AbstractFunctionSpaceElement> df((basis[0])->clone());
  linN.evald(*df);

  for(int i=0; i<sz;++i)
  {
    gradientcoeff[i]=0.0;
    for(int j=0; j<basis[i]->nComponents();++j)
      if(basis[i]->getRole(j) != "dual")
        gradientcoeff[i]+=df->applyAsDualTo(*(basis[i]),j);
  }
  for(int i=0; i<sz;++i)
  {
    *df *= 0.0;
//    linT.ddxpy(*df,*(basis[i]),0,dimx,0,dimx);

    ddxpyprimal(linT,*df,*(basis[i]));
    for(int j=0; j<sz;++j)
    {
      hessianblf(i,j)=df->applyAsDualTo(*(basis[j]));
      normLblf(i,j)=scL(*(basis[i]),*(basis[j]));                                  
      normCblf(i,j)=scC(*(basis[i]),*(basis[j]));                                  
    }
  }
}


// subspace minimization, create cubic model wrt. basis
CUBThetaModelFunction::CUBThetaModelFunction(AbstractLinearization const& linT, AbstractLinearization const& linN,
                                             std::vector<AbstractFunctionSpaceElement *>& basis, 
                                             AbstractScalarProduct& scL, AbstractScalarProduct& scC, 
                                             AbstractTangentialSolver const* tS,          
                                             SearchSpaceCreator& sp, 
                                             bool useHessianAsNorm)
  :
  sz(sp.getDimension()),
  nidx(sz-sp.hasEqualityConstraints()),
  tidx0(0),
  tidx1(nidx),
  hessianblf(sz, sz),
  normHblf(sz,sz),
  ablf(sz,sz),
  normLblf(sz, sz),
  normCblf(sz, sz),
  normMblf(sz, sz),
  gradientcoeff(sz),
  fixednu(false),
  normalstep(sp.hasEqualityConstraints()),
  tangentialstep(sp.hasNontrivialFunctional())
{
  std::unique_ptr<AbstractFunctionSpaceElement> df((basis[0])->clone());
  // get derivative of normal functional wrt. to FE basis
  linN.evald(*df);

  // get derivative wrt to basis vectors
  for(int i=0; i<sz;++i)
  {
    gradientcoeff[i]=0.0;
    for(int j=0; j<basis[i]->nComponents();++j)
      if(basis[i]->getRole(j) != "dual")
        gradientcoeff[i]+=df->applyAsDualTo(*(basis[i]),j);
  }


  // TCG did get  norm of tangential step don't calculate it again
  std::vector<double> m,a;
  if(tS->getNormInfo(m,a)==true)
  {

    normLblf(0,0)=m[0];  
    ablf(0,0)=a[0];
    if(tidx1>0)
    {
      normLblf(1,1)=m[1];
      normLblf(1,0)=m[2];
      normLblf(0,1)=m[2];
      ablf(1,1)=a[1];
      ablf(0,1)=a[2];
      ablf(1,0)=a[2];
    }
    for(int i=0; i<sz;++i)
    {
      *df *= 0.0;
      //      linT.ddxpy(*df,*(basis[i]),0,dimx,0,dimx);
      ddxpyprimal(linT,*df,*(basis[i]));

      for(int j=0; j<sz;++j)
      {

        normCblf(i,j)=normLblf(i,j);
        hessianblf(i,j)=df->applyAsDualTo(*(basis[j]));
        ablf(i,j)=hessianblf(i,j);
      }
    }
  }

  // standard case
  else
  {
    for(int i=0; i<sz;++i)
    {
      *df *= 0.0;
//    linT.ddxpy(*df,*(basis[i]),0,dimx,0,dimx);
      ddxpyprimal(linT,*df,*(basis[i]));
      for(int j=0; j<sz;++j)
      {
        // hessian bzgl. basis
        hessianblf(i,j)=df->applyAsDualTo(*(basis[j]));
        // get norm
        normLblf(i,j)=scL(*(basis[i]),*(basis[j]));
        normCblf(i,j)=scC(*(basis[i]),*(basis[j]));                                  
      }
    }
  }

  // not interesting
  normHblf(0,0)=fabs(hessianblf(0,0));
  normHblf(1,1)=fabs(hessianblf(1,1));
  normHblf(1,0)=0.0;
  normHblf(0,1)=0.0;

  if(useHessianAsNorm)
  {
    for(int i=0; i<sz;++i)
    {
      for(int j=0; j<sz;++j)
      {
        normMblf(i,j)=  normHblf(i,j);
      }
    }
  } 
  // used norm
  else
  {
    for(int i=0; i<sz;++i)
    {
      for(int j=0; j<sz;++j)
      {
        normMblf(i,j)=  normLblf(i,j);
      }
    }
  }
}

void CUBThetaModelFunction::plausibility()
{
//   for(int i=0; i<sz; ++i)
//   {
//     for(int j=0; j<sz; ++j)
//       if(hessianblf(i,j)-ablf(i,j) > 0.5*(std::fabs(hessianblf(i,i))+std::fabs(ablf(i,i))))
//         std::cout << "Plausibility:" << (hessianblf(i,j)-ablf(i,j)) << ", ";
//     std::cout << std::endl;
//   }
}



void CUBThetaModelFunction::debugInfo()
{
  std::cout << "------------------Debug Info:------------------" << std::endl;
  std::cout << "Gradients:" << std::endl;
  std::cout << "[ ";
  for(int i=0; i<sz; ++i)
    std::cout << gradientcoeff[i] << ", ";
  std::cout << " ]" << std::endl;
  std::cout << "Hessian:" << std::endl;
  std::cout << "[ ";
  for(int i=0; i<sz; ++i)
  {
    for(int j=0; j<sz; ++j)
      std::cout << hessianblf(i,j) << ", ";
    std::cout << std::endl;
  }
  std::cout << " ]" << std::endl;

  std::cout << "A:" << std::endl;
  std::cout << "[ ";
  for(int i=0; i<sz; ++i)
  {
    for(int j=0; j<sz; ++j)
      std::cout << ablf(i,j) << ", ";
    std::cout << std::endl;
  }
  std::cout << " ]" << std::endl;
  std::cout << "NormC:" << std::endl;
  std::cout << "[ ";
  for(int i=0; i<sz; ++i)
  {
    for(int j=0; j<sz; ++j)
      std::cout << normCblf(i,j) << ", ";
    std::cout << std::endl;
  }
  std::cout << " ]" << std::endl;
  std::cout << "NormL:" << std::endl;
  std::cout << "[ ";
  for(int i=0; i<sz; ++i)
  {
    for(int j=0; j<sz; ++j)
      std::cout << normLblf(i,j) << ", ";
    std::cout << std::endl;
  }
  std::cout << " ]" << std::endl;
    std::vector<double> eigs(sz);
    SymmetricEigenvalues(normCblf,eigs);
    std::cout << "NormCEigenvalues:";
    for(int l=0; l< sz; ++l)
      std::cout << eigs[l] << " ";
    std::cout << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
}

void CUBThetaModelFunction::getMinimalCoefficients(std::vector<double>& coeff)
{	
  gamma = 1e-20;
  fixednu=true;
  scalef=0.0;
  double typical_size[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  for(int i=0; i<sz;++i)
  {
    scalef += 0.5*std::fabs(gradientcoeff[i]);
    for(int j=0; j<sz;++j)
    {
      scalef+=std::sqrt(std::abs(hessianblf(i,j))+std::abs(normMblf(i,j)));
    }
  }
  double fpls;
  coeff.resize(0);
  coeff.resize(3,0.0);
  coeff.resize(sz,0.0);
  Uncmin<dvec, dmat, CUBThetaModelFunction> min(this); // create Uncmin object
  dvec xpls(tidx1-tidx0),xpls2(tidx1-tidx0);
    
  min.SetScaleFunc(scalef);
  min.SetNumDigits(80);
  min.SetMaxIter(1000);
  min.SetChecks(0);
	
  double normn=std::sqrt(std::fabs(normCblf(nidx,nidx)));
  double h = omegaC*normn;
  if(normn+1.0==1.0) 
    lopt =1.0; // nu = lopt
  else
    lopt = std::min(1.0,2.0*normalContraction/h);
//  std::cout << "Scaling:" << scalef << std::endl;

  min.SetFunction(this);
  // starting values

  dvec gpls(tidx1-tidx0);
  dvec typ(typical_size, typical_size + tidx1-tidx0);
  min.SetScaleArg(typ);

  dvec start(tidx1-tidx0);

  for(int i=tidx0; i<tidx1;++i)
    start.data[i-tidx0]=1.0;

  int k(0),msg;

// Try optimistic initial values

  do {
    ++k;
    min.Minimize(start, xpls2,  fpls,  gpls);
    msg = min.GetMessage();
     std::cout<< "Message:" << msg << std::endl;
    if(msg != 0 && msg != 1 && msg != 2  && msg != 3)
    {
      if(msg > -20) 
      {
        gamma *= 5;
        gamma = std::min(gamma,1.0);
      }
      else 
      {
        for(int i=tidx0;i<tidx1;++i)
        {
          start.data[i-tidx0] *= 0.25;
        }
      }
    }
      
  } while(msg != 3 && msg != 2 && msg != 1 && msg != 0 && k<20);


  if(k==20)
  {
    int sz0=sz;
    sz = 1+(fixednu==true && normalstep);
    Uncmin<dvec, dmat, CUBThetaModelFunction> min2(this); // create Uncmin object
    dvec typ2(typical_size, typical_size + 1);
    min2.SetScaleArg(typ2);
    min2.SetScaleFunc(scalef);
    min2.SetNumDigits(80);
    min2.SetMaxIter(1000);
    min2.SetChecks(1);
    min2.SetFunction(this);
    std::cout << "Subspace minimization failed. Trying one-dimensional minimization" << std::endl;
    dvec start2(1), xpls3(1),gpls2(1);
    start2.data[0]=1.0;
    do {
      ++k;
      min2.Minimize(start2, xpls3,  fpls,  gpls2);
      msg = min2.GetMessage();
     std::cout<< "Message:" << msg << std::endl;
     if(msg != 0 && msg != 1 && msg != 2)
      {
        if(msg > -20) 
        {
          gamma *= 5;
          gamma = std::min(gamma,1.0);
        }
        else 
          start2.data[0] *= 0.25;
      }
      
    } while(msg != 2 && msg != 1 && msg != 0 && k<40);
    for(int i=tidx0; i<tidx1; ++i)
      xpls2.data[i-tidx0]=0.0;
    xpls2.data[0]=xpls3.data[0];
    sz=sz0;
  }

  if(k>21) {
    std::cout << "Warning: Minimizer terminates after " << k << " attempts." << std::endl;
  }

    for(int i=tidx0; i<tidx1; ++i)
    {
//      coeff[i]=std::min(std::max(xpls2.data[i-tidx0],-10.0),10.0);
      coeff[i]=xpls2.data[i-tidx0];
      if(std::fabs(coeff[i]-1.0)<0.005 && tidx1-tidx0==1) coeff[i]=1.0;
    }

  if(msg==-20 || msg==-7)
    for(int i=tidx0; i<tidx1; ++i)
      coeff[i]=start.data[i-tidx0];

  if(k==40)
    for(int i=tidx0; i<tidx1; ++i)
      coeff[i]=0.0;

  if(msg != 0 && msg != 1) std::printf("Minimizer terminates with: %d\n", msg);
  if(normalstep)
    if(fixednu) coeff[nidx]=lopt;    
    else coeff[nidx]=xpls2.data[nidx];
  else
    coeff[nidx]=0.0/0.0;
  std::cout << "ModelValue:" << evalRegModel(coeff,omegaL) << std::endl;
}

}  // namespace Kaskade
