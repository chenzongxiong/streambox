#include "hyp_ip.hh"
#include <cmath>
#include <fstream>

namespace Kaskade
{

void HypIP::terminationMessage(int flag)
{
  std::cout << "Composite Step: "; 
  switch((int)flag)
  {
  case 1  : std::cout << "Desired accuracy reached!" << std::endl; break;
  case -1 : std::cout << "Regularity test failed!" << std::endl;; break;
  case -2 : std::cout << "Maximum number of iterations reached!" << std::endl; break;
  default : Algorithm::terminationMessage(flag); 
  }
}

  /// Return true, if convergence is detected, false otherwise
  Convergence HypIP::convergenceTest
  (AbstractFunctionSpaceElement const& correction, AbstractFunctionSpaceElement const& iterate, std::vector<double>& coeff)
  {
    double nrm=normC(correction);
    std::cout << "|dx|" << nrm << std::endl;
    for(int i=0; i< coeff.size(); ++i)
      if(coeff[i] < 0.5) return Convergence::Missed;
    if(nrm < p.desiredAccuracy) return Convergence::Achieved;
    
    return Convergence::Missed;
  }


  AcceptanceTest HypIP::evaluateCorrection(
    AbstractFunctionSpaceElement& correction, 
    AbstractLinearization& lin,
    CUBThetaModelFunction& mF,
    double f0,
    std::vector<double>& coeff)
{ 
  double ThetaMax = p.ThetaMax;
  double plausiblefactor=10.0;
  bool passed(true);
  double model2=mF.evalModel(coeff);
  
  double normCcorr = normC(correction);
  double normLcorr = normL(correction);

// Acceptance test w.r.t the equality constraints

  chart.addPerturbation(*trialIterate,correction,lin,basisVectors);

  int nidx = searchSpace.getDimension()-1;
  
  std::vector<double> coeffnormal(coeff.size(),0.0);
    
  coeffnormal[nidx]=coeff[nidx];

  double modelnormal=mF.evalModel(coeffnormal);


  if(normalSolver)
  {    

    std::unique_ptr<AbstractLinearization> sLin=functionalN->getLinearization(*trialIterate);

    normCcorr = normC(correction);

    normalSolver->resolve(*scorrection,*sLin);

    if(report >= 2) { scorrection->print("  dx_");}

    double normCscorr = normC(*scorrection);

    scorrection->axpy(coeff[nidx]-1.0,*basisVectors[nidx]);

    if(report >= 2) { scorrection->print("  ds:");}
    
    snormalstepnorm=normC(*scorrection);

    if(report >= 1) {std::cout << "|dn|:" << normalstepnorm << " |dx|:" << normCcorr << " |dx_|:" << normCscorr << " |ds|:" << snormalstepnorm << std::endl;}
    if(normCcorr > 1e-15)
      Theta=snormalstepnorm/normCcorr;
    else
      Theta = 0.0;
    double oCnew=2*snormalstepnorm/normCcorr/normCcorr;
    std::cout << "oC_2:" << oCnew << " Theta:" << Theta << std::endl;
    if(Theta < 5e-2 && oCnew > omegaC)
    {
      std::cout << "No oC_1 update: " << oCnew << " round-off errors?" << std::endl;
    }
    else
    {
      if(oCnew-omegaC > p.sigma*omegaC)
      {
        std::cout << " significant update from oC_1 = " << omegaC << std::endl;
        passed=false;
      }  
      else 
      {
        std::cout << std::endl;
      }
      omegaC=std::min(oCnew,plausiblefactor*omegaC);
    }
      std::cout << "|dx_|:" << normCscorr << " |dn|:" << normalstepnorm << std::endl;
    if(Theta > ThetaMax)
    {
      std::cout << "Failed by contraction" << std::endl;
      passed=false;
      if(Theta >= 5*ThetaMax)
        return AcceptanceTest::Failed;
    }
  
    scorrection->axpy(1.0,correction);
    chart.addPerturbation(*trialIterate,*scorrection,lin,basisVectors);
  }
// Acceptance Test due to nonlinearity of the constraints passed!
// Continuing with acceptance test for functional descent


  std::unique_ptr<AbstractLinearization> barLin=functionalN->getLinearization(*trialIterate);
  double fval=barLin->eval();
  double modelerror=fval-f0-model2;
  double oLnew=std::min(omegaL*plausiblefactor,std::fabs(modelerror)/std::pow(normLcorr,3)*6.0);

  if(report>=2) {std::cout << "f0:" << f0 << " f:" << fval <<  " m:" << f0+model2 << " f-f0:" << fval-f0 << " f-m:" << fval-f0-model2 << " mn:" << f0+modelnormal << std::endl;}
  if((std::fabs(fval)+1.0)*1e-14 < std::fabs(modelerror) || oLnew < omegaL)
  {
    std::cout << " oL_2 = " << oLnew << " oL_1 = "<< oLnew*normLcorr << " oL_0 = "<< oLnew*normLcorr*normLcorr;
    if(oLnew-omegaL> p.sigma*omegaL)
    {
      std::cout << " significant update from oL_2 = " << omegaL;
           if(fval > f0  || oLnew > plausiblefactor*omegaL)// && std::fabs(coeff[nidx]-1.0)<= 1e-5) 
      {
        std::cout << " and tendency to ascent";
        passed=false;
      }
    }
    std::cout << std::endl;
    if(passed) omegaL=std::max(oLnew,1.0/plausiblefactor*omegaL); else omegaL = std::max(oLnew,omegaL);
  }
  else
  {
    std::cout << "No oL_2 update: " << oLnew << " round-off errors?" << std::endl;
  }
  
  if(passed)
  {
    RemainderTerm = snormalstepnorm;
    return AcceptanceTest::Passed;
  }
  else
    return AcceptanceTest::Failed;
}

//   AcceptanceTest HypIP::evaluateCorrectionOld(
//     AbstractFunctionSpaceElement& correction, 
//     AbstractLinearization& lin,
//     CUBModelFunction& mF,
//     double f0,
//     std::vector<double>& coeff)
// { 
//   double plausiblefactor=10.0;
//   bool passed(true);
//   double model2=mF.evalModel(coeff);
//   double normCcorr = normC(correction);
//   double normLcorr = normL(correction);

// // Acceptance test w.r.t the equality constraints

//   chart.addPerturbation(*trialIterate,correction,lin,basisVectors);

//   if(normalSolver)
//   {    

//     std::unique_ptr<AbstractLinearization> sLin=functionalN->getLinearization(*trialIterate);

//     normalSolver->resolve(*scorrection,*sLin);

//     if(report >= 2) { scorrection->print("  dx_");}

//     double normscorr = normC(*scorrection);
//     int nidx = searchSpace.getDimension()-1;
//     scorrection->axpy(coeff[nidx]-1.0,*basisVectors[nidx]);

//     if(report >= 2) { scorrection->print("  ds:");}
    
//     double numerator=normC(*scorrection);
//     double normalstepnorm= normC(*basisVectors[nidx]);
//     if(report >= 1) {std::cout << "|dn|:" << normalstepnorm << " |dx|:" << normCcorr << " |dx_|:" << normscorr << " |ds|:" << numerator << std::endl;}
//     double Theta=numerator/normCcorr;
//     double oCnew=2*numerator/normCcorr/normCcorr;
//     std::cout << "oC_2:" << oCnew << " Theta:" << Theta << std::endl;
//     if(((normCcorr<1e-8 && normscorr < 1e-2*normCcorr) || normCcorr < 1e-4) && oCnew > omegaC)
//     {
//       std::cout << "No oC_1 update: " << oCnew << " round-off errors?" << std::endl;
//     }
//     else
//     {
//       if(oCnew-omegaC > p.sigma*omegaC)
//       {
//         std::cout << " significant update from oC_1 = " << omegaC << std::endl;
// //        passed=false;
//       }  
//       else 
//       {
//         std::cout << std::endl;
//       }
//       omegaC=std::min(oCnew,plausiblefactor*omegaC);
//     }
//       std::cout << "|dx_|:" << normscorr << " |ds|*oC*|dx|:" << normscorr*omegaC*normCcorr << std::endl;
//       std::cout << "q:" << normscorr*std::min(omegaC*normCcorr,1.0)  << " |dn|:" << normalstepnorm << " 2*Td/oC:" << 2*p.Thetad/omegaC << std::endl;
//     if(normscorr*std::min(omegaC*normCcorr,1.0) >= std::max(normalstepnorm,2*p.Thetad/omegaC))
//     {
//       std::cout << "Failed by contraction" << std::endl;
//       passed=false;
//       if(normscorr*omegaC*normCcorr >= 5*std::max(normalstepnorm,2*p.Thetad/omegaC))
//         return AcceptanceTest::Failed;
//     }
  
//     scorrection->axpy(1.0,correction);
//     chart.addPerturbation(*trialIterate,*scorrection,lin,basisVectors);
//   }
// // Acceptance Test due to nonlinearity of the constraints passed!
// // Continuing with acceptance test for functional descent


//   std::unique_ptr<AbstractLinearization> barLin=functionalN->getLinearization(*trialIterate);
//   double fval=barLin->eval();
//   double modelerror=fval-f0-model2;
//   double oLnew=std::min(std::fabs(modelerror)/std::pow(normLcorr,3)*6.0,plausiblefactor*omegaL);

//   if(report>=2) {std::cout << "fo:" << f0 << " fn:" << fval <<  "fm:" << f0+model2 << " fn-fo:" << fval-f0 << " fn-fm:" << fval-f0-model2 << std::endl;}
//   if(std::fabs(fval)*1e-8 < std::fabs(modelerror) || oLnew < omegaL)
//   {
//     std::cout << " oL_2 = " << oLnew << " oL_1 = "<< oLnew*normLcorr << " oL_0 = "<< oLnew*normLcorr*normLcorr;
//     if(oLnew-omegaL> p.sigma*omegaL)
//     {
//       std::cout << " significant update from oL_2 = " << omegaL;
//       if(modelerror > 0) 
//       {
//         std::cout << " and tendency to ascent";
//         passed=false;
//       }
//     }
//     std::cout << std::endl;
//     if(passed) omegaL=oLnew; else omegaL = std::max(oLnew,omegaL);
//   }
//   else
//   {
//     std::cout << "No oL_2 update: " << oLnew << " round-off errors?" << std::endl;
//   }
  
//   if(passed)
//   {
//     return AcceptanceTest::Passed;
//   }
//   else
//     return AcceptanceTest::Failed;
// }

void HypIP::dampingLoop()
{
  double nu0last=nu0;
  
    normalLinearization = functionalN->getLinearization(*iterate);
    tangentialLinearization = functionalT->getLinearization(*iterate);


    normalLinearization->precompute();

    if(searchSpace.normalFactorizationPresent())
    {
      searchSpace.getNormalSolver()->resolve(*scorrection,*normalLinearization);
      sbarstepnorm=normC(*scorrection);
      std::cout << "Test: " << sbarstepnorm << std::endl;
      theta2=sbarstepnorm/snormalstepnorm;
    }

    if(report>=2) { iterate->print("   x:"); }
    normL.setOrigin(*normalLinearization);
    normC.setOrigin(*normalLinearization);


    searchSpace.computeBasisVectors(basisVectors,
                                    *iterate,
                                    *normalLinearization,
                                    *tangentialLinearization,
                                    normC,p.ThetaAim,omegaC,omegaL,omegaL,report,nu0,normalstepnorm);

    RelNormalSteps = normalstepnorm/RemainderTerm;

    if(nu0last==1.0 && RelNormalSteps > 1.0) {
      p.safetyFactor *= 2.0;
      std::cout << "Issue with feasibility: " << RelNormalSteps << std::endl;
    }
    else 
    {
      p.safetyFactor = 1.0;
    }

  for(int i=0; i<searchSpace.getDimension(); ++i)
    std::cout << "|x["<<i <<"]|: " << normL(*(basisVectors[i])) << "  ";
  std::cout << std::endl;
//  iterate->writeToFile("iterate.log",true);
  double f0=normalLinearization->eval();
  CUBThetaModelFunction modelFunction(*tangentialLinearization,*normalLinearization, basisVectors, normL,normC, searchSpace,dimx); 
  do
  {
    modelFunction.reparametrize(omegaL,omegaC,omegaL*p.gamma,p); 
    modelFunction.getMinimalCoefficients(coeff);
    if(report>=1) 
    {
      std::cout << "Coeff:[ ";
      for(int i=0; i<searchSpace.getDimension(); ++i)
        std::cout << coeff[i] << " ";
      std::cout << "]" <<std::endl;
    }
    searchSpace.getLinearCombination(basisVectors, coeff, *correction);
    if(report>=2) { correction->print("   dx:"); }
      

  }
  while(evaluateCorrection(*correction ,  
                           *tangentialLinearization, modelFunction, f0, coeff)==AcceptanceTest::Failed);

}

int HypIP::refinementLoop()
{
  int maxRef = 10;
  int maxPatches = 500000;
  double relerror = 0.25;
  for(int i=0; i<maxRef; ++i)
  {
    
    grid->adapt();
    if(report) std::cout << "------- Refinement Step: " << i+1 << " ------ #Patches: " << grid->size() << " ---- " << std::endl;

    dampingLoop();
    estimate.reset(errorEstimator->createEstimate(*trialIterate,basisVectors,coeff,*normalLinearization,*tangentialLinearization).release());
    double err=estimate->absoluteError();
    double nc=normC(*scorrection);
    double requirederror=relerror*nc/coeff[searchSpace.getDimension()-1];
    bool accuracyReached = err<requirederror || err+nc < p.desiredAccuracy;
    double desacc=std::max(requirederror,p.desiredAccuracy-err);
    std::cout << "Error: " << err << " Correction: " << requirederror << " Relation: " << err/requirederror << std::endl;
    if(estimate.get() && grid->size() < maxPatches && maxRef!=0)
    {
      double cbulk = std::sqrt(1-std::pow(desacc/err*0.8,2));
      double bulk = std::max(0.5, std::min(cbulk,0.9));
      std::cout << "Bulk-Parameter: " << bulk << std::endl;
      grid->mark(*estimate,bulk);
    }
    
    if(accuracyReached || grid->size() >= maxPatches) return 1;
  }
  return 2;
 }



int HypIP::runAlgorithm() 
{
  std::_Ios_Openmode mode= std::ios::out;
  if(!p.firstmu) mode= std::ios::app;
  std::ofstream omegafile("omega.log",mode);
  std::ofstream diagfile("diagnostics.log",mode);
  std::ofstream coefffile("coeff.log",mode);
  std::ofstream values("values.log",mode);
  std::ofstream normcorr("normcorr.log",mode);
  std::ofstream steps("steps.log",mode);
  bV.resize(0);
  basisVectors.resize(0);
  normalSolver=searchSpace.getNormalSolver();
  for(int i=0; i<searchSpace.getMaxDimension(); ++i)
  {
    bV.push_back(std::tr1::shared_ptr<AbstractFunctionSpaceElement>(iterate->clone().release()));
    basisVectors.push_back(bV[i].get());
  }
  if(p.firstmu)
  {
    omegaL=0.3;
    omegaC=1.0;
//    iterate->writeToFile("iterate.log",false);
//    iterate->writeToFile("normalstep.log",false);
//    iterate->writeToFile("tangentstep.log",false);
  } else
  {
    std::cout << "Starting values for omegas: " << omegaL << "  " << omegaC << std::endl;
  }
  // Main loop
  for(step=1; step <= maxSteps(); step++)
  {     
    if(report) std::cout <<  "\n --------------- Optimization-Step: " << step << " ----------------------" << std::endl;
    if(grid==0 || errorEstimator == 0)
      dampingLoop();
    else
    {
      refinementLoop();
    }
    iterate->swap(*trialIterate);

//    (basisVectors[0])->writeToFile("tangentstep.log",true);
//    (basisVectors[1])->writeToFile("normalstep.log",true);
    normalLinearization.reset();

    omegafile << omegaL << " " << omegaC << " " << std::endl;
    diagfile << omegaL << " " << omegaC << " | " << Theta << " " << theta2 << " " << RelNormalSteps << " | " << coeff[0] << " " << coeff[1] << " " << coeff[2] << std::endl;
    coefffile << coeff[0] << " " << coeff[1] << " " << coeff[2] << std::endl;
    normcorr << normC(*scorrection) << std::endl;

    if(convergenceTest(*scorrection, *iterate,coeff)==Convergence::Achieved)
    {
//      iterate->writeToFile("iterate.log",true);
      normalLinearization = functionalN->getLinearization(*iterate);
      values << normalLinearization->eval() << std::endl;
      steps << step << std::endl;
      return 1;
    }
  }
  return -2;
}

void HypIP::solve(AbstractFunctional* fN, AbstractFunctional* fT,AbstractFunctionSpaceElement& x)
{
  functionalN=fN;
  functionalT=fT;
  iterate = x.clone();
  trialIterate=x.clone(); 
  correction=x.clone();
  scorrection=x.clone(); 
  p.termination=algorithmWrapper();
  if(p.termination > 0 || p.termination == -2) x=*iterate;
}

} // namespace Kaskade
