/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/lagrangespace.hh"
#include "fem/hierarchicspace.hh"
#include "linalg/direct.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare, createUnitCube
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"

#include "mg/pcg.hh"
#include "mg/apcg.hh"

#include "utilities/kaskopt.hh"

using namespace Kaskade;
#include "poisson.hh"


bool compareAbs(const double x1, const double  x2)
{
  return fabs(x1)>fabs(x2);
}

int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer tutorial program using HB error estimation (SpaceDimension=2)"
            << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true;
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int refinements = getParameter(pt, "refinements", 5),
      order =  getParameter(pt, "order", 1),
      maxAdaptSteps =  getParameter(pt, "maxAdaptSteps", 10),
      verbosity = getParameter(pt, "verbosity", 1);
  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "max. adaptive refine steps     : " << maxAdaptSteps << std::endl;
  std::cout << "output level (verbosity)       : " << verbosity << std::endl << std::endl;
      
  int  direct, onlyLowerTriangle = false;
  
  DirectType directType;
  IterateType iterateType = IterateType::CG;
  PrecondType precondType = PrecondType::NONE;
  MatrixProperties property;
  std::string empty;

  // a typical input from console :
  // ht --maxAdaptSteps 20 --solver.type iterate --solver.iterate IterateType::PCG   oder
  // ht --maxAdaptSteps 20 --solver.type direct --solver.direct PARADISO

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);

  std::cout << "Solver type is " << (direct ? "direct" : "iterative") << std::endl;

  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  directType = static_cast<DirectType>(getParameter(pt, s, 2));

  s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
  iterateType = static_cast<IterateType>(getParameter(pt, s, 0));

  property = MatrixProperties::SYMMETRIC;    //MatrixProperties::POSITIVEDEFINITE;

  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || ( (precondType == PrecondType::ICC) && !direct ) )
  {
    onlyLowerTriangle = true;
    std::cout << 
      "Note: direct solver MUMPS/PARADISO or ICC preconditioner ===> onlyLowerTriangle is set to true!" 
      << std::endl << std::endl;
  }

  //   two-dimensional space: dim=2
  constexpr int dim=2;        
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,1,0> >;
  using VSD = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = PoissonFunctional<double,VSD>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using CoefficientVectors = VSD::CoefficientVectorRepresentation<0,neq>::type;
  using H1ExSpace = FEFunctionSpace<ContinuousHierarchicExtensionMapper<double,LeafView> >;
  using H1ExSpaces = boost::fusion::vector<H1Space const*,H1ExSpace const*>;
  using ExVariableDescriptions = boost::fusion::vector<VariableDescription<1,1,0> >;
  using ExVSD = VariableSetDescription<H1ExSpaces,ExVariableDescriptions>;
  using ErrorEstimator = HierarchicErrorEstimator<LinearizationAt<Functional>,ExVSD>;
  using EstAssembler = VariationalFunctionalAssembler<ErrorEstimator>;
  using ExCoefficientVectors = ExVSD::CoefficientVectorRepresentation<0,1>::type;
  using LinearSpace = VSD::CoefficientVectorRepresentation<0,neq>::type;
  using AssOperator = AssembledGalerkinOperator<Assembler,0,1,0,1>;
  using AssEstOperator = AssembledGalerkinOperator<EstAssembler>;

  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " triangles, " << std::endl;
  std::cout << "      " << gridManager.grid().size(dim-1) << " edges, " << std::endl;
  std::cout << "      " << gridManager.grid().size(dim) << " points" << std::endl;
  // a gridmanager is constructed
  // as connector between geometric and algebraic information
  gridManager.setVerbosity(verbosity);


  // construction of finite element space for the scalar solution T
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  std::string varNames[1] = { "T" };
  VSD variableSetDescription(spaces,varNames);

  Functional F;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  std::cout << "no of variables = " << nvars << std::endl;
  std::cout << "no of equations = " << neq   << std::endl << std::endl;

  Assembler assembler(gridManager,spaces);

  H1ExSpace spaceEx(gridManager,gridManager.grid().leafGridView(), order+1);

  H1ExSpaces exSpaces(&temperatureSpace,&spaceEx);

  std::string exVarNames[2] = { "l", "e"};
  ExVSD exVSD(exSpaces, exVarNames);


  EstAssembler estAssembler(gridManager,exSpaces);
  auto const& is = gridManager.grid().leafIndexSet();

  VSD::VariableSet x(variableSetDescription), dx(variableSetDescription);

  double rTolX = 1.0e-3, aTolX = 1.0e-3, minRefine = 0.2;

  std::vector<std::pair<double,double> > tolX(variableSetDescription.noOfVariables);
  std::vector<std::pair<double,double> > tolXC(variableSetDescription.noOfVariables);
  for (int i=0; i<tolX.size(); ++i)
  {
    tolX[i] = std::make_pair(aTolX,rTolX);
    tolXC[i] = std::make_pair(aTolX/100,rTolX/100);
  }
  std::cout << "Accuracy: rTol = " << rTolX << ", aTol = " << aTolX << std::endl;
  std::cout << "Minimal number of cells to refine in an adaptive step: " << minRefine *100
            << " % of all" << std::endl << std::endl;


  int refSteps = -1;
  bool accurate = false;
  int iteSteps = getParameter(pt, "solver.iteMax", 1000);
  double iteEps = getParameter(pt, "solver.iteEps", 1.0e-6);
  if (!direct)
  {
  std::cout << "max. number of iterations in iterative solver = " << iteSteps << std::endl;
  std::cout << "requested accuracy in iterate solver = " << iteEps << std::endl << std::endl;
  };
  Dune::InverseOperatorResult res;
  double errNorm = 0;

  size_t size = variableSetDescription.degreesOfFreedom(0,1);
  double gamma = 1.0, d = 2;
  double beta = 1.0/sqrt(d*gamma);
  double alpha = (d*gamma-1.0)/(d*(1.0+gamma));
  double qk = 1.0, dNk = size, zk = 0.0;
  double requested = sqrt(1-beta*beta)*tolX[0].first;
  double safety = 1.0;
  double yk = pow(dNk,alpha);
  // printf("gamma=%e, d=%e, beta=%e, alpha=%e, requested=%e, safety=%e\n",
  //         gamma,d,beta,alpha,requested,safety);
      
  do
  {
    refSteps++;

    CoefficientVectors solution(VSD::CoefficientVectorRepresentation<0,neq>::init(spaces));
    solution = 0;
    CoefficientVectors hilfe(VSD::CoefficientVectorRepresentation<0,neq>::init(spaces));
    hilfe = 0;
    assembler.assemble(linearization(F,x));
    CoefficientVectors rhs(assembler.rhs());
    AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, onlyLowerTriangle);

    if (direct)
    {
      directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
    }
    else
    {
      switch (iterateType)
      {
      case IterateType::CG:
      {
        if ( verbosity>0 ) std::cout << "preconditioned cg solver is used" << std::endl;
        JacobiPreconditioner<AssOperator> jacobi(A,1.0);
        const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
        StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
        int lookAhead = getParameter(pt, "solver.lookAhead", 3);
        termination.setLookAhead(lookAhead);
        CG<LinearSpace,LinearSpace> cg(A,jacobi,defaultScalarProduct,termination,verbosity-1);
        cg.apply(hilfe,rhs,res);
      }
      break;
      case IterateType::PCG:
      {
        if ( verbosity>0) std::cout << "   preconditioned cascadic multigrid solver I is used" << std::endl;
        JacobiPreconditioner<AssOperator> jacobi(A,1.0);
        NMIIIPCGSolver<LinearSpace> pcg(A,jacobi,iteEps,iteSteps,verbosity-1);
        pcg.apply(hilfe,rhs,res);
      }
      break;
      case IterateType::APCG:
      {
        if ( verbosity>0) std::cout << "   preconditioned cascadic multigrid solver II is used" << std::endl;
        int addedIterations = getParameter(pt, "solver.APCG.addedIterations", 10);
        JacobiPreconditioner<AssOperator> jacobiPCG(A,1.0);
        NMIIIAPCGSolver<LinearSpace> apcg(A,jacobiPCG,iteEps,iteSteps,verbosity-1,addedIterations);
        apcg.apply(hilfe,rhs,res);
      }
      break;
      default:
        std::cout << "Solver not available" << std::endl;
        throw -111;
      }
      solution.axpy(-1,hilfe);
    }
    
    dx.data = solution.data;

    // Do hierarchical error estimation. Remember to provide the very same underlying problem to the
    // error estimator functional as has been used to compute dx (do not modify x!).
    std::vector<double> errorDistribution(is.size(0),0.0);
    double maxErr = 0.0;
    double errLevel = 0.0;
      
    if (!tolX.empty())
    {
      estAssembler.assemble(ErrorEstimator(LinearizationAt<Functional>(F,x),dx));
      constexpr int estNvars = ErrorEstimator::AnsatzVars::noOfVariables;
      constexpr int estNeq = ErrorEstimator::TestVars::noOfVariables;
      size_t  estNnz = estAssembler.nnz(0,estNeq,0,estNvars,false);
      size_t  estSize = exVSD.degreesOfFreedom(0,estNvars);

      std::vector<double> estData(estNnz), estRhs(estSize), estSolVec(estSize);
      estAssembler.toSequence(0,estNeq,estRhs.begin());

      // iterative solution of error estimator

      AssEstOperator agro(estAssembler);
      Dune::InverseOperatorResult estRes;
      ExCoefficientVectors estRhside(estAssembler.rhs());
      ExCoefficientVectors estSol(ExVSD::CoefficientVectorRepresentation<0,1>::init(exSpaces));
      estSol = 1.0 ;
      JacobiPreconditioner<AssEstOperator> jprec(agro, 1.0);
      jprec.apply(estSol,estRhside); //single Jacobi iteration
      estSol.write(estSolVec.begin());

      // Transfer error indicators to cells.
      for (auto ci=variableSetDescription.gridView.begin<0>(); ci!=variableSetDescription.gridView.end<0>(); ++ci)
      {
        double err = 0;
        auto gix = spaceEx.mapper().globalIndices(*ci);
        for (auto j=gix.begin(); j!=gix.end(); ++j)
          err += fabs(component<0>(estSol)[*j]);
        errorDistribution[is.index(*ci)] = err;
        if (fabs(err)>maxErr) maxErr = fabs(err);
      }

      errLevel = 0.5*maxErr;
      if (minRefine>0.0)
      {
        std::vector<double> eSort(errorDistribution);
        std::sort(eSort.begin(),eSort.end(),compareAbs);
        int minRefineIndex = minRefine*(eSort.size()-1);
        double minErrLevel = fabs(eSort[minRefineIndex])+1.0e-15;
        if (minErrLevel<errLevel)
          errLevel = minErrLevel;
      }

      errNorm = 0 ;
      for (size_t k=0; k < estRhs.size() ; k++ )
        errNorm +=  estRhs[k]*estSolVec[k] ;
    }
    errNorm = sqrt(errNorm) ;
    // apply the Newton correction here
    x += dx;

    if ( verbosity>0) 
      std::cout << "step = " << refSteps << ": " << size << " points,   ||estim. err|| = " << std::setprecision(3) << std::scientific << errNorm << std::resetiosflags( ::std::ios::scientific ) << std::setprecision(6) << std::endl;


	// graphical output of solution, mesh is not yet refined
    std::ostringstream fn;
    fn << "graph/peak-grid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();
    writeVTKFile(x,fn.str(),IoOptions().setOrder(order).setPrecision(7));
    std::cout << "   output written to file " << fn.str() << std::endl;

    if (refSteps>maxAdaptSteps) 
    {
      std::cout << "max. number of refinement steps is reached" << std::endl;
      break;
    }

	// Evaluation of (global/local) error estimator information 
    if (!tolX.empty())
    {
      if (errNorm<requested)
      {
        accurate = true ;
        std::cout << "||estim. error|| is smaller than requested" << std::endl;
      }
      else
      {
        // Refine mesh.
        auto cend = variableSetDescription.gridView.end<0>();
        for (auto ci=variableSetDescription.gridView.begin<0>(); ci!=cend; ++ci)
          if (fabs(errorDistribution[is.index(*ci)]) >= errLevel) gridManager.mark(1,*ci);

        accurate = !gridManager.adaptAtOnce();
      }
    }
    
    size = variableSetDescription.degreesOfFreedom(0,1);
    qk = size/dNk;
    dNk = size;
    yk += pow(dNk,alpha);
    zk = pow(dNk,alpha)*(pow(errNorm/requested,d*alpha)-pow(qk,alpha))/(pow(qk,alpha)-1.0);
    iteEps = safety*beta*errNorm*yk/(yk+zk);
    
  } while (!accurate);

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End heat transfer (peak source) tutorial program" << std::endl;
}
