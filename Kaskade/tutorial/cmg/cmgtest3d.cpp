/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <complex>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <boost/timer.hpp>

#include "dune/common/stdstreams.hh"
#include "dune/grid/sgrid.hh"
#include "dune/grid/uggrid.hh"
#include "dune/grid/common/gridinfo.hh"

#include "fem/assemble.hh"
#include "fem/embedded_errorest.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"
#include "fem/iterate_grid.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/mumps_solve.hh"

#include "io/vtk.hh"
#include "io/amira.hh"

#include "mg/pcg.hh"
#include "mg/apcg.hh"

#include "utilities/kaskopt.hh"

#include "poisson.hh"
#include "amiramesh.hh"

int problemNo = 1;


bool compareAbs(const double x1, const double  x2)
  {
    return fabs(x1)>fabs(x2);
  }

int main(int argc, char *argv[])
  {
    int verbosity = 1;
    bool dump = false;
    boost::property_tree::ptree *pt = GetKaskadeOptions(argc, argv, verbosity, dump);

    std::cout << "Start cascadic multigrid test program" << std::endl;

    int  direct, refs, order, onlyLowerTriangle = false, maxAdaptSteps;
    SolverType directType;
    IterateType iterateType = IterateType::PCG;
    MatrixProperties property;
    std::string empty, problem, functional, geoFile;

	problem = GetParameter(pt, "problem", empty);
	refs = GetParameter(pt, problem+".refs", 0),
    order =  GetParameter(pt, problem+".order", 1),
    maxAdaptSteps = GetParameter(pt, problem+".maxAdaptSteps", 10);
    geoFile = GetParameter(pt, problem+".geoFile", empty);
    functional = GetParameter(pt, problem+".functional", empty);
    problemNo = GetParameter(pt, "names.functional."+functional, 1);
    std::cout << "selected problem " << problem << ", functional=" <<
                 functional << "(" << problemNo <<  "), geoFile=" <<
                 geoFile << std::endl;

    std::string s("names.type.");
    s += GetParameter(pt, "solver.type", empty);
    direct = GetParameter(pt, s, 0);

    s = "names.direct." + GetParameter(pt, "solver.direct", empty);
    directType = static_cast<SolverType>(GetParameter(pt, s, 2));

	s = "names.iterate." + GetParameter(pt, "solver.iterate", empty);
	iterateType = static_cast<IterateType>(GetParameter(pt, s, 0));

    property = MatrixProperties::POSITIVEDEFINITE;

    if (((property == MatrixProperties::SYMMETRIC)||(property == MatrixProperties::POSITIVEDEFINITE))&&
        ((directType == DirectType::MUMPS)||(directType == DirectType::PARDISO)))
      {
        onlyLowerTriangle = true;
      }

    int const dim = 3;
    int k;
    typedef Dune::UGGrid<dim> Grid;
    Grid *grid = new Grid(20000);
    Dune::GridFactory<Grid> factory(grid);

    DuneAmiraMesh<dim,dim,Grid> mesh(geoFile.c_str());
    mesh.InsertUGGrid(factory);
	    
//    std::auto_ptr<Grid> grid( factory.createGrid() );
    grid = factory.createGrid();
    // the coarse grid will be refined refs times
    for (k=0; k<refs; k++)
      {
        grid->globalRefine(1);
	  }
    // some information on the refined mesh
    std::cout << "Grid: " << grid->size(0) << " tetrahedra, " << std::endl;
    std::cout << "      " << grid->size(1) << " triangles, " << std::endl;
    std::cout << "      " << grid->size(1) << " edges, " << std::endl;
    std::cout << "      " << grid->size(dim) << " points" << std::endl;
    // a gridmanager is constructed 
    // as connector between geometric and algebraic information
    GridManager<Grid> gridManager(grid);    

    typedef Grid::LeafGridView LeafView;
    // construction of finite element space for the scalar solution T
    typedef FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> > H1Space;

    H1Space temperatureSpace(gridManager,gridManager.grid().leafView(),
                             order);
    typedef boost::fusion::vector<H1Space const*> Spaces;
    Spaces spaces(&temperatureSpace);
    // VariableDescription<int spaceId, int components, int Id>
    // spaceId: number of associated FEFunctionSpace
    // components: number of components in this variable
    // Id: number of this variable
    typedef boost::fusion::vector<VariableDescription<0,1,0> >
      VariableDescriptions;
    std::string varNames[1] = { "T" };
    typedef VariableSetDescription<Spaces,VariableDescriptions> VariableSet;
    VariableSet variableSet(spaces,varNames);

    typedef PoissonFunctional<double,VariableSet> Functional;
    Functional F;
    typedef VariationalFunctionalAssembler<LinearizationAt<Functional> > GOP;
    typedef GOP::TestVariableRepresentation<>::type Rhs;
    GOP gop(gridManager.signals,spaces);

    typedef FEFunctionSpace<ContinuousHierarchicExtensionMapper<double,LeafView> > H1ExSpace;
    H1ExSpace spaceEx(gridManager,gridManager.grid().leafView(), order+1);

    typedef boost::fusion::vector<H1Space const*,H1ExSpace const*> H1ExSpaces;
    H1ExSpaces exSpaces(&temperatureSpace,&spaceEx);

    typedef boost::fusion::vector<VariableDescription<1,1,0> > ExVariableDescriptions;
    typedef VariableSetDescription<H1ExSpaces,ExVariableDescriptions> ExVariableSet;
    std::string exVarNames[2] = { "l", "e"};
    ExVariableSet exVariableSet(exSpaces, exVarNames);

    typedef HierarchicErrorEstimator<LinearizationAt<Functional>,ExVariableSet> ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> EstGOP;

    EstGOP estGop(gridManager.signals,exSpaces);
	typedef EstGOP::AnsatzVariableRepresentation<> Ansatz;
	typedef EstGOP::TestVariableRepresentation<> Test;

	typedef VariableSet::Grid::Traits::LeafIndexSet IS ;
	IS const& is = gridManager.grid().leafIndexSet();

    VariableSet::VariableSet x(variableSet), dx(variableSet), tmp(variableSet);

    double rTolX = GetParameter(pt, problem+".rTolX", 1.0e-4),
           aTolX = GetParameter(pt, problem+".aTolX", 1.0-2),
           minRefine = GetParameter(pt, problem+".minRefine", 0.0);
    std::vector<std::pair<double,double> > tolX(variableSet.noOfVariables);
    std::vector<std::pair<double,double> > tolXC(variableSet.noOfVariables);
    for (int i=0; i<tolX.size(); ++i)
      {
        tolX[i] = std::make_pair(aTolX,rTolX);
        tolXC[i] = std::make_pair(aTolX/100,rTolX/100);
     }
	std::cout << "rTolX = " << rTolX << ", aTolX = " << aTolX <<
	             ", minRefine = " << minRefine << std::endl ;

	IoOptions options;
	options.outputType = IoOptions::ascii;

    int refSteps = 0;
    bool accurate = false;
	int iteSteps = GetParameter(pt, "solver.iteMax", 1000);
	int verbose = GetParameter(pt, "solver.verbose", 1);
	double iteEps = GetParameter(pt, "solver.iteEps", 1.0e-6);
	typedef GOP::TestVariableRepresentation<>::type LinearSpace;
	Dune::InverseOperatorResult res;
	double errNorm;

	size_t size = variableSet.degreesOfFreedom(0,1);
	double gamma = 1.0, d = dim;
	double beta = 1.0/sqrt(d*gamma);
	double alpha = (d*gamma-1.0)/(d*(1.0+gamma));
	double qk = 1.0, dNk = size, zk = 0.0;
	double requested = sqrt(1-beta*beta)*tolX[0].first;
	double safety = 1.0;
	double yk = pow(dNk,alpha);
printf("gamma=%e, d=%e, beta=%e, alpha=%e, requested=%e, safety=%e\n",
       gamma,d,beta,alpha,requested,safety);
//printf("    %10ld %e   %e   %e   %e\n", size, qk, dNk, yk, zk);
    do
      {
		typedef GOP::AnsatzVariableRepresentation<>::type Sol;
		Sol solution(GOP::AnsatzVariableRepresentation<>::init(gop));
		Sol hilfe(GOP::AnsatzVariableRepresentation<>::init(gop));
		solution = 0;
		hilfe = 0;
		gop.assemble(linearization(F,x));
		Rhs rhs(GOP::TestVariableRepresentation<>::rhs(gop));
		AssembledGalerkinOperator<GOP,0,1,0,1> A(gop, onlyLowerTriangle);
		AssembledGalerkinOperator<GOP,0,1,0,1>::matrix_type tri(A.getmat());

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
				   JacobiPreconditioner<GOP,0,1,0> jacobi(gop,1.0);
				   Dune::CGSolver<LinearSpace> cg(A,jacobi,iteEps,iteSteps,verbose);
				   cg.apply(hilfe,rhs,res);
			     }
			     break;
			   case IterateType::PCG:
			     {
			       JacobiPreconditioner<GOP,0,1,0> jacobi(gop,1.0);
			       Kaskade7::NMIIIPCGSolver<LinearSpace> pcg(A,jacobi,iteEps,iteSteps,verbose);
			       pcg.apply(hilfe,rhs,res);
			     }
			     break;
			   case IterateType::APCG:
			     {
    			   int addedIterations = GetParameter(pt, "solver.APCG.addedIterations", 10);
			       JacobiPreconditioner<GOP,0,1,0> jacobiPCG(gop,1.0);
			       Kaskade7::NMIIIAPCGSolver<LinearSpace> apcg(A,jacobiPCG,iteEps,iteSteps,verbose,addedIterations);
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
printf("Solved!\n");

		std::ostringstream fn;
		fn << "graph3d/cmg-grid";
		fn.width(3);
		fn.fill('0');
		fn.setf(std::ios_base::right,std::ios_base::adjustfield);
		fn << refSteps;
		fn.flush();
	    LeafView leafView = gridManager.grid().leafView();
	    writeVTKFile(leafView,variableSet,x,fn.str(),options,order);

        // Do hierarchical error estimation. Remember to provide the very same underlying problem to the
        // error estimator functional as has been used to compute dx (do not modify x!).
        if (!tolX.empty())
          {
	        tmp *= 0 ;
	        estGop.assemble(ErrorEstimator(LinearizationAt<Functional>(F,x),dx));
		    int const estNvars = ErrorEstimator::AnsatzVars::noOfVariables;
		    int const estNeq = ErrorEstimator::TestVars::noOfVariables;
		    size_t  estNnz = estGop.nnz(0,estNeq,0,estNvars,false);
		    size_t  estSize = exVariableSet.degreesOfFreedom(0,estNvars);

		    std::vector<int> estRidx(estNnz), estCidx(estNnz);
		    std::vector<double> estData(estNnz), estRhs(estSize), estSolVec(estSize);
		    estGop.toSequence(0,estNeq,estRhs.begin());

		  // iterative solution of error estimator

		    AssembledGalerkinOperator<EstGOP> E(estGop);
		    Dune::InverseOperatorResult estRes;
		    Test::type estRhside(Test::rhs(estGop) ) ;
		    Ansatz::type estSol(Ansatz::init(estGop) ) ;
		    estSol = 1.0 ;
		    JacobiPreconditioner<EstGOP> jprec(estGop, 1.0);
		    jprec.apply(estSol,estRhside); //single Jacobi iteration
		    estSol.write(estSolVec.begin());
	  
	  // Transfer error indicators to cells.
		    std::vector<double> errorDistribution(is.size(0),0.0);
		    typedef VariableSet::GridView::Codim<0>::Iterator CellIterator ;
		    double maxErr = 0.0;
		    for (CellIterator ci=variableSet.gridView.begin<0>(); ci!=variableSet.gridView.end<0>(); ++ci)
		      {
				typedef H1ExSpace::Mapper::GlobalIndexRange GIR;
				double err = 0;
				GIR gix = spaceEx.mapper().globalIndices(*ci);
				for (GIR::iterator j=gix.begin(); j!=gix.end(); ++j)
				  err += fabs(boost::fusion::at_c<0>(estSol.data)[*j]);
				errorDistribution[is.index(*ci)] = err;
				if (fabs(err)>maxErr) maxErr = fabs(err);
		      }
		  
		    double errLevel = 0.5*maxErr;

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
		    for (k=0; k < estRhs.size() ; k++ )
		      errNorm +=  fabs(estRhs[k]*estSolVec[k]) ;

	        if (errNorm<requested)
	          accurate = true ;
	        else
	          {

	  // Refine mesh.
	
				int noToRefine = 0;
				double alphaSave = 1.0 ;
				std::vector<bool> toRefine( is.size(0), false ) ; //for adaptivity in compression
				std::vector< std::vector<bool> > refinements;
				for (CellIterator ci=variableSet.gridView.begin<0>(); ci!=variableSet.gridView.end<0>(); ++ci)
				  if (fabs(errorDistribution[is.index(*ci)]) >= alphaSave*errLevel)
					{
					  noToRefine++;
					  toRefine[is.index(*ci)] = true ;
					  gridManager.mark(1,*ci);
					}
	
				refinements.push_back(toRefine);
				accurate = !gridManager.adaptAtOnce(); 
	          }
        }
        // apply the Newton correction here
       x += dx;
        
        
//           { // Compute spatial error estimate
//     		   std::vector<double> norm2, error2;
//     		   norm2.resize(1);
//     		   error2.resize(1);
//             tmp = x;
//             projectHierarchically(variableSet,tmp);
//             tmp -= x;
// 
//           // perform mesh adaptation
//             accurate = embeddedErrorEstimator(variableSet, tmp, x, scal, tolX, gridManager, norm2, error2);
//             printf(" %e %e\n", norm2[0], error2[0]);
//             if (!accurate)
//               {
//                 nnz = gop.nnz(0,1,0,1,onlyLowerTriangle);
//                 size = variableSet.degreesOfFreedom(0,1);
// 				ridx.resize(nnz);
// 				cidx.resize(nnz);
// 				data.resize(nnz);
// 				rhs.resize(size);
// 				sol.resize(size);
//               }
//           }
        refSteps++;

        printf("%3d %10ld %5d   %e   %e\n", refSteps, size, res.iterations, errNorm,
               iteEps);

        size = variableSet.degreesOfFreedom(0,1);
	    qk = size/dNk;
		dNk = size;
	    yk += pow(dNk,alpha);
	    zk = pow(dNk,alpha)*(pow(errNorm/requested,d*alpha)-pow(qk,alpha))/(pow(qk,alpha)-1.0);
        iteEps = safety*beta*errNorm*yk/(yk+zk);

//        printf("  %10ld %e   %e   %e\n", size, yk, zk, yk+zk);

        if (refSteps>maxAdaptSteps) break;
// 
// 		std::ostringstream fn;
// 		fn << "graph3d/cmg-grid";
// 		fn.width(3);
// 		fn.fill('0');
// 		fn.setf(std::ios_base::right,std::ios_base::adjustfield);
// 		fn << refSteps;
// 		fn.flush();
// 	    LeafView leafView = gridManager.grid().leafView();
// 	    writeVTKFile(leafView,variableSet,x,fn.str(),options,order);

      } while (!accurate);     

    std::cout << "End cmgtest" << std::endl;
  }
