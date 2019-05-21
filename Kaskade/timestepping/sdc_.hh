/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SDC__HH
#define SDC__HH

#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <boost/timer/timer.hpp>

#include "dune/istl/solvers.hh"
#include "dune/istl/preconditioners.hh"
#include "dune/istl/bcrsmatrix.hh"

#include "fem/embedded_errorest.hh"
#include "fem/iterate_grid.hh"
#include "timestepping/extrapolation.hh"
#include "timestepping/semieuler.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"
#include "fem/norms.hh"

#include "linalg/factorization.hh"
#include "linalg/umfpack_solve.hh"
#include "linalg/mumps_solve.hh"
#include "linalg/superlu_solve.hh"

#include "linalg/trivialpreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/iluprecond.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/hyprecond.hh"

#include "utilities/enums.hh"

#include "sdc.hh"


#include <iostream>
#include "dune/common/dynmatrix.hh"
#include "dune/common/dynvector.hh"

#include "utilities/enums.hh"


namespace Kaskade
{
  
template <class Eq>
class Sdc_
{
  public:
    typedef Eq                                                     EvolutionEquation;
    typedef typename EvolutionEquation::AnsatzVars::VariableSet State;

  private:
    typedef SemiLinearizationAt<SemiImplicitEulerStep<EvolutionEquation> > Linearization;
    typedef VariationalFunctionalAssembler<Linearization> Assembler;
    typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,EvolutionEquation::TestVars::noOfVariables>::type CoefficientVectors;

  public:
    /**
     * Constructs an ODE integrator. The arguments eq and ansatzVars
     * have to exist during the lifetime of the integrator.
     */
  Sdc_(GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager_,
          EvolutionEquation& eq_, typename EvolutionEquation::AnsatzVars const& ansatzVars_,
          typename EvolutionEquation::AnsatzVars::VariableSet x_,
          DirectType st=DirectType::MUMPS, PrecondType precondType_ = PrecondType::ILUK,
          int verbosity_ = 1):
          gridManager(gridManager_), ansatzVars(ansatzVars_), eq(&eq_,0),
          assembler(gridManager,ansatzVars.spaces) , u(x_),dx(x_),
          rhsAssemblyTime(0.0), matrixAssemblyTime(0.0), precAssemblyTime(0.0), factorizationTime(0.0), solutionTime(0.0),
          initSolutionTime(0.0), estimateTime(0), directType(st), precondType(precondType_), verbosity(verbosity_)
          
    {}

 
  State const& stepCardio(State const& x, double const dt, std::vector<std::pair<double,double> > const& tolX)
  {		  		 
    boost::timer::cpu_timer timer;

    int const n = 3, kEnd = 3;   
    //std::cout << "SDC method: n = " << n << ",  kEnd = " << kEnd << std::endl;

    //LobattoTimeGrid grid(n,eq.time(),eq.time()+dt);
    RadauTimeGrid grid(n,eq.time(),eq.time()+dt);
    std::cout << "SDC-Matrix method, reaction terms evaluated only in nodes: n = " << n << ",  kEnd = " << kEnd << std::endl;

    Dune::DynamicVector<double> const& pts = grid.points();

    // S is provided an used in sdcIterationStep()
    //Dune::DynamicMatrix<double> const& S = grid.integrationMatrix();
     
    // Shat is computed in sdcIterationStep(...) called below
    //SDCTimeGrid::RealMatrix Shat;
    //eulerIntegrationMatrix(grid, Shat);
    //luIntegrationMatrix(grid, Shat);
    //std::cout << "Shat = \n" << Shat << std::endl;
	
    double tau1=0;
							
    int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
    int const neq = EvolutionEquation::TestVars::noOfVariables;
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    size_t  size = ansatzVars.degreesOfFreedom(0,nvars);
    
    std::cout << "nvars = " << nvars << std::endl;
    std::cout << "  neq = " <<   neq << std::endl;
    std::cout << "  dof = " <<  size << std::endl;

    typedef  AssembledGalerkinOperator<Assembler,0,neq,0,neq> Op;

    double const t = eq.time();
    eq.temporalEvaluationRange(t,t+dt);

    double const tau = dt;
    eq.setTau(tau);

    // mesh adaptation loop. Just once if i>0.
    gridManager.setVerbosity(verbosity);
    bool accurate = false;

    bool const iterative = true;       // false: use a direct solver,      true: use an iterative solver
    std::vector<double> rhs(size), sol(size);

    // Evaluate and factorize matrix B(t)-tau*J

    std::vector<State> currSol(n+1,x);
    std::vector<State> lastSol(n+1,x);
				
      typedef typename Op::field_type field_type;
      typedef Dune::FieldMatrix<field_type,1,1> Field_Matrix;
      typedef Dune::BCRSMatrix<Field_Matrix> BCRS_Matrix;
      BCRS_Matrix  matM;
      BCRS_Matrix  matStiff;

		
    //for (int kk=0; kk<kEnd; kk++)
    int kk=-1;
    do
    {
      kk++;
      size_t  dofs = ansatzVars.degreesOfFreedom(0,nvars);
      std::cout << "sweep: k = " << kk << ", degree of freedom = " << dofs << std::endl;

      
      eq.time(t);
      
      CoefficientVectors rU(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      
      //for (int k=0; k< boost::fusion::at_c<0>(rhs1[0].data).size(); k++)
      //std::cout << " size of rU(0) = " << boost::fusion::at_c<0>(rU.data).size() << std::endl;
      //std::cout << " size of rU(1) = " << boost::fusion::at_c<1>(rU.data).size() << std::endl;
      
      rU = 0;
      std::vector<CoefficientVectors> rUiMat(n+1,rU);
      std::vector<CoefficientVectors> rJiMat(n+1,rU);
      std::vector<CoefficientVectors> integMatS(n,rU);
      std::vector<CoefficientVectors> dUiMat(n+1,rU);
      CoefficientVectors sHatDu(rU);
      CoefficientVectors currsol(rU);
      CoefficientVectors lastsol(rU);						

      CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      solution=0;
      dUiMat[0] = solution;





//       typedef typename Op::field_type field_type;
//       typedef Dune::FieldMatrix<field_type,1,1> Field_Matrix;
//       typedef Dune::BCRSMatrix<Field_Matrix> BCRS_Matrix;
      
      if (!accurate)
      {
      // mass matrix M 
        boost::timer::cpu_timer massMatrixTimer;
        timer.start();
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(false);
        eq.setOnlyStiffnessMatrix(false);
        eq.setOnlyF_uMatrix(false);
        eq.setOnlyMassMatrix(true);
        assembler.assemble(Linearization(eq,lastSol[0],lastSol[0],dx),Assembler::MATRIX);               
        eq.setOnlyMassMatrix(false);
        Op M(assembler);      
        //BCRS_Matrix  matM(M.template get<BCRS_Matrix>());    //  wenn matM lokal in loop
        matM = M.template get<BCRS_Matrix>();
        std::cout << "computing time for mass matrix: " << boost::timer::format(massMatrixTimer.elapsed()) << "\n";
        matrixAssemblyTime += (double)(timer.elapsed().wall)/1e9;
        
        // stiffness matrix A 
        boost::timer::cpu_timer stiffMatrixTimer;
        timer.start();
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(false);
        eq.setOnlyMassMatrix(false);
        eq.setOnlyF_uMatrix(false);
        eq.setOnlyStiffnessMatrix(true);
        assembler.assemble(Linearization(eq,lastSol[0],lastSol[0],dx),Assembler::MATRIX);               
        eq.setOnlyStiffnessMatrix(false);
        Op Stiff(assembler);
        //BCRS_Matrix  matStiff(Stiff.template get<BCRS_Matrix>());
        matStiff = Stiff.template get<BCRS_Matrix>();
        //std::cout << "Stiff = \n" << matStiff << std::endl;
        matrixAssemblyTime += (double)(timer.elapsed().wall)/1e9;
        std::cout << "computing time for stiffness matrix: " << boost::timer::format(stiffMatrixTimer.elapsed()) << "\n";
      }




      // pointwise evaluation of RHS  rUi = f(u)
      timer.start();
      boost::timer::cpu_timer reactionTimer;
      for (int ii=0; ii<n+1; ii++)
      {
        if (ii<n) tau1=pts[ii+1]-pts[ii]; 
        double tt = eq.time();
        eq.setTau(tau1);				
		CoefficientVectors rU(solution);
        int dofsV = boost::fusion::at_c<1>(u.data).space().degreesOfFreedom();
        //std::cout << "\n dofs(v) = " << dofsV << std::endl;
        std::vector<double> uVec(dofsV), vVec(dofsV), fVec(dofsV), gVec(dofsV);
        for (int j=0; j<dofsV; ++j) 
        {
          uVec[j] = (boost::fusion::at_c<0>(lastSol[ii].data)).coefficients()[j];
          vVec[j] = (boost::fusion::at_c<1>(lastSol[ii].data)).coefficients()[j];
        }
        for (int k=0; k<dofsV; ++k) 
        {
          fVec[k] = eq.f(uVec[k],vVec[k]);  
          gVec[k] = eq.g(uVec[k],vVec[k]);  
        }
        CoefficientVectors rhs(rU);
        for (int k=0; k<dofsV; ++k) boost::fusion::at_c<0>(rhs.data)[k] = fVec[k];
        for (int k=0; k<dofsV; ++k) boost::fusion::at_c<1>(rhs.data)[k] = gVec[k];
        rUiMat[ii]=rhs; 

        eq.time(eq.time()+tau1);
      }
      std::cout << "computing time for reaction terms: " << boost::timer::format(reactionTimer.elapsed()) << "\n";
      rhsAssemblyTime += (double)(timer.elapsed().wall)/1e9;



      std::vector<BCRS_Matrix> allMatF_u(n+1,matStiff);

      // pointwise evaluation of diagonals in Jacobian matrices      
      timer.start();
      boost::timer::cpu_timer jacobianTimer;
      for (int ii=0; ii<n+1; ii++)
      {
        if (ii<n) tau1=pts[ii+1]-pts[ii]; 
        double tt = eq.time();
        eq.setTau(tau1);				
		CoefficientVectors rU(solution);
        int dofsV = boost::fusion::at_c<1>(u.data).space().degreesOfFreedom();
        std::vector<double> uVec(dofsV), vVec(dofsV), fuVec(dofsV), gvVec(dofsV);
        for (int j=0; j<dofsV; ++j) 
        {
          uVec[j] = (boost::fusion::at_c<0>(lastSol[ii].data)).coefficients()[j];
          vVec[j] = (boost::fusion::at_c<1>(lastSol[ii].data)).coefficients()[j];
        }
        for (int k=0; k<dofsV; ++k) 
        {
          fuVec[k] = eq.fu(uVec[k],vVec[k]);  
          gvVec[k] = eq.gv(uVec[k],vVec[k]);  
        }
        CoefficientVectors rhs(rU);
        for (int k=0; k<dofsV; ++k) boost::fusion::at_c<0>(rhs.data)[k] = fuVec[k];
        for (int k=0; k<dofsV; ++k) boost::fusion::at_c<1>(rhs.data)[k] = gvVec[k];
        rJiMat[ii]=rhs; 
        
        BCRS_Matrix  matF_u(matM);
        matF_u *= 0;
        for (int j=0; j<dofsV; ++j) matF_u[j][j]             = matM[j][j] * fuVec[j];
        for (int j=0; j<dofsV; ++j) matF_u[j+dofsV][j+dofsV] = matM[j+dofsV][j+dofsV] * gvVec[j];
        allMatF_u[ii] = matF_u;  

        eq.time(eq.time()+tau1);
      }
      std::cout << "computing time for diagonals of Jacobian matrices: " << boost::timer::format(jacobianTimer.elapsed()) << "\n";
      matrixAssemblyTime += (double)(timer.elapsed().wall)/1e9;
      



      // Reaction matrix f_u for u=Sol[ii-1]
      /*
      boost::timer::cpu_timer jacobianTimer;
      for (int ii=0; ii<n+1; ii++)
      {
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(false);
        eq.setOnlyMassMatrix(false);
        eq.setOnlyStiffnessMatrix(false);
        eq.setOnlyF_uMatrix(true);
        assembler.assemble(Linearization(eq,lastSol[ii],lastSol[ii],dx),Assembler::MATRIX);               
        eq.setOnlyF_uMatrix(false);
        Op F_u(assembler);
        //BCRS_Matrix  matF_u(F_u.template get<BCRS_Matrix>());
        allMatF_u[ii] = F_u.template get<BCRS_Matrix>();  
        //std::cout << "ii = " << ii << ", F_u[ii] = \n" << allMatF_u[ii] << std::endl;

      }
      std::cout << "computing time for Jacobian matrices: " << boost::timer::format(jacobianTimer.elapsed()) << "\n";
      */

      const int blocksize = 1;
      typedef Dune::FieldVector<double,blocksize> BlockType;
      typedef typename Dune::BlockVector<BlockType> B_Vector;

      CoefficientVectors u(rU);
      u = lastSol[0];
      //size_t  dofs = ansatzVars.degreesOfFreedom(0,nvars);
      B_Vector uVec(dofs);
      std::vector<double> uStdVec(dofs);
      IstlInterfaceDetail::toVector(u, uStdVec);
      for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
 
 
 
      std::vector<B_Vector> rUi(n+1,uVec);
      for (int ii=0; ii<n+1; ii++) 
      {
        u = rUiMat[ii];
        IstlInterfaceDetail::toVector(u, uStdVec);
        for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
        rUi[ii] = uVec;
      }

      std::vector<B_Vector> solLast(n+1,uVec);
      for (int ii=0; ii<n+1; ii++) 
      {
        u = lastSol[ii];
        IstlInterfaceDetail::toVector(u, uStdVec);
        for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
        solLast[ii] = uVec;
      }
      
      B_Vector dxVec(dofs);
      dxVec = 0;

      timer.start();
      sdcIterationStepCardio(grid, matM, matStiff, rUi, allMatF_u, solLast, dxVec);
      solutionTime += (double)(timer.elapsed().wall)/1e9;
      
      

      eq.time(t);


      //currSol <--- solLast  :  after calling sdcIterationStep we have the corrected solution in solLast
      for (int ii=0; ii<n+1; ii++) 
      {
        for(int j=0;j<dofs;j++) uStdVec[j] = solLast[ii][j];
        IstlInterfaceDetail::fromVector(uStdVec, u);
        currSol[ii] = u;
      }
      lastSol=currSol;	
      
      eq.time(eq.time()+dt);

      
      for(int j=0;j<dofs;j++) uStdVec[j] = dxVec[j];
      IstlInterfaceDetail::fromVector(uStdVec, dx);
      
      // estimate spatial error
      // dx is the defect in the right Lobatto/Radau point
      //if (!tolX.empty() && kk==0) 
      //accurate = true;
      if (!tolX.empty() && !accurate) 
      {
        State tmp(dx);
        projectHierarchically(ansatzVars,tmp);
        tmp -= dx;

        // perform mesh adaptation
        //accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,verbosity);
        accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,1);
        if (!accurate) {
          nnz = assembler.nnz(0,neq,0,nvars,false);
          size = ansatzVars.degreesOfFreedom(0,nvars);
          rhs.resize(size);
          sol.resize(size);
          std::cout << "accurate = " << accurate << ",   dofs after mesh refinement: " << size << std::endl;
        }
      }
    }
    while ( ( kk<(kEnd-1) ) || !accurate );		
		
    u = currSol[n];	
    return u;
	
  }




  State const& stepMOL(State const& x, double const dt, std::vector<std::pair<double,double> > const& tolX)
  {		  		 
    int const n = 6, kEnd = 6;   
    //std::cout << "SDC method: n = " << n << ",  kEnd = " << kEnd << std::endl;

    //LobattoTimeGrid grid(n,eq.time(),eq.time()+dt);
    RadauTimeGrid grid(n,eq.time(),eq.time()+dt);
    std::cout << "SDC-Matrix method: n = " << n << ",  kEnd = " << kEnd << std::endl;

    Dune::DynamicVector<double> const& pts = grid.points();

    // S is provided an used in sdcIterationStep()
    //Dune::DynamicMatrix<double> const& S = grid.integrationMatrix();
     
    // Shat is computed in sdcIterationStep(...) called below
    //SDCTimeGrid::RealMatrix Shat;
    //eulerIntegrationMatrix(grid, Shat);
    //luIntegrationMatrix(grid, Shat);
    //std::cout << "Shat = \n" << Shat << std::endl;
	
    double tau1=0;
							
    int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
    int const neq = EvolutionEquation::TestVars::noOfVariables;
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    size_t  size = ansatzVars.degreesOfFreedom(0,nvars);
    
    std::cout << "nvars = " << nvars << std::endl;
    std::cout << "  neq = " <<   neq << std::endl;
    std::cout << "  dof = " <<  size << std::endl;

    typedef  AssembledGalerkinOperator<Assembler,0,neq,0,neq> Op;

    double const t = eq.time();
    eq.temporalEvaluationRange(t,t+dt);

    double const tau = dt;
    eq.setTau(tau);

    // mesh adaptation loop. Just once if i>0.
    gridManager.setVerbosity(verbosity);
    bool accurate = false;

    bool const iterative = true;       // false: use a direct solver,      true: use an iterative solver
    std::vector<double> rhs(size), sol(size);

    // Evaluate and factorize matrix B(t)-tau*J

    std::vector<State> currSol(n+1,x);
    std::vector<State> lastSol(n+1,x);
				
			
    //for (int kk=0; kk<kEnd; kk++)
    int kk=-1;
    do
    {
      kk++;
      std::cout << "sweep: k = " << kk << std::endl;
      
      eq.time(t);
      
      CoefficientVectors rU(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      
      //for (int k=0; k< boost::fusion::at_c<0>(rhs1[0].data).size(); k++)
      //std::cout << " size of rU(0) = " << boost::fusion::at_c<0>(rU.data).size() << std::endl;
      //std::cout << " size of rU(1) = " << boost::fusion::at_c<1>(rU.data).size() << std::endl;
      
      rU = 0;
      std::vector<CoefficientVectors> rUiMat(n+1,rU);
      std::vector<CoefficientVectors> integMatS(n,rU);
      std::vector<CoefficientVectors> dUiMat(n+1,rU);
      CoefficientVectors sHatDu(rU);
      CoefficientVectors currsol(rU);
      CoefficientVectors lastsol(rU);						

      CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      solution=0;
      dUiMat[0] = solution;

      for (int ii=0; ii<n+1; ii++)
      {
        if (ii<n) tau1=pts[ii+1]-pts[ii];
        //tau1=pts[ii]-pts[ii-1];
        //eq.time(eq.time()+tau1);
        double tt = eq.time();
        //eq.temporalEvaluationRange(tt,tt+tau1);
        eq.setTau(tau1);				

        dx *= 0;              
        eq.setSecondRHS(false);
        // righthand-side r of system  [B(t)-tau*J] u_{i+1} = r(u_{i}) 
        //std::cout << "init RHS: ii = " << ii << ", time = " << eq.time() << std::endl;
        assembler.assemble(Linearization(eq,lastSol[ii],lastSol[ii],dx),Assembler::RHS);  

		CoefficientVectors rU(assembler.rhs());
        CoefficientVectors rhs(rU);
        rhs = 0;
        rhs.axpy(1.0/tau1, rU);
        
        //std::cout << "righthand-side ii = " << ii << std::endl;
        //for (int k=0; k<size; ++k) std::cout << " rhs = " << boost::fusion::at_c<0>(rhs.data)[k] << std::endl;
        
        rUiMat[ii]=rhs; 
        
        //if (kk>0) eq.time(eq.time()+tau1);
        eq.time(eq.time()+tau1);
      }


      typedef typename Op::field_type field_type;
      typedef Dune::FieldMatrix<field_type,1,1> Field_Matrix;
      typedef Dune::BCRSMatrix<Field_Matrix> BCRS_Matrix;
      
      // mass matrix M 
      dx *= 0;               
      eq.setSecondRHS(false);
      eq.setOnlyJacobian(false);
      eq.setOnlyStiffnessMatrix(false);
      eq.setOnlyF_uMatrix(false);
      eq.setOnlyMassMatrix(true);
      assembler.assemble(Linearization(eq,lastSol[0],lastSol[0],dx),Assembler::MATRIX);               
      eq.setOnlyMassMatrix(false);
      Op M(assembler);      
      BCRS_Matrix  matM(M.template get<BCRS_Matrix>());

      // stiffness matrix A 
      dx *= 0;               
      eq.setSecondRHS(false);
      eq.setOnlyJacobian(false);
      eq.setOnlyMassMatrix(false);
      eq.setOnlyF_uMatrix(false);
      eq.setOnlyStiffnessMatrix(true);
      assembler.assemble(Linearization(eq,lastSol[0],lastSol[0],dx),Assembler::MATRIX);               
      eq.setOnlyStiffnessMatrix(false);
      Op Stiff(assembler);
      BCRS_Matrix  matStiff(Stiff.template get<BCRS_Matrix>());
      //std::cout << "Stiff = \n" << matStiff << std::endl;

      // Reaction matrix f_u for u=Sol[ii-1]
      std::vector<BCRS_Matrix> allMatF_u(n+1,matStiff);
      for (int ii=0; ii<n+1; ii++)
      {
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(false);
        eq.setOnlyMassMatrix(false);
        eq.setOnlyStiffnessMatrix(false);
        eq.setOnlyF_uMatrix(true);
        assembler.assemble(Linearization(eq,lastSol[ii],lastSol[ii],dx),Assembler::MATRIX);               
        eq.setOnlyF_uMatrix(false);
        Op F_u(assembler);
        //BCRS_Matrix  matF_u(F_u.template get<BCRS_Matrix>());
        allMatF_u[ii] = F_u.template get<BCRS_Matrix>();  
        //std::cout << "ii = " << ii << ", F_u[ii] = \n" << allMatF_u[ii] << std::endl;

      }


      const int blocksize = 1;
      typedef Dune::FieldVector<double,blocksize> BlockType;
      typedef typename Dune::BlockVector<BlockType> B_Vector;

      CoefficientVectors u(rU);
      u = lastSol[0];
      size_t  dofs = ansatzVars.degreesOfFreedom(0,nvars);
      B_Vector uVec(dofs);
      std::vector<double> uStdVec(dofs);
      IstlInterfaceDetail::toVector(u, uStdVec);
      for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
 
 
 
      std::vector<B_Vector> rUi(n+1,uVec);
      for (int ii=0; ii<n+1; ii++) 
      {
        u = rUiMat[ii];
        IstlInterfaceDetail::toVector(u, uStdVec);
        for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
        rUi[ii] = uVec;
      }

      std::vector<B_Vector> solLast(n+1,uVec);
      for (int ii=0; ii<n+1; ii++) 
      {
        u = lastSol[ii];
        IstlInterfaceDetail::toVector(u, uStdVec);
        for(int j=0;j<dofs;j++) uVec[j] = uStdVec[j];
        solLast[ii] = uVec;
      }
      
      B_Vector dxVec(dofs);
      dxVec = 0;

      sdcIterationStep(grid, matM, matStiff, rUi, allMatF_u, solLast, dxVec);
      
      
      // sum in the rhs
//      int matCalc=summation(grid,pts,rUiMat,integMatS);	
	    	


      eq.time(t);


      //currSol <--- solLast  :  after calling sdcIterationStep we have the corrected solution in solLast
      for (int ii=0; ii<n+1; ii++) 
      {
        for(int j=0;j<dofs;j++) uStdVec[j] = solLast[ii][j];
        IstlInterfaceDetail::fromVector(uStdVec, u);
        currSol[ii] = u;
      }
      lastSol=currSol;	
      
      eq.time(eq.time()+dt);

      
      for(int j=0;j<dofs;j++) uStdVec[j] = dxVec[j];
      IstlInterfaceDetail::fromVector(uStdVec, dx);
      
      // estimate spatial error
      // dx is the defect in the right Lobatto/Radau point
      //if (!tolX.empty() && kk==0) 
      //accurate = true;
      if (!tolX.empty() && !accurate) 
      {
        State tmp(dx);
        projectHierarchically(ansatzVars,tmp);
        tmp -= dx;

        // perform mesh adaptation
        //accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,verbosity);
        accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,1);
        if (!accurate) {
          nnz = assembler.nnz(0,neq,0,nvars,false);
          size = ansatzVars.degreesOfFreedom(0,nvars);
          rhs.resize(size);
          sol.resize(size);
          std::cout << "accurate = " << accurate << ",   dofs after mesh refinement: " << size << std::endl;
        }
      }
    }
    while ( ( kk<(kEnd-1) ) || !accurate );		
		
    u = currSol[n];	
    return u;
	
  }





  State const& step(State const& x, double const dt, std::vector<std::pair<double,double> > const& tolX)
  {		  		 
    int const n = 4, kEnd = 4;  
    //std::cout << "SDC method: n = " << n << ",  kEnd = " << kEnd << std::endl;

    //LobattoTimeGrid grid(n,eq.time(),eq.time()+dt);
    RadauTimeGrid grid(n,eq.time(),eq.time()+dt);
    std::cout << "SDC method: n = " << n << ",  kEnd = " << kEnd << std::endl;

    Dune::DynamicVector<double> const& pts = grid.points();
//     std::cout << "Lobatto/Radau - Points:\n";
//     for (int j=0; j<=n; j++) {
//       std::cout << pts[j] << std::endl;
//     }

    Dune::DynamicMatrix<double> const& S = grid.integrationMatrix();
//     std::cout << "Integration-matrix S:\n";
//     for (int i=0; i<n; i++)
//     {	
//       for (int j=0; j<=n; j++)
//       {
//           std::cout << "S[" << i << "][" << j << "] = " << S[i][j] << "   ";
//       }   
//       std::cout << std::endl; 
//      }

//     std::cout << "transposed Integration-matrix S:\n";
//     for (int j=0; j<=n; j++)
//     {	
//       for (int i=0; i<n; i++)
//       {
//           std::cout << "S'[" << j << "][" << i << "] = " << S[i][j] << "   ";
//       }   
//       std::cout << std::endl; 
//      }
     
    SDCTimeGrid::RealMatrix Shat;
    //eulerIntegrationMatrix(grid, Shat);
    luIntegrationMatrix(grid, Shat);
    // std::cout << "Shat = \n" << Shat << std::endl;
    
//     std::cout << "Shat-matrix:\n";
//     for (int i=0; i<n; i++)
//     {	
//       for (int j=0; j<=n; j++)
//       {
//           std::cout << "Shat[" << i << "][" << j << "] = " << Shat[i][j] << "   ";
//       }   
//       std::cout << std::endl; 
//      }
	
    double tau1=0;
							
    int const nvars = EvolutionEquation::AnsatzVars::noOfVariables;
    int const neq = EvolutionEquation::TestVars::noOfVariables;
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    size_t  size = ansatzVars.degreesOfFreedom(0,nvars);
    
    std::cout << "nvars = " << nvars << std::endl;
    std::cout << "  neq = " <<   neq << std::endl;
    std::cout << "  dof = " <<  size << std::endl;

    typedef  AssembledGalerkinOperator<Assembler,0,neq,0,neq> Op;

    double const t = eq.time();
    eq.temporalEvaluationRange(t,t+dt);

    double const tau = dt;
    eq.setTau(tau);

    // mesh adaptation loop. Just once if i>0.
    gridManager.setVerbosity(verbosity);
    bool accurate = false;

    bool const iterative = true;       // false: use a direct solver,      true: use an iterative solver
    std::vector<double> rhs(size), sol(size);

    // Evaluate and factorize matrix B(t)-tau*J

    std::vector<State> currSol(n+1,x);
    std::vector<State> lastSol(n+1,x);
				
			
    //for (int kk=0; kk<kEnd; kk++)
    int kk=-1;
    do
    {
      kk++;
      std::cout << "sweep: k = " << kk << std::endl;
      
      eq.time(t);
      
      CoefficientVectors rU(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      
      //for (int k=0; k< boost::fusion::at_c<0>(rhs1[0].data).size(); k++)
      //std::cout << " size of rU(0) = " << boost::fusion::at_c<0>(rU.data).size() << std::endl;
      //std::cout << " size of rU(1) = " << boost::fusion::at_c<1>(rU.data).size() << std::endl;
      
      rU = 0;
      std::vector<CoefficientVectors> rUiMat(n+1,rU);
      std::vector<CoefficientVectors> integMatS(n,rU);
      std::vector<CoefficientVectors> dUiMat(n+1,rU);
      CoefficientVectors sHatDu(rU);
      CoefficientVectors currsol(rU);
      CoefficientVectors lastsol(rU);						

      CoefficientVectors solution(EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::init(ansatzVars));
      solution=0;
      dUiMat[0] = solution;

      for (int ii=0; ii<n+1; ii++)
      {
        if (ii<n) tau1=pts[ii+1]-pts[ii];
        //tau1=pts[ii]-pts[ii-1];
        //eq.time(eq.time()+tau1);
        double tt = eq.time();
        //eq.temporalEvaluationRange(tt,tt+tau1);
        eq.setTau(tau1);				

        dx *= 0;              
        eq.setSecondRHS(false);
        // righthand-side r of system  [B(t)-tau*J] u_{i+1} = r(u_{i}) 
        //std::cout << "init RHS: ii = " << ii << ", time = " << eq.time() << std::endl;
        assembler.assemble(Linearization(eq,lastSol[ii],lastSol[ii],dx),Assembler::RHS);  

		CoefficientVectors rU(assembler.rhs());
        CoefficientVectors rhs(rU);
        rhs = 0;
        rhs.axpy(1.0/tau1, rU);
        
        //std::cout << "righthand-side ii = " << ii << std::endl;
        //for (int k=0; k<size; ++k) std::cout << " rhs = " << boost::fusion::at_c<0>(rhs.data)[k] << std::endl;
        
        rUiMat[ii]=rhs; 
        
        //if (kk>0) eq.time(eq.time()+tau1);
        eq.time(eq.time()+tau1);
      }


      // sum in the rhs
      int matCalc=summation(grid,pts,rUiMat,integMatS);	
	    	

      eq.time(t);


      for (int ii=1; ii<n+1; ii++)
      {
        tau1=pts[ii]-pts[ii-1];
        eq.time(eq.time()+tau1);
        double tt = eq.time();
        //eq.temporalEvaluationRange(tt,tt+tau1);

        //std::cout << "ii = " << ii << ": tau = " << tau1 << ", Shat[ii-1][ii] = " << Shat[ii-1][ii] << std::endl;
        tau1 = Shat[ii-1][ii];           // neu !!!
        eq.setTau(tau1);
        
			
        // Jacobian r'(u) 
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(true);
        assembler.assemble(Linearization(eq,lastSol[ii-1],lastSol[ii-1],dx),Assembler::MATRIX);               
        Op Jac(assembler);
        

        CoefficientVectors du(dUiMat[ii-1]);
        dUiMat[ii-1] = 0;
        Jac.applyscaleadd(+1,du,dUiMat[ii-1]);		// d.h. dUiMat[ii-1] += Jac*du
 
        // matrix B(t)-tau*J 
        dx *= 0;               
        eq.setSecondRHS(false);
        eq.setOnlyJacobian(false);
        assembler.assemble(Linearization(eq,lastSol[ii],lastSol[ii],dx),Assembler::MATRIX);               
        Op A(assembler);
					
        // get rhs for linear system
			
        State nn(currSol[ii-1]);                 
        currsol = nn;                  
        lastsol=lastSol[ii];
        eq.setR(currsol,lastsol);  
        
        dx *= 0;              
        eq.setSecondRHS(true);
        assembler.assemble(Linearization(eq,currSol[ii-1],currSol[ii-1],dx),Assembler::RHS);  
        
        CoefficientVectors rU(assembler.rhs());	 
        CoefficientVectors newRHS(rU);
        newRHS = 0;
        newRHS.axpy(1.0/tau1, rU);
    
        newRHS += integMatS[ii-1];
        
        // start computation Shat * dU
        sHatDu = 0;
        // j basis function	 
        for (int j=0; j<ii; j++)
        {
            sHatDu.axpy(Shat[ii-1][j],dUiMat[j]);
        }       
        // end computation Shat * dU

        newRHS+= sHatDu;

						
        // solve linear system				   
        typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,neq>::type LinearSpaceX;

        if(iterative) 
        {
          //solution=0;
          //JacobiPreconditioner<Op> jacobi(A,1.0);
          int fill_lev = 0;
          ILUKPreconditioner<Op> iluk(A,fill_lev,verbosity);
          int iteSteps = 5000;
          double iteEps = 1.0e-12;
          //Dune::BiCGSTABSolver<LinearSpaceX> cg(A,jacobi,iteEps,iteSteps,1);
          Dune::BiCGSTABSolver<LinearSpaceX> cg(A,iluk,iteEps,iteSteps,1);
          Dune::InverseOperatorResult res;
          cg.apply(solution,newRHS,res);
          dx=solution;
        }
        else
        {
            MatrixAsTriplet<double> triplet(nnz);
            triplet = A.template get<MatrixAsTriplet<double> >();

            Factorization<double> *matrix = 0;
            std::cout << "direct type = " << directType << std::endl;
            switch (directType) {
            case DirectType::UMFPACK:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::UMFPACK3264:
              matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::MUMPS:
              matrix = new MUMPSFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            case DirectType::SUPERLU:
              matrix = new SUPERLUFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
              break;
            default:
              throw -321;
              break;
            }
            matrix->setVerbose(verbosity);
            //factorizationTime += (double)(timer.elapsed().wall)/1e9;

            // First right hand side (j=0) has been assembled together with matrix.
            //timer.start();
            IstlInterfaceDetail::toVector(newRHS, rhs);
            //newRHS.write(begin(rhs));
            //for (int k=0; k<rhs.size(); ++k) std::cout << "RHS = " << rhs[k] << std::endl;
            matrix->solve(rhs,sol);
            delete matrix;
            
            for (int k=0; k<sol.size(); ++k) assert(std::isfinite(sol[k]));
            dx.read(sol.begin()); 
            //for (int k=0; k<sol.size(); ++k) at_c<0>(solution.data)[k] = sol[k];
            solution.read(sol.begin());    //   ???
        }
        
        dUiMat[ii] = solution;

        L2Norm l2;
        //State defect(nn);
        double nrm2 = l2( boost::fusion::at_c<0>(dx.data) );
        //double nrm2 = l2( dx.data );
        std::cout << " sdc: k = " << kk << ",  ||dx|| = " << nrm2 << std::endl;

        
        currSol[ii]+=solution;                           
//        eq.time(eq.time()+tau1);
        //std::cout << " time in ii-loop  at end= " << eq.time() << std::endl;
      }
      
      lastSol=currSol;	
      
      // estimate spatial error
      // dx is the defect in the right Lobatto/Radau point
      //if (!tolX.empty() && kk==0) 
      if (!tolX.empty() && !accurate) 
      {
        State tmp(dx);
        projectHierarchically(ansatzVars,tmp);
        tmp -= dx;

        // perform mesh adaptation
        //accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,verbosity);
        accurate = embeddedErrorEstimator(ansatzVars,tmp,x,eq.scaling(),tolX,gridManager,1);
        if (!accurate) {
          nnz = assembler.nnz(0,neq,0,nvars,false);
          size = ansatzVars.degreesOfFreedom(0,nvars);
          rhs.resize(size);
          sol.resize(size);
          std::cout << "accurate = " << accurate << ",   dofs after mesh refinement: " << size << std::endl;
        }
      }
    }
    while ( ( kk<(kEnd-1) ) || !accurate );		
		
    u = currSol[n];	
    return u;
	
  }


									
  int summation(LobattoTimeGrid const& grid,
                Dune::DynamicVector<double>  pts, 
                std::vector<CoefficientVectors> const& rhs1,
                std::vector<CoefficientVectors>& integMatS )
  {
    Dune::DynamicMatrix<double> const& A = grid.integrationMatrix();
    
    int n=pts.size()-1;

// 	double sumWeights = 0;
//     for (int i=0; i<n; i++)
//     {	
//       for (int j=0; j<=n; j++)
//       {
//           sumWeights += A[i][j];
//           std::cout << "S[" << i << "][" << j << "] = " << A[i][j] << "   ";
//       }   
//       std::cout << std::endl; 
//      }
//      std::cout << "sum of weights = " << sumWeights << std::endl;
	
    //  Ensure that integMatS is initialized by Zero         
    for (int i=0; i<n; i++) integMatS[i] = 0;
    	
	// j basis function	 
    for (int j=0; j<=n; j++)
    {
      // i for integration interval [t_i,t_i+1]
      for (int i=0; i<n; i++)
      {	
        // k for position x_k
        integMatS[i].axpy(A[i][j],rhs1[j]);
      }                    
    }
	return true;
  }



  int summation(RadauTimeGrid const& grid,
                Dune::DynamicVector<double>  pts, 
                std::vector<CoefficientVectors> const& rhs1,
                std::vector<CoefficientVectors>& integMatS )
  {
    Dune::DynamicMatrix<double> const& A = grid.integrationMatrix();
	
    int n=pts.size()-1;
    
// 	double sumWeights = 0;
//     for (int i=0; i<n; i++)
//     {	
//       for (int j=0; j<=n; j++)
//       {
//           sumWeights += A[i][j];
//           std::cout << "S[" << i << "][" << j << "] = " << A[i][j] << "   ";
//       }   
//       std::cout << std::endl; 
//      }
//      std::cout << "sum of weights = " << sumWeights << std::endl;

	
    //  Ensure that integMatS is initialized by Zero         
    for (int i=0; i<n; i++) integMatS[i] = 0;
    	
	// j basis function	 
    for (int j=0; j<=n; j++)
    {
      // i for integration interval [t_i,t_i+1]
      for (int i=0; i<n; i++)
      {	
        // k for position x_k
          integMatS[i].axpy(A[i][j],rhs1[j]);
       }                    
     }
	return true;
  }


    /**
     * Estimates the time discretization error of the previously
     * computed step by taking the difference between the diagonal and
     * subdiagonal extrapolation values of maximal order. This requires
     * that order>1 has been given for the last step.
     */
  std::vector<std::pair<double,double> > estimateError(State const& x,int i, int j) const
  {
    std::vector<std::pair<double,double> > e(ansatzVars.noOfVariables);
    return e;
  }

  template <class OutStream>
  void reportTime(OutStream& out) const 
  {
      out << "sdc time: " << matrixAssemblyTime << "s matrix assembly\n"
          << "            " << rhsAssemblyTime << "s rhs assembly\n"
          << "            " << precAssemblyTime << "s preconditioner assembly\n"
          << "            " << factorizationTime << "s factorization\n"
          << "            " << initSolutionTime << "s init solution\n"
          << "            " << solutionTime << "s solution\n"
          << "            " << estimateTime << "s estimate\n";
  }

  void advanceTime(double dt) { eq.time(eq.time()+dt); }


  private:
    GridManager<typename EvolutionEquation::AnsatzVars::Grid>& gridManager;
    typename EvolutionEquation::AnsatzVars const& ansatzVars;
    SemiImplicitEulerStep<EvolutionEquation>      eq;
    Assembler                                     assembler;
    

  public:
    State u;
    State dx;    
    double rhsAssemblyTime, matrixAssemblyTime, precAssemblyTime, factorizationTime, solutionTime;
    double elapsedTimeSinceReset, initSolutionTime, estimateTime;
    DirectType directType;
    PrecondType precondType;
    int verbosity;
  };
} // namespace Kaskade
#endif
