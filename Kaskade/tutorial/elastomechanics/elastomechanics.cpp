/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef Cygwin
// availability of std::to_string
#define _GLIBCXX_USE_C99 1
#endif

#include <iostream>

// #include <boost/timer/timer.hpp>


#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include <dune/istl/operators.hh>

#include "utilities/enums.hh"
#include "fem/assemble.hh"
#include "fem/spaces.hh"
#include "io/vtk.hh"
#include "io/matlab.hh"
#include "linalg/direct.hh"
#include "linalg/apcg.hh"
#include "linalg/threadedMatrix.hh"
#include "linalg/triplet.hh"
#include "mg/multigrid.hh"
#include "mg/additiveMultigrid.hh"
#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"
#include "utilities/memory.hh"

using namespace Kaskade;
#include "elastomechanics.hh"

int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start elastomechanics tutorial program" << std::endl;
  
  int coarseGridSize, maxit, order, refinements, solver, verbose;
  bool additive, direct, vtk;
  double atol;
  std::string material;
  
  if (getKaskadeOptions(argc,argv,Options
    ("refinements",      refinements,     3,          "number of uniform grid refinements")
    ("coarse",           coarseGridSize,  1,          "number of coarse grid elements along each cube edge")
    ("order",            order,           1,          "finite element ansatz order")
    ("material",         material,        "steel",    "type of material")
    ("direct",           direct,          true,       "if true, use a direct solver")
    ("solver",           solver,          2,          "0=UMFPACK, 1=PARDISO 2=MUMPS 3=SUPERLU 4=UMFPACK32/64 5=UMFPACK64")
    ("additive",         additive,        false,      "use additive multigrid")
    ("verbosity",        verbose,         0,          "amount of reported details")
    ("vtk",              vtk,             false,      "write solution to VTK file")
    ("atol",             atol,            1e-8,       "absolute energy error tolerance for iterative solver")
    ("maxit",            maxit,           100,        "maximum number of iterations")))
    return 1;

  boost::timer::cpu_timer totalTimer;

  std::cout << "refinements of original mesh : " << refinements << std::endl;
  std::cout << "discretization order         : " << order << std::endl;


  constexpr int DIM = 3;
  using Grid = Dune::UGGrid<DIM>;
  using Spaces = boost::fusion::vector<H1Space<Grid> const*>;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,DIM,0> >;
  using VarSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = ElasticityFunctional<VarSetDesc>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using CoefficientVectors = VarSetDesc::CoefficientVectorRepresentation<0,1>::type;

  Dune::FieldVector<double,DIM> c0(0.0), dc(1.0);
  GridManager<Grid> gridManager( createCuboid<Grid>(c0,dc,1.0/coarseGridSize,true) );
  gridManager.globalRefine(refinements);
  gridManager.enforceConcurrentReads(true);
  std::cout << "grid creation & refinement time: " << totalTimer.format();

	
  // construction of finite element space for the scalar solution T.
  H1Space<Grid> h1Space(gridManager,gridManager.grid().leafGridView(),order);
	
  Spaces spaces(&h1Space);
	
  std::string varNames[1] = { "u" };
	
  VarSetDesc varSetDesc(spaces,varNames);

  // Create the variational functional.
  Functional F(ElasticModulus::material(material));
	
  // construct Galerkin representation
  Assembler assembler(gridManager,spaces);
  VarSetDesc::VariableSet x(varSetDesc);
  VarSetDesc::VariableSet dx(varSetDesc);

	
	
  boost::timer::cpu_timer assembTimer;
  assembler.assemble(linearization(F,x));
  std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed()) << "\n";
  
  CoefficientVectors rhs(assembler.rhs());
  CoefficientVectors solution(VarSetDesc::CoefficientVectorRepresentation<0,1>::init(spaces));
  boost::timer::cpu_timer solveTimer;
  if (direct)
  {
    DirectType directType = static_cast<DirectType>(solver);
    AssembledGalerkinOperator<Assembler> A(assembler, directType == DirectType::MUMPS || directType == DirectType::PARDISO);
    directInverseOperator(A,directType,MatrixProperties::POSITIVEDEFINITE).applyscaleadd(-1.0,rhs,solution);
    x.data = solution.data;
  }
  else
  {
    using X = Dune::BlockVector<Dune::FieldVector<double,DIM>>;
    DefaultDualPairing<X,X> dp;
    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,DIM,DIM>>;
    using LinOp = Dune::MatrixAdapter<Matrix,X,X>;
    Matrix Amat(assembler.get<0,0>(),true);
    LinOp A(Amat);
    SymmetricLinearOperatorWrapper<X,X> sa(A,dp);
    PCGEnergyErrorTerminationCriterion<double> term(atol,maxit);
    
    
    Dune::InverseOperatorResult res;
    X xi(component<0>(rhs).N());

    std::unique_ptr<SymmetricPreconditioner<X,X>> mg;
    if (additive)
    {
      if (order==1)
        mg = moveUnique(makeBPX(Amat,gridManager));
      else
      {
        H1Space<Grid> p1Space(gridManager,gridManager.grid().leafGridView(),1);
        mg = moveUnique(makePBPX(Amat,h1Space,p1Space));
      }
    }
    else
      mg = makeMultigrid(Amat,h1Space);


    Pcg<X,X> pcg(sa,*mg,term,verbose);
    pcg.apply(xi,component<0>(rhs),res);
    xi *= -1;
    component<0>(x) = xi;
  }
  std::cout << "computing time for solve: " << solveTimer.format();

  // output of solution in VTK format for visualization,
  // the data are written as ascii stream into file elasto.vtu,
  // possible is also binary
  if (vtk)
    writeVTKFile(x,"elasto",IoOptions().setOrder(order).setPrecision(7));

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End elastomechanics tutorial program" << std::endl;
  
  return 0;
}
