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

#include <iostream>
#include <utility> // std::move

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"

#include "dune/grid/onedgrid.hh"
#include "dune/grid/onedgrid/onedgridfactory.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "fem/istlinterface.hh"
#include "linalg/umfpack_solve.hh"
#include "linalg/triplet.hh"
#include "io/gnuplot.hh"
#include "utilities/kaskopt.hh"

using namespace Kaskade;
#include "atp.hh"

int main(int argc, char *argv[])
  {
  using namespace boost::fusion;

  std::cout << "Start atp transfer tutorial program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosity = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosity, dump);

  int  direct = 1,
       refinements = getParameter(pt, "refinement", 5),
       order =  getParameter(pt, "order", 2);
  double fsign = static_cast<double>(getParameter(pt, "sign", -1.0));
  int graphicalOutput= getParameter(pt, "graphicalOutput", 1);

  std::cout << "refinements of original mesh : " << refinements << std::endl;
  std::cout << "discretization order         : " << order << std::endl;

  std::cout << "fsign = " << fsign << "\n";


  //   one-dimensional space: dim=1
  int const dim=1;    
  using Grid = Dune::OneDGrid;
  using LeafView = Grid::LeafGridView;
  // construction of finite element space for the scalar solution T
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  // using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,1,0> >;
  using VarSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = ATPFunctional<double,VarSetDesc>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  using CoefficientVectors = VarSetDesc::CoefficientVectorRepresentation<0,1>::type;
  constexpr int nvars = ATPFunctional<double,VarSetDesc>::AnsatzVars::noOfVariables;
  constexpr int neq = ATPFunctional<double,VarSetDesc>::TestVars::noOfVariables;

  Dune::GridFactory<Grid> factory;

  // point (in case of dimension>1: vertex) coordinates v[0]
  Dune::FieldVector<double,dim> v; 
  v[0]=-3; factory.insertVertex(v);
  v[0]=0; factory.insertVertex(v);
  v[0]=3; factory.insertVertex(v);
  std::vector<unsigned int> vid(2);
  Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
  vid[0]=0; vid[1]=1; factory.insertElement(gt,vid);
  vid[0]=1; vid[1]=2; factory.insertElement(gt,vid);
  // interval defined by 2 point indices
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  // the coarse grid will be refined refinements times
  grid->globalRefine(refinements);
  // some information on the refined mesh
  std::cout << "Grid: " << grid->size(1) << " points " << std::endl;
  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));   

  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),
               order);
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  VarSetDesc varSetDesc(spaces,{ "T" });
  Functional F(fsign);
  Assembler assembler(gridManager,spaces);
  VarSetDesc::VariableSet x(varSetDesc);
  VarSetDesc::VariableSet dx(varSetDesc);

  // set nnz to the number of structural nonzero elements of the matrix to be assembled below
  size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
  size_t  size = varSetDesc.degreesOfFreedom(0,nvars);
  AssembledGalerkinOperator<Assembler,0,neq,0,nvars> A(assembler);
  MatrixAsTriplet<double> triplet(nnz);
      
  std::vector<double> rhs(size), sol(size);
  
  int k=0;
  L2Norm l2Norm;
  double norm_dx, norm_rhs;
  x=0;
  
  std::cout << std::endl << "Newton iteration starts:" << std::endl <<
            "iter   ||correction||            ||F||  assemble time  linsolve time"
            << std::endl;

// begin of ordinary Newton iteration loop
  do 
  {
    boost::timer::cpu_timer assembTimer;
    assembler.assemble(linearization(F,x));
	double assembleTime = (double)(assembTimer.elapsed().user)/1e9;
    triplet = A.get<MatrixAsTriplet<double> >();
    //     for (k=0; k< nnz; k++)
    //       {
    //         printf("%3d %3d %e\n", triplet.ridx[k], triplet.cidx[k], triplet.data[k]);
    //       }
    boost::timer::cpu_timer directTimer;
    Factorization<double> *matrix = 0;
    matrix = new UMFFactorization<double>(size,0,triplet.ridx,triplet.cidx,triplet.data);
    assembler.toSequence(0,neq,rhs.begin());
    for (int l=0; l<rhs.size(); ++l) assert(std::isfinite(rhs[l]));
    matrix->solve(rhs,sol);
    double solveTime = (double)(directTimer.elapsed().user)/1e9;
    delete matrix;
    for (int l=0; l<sol.size(); ++l) assert(std::isfinite(sol[l]));
    dx.read(sol.begin());
    dx *= -1;
    x += dx;
    norm_dx=l2Norm(component<0>(dx));
    VarSetDesc::VariableSet tmp(dx);
    tmp.read(rhs.begin());
    norm_rhs=l2Norm(component<0>(tmp));
    std::cout << std::setw(4) << k+1 << "  " 
              << std::setw(15) << std::setprecision(5) << std::scientific << norm_dx << "  "  
              << std::setw(15) << norm_rhs << "  " << std::setprecision(3);
    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::fixed << std::setw(6) << "       " << std::setprecision(3)
              << assembleTime << "s        " << std::setw(6) << solveTime << "s " << std::endl;
    // output of solution in VTK format for visualization,
    // the data are written as ascii stream into file temperature.vtu,
    // possible is also binary
    // char fname[6];
    // IoOptions options;
    // options.outputType = IoOptions::ascii;
    // LeafView leafView = gridManager.grid().leafGridView();
    // sprintf(fname,"atp_%#02d",k);
    // writeVTKFile(x,fname,options,order);
    if ( graphicalOutput!=0 )
    {
      IoOptions gnuplotOptions{};
      std::string empty;
      std::string s = "names.gnuplotinfo." + getParameter(pt, "gnuplotinfo", empty);
      gnuplotOptions.info = static_cast<IoOptions::Info>(getParameter(pt,s,0));
      writeGnuplotFile(x,"function",gnuplotOptions);
    };
    // output of solution for Amira visualization,
    // the data are written in binary format into file temperature.am,
    // possible is also ascii
    //  IoOptions options;
    //  options.outputType = IoOptions::ascii;
    //  LeafView leafView = gridManager.grid().leafGridView();
    //  writeAMIRAFile(leafView,varSetDesc,x,"atp",options);
    //StopSnippet6
    k ++;
  }
  while ( norm_dx > 1.0e-5 );

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End atp transfer tutorial program" << std::endl;
}
