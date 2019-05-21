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

#define FUSION_MAX_VECTOR_SIZE 25  // increase maximal vector size for boost/fusion vectors
#include <iostream>
#include <algorithm>  // min, max

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/istl/solvers.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "fem/norms.hh"
#include "linalg/iluprecond.hh"
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

using namespace Kaskade;
#include "sst.hh"

struct InitialValue 
{
  using Scalar = double;
  static constexpr int components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  InitialValue(int c): component(c) {}
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
          Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> const& localCoordinate) const 
  {
  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
  if (component==0) 
    return 1.306028e6;
  else if (component==1) 
    return 1.076508e12;
  else if (component==2) 
    return 6.457715e10;
  else if (component==3) 
    return 3.542285e10;
  else
    assert("wrong index!\n"==0);
  return 0;
  
  }

private:
  int component;
};

int main(int argc, char *argv[])
{
  using Scalar = double;
  using namespace boost::fusion;

  std::cout << "Start sst transfer tutorial program with ordinary Newton iteration" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int  direct,
     refinements = getParameter(pt, "refinement", 5),
     order =  getParameter(pt, "order", 2),
     verbosity   = getParameter(pt, "verbosity", 0);
  Scalar tol = getParameter(pt, "tolerance", 1.0e-10);
  //   IterateType iterateType = IterateType::CG;
  //   PrecondType precondType = PrecondType::NONE;
  Scalar fsign = static_cast<Scalar>(getParameter(pt, "sign", 0.0));
  std::string empty;

  std::cout << "refinements of original mesh   : " << refinements << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "tolerance for Newton iteration : " << tol << std::endl;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);

  //   s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
  //   iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  //   s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  //   precondType = static_cast<PrecondType>(getParameter(pt, s, 0));

  constexpr int dim=2;    
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<Scalar,LeafView> >;
// using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<Scalar,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<1> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<2> >,
                               Variable<SpaceIndex<0>,Components<1>,VariableId<3> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = SSTFunctional<Scalar,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = SSTFunctional<Scalar,VariableSet>::TestVars::noOfVariables;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;

  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);
  std::cout << std::endl << "Grid: " << gridManager.grid().size(0) << " triangles, " << std::endl;
  std::cout << "      " << gridManager.grid().size(1) << " edges, " << std::endl;
  std::cout << "      " << gridManager.grid().size(2) << " points" << std::endl;

  // construction of finite element space for the scalar solution T
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),
               order);
  Spaces spaces(&temperatureSpace);
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
  std::string varNames[4] = { "u0", "u1", "u2", "u3" };
  VariableSet variableSet(spaces,varNames);

  Functional F;
  Assembler assembler(gridManager,spaces);
  VariableSet::VariableSet x(variableSet);
  VariableSet::VariableSet newtonCorr(variableSet);
  VariableSet::VariableSet tmp(variableSet);

  constexpr int nvars = SSTFunctional<Scalar,VariableSet>::AnsatzVars::noOfVariables;
  size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
  size_t  size = variableSet.degreesOfFreedom(0,nvars);
  std::cout << "nvars=" << nvars << ", neq=" << neq << ", size=" << size << ", nnz=" << nnz << "\n";
    
  std::vector<Scalar> xdata(size), scal(size), dxdata(size);
  
  int k=0;
  L2Norm l2Norm;
  Scalar norm_dx, norm_rhs;
  std::vector<Scalar> norm_dx_comp(4), norm_rhs_comp(4);
  
  F.scaleInitialValue<0>(InitialValue(0),x);
  F.scaleInitialValue<1>(InitialValue(1),x);
  F.scaleInitialValue<2>(InitialValue(2),x);
  F.scaleInitialValue<3>(InitialValue(3),x);
  x.write(xdata.begin());
  for (int i=0;i<size;i++) scal[i]=std::max(std::abs(xdata[i]),1.0);
  
  CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,neq>::init(spaces));
  
  writeVTKFile(x,"graph/sst_start",IoOptions().setOrder(std::min(order,2)).setPrecision(7));
  gridManager.enforceConcurrentReads(false);
  Dune::InverseOperatorResult res;
  
  std::cout << std::endl << "Newton iteration starts:" << std::endl <<
            "iter  scaled ||corr||            ||F||    itsol: eps     result  steps      rate        time"
            << std::endl;

// begin of ordinary Newton iteration loop
  do {
    boost::timer::cpu_timer assembTimer;
    assembler.assemble(linearization(F,x));
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    AssembledGalerkinOperator<Assembler,0,neq,0,neq> A(assembler);
    CoefficientVectors rhs(assembler.rhs());
  
    tmp.data = rhs.data;
    norm_rhs_comp[0]=l2Norm(boost::fusion::at_c<0>(tmp.data));
    norm_rhs_comp[1]=l2Norm(boost::fusion::at_c<1>(tmp.data));
    norm_rhs_comp[2]=l2Norm(boost::fusion::at_c<2>(tmp.data));
    norm_rhs_comp[3]=l2Norm(boost::fusion::at_c<3>(tmp.data));
    norm_rhs = sqrt(norm_rhs_comp[0]*norm_rhs_comp[0]+norm_rhs_comp[1]*norm_rhs_comp[1]+
                    norm_rhs_comp[2]*norm_rhs_comp[2]+norm_rhs_comp[3]*norm_rhs_comp[3]);
  
    boost::timer::cpu_timer iteTimer;
    int iteSteps = getParameter(pt, "solver.iteMax", 1000);
    Scalar iteEps = getParameter(pt, "solver.iteEps", 1.0e-12);
    int fill_lev = getParameter(pt, "solver.ILUK.fill_lev", 0);
    ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,fill_lev,verbosity);
    Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity);
    solution = 0;
    cg.apply(solution,rhs,res);
    newtonCorr.data = solution.data;
    newtonCorr *= -1;
    // add Newton-correction
    x += newtonCorr;
  
    newtonCorr.write(dxdata.begin());
    for (int i=0;i<size;i++) dxdata[i] /= scal[i];
    tmp.read(dxdata.begin());
    norm_dx_comp[0]=l2Norm(boost::fusion::at_c<0>(tmp.data));
    norm_dx_comp[1]=l2Norm(boost::fusion::at_c<1>(tmp.data));
    norm_dx_comp[2]=l2Norm(boost::fusion::at_c<2>(tmp.data));
    norm_dx_comp[3]=l2Norm(boost::fusion::at_c<3>(tmp.data));
    norm_dx = sqrt(norm_dx_comp[0]*norm_dx_comp[0]+norm_dx_comp[1]*norm_dx_comp[1]+
                   norm_dx_comp[2]*norm_dx_comp[2]+norm_dx_comp[3]*norm_dx_comp[3]);
    
    x.write(xdata.begin());
    for (int i=0;i<size;i++) scal[i]=std::max(std::abs(xdata[i]),1.0);
  
    std::cout << std::setw(4) << k+1 << "  " 
              << std::setw(15) << std::setprecision(5) << std::scientific << norm_dx << "  "  
              << std::setw(15) << norm_rhs << "  " << std::setprecision(6);
    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setw(12) << iteEps << "  " 
              << (res.converged?"converged":"failed  ") << "  "
              << std::setw(5) << res.iterations << "  " 
              << std::setw(8) << res.conv_rate << "  "
              << std::setw(9) << (Scalar)(iteTimer.elapsed().user)/1e9 << "s " << std::endl;
  
    // output of solution in VTK format for visualization,
    // the data are written as ascii stream into file temperature.vtu,
    // possible is also binary
    std::ostringstream fname;
    fname << "graph/sst_";
    fname.width(2);
    fname.fill('0');
    fname.setf(std::ios_base::right,std::ios_base::adjustfield);
    fname << k;
    fname.flush();
    writeVTKFile(x,fname.str(),IoOptions().setOrder(std::min(order,2)).setPrecision(7));
    
    // output of solution for Amira visualization,
    // the data are written in binary format into file temperature.am,
    // possible is also ascii
    //IoOptions options;
    //options.outputType = IoOptions::ascii;
    //LeafView leafGridView = gridManager.grid().leafGridView();
    //writeAMIRAFile(leafGridView,x,"sst",options);
  k ++;
  }
  while ( norm_dx > tol );
 // end of ordinary Newton iteration loop
 
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End sst transfer tutorial program" << std::endl;
}
