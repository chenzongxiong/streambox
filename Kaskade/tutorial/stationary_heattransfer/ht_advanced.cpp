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

#include <iostream>

#ifndef GRIDTYPE
#define GRIDTYPE 1 // UGGrid
#endif

#if ( GRIDTYPE < 1 ) || ( GRIDTYPE > 4 )
#error "GRIDTYPE must be 1 (UGGrid), 2 (ALUSImplexGrid), 3 (ALUConformGrid) or 4 (AlbertaGrid)"
#endif

#if GRIDTYPE==1
#define CHARGRIDTYPE "UGGrid"
#elif  GRIDTYPE==2
#define CHARGRIDTYPE "ALUSimplexGrid"
#elif  GRIDTYPE==3
#define CHARGRIDTYPE "ALUConformGrid"
#elif  GRIDTYPE==4
#define CHARGRIDTYPE "AlbertaGrid"
#endif

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#ifndef SPACEDIM
#define SPACEDIM 2
#endif

#if ( SPACEDIM < 2 ) || ( SPACEDIM > 3 )
#error 'Dimension SPACEDIM must be 2 or 3'
#endif

#if ( GRIDTYPE==2 ) || ( GRIDTYPE==3 )
#define HAVE_DUNE_ALUGRID 1
#define ENABLE_DUNE_ALUGRID 1
#if SPACEDIM==2
#include "dune/alugrid/2d/alu2dinclude.hh"
#include "dune/alugrid/2d/alugrid.hh"
#include "dune/alugrid/2d/gridfactory.hh"
#else
#include "dune/alugrid/3d/alu3dinclude.hh"
#include "dune/alugrid/3d/alugrid.hh"
#include "dune/alugrid/3d/gridfactory.hh"
#endif
#endif

#if GRIDTYPE==4
#define ENABLE_ALBERTA 1
#define ALBERTA_DIM SPACEDIM
#include "dune/grid/albertagrid/agrid.hh"
#include "dune/grid/albertagrid/gridfactory.hh"
#endif

#include "utilities/enums.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "linalg/direct.hh"
#include "linalg/trivialpreconditioner.hh"
//#include "linalg/partialDirectPreconditioner.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/iluprecond.hh"      // ILUT, ILUK, ARMS
#include "linalg/iccprecond.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG, Euclid
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"

using namespace Kaskade;
#include "ht.hh"

#if SPACEDIM==3
#include "cubus.hh"
#endif

#if GRIDTYPE==4
#define DEFAULT_REFINEMENTS 14
#else
#if GRIDTYPE==3
#define DEFAULT_REFINEMENTS 7
#else
#if SPACEDIM==2
#define DEFAULT_REFINEMENTS 5
#endif
#if SPACEDIM==3
#define DEFAULT_REFINEMENTS 3
#endif
#endif
#endif

int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  std::cout << "Start heat transfer tutorial program (with GridType=" << CHARGRIDTYPE <<
              " and SpaceDimension=" << SPACEDIM << ")" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = true; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int  refinements = getParameter(pt, "refinements", DEFAULT_REFINEMENTS),
       order       =  getParameter(pt, "order", 2),
       verbosity   = getParameter(pt, "verbosity", 1);
  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "output level (verbosity)       : " << verbosity << std::endl;

  int  direct, onlyLowerTriangle = false;
    
  DirectType directType;
//  IterateType iterateType = IterateType::CG;
  MatrixProperties property = MatrixProperties::SYMMETRIC;
  PrecondType precondType = PrecondType::NONE;
  std::string empty;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);
    
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  directType = static_cast<DirectType>(getParameter(pt, s, 0));

//  s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
//  iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  precondType = static_cast<PrecondType>(getParameter(pt, s, 0));
  
  int blocks = getParameter(pt,"blocks",40);
  int nthreads = getParameter(pt,"threads",4);
  double rowBlockFactor = getParameter(pt,"rowBlockFactor",2.0);

  property = MatrixProperties::SYMMETRIC;

  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || ( (precondType == PrecondType::ICC) && !direct ) )
  {
    onlyLowerTriangle = true;
    std::cout << 
      "Note: direct solver MUMPS/PARADISO or PrecondType::ICC preconditioner ===> onlyLowerTriangle is set to true!" 
      << std::endl;
  }

  boost::timer::cpu_timer gridTimer;
#if SPACEDIM==2
  //   two-dimensional space: dim=2
  constexpr int dim=2;        
#if GRIDTYPE==1
  using Grid = Dune::UGGrid<dim>;
#endif
  // There are alternatives to UGGrid: ALUSimplexGrid (red refinement), ALUConformGrid (bisection)
  // and AlbertaGrid
#if GRIDTYPE==2
  using Grid = Dune::ALUGrid<dim,dim,Dune::ALUGridElementType::simplex,Dune::ALUGridRefinementType::nonconforming>;
#endif
#if GRIDTYPE==3
  using Grid = Dune::ALUGrid<dim,dim,Dune::ALUGridElementType::simplex,Dune::ALUGridRefinementType::conforming>;
#endif
#if GRIDTYPE==4
  using Grid = Dune::AlbertaGrid<dim,dim>;
#endif
  Dune::GridFactory<Grid> factory;

  // vertex coordinates v[0], v[1]
  Dune::FieldVector<double,dim> v;    
  v[0]=0; v[1]=0; factory.insertVertex(v);
  v[0]=1; v[1]=0; factory.insertVertex(v);
  v[0]=1; v[1]=1; factory.insertVertex(v);
  v[0]=0; v[1]=1; factory.insertVertex(v);
  // triangle defined by 3 vertex indices
  std::vector<unsigned int> vid(3);
  Dune::GeometryType gt(Dune::GeometryType::simplex,2);
  vid[0]=0; vid[1]=1; vid[2]=2; factory.insertElement(gt,vid);
  vid[0]=0; vid[1]=2; vid[2]=3; factory.insertElement(gt,vid);
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  // the coarse grid will be refined three times
  grid->globalRefine(refinements);
  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(2) << " points" << std::endl;

  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
#else
  //  three-dimensional space: dim=3
  //  we offer 2 geometries:
  //  - very simple:  1 tetrahedron defined by 4 vertices
  //  - more complex: 1 cube defined by 48 tetrahedra, provided in cubus.hh
  constexpr int dim=3;
    
  //    
  //definition of 1 tetrahedron by 4 vertices
#if GRIDTYPE==4
  using Grid = Dune::AlbertaGrid<dim,dim>;
  Dune::GridFactory<Grid> factory;
  // vertex coordinates v[0], v[1]
  Dune::FieldVector<double,dim> v;    
  v[0]=0; v[1]=0; v[2]=0; factory.insertVertex(v);
  v[0]=1; v[1]=0; v[2]=0; factory.insertVertex(v);
  v[0]=0; v[1]=1; v[2]=0; factory.insertVertex(v);
  v[0]=0; v[1]=0; v[2]=1; factory.insertVertex(v);
  // tetrahedron defined by 4 vertex indices
  std::vector<unsigned int> vid(4);
  Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
  vid[0]=0; vid[1]=1; vid[2]=2; vid[3]=3; factory.insertElement(gt,vid);
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  // the coarse grid will be refined three times
  grid->globalRefine(refinements);
  //
#else 
  // definition of a more complex mesh using cubus.hh
  // note: trying to use the following code with AlbertaGrid will lead to a crash
  // during runtime due to a bug in the AlbertaGrid refinement routine
  int heapSize=1024;
#if GRIDTYPE==1
  using Grid = Dune::UGGrid<dim>;
#endif
  // There are alternatives to UGGrid: ALUSimplexGrid (red refinement)
#if GRIDTYPE==2
  using Grid = Dune::ALUGrid<dim,dim,Dune::ALUGridElementType::simplex,Dune::ALUGridRefinementType::nonconforming>;
#endif
#if GRIDTYPE==3
#error ALUCONFORM GridType not supported by DUNE for Dimension=3
#endif
  std::unique_ptr<Grid> grid( RefineGrid<Grid>(refinements, heapSize) );
#endif

  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " tetrahedra, " << std::endl;
  std::cout << "      " << grid->size(1) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(dim-1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(dim) << " points" << std::endl;
  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
#endif
  std::cout << "computing time for generation of initial mesh: " << boost::timer::format(gridTimer.elapsed()) << "\n";
    
  using LeafView = Grid::LeafGridView;
    
  // construction of finite element space for the scalar solution T.
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  // using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = HeatFunctional<double,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = Functional::TestVars::noOfVariables;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  // avoid collision of reference to the Kaskade CG with an equal named enum value in the Alberta headers
  using CG = Kaskade::CG<LinearSpace,LinearSpace>;
    
  // construction of finite element space for the scalar solution T.
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
    
  Spaces spaces(&temperatureSpace);
    
  // construct variable list.
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
        
  std::string varNames[1] = { "u" };
    
  VariableSet variableSet(spaces,varNames);

  // construct variational functional
    
  double kappa = 1.0;
  double q = 1.0;
  Functional F(kappa,q);
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  std::cout << std::endl << "no of variables = " << nvars << std::endl;
  std::cout << "no of equations = " << neq   << std::endl;
  size_t dofs = variableSet.degreesOfFreedom(0,nvars);
  std::cout << "number of degrees of freedom = " << dofs   << std::endl;

  //construct Galerkin representation
  Assembler assembler(gridManager,spaces);
  VariableSet::VariableSet u(variableSet);
  VariableSet::VariableSet du(variableSet);

  size_t nnz = assembler.nnz(0,neq,0,nvars,onlyLowerTriangle);
  std::cout << "number of nonzero elements in the stiffness matrix: " << nnz << std::endl << std::endl;
  
  boost::timer::cpu_timer assembTimer;
  CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,neq>::init(spaces));
  solution = 0;
  
  // UG seems to admit concurrent reads while claiming not to be thread safe. In this case we enforce multithreading during assembly.
  gridManager.enforceConcurrentReads(std::is_same<Grid,Dune::UGGrid<dim> >::value);
  assembler.setNSimultaneousBlocks(blocks);
  assembler.setRowBlockFactor(rowBlockFactor);
  assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
  std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed()) << "\n";
  
  CoefficientVectors rhs(assembler.rhs());
  AssembledGalerkinOperator<Assembler,0,neq,0,nvars> A(assembler, onlyLowerTriangle);
  
  // matrix may be used in triplet format, e.g.,
  // MatrixAsTriplet<double> tri = A.get<MatrixAsTriplet<double> >();
  //     for (k=0; k< nnz; k++)
  //       {
  //         printf("%3d %3d %e\n", tri.ridx[k], tri.cidx[k], tri.data[k]);
  //       }

  if (direct)
  {
    boost::timer::cpu_timer directTimer;
    directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
    u.data = solution.data;
    std::cout << "computing time for direct solve: " << boost::timer::format(directTimer.elapsed()) << "\n";
  }
  else
  {
    boost::timer::cpu_timer iteTimer;
    Dune::InverseOperatorResult res;
    const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
    int iteSteps = getParameter(pt, "solver.iteMax", 2000);
    double iteEps = getParameter(pt, "solver.iteEps", 1.0e-10);
    StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
    int lookAhead;
    switch (precondType)
    {
      case PrecondType::NONE:
      case PrecondType::HB:   lookAhead=50; break;
      default:                lookAhead=3; break;
    }
    lookAhead = getParameter(pt, "solver.lookAhead", lookAhead);
    termination.setLookAhead(lookAhead);

    switch (precondType)
    {
      case PrecondType::NONE:
      {
        std::cout << "selected preconditioner: NONE" << std::endl;
        TrivialPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > trivial;
        CG cg(A,trivial,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ADDITIVESCHWARZ:
      {
        std::cout << "selected preconditioner: ADDITIVESCHWARZ" << std::endl;
        std::pair<size_t,size_t> idx = temperatureSpace.mapper().globalIndexRange(gridManager.grid().leafIndexSet().geomTypes(dim)[0]);
        AdditiveSchwarzPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > addschwarz(A,idx.first,idx.second,verbosity);
        CG cg(A,addschwarz,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUT:
      {
        std::cout << "selected preconditioner: ILUT" << std::endl;
//        std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int lfil = getParameter(pt, "solver.ILUT.lfil", 140);
        double dropTol = getParameter(pt, "solver.ILUT.dropTol", 0.01);
        ILUTPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > ilut(A,lfil,dropTol,verbosity);
        CG cg(A,ilut,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::BiCGSTABSolver<LinearSpace> cg(A,ilut,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUK:
      {
        std::cout << "selected preconditioner: ILUK" << std::endl;
//        std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int fill_lev = getParameter(pt, "solver.ILUK.fill_lev", 3);
        ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,fill_lev,verbosity);
        CG cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ARMS:
      {
        int lfil = getParameter(pt, "solver.ARMS.lfil", 140);
        int lev_reord = getParameter(pt, "solver.ARMS.lev_reord", 1);
        double dropTol = getParameter(pt, "solver.ARMS.dropTol", 0.01);
        double tolind = getParameter(pt, "solver.ARMS.tolind", 0.2);
        ARMSPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,lfil,dropTol,lev_reord,tolind,verbosity);
        CG cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ICC:
      {
        std::cout << "selected preconditioner: ICC" << std::endl;
        if (property != MatrixProperties::SYMMETRIC) 
        {
          std::cout << "PrecondType::ICC preconditioner of TAUCS lib has to be used with matrix.property==MatrixProperties::SYMMETRIC\n";
          std::cout << "i.e., call the executable with option --solver.property MatrixProperties::SYMMETRIC\n\n";
        }
        double dropTol = getParameter(pt, "solver.ICC.dropTol", 0.01);;
        ICCPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > icc(A,dropTol);
        CG cg(A,icc,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ICC0:
      {
        std::cout << "selected preconditioner: ICC0" << std::endl;
        ICC_0Preconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > icc0(A);
        CG cg(A,icc0,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::HB:
      {
        std::cout << "selected preconditioner: HB" << std::endl;
        HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type, AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type > hb(gridManager.grid());
        CG cg(A,hb,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::BOOMERAMG:
      {
        int steps = getParameter(pt, "solver.BOOMERAMG.steps", iteSteps);
        int coarsentype = getParameter(pt, "solver.BOOMERAMG.coarsentype", 21);
        int interpoltype = getParameter(pt, "solver.BOOMERAMG.interpoltype", 0);
        int cycleType = getParameter(pt, "solver.BOOMERAMG.cycleType", 1);
        int relaxType = getParameter(pt, "solver.BOOMERAMG.relaxType", 3);
        int variant = getParameter(pt, "solver.BOOMERAMG.variant", 0);
        int overlap = getParameter(pt, "solver.BOOMERAMG.overlap", 1);
        double tol = getParameter(pt, "solver.BOOMERAMG.tol", iteEps);
        double strongThreshold = getParameter(pt, "solver.BOOMERAMG.strongThreshold", (dim==2)?0.25:0.6);
        BoomerAMG<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> >
                  boomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
                  strongThreshold,variant,overlap,1,verbosity);
        CG cg(A,boomerAMGPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::LoopSolver<LinearSpace> cg(A,boomerAMGPrecon,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::EUCLID:
      {
        std::cout << "selected preconditioner: EUCLID" << std::endl;
        int level      = getParameter(pt, "solver.EUCLID.level",1);
        double droptol = getParameter(pt, "solver.EUCLID.droptol",0.01);
        int printlevel = 0;
        if (verbosity>2) printlevel=verbosity-2;
        printlevel = getParameter(pt,"solver.EUCLID.printlevel",printlevel);
        int bj = getParameter(pt, "solver.EUCLID.bj",0);
        Euclid<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > EuclidPrecon(A,level,droptol,printlevel,bj,verbosity);
        CG cg(A,EuclidPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
     }
      break;
      case PrecondType::JACOBI:
      default:
      {
        std::cout << "selected preconditioner: JACOBI" << std::endl;
        JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobi(A,1.0);
        CG cg(A,jacobi,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
    }
    solution *= -1.0;
    u.data = solution.data;
    
    std::cout << "iterative solve eps= " << iteEps << ": " 
              << (res.converged?"converged":"failed") << " after "
              << res.iterations << " steps, rate="
              << res.conv_rate << ", computing time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
  }
  
    
  
  // compute L2 norm of the solution
  boost::timer::cpu_timer outputTimer;
  L2Norm l2Norm;
  std::cout << "L2norm(solution) = " << l2Norm(boost::fusion::at_c<0>(u.data)) << std::endl;
    
    

  // output of solution in VTK format for visualization,
  // the data are written as ascii stream into file temperature.vtu,
  // possible is also binary
  writeVTKFile(u,"temperature",IoOptions().setOrder(order).setPrecision(7));

  std::cout << "graphical output finished, data in VTK format is written into file temperature.vtu \n";
    
  // output of solution for Amira visualization,
  // the data are written in binary format into file temperature.am,
  // possible is also ascii
  // IoOptions options;
  // options.outputType = IoOptions::ascii;
  // LeafView leafGridView = gridManager.grid().leafGridView();
  // writeAMIRAFile(leafGridView,variableSet,u,"temperature",options);

  std::cout << "computing time for output: " << boost::timer::format(outputTimer.elapsed()) << "\n";

  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End heat transfer tutorial program" << std::endl;
}
