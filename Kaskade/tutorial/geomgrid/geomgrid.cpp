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
#include "dune/grid/geometrygrid.hh"
#include "dune/grid/common/gridinfo.hh"
#include "dune/istl/solvers.hh"

#include "fem/assemble.hh"
#include "fem/embedded_errorest.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"
#include "linalg/trivialpreconditioner.hh"
#include "linalg/direct.hh"
#include "linalg/triplet.hh"
#include "linalg/iluprecond.hh"
#include "linalg/additiveschwarz.hh"
#include "linalg/hyprecond.hh"
#include "io/vtk.hh"
#include "io/amira.hh"

#include "utilities/kaskopt.hh"

#include "ht.hh"

 namespace Dune
  {

    namespace Capabilities
    {

      template< class Grid >
      struct hasHierarchicIndexSet
      {
        static const bool v = false;
      };

      template< class Grid >
      struct hasHierarchicIndexSet< const Grid >
      {
        static const bool v = hasHierarchicIndexSet< Grid >::v;
      };

    }

  }


class SquareToCircle
  : public Dune::AnalyticalCoordFunction< double, 2, 2, SquareToCircle >
  {
    typedef SquareToCircle This;
    typedef Dune::AnalyticalCoordFunction< double, 2, 2, This > Base;

  public:
    typedef Base::DomainVector DomainVector;
    typedef Base::RangeVector RangeVector;

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      double enorm = sqrt(x[0]*x[0]+x[1]*x[1]);

      if(enorm > 0.00001)
      {

        double radius = sqrt(2.0)*std::max(std::fabs(x[0]),std::fabs(x[1]));

      

      double scaling = radius/enorm;


      y[ 0 ] = x[ 0 ]*scaling;
      y[ 1 ] = x[ 1 ]*scaling;
      } else
      {
      y[ 0 ] = x[ 0 ];
      y[ 1 ] = x[ 1 ];
      }
        
    }
  };


int main(int argc, char *argv[])
  {

	//   two-dimensional space: dim=2
	int const dim=2; 		
	typedef Dune::UGGrid<dim> Grid;

        typedef Dune::GeometryGrid<Grid,SquareToCircle> GGrid;

	Dune::GridFactory<Grid> factory;

	// vertex coordinates v[0], v[1]
	Dune::FieldVector<double,dim> v; 	
	v[0]=-1; v[1]=-1; factory.insertVertex(v);
	v[0]=1; v[1]=-1; factory.insertVertex(v);
	v[0]=1; v[1]=1; factory.insertVertex(v);
	v[0]=-1; v[1]=1; factory.insertVertex(v);
	v[0]=0; v[1]=0; factory.insertVertex(v);
	// triangle defined by 3 vertex indices
	std::vector<unsigned int> vid(3);
	Dune::GeometryType gt(Dune::GeometryType::simplex,2);
	vid[0]=0; vid[1]=1; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=1; vid[1]=2; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=2; vid[1]=3; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=3; vid[1]=0; vid[2]=4; factory.insertElement(gt,vid);
	std::auto_ptr<Grid> grid( factory.createGrid() ) ;

        

	// the coarse grid will be refined three times
	// some information on the refined mesh
	std::cout << "Grid: " << grid->size(0) << " triangles, " << std::endl;
	std::cout << "      " << grid->size(1) << " edges, " << std::endl;
	std::cout << "      " << grid->size(2) << " points" << std::endl;
	// a gridmanager is constructed 
	// as connector between geometric and algebraic information
        
        SquareToCircle def;

        std::auto_ptr<GGrid> ggrid(new GGrid(*grid,def));

	ggrid->globalRefine(6);

	GridManager<GGrid> gridManager(ggrid);    


//StartSnippet2
	// construction of finite element space for the scalar solution T
	typedef GGrid::LeafGridView GridView;
	typedef FEFunctionSpace<ContinuousLagrangeMapper<double,GridView> > H1Space;
	H1Space temperatureSpace(gridManager,gridManager.grid().leafView(),
							 1);
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
//StopSnippet2

//StartSnippet3
	typedef HeatFunctional<double,VariableSet> Functional;
	Functional F;
	typedef VariationalFunctionalAssembler<LinearizationAt<Functional> > GOP;
	typedef GOP::TestVariableRepresentation<>::type Rhs;
	GOP gop(gridManager.signals,spaces);
	VariableSet::VariableSet x(variableSet);
//StopSnippet3

//StartSnippet4

	typedef GOP::AnsatzVariableRepresentation<>::type Sol;
	Sol solution(GOP::AnsatzVariableRepresentation<>::init(gop));
	solution = 0;
	gop.assemble(linearization(F,x));


	Rhs rhs(GOP::TestVariableRepresentation<>::rhs(gop));
	AssembledGalerkinOperator<GOP,0,1,0,1> A(gop, false);
	AssembledGalerkinOperator<GOP,0,1,0,1>::matrix_type tri(A.getmat());
//	A.getmat().print();


//StartSnippet5
    boost::timer directTimer;
    directInverseOperator(A,DirectType::UMFPACK,MatrixProperties::GENERAL).applyscaleadd(-1.0,rhs,solution);
    std::cout << "direct solve: " << directTimer.elapsed() << "s\n";
    x.data = solution.data;
//StopSnippet5

//StartSnippet6
	// output of solution in VTK format for visualization,
	// the data are written as ascii stream into file temperature.vtu,
	// possible is also binary
 	IoOptions options;
 	options.outputType = IoOptions::ascii;
 
    typedef GGrid::LeafGridView LeafGridView;
    LeafGridView leafView = gridManager.grid().leafView();
    writeVTKFile(leafView,variableSet,x,"temperature",options);

  }
