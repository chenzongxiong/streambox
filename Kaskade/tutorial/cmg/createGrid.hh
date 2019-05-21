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

#ifndef CREATEGRID_HH
#define CREATEGRID_HH

template <class Grid, const int dim>
void createUnitSquare(Dune::GridFactory<Grid> &factory)
  {
	Dune::FieldVector<double,dim> v; 	
	std::vector<unsigned int> vid(dim+1);
	Dune::GeometryType gt(Dune::GeometryType::simplex,dim);

	v[0]=0; v[1]=0; factory.insertVertex(v);
	v[0]=1; v[1]=0; factory.insertVertex(v);
	v[0]=1; v[1]=1; factory.insertVertex(v);
	v[0]=0; v[1]=1; factory.insertVertex(v);

	vid[0]=0; vid[1]=1; vid[2]=2; factory.insertElement(gt,vid);
	vid[0]=2; vid[1]=3; vid[2]=0; factory.insertElement(gt,vid);
  }

template <class Grid, const int dim>
void createUnitCrack(Dune::GridFactory<Grid> &factory, double width = 0.03)
  {
	Dune::FieldVector<double,dim> v; 	
	v[0]=0; v[1]=0; factory.insertVertex(v);
	v[0]=1; v[1]=0; factory.insertVertex(v);
	v[0]=1; v[1]=1; factory.insertVertex(v);
	v[0]=0; v[1]=1; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.5; factory.insertVertex(v);
	v[0]=1; v[1]=0.49; factory.insertVertex(v);
	v[0]=1; v[1]=0.51; factory.insertVertex(v);

	std::vector<unsigned int> vid(dim+1);
	Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
	vid[0]=0; vid[1]=1; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=1; vid[1]=5; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=6; vid[1]=2; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=4; vid[1]=2; vid[2]=3; factory.insertElement(gt,vid);
	vid[0]=0; vid[1]=4; vid[2]=3; factory.insertElement(gt,vid);
  }

template <class Grid, const int dim>
void createOriginalBoDD(Dune::GridFactory<Grid> &factory, double width = 0.03)
  {
	Dune::FieldVector<double,dim> v; 	
	v[0]=-1; v[1]=-1; factory.insertVertex(v);
	v[0]=1; v[1]=-1; factory.insertVertex(v);
	v[0]=1; v[1]=1; factory.insertVertex(v);
	v[0]=-1; v[1]=1; factory.insertVertex(v);
	v[0]=0.0; v[1]=0.0; factory.insertVertex(v);
	v[0]=1; v[1]=-width; factory.insertVertex(v);
	v[0]=1; v[1]=width; factory.insertVertex(v);

	std::vector<unsigned int> vid(dim+1);
	Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
	vid[0]=0; vid[1]=1; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=1; vid[1]=5; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=6; vid[1]=2; vid[2]=4; factory.insertElement(gt,vid);
	vid[0]=4; vid[1]=2; vid[2]=3; factory.insertElement(gt,vid);
	vid[0]=0; vid[1]=4; vid[2]=3; factory.insertElement(gt,vid);
  }

template <class Grid, const int dim>
void createAreaGrid(Dune::GridFactory<Grid> &factory)
  {
	Dune::FieldVector<double,dim> v; 	
	v[0]=0.0; v[1]=0.0; factory.insertVertex(v);
	v[0]=0.3; v[1]=0.0; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.0; factory.insertVertex(v);
	v[0]=0.7; v[1]=0.0; factory.insertVertex(v);
	v[0]=1.0; v[1]=0.0; factory.insertVertex(v);
	v[0]=0.0; v[1]=0.3; factory.insertVertex(v);
	v[0]=0.3; v[1]=0.3; factory.insertVertex(v);
	v[0]=0.7; v[1]=0.3; factory.insertVertex(v);
	v[0]=1.0; v[1]=0.3; factory.insertVertex(v);
	v[0]=0.0; v[1]=0.5; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.5; factory.insertVertex(v);
	v[0]=1.0; v[1]=0.5; factory.insertVertex(v);
	v[0]=0.0; v[1]=0.7; factory.insertVertex(v);
	v[0]=0.3; v[1]=0.7; factory.insertVertex(v);
	v[0]=0.7; v[1]=0.7; factory.insertVertex(v);
	v[0]=1.0; v[1]=0.7; factory.insertVertex(v);
	v[0]=0.0; v[1]=1.0; factory.insertVertex(v);
	v[0]=0.3; v[1]=1.0; factory.insertVertex(v);
	v[0]=0.5; v[1]=1.0; factory.insertVertex(v);
	v[0]=0.7; v[1]=1.0; factory.insertVertex(v);
	v[0]=1.0; v[1]=1.0; factory.insertVertex(v);

	std::vector<unsigned int> vid(dim+1);
	Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
	vid[0]= 0; vid[1]= 1; vid[2]= 6; factory.insertElement(gt,vid);
	vid[0]= 0; vid[1]= 6; vid[2]= 5; factory.insertElement(gt,vid);
	vid[0]= 1; vid[1]= 2; vid[2]= 6; factory.insertElement(gt,vid);
	vid[0]= 6; vid[1]= 2; vid[2]= 7; factory.insertElement(gt,vid);
	vid[0]= 2; vid[1]= 3; vid[2]= 7; factory.insertElement(gt,vid);
	vid[0]= 3; vid[1]= 4; vid[2]= 7; factory.insertElement(gt,vid);
	vid[0]= 7; vid[1]= 4; vid[2]= 8; factory.insertElement(gt,vid);

	vid[0]= 5; vid[1]= 6; vid[2]= 9; factory.insertElement(gt,vid);
	vid[0]=12; vid[1]= 9; vid[2]=13; factory.insertElement(gt,vid);
	vid[0]= 9; vid[1]= 6; vid[2]=13; factory.insertElement(gt,vid);
	vid[0]=13; vid[1]= 6; vid[2]=10; factory.insertElement(gt,vid);
	vid[0]= 6; vid[1]= 7; vid[2]=10; factory.insertElement(gt,vid);
	vid[0]=13; vid[1]=10; vid[2]=14; factory.insertElement(gt,vid);
	vid[0]=10; vid[1]= 7; vid[2]=14; factory.insertElement(gt,vid);
	vid[0]=14; vid[1]= 7; vid[2]=11; factory.insertElement(gt,vid);
	vid[0]= 7; vid[1]= 8; vid[2]=11; factory.insertElement(gt,vid);
	vid[0]=14; vid[1]=11; vid[2]=15; factory.insertElement(gt,vid);

	vid[0]=12; vid[1]=13; vid[2]=16; factory.insertElement(gt,vid);
	vid[0]=16; vid[1]=13; vid[2]=17; factory.insertElement(gt,vid);
	vid[0]=17; vid[1]=13; vid[2]=18; factory.insertElement(gt,vid);
	vid[0]=13; vid[1]=14; vid[2]=18; factory.insertElement(gt,vid);
	vid[0]=18; vid[1]=14; vid[2]=19; factory.insertElement(gt,vid);
	vid[0]=19; vid[1]=14; vid[2]=20; factory.insertElement(gt,vid);
	vid[0]=14; vid[1]=15; vid[2]=20; factory.insertElement(gt,vid);
  }

#endif
