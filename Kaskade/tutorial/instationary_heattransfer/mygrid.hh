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

#ifndef MYGRID_HH
#define MYGRID_HH

#include <memory>

template <class Grid>
std::unique_ptr<Grid> RefineGrid(int refinements, int heapSize=500)
{
  // Grid::setDefaultHeapSize(heapSize);
  Dune::GridFactory<Grid> factory;
  Dune::FieldVector<double,2> v;
  
  
  /*
  // vertices
  v[0]=0; v[1]=0;
  factory.insertVertex(v);
  v[0]=1; v[1]=0;
  factory.insertVertex(v);
  v[0]=1; v[1]=1;
  factory.insertVertex(v);
  v[0]=0; v[1]=1;
  factory.insertVertex(v);
  // elements
  std::vector<unsigned int> vid(3);
  vid[0]=0; vid[1]=1; vid[2]=2;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=2; vid[2]=3;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  */

  // vertices
  //   (0,0)...(1,1)
  v[0]=0.5; v[1]=0.5;
  factory.insertVertex(v);
  v[0]=1; v[1]=0.5;
  factory.insertVertex(v);
  v[0]=1; v[1]=1;
  factory.insertVertex(v);
  v[0]=0.5; v[1]=1;
  factory.insertVertex(v);
  v[0]=0.75; v[1]=0.75;
  factory.insertVertex(v);

  v[0]=0; v[1]=1;
  factory.insertVertex(v);
  v[0]=0; v[1]=0.5;
  factory.insertVertex(v);
  v[0]=0.25; v[1]=0.75;
  factory.insertVertex(v);

  v[0]=0; v[1]=0;
  factory.insertVertex(v);
  v[0]=0.5; v[1]=0;
  factory.insertVertex(v);
  v[0]=0.25; v[1]=0.25;
  factory.insertVertex(v);

  v[0]=1; v[1]=0;
  factory.insertVertex(v);
  v[0]=0.75; v[1]=0.25;
  factory.insertVertex(v);
	
	
  // elements
  std::vector<unsigned int> vid(3);
  vid[0]=0; vid[1]=1; vid[2]=4;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=4; vid[2]=3;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=2; vid[1]=3; vid[2]=4;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=2; vid[1]=4; vid[2]=1;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);

  vid[0]=3; vid[1]=5; vid[2]=7;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=3; vid[2]=7;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=5; vid[1]=6; vid[2]=7;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=7; vid[2]=6;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);

  vid[0]=0; vid[1]=6; vid[2]=10;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=10; vid[2]=9;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=6; vid[1]=8; vid[2]=10;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=8; vid[1]=9; vid[2]=10;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);

  vid[0]=0; vid[1]=9; vid[2]=12;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=12; vid[2]=1;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=1; vid[1]=12; vid[2]=11;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=9; vid[1]=11; vid[2]=12;
  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  
	std::unique_ptr<Grid> grid( factory.createGrid() );
	grid->globalRefine(refinements);  
	return grid;
}

#endif
