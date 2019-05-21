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

#ifndef CUBUS_HH
#define CUBUS_HH

template <class GRID>
std::unique_ptr<GRID, std::default_delete<GRID> > create(int heapSize=500)
{
//  GRID::setDefaultHeapSize(heapSize);
  Dune::GridFactory<GRID> factory;
	
	// vertices
	Dune::FieldVector<double,3> v;
	v[0]=0; v[1]=0; v[2]=1; factory.insertVertex(v);
	v[0]=0; v[1]=1; v[2]=1; factory.insertVertex(v);
	v[0]=1; v[1]=1; v[2]=1; factory.insertVertex(v);
	v[0]=1; v[1]=0; v[2]=1; factory.insertVertex(v);
	v[0]=0; v[1]=0; v[2]=0; factory.insertVertex(v);
	v[0]=0; v[1]=1; v[2]=0; factory.insertVertex(v);
	v[0]=1; v[1]=1; v[2]=0; factory.insertVertex(v);
	v[0]=1; v[1]=0; v[2]=0; factory.insertVertex(v);
	v[0]=0.5; v[1]=1; v[2]=1; factory.insertVertex(v);
	v[0]=0.5; v[1]=1; v[2]=0.5; factory.insertVertex(v);
	v[0]=0; v[1]=1; v[2]=0.5; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.5; v[2]=1; factory.insertVertex(v);
	v[0]=0; v[1]=0.5; v[2]=0.5; factory.insertVertex(v);
	v[0]=0; v[1]=0.5; v[2]=1; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.5; v[2]=0.5; factory.insertVertex(v);
	v[0]=0; v[1]=0.5; v[2]=0; factory.insertVertex(v);
	v[0]=0; v[1]=0; v[2]=0.5; factory.insertVertex(v);
	v[0]=0.5; v[1]=1; v[2]=0; factory.insertVertex(v);
	v[0]=0.5; v[1]=0.5; v[2]=0; factory.insertVertex(v);
	v[0]=1; v[1]=1; v[2]=0.5; factory.insertVertex(v);
	v[0]=1; v[1]=0.5; v[2]=1; factory.insertVertex(v);
	v[0]=0.5; v[1]=0; v[2]=0.5; factory.insertVertex(v);
	v[0]=0.5; v[1]=0; v[2]=1; factory.insertVertex(v);
	v[0]=1; v[1]=0.5; v[2]=0; factory.insertVertex(v);
	v[0]=0.5; v[1]=0; v[2]=0; factory.insertVertex(v);
	v[0]=1; v[1]=0.5; v[2]=0.5; factory.insertVertex(v);
	v[0]=1; v[1]=0; v[2]=0.5; factory.insertVertex(v);

	// tetrahedron defined by 4 vertex indices
	std::vector<unsigned int> vid(4);
	Dune::GeometryType gt(Dune::GeometryType::simplex,3);
   	vid[0]=13;   vid[1]=12;   vid[2]=11;   vid[3]= 0; factory.insertElement(gt,vid);
   	vid[0]=13;   vid[1]= 8;   vid[2]=10;   vid[3]= 1; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]= 9;   vid[2]= 8;   vid[3]= 2; factory.insertElement(gt,vid);
   	vid[0]=12;   vid[1]=10;   vid[2]= 9;   vid[3]= 5; factory.insertElement(gt,vid);
   	vid[0]= 9;   vid[1]=12;   vid[2]=13;   vid[3]=10; factory.insertElement(gt,vid);
   	vid[0]=13;   vid[1]= 9;   vid[2]=10;   vid[3]= 8; factory.insertElement(gt,vid);
   	vid[0]=13;   vid[1]=11;   vid[2]=12;   vid[3]= 9; factory.insertElement(gt,vid);
   	vid[0]= 9;   vid[1]= 8;   vid[2]=13;   vid[3]=11; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=12;   vid[2]=16;   vid[3]= 0; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=14;   vid[2]= 9;   vid[3]= 2; factory.insertElement(gt,vid);
   	vid[0]=16;   vid[1]=15;   vid[2]=14;   vid[3]= 4; factory.insertElement(gt,vid);
   	vid[0]=12;   vid[1]= 9;   vid[2]=15;   vid[3]= 5; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=12;   vid[2]=11;   vid[3]= 9; factory.insertElement(gt,vid);
   	vid[0]=12;   vid[1]= 9;   vid[2]=14;   vid[3]=15; factory.insertElement(gt,vid);
   	vid[0]=12;   vid[1]=11;   vid[2]=16;   vid[3]=14; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=12;   vid[2]=15;   vid[3]=16; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=19;   vid[2]= 9;   vid[3]= 2; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=15;   vid[2]=18;   vid[3]= 4; factory.insertElement(gt,vid);
   	vid[0]= 9;   vid[1]=17;   vid[2]=15;   vid[3]= 5; factory.insertElement(gt,vid);
   	vid[0]=19;   vid[1]=18;   vid[2]=17;   vid[3]= 6; factory.insertElement(gt,vid);
   	vid[0]=17;   vid[1]=19;   vid[2]=14;   vid[3]=18; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=17;   vid[2]=18;   vid[3]=15; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]= 9;   vid[2]=19;   vid[3]=17; factory.insertElement(gt,vid);
   	vid[0]=17;   vid[1]=15;   vid[2]=14;   vid[3]= 9; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=16;   vid[2]=22;   vid[3]= 0; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=20;   vid[2]=14;   vid[3]= 2; factory.insertElement(gt,vid);
   	vid[0]=22;   vid[1]=21;   vid[2]=20;   vid[3]= 3; factory.insertElement(gt,vid);
   	vid[0]=16;   vid[1]=14;   vid[2]=21;   vid[3]= 4; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=16;   vid[2]=11;   vid[3]=14; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=21;   vid[2]=14;   vid[3]=20; factory.insertElement(gt,vid);
   	vid[0]=11;   vid[1]=22;   vid[2]=16;   vid[3]=21; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=20;   vid[2]=11;   vid[3]=22; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=26;   vid[2]=25;   vid[3]= 3; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=18;   vid[2]=24;   vid[3]= 4; factory.insertElement(gt,vid);
   	vid[0]=25;   vid[1]=23;   vid[2]=18;   vid[3]= 6; factory.insertElement(gt,vid);
   	vid[0]=26;   vid[1]=24;   vid[2]=23;   vid[3]= 7; factory.insertElement(gt,vid);
   	vid[0]=23;   vid[1]=26;   vid[2]=21;   vid[3]=24; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=23;   vid[2]=24;   vid[3]=18; factory.insertElement(gt,vid);
   	vid[0]=21;   vid[1]=25;   vid[2]=26;   vid[3]=23; factory.insertElement(gt,vid);
   	vid[0]=23;   vid[1]=18;   vid[2]=21;   vid[3]=25; factory.insertElement(gt,vid);
   	vid[0]=20;   vid[1]=19;   vid[2]=14;   vid[3]= 2; factory.insertElement(gt,vid);
   	vid[0]=20;   vid[1]=21;   vid[2]=25;   vid[3]= 3; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=18;   vid[2]=21;   vid[3]= 4; factory.insertElement(gt,vid);
   	vid[0]=19;   vid[1]=25;   vid[2]=18;   vid[3]= 6; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=25;   vid[2]=19;   vid[3]=20; factory.insertElement(gt,vid);
   	vid[0]=25;   vid[1]=19;   vid[2]=18;   vid[3]=14; factory.insertElement(gt,vid);
   	vid[0]=25;   vid[1]=14;   vid[2]=21;   vid[3]=20; factory.insertElement(gt,vid);
   	vid[0]=14;   vid[1]=25;   vid[2]=21;   vid[3]=18; factory.insertElement(gt,vid);

	return static_cast<std::unique_ptr<GRID, std::default_delete<GRID> > >(factory.createGrid());
}
#endif
