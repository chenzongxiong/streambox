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

#include <cstdlib>
#include <iostream>
#include <boost/timer.hpp>
#include "dune/config.h"
#include "dune/common/stdstreams.hh"
#include "dune/grid/uggrid.hh"
#include "lagrangespace.hh"
#include "fetransfer.hh"

bool check(int order)
{
  int const dim = 2;
  int const level = 3;
  
  typedef Dune::UGGrid<dim> Grid;
  typedef Grid::Traits::LeafIndexSet IS;
  
  std::auto_ptr<Grid> grid(new Grid);
  grid->createBegin();
  Dune::FieldVector<double,dim> v;
  // vertices
  v[0]=0; v[1]=0;
  grid->insertVertex(v);
  v[0]=1; v[1]=0;
  grid->insertVertex(v);
  v[0]=1; v[1]=1;
  grid->insertVertex(v);
  v[0]=0; v[1]=1;
  grid->insertVertex(v);
  // elements
  std::vector<unsigned int> vid(3);
  vid[0]=0; vid[1]=1; vid[2]=2;
  grid->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  vid[0]=0; vid[1]=2; vid[2]=3;
  grid->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),vid);
  grid->createEnd();  
  grid->globalRefine(level);

  GridManager<Grid> gridManager(grid);

  bool pass = true;
  

  typedef FEFunctionSpace<ContinuousLagrangeMapper<double,Grid> > H1Space;
  H1Space spaceH1(gridManager,gridManager.grid().leafIndexSet(),order);

  typedef FEFunctionSpace<DiscontinuousLagrangeMapper<double,Grid> > L2Space;
  L2Space spaceL2(gridManager,gridManager.grid().leafIndexSet(),order-1);
  
  
  
  H1Space::Element<2>::type f(spaceH1);
  for (int i=0; i<spaceH1.dimension(); ++i)
    (*f)[i] = (rand()%1000)/1000.0;

  L2Space::Element<3>::type g(spaceL2);
  for (int i=0; i<spaceH1.dimension(); ++i)
    (*g)[i] = (rand()%1000)/1000.0;


  std::vector<Dune::FieldVector<double,dim> > x(1000);
  std::vector<H1Space::Element<2>::type::ValueType> fVal(x.size());
  std::vector<L2Space::Element<3>::type::ValueType> gVal(x.size());
  for (int i=0; i<x.size(); ++i) {
    for (int j=0; j<dim; ++j)
      x[i][j] = (rand()%100002)/100001.0; // crude numbers to prevent test points to lie on cell boundaries.
    fVal[i] = f.value(x[i]);
    gVal[i] = g.value(x[i]);
  }

  for (int l=0; l<7; ++l) {
    for (Grid::Codim<0>::Partition<Dune::All_Partition>::LeafIterator ci=gridManager.grid().leafbegin<0>();
         ci!=gridManager.grid().leafend<0>(); ++ci) {
      int r = rand()%100;
      if (r < 20)
        gridManager.mark((l>3 && ci->level()>level)? -1: 1,ci);
    }
    gridManager.adaptAtOnce();

    double ferr = 0, gerr = 0, favg = 0, gavg = 0;
    for (int i=0; i<x.size(); ++i) {
      double fd = (fVal[i]-f.value(x[i])).two_norm();
      double gd = (gVal[i]-g.value(x[i])).two_norm();
      ferr = std::max(ferr,fd);
      gerr = std::max(gerr,gd);
      favg += fd/x.size();
      gavg += gd/x.size();
    }
    
    if (ferr > 1e-13) {
      std::cerr << "Continuous mapper of order " << order << " failed in step " << l
                << ": maximal error " << ferr << " average error " << favg << '\n';
      pass = false;
    }
    
    if (gerr > 1e-13) {
      std::cerr << "Disontinuous mapper of order " << order << " failed in step " << l
                << ": maximal error " << gerr << " average error " << gavg << '\n';
      pass = false;
    }
  }
  return pass;
}

int main(void) 
{
  bool pass = true;
  for (int i=1; i<5; ++i)
    pass &= check(i);

  if (pass)
    std::cout << "passed.\n";

  return 0;
}
