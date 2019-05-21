/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */





#ifdef UNITTEST

#include <iostream>
#include <cmath>

#include "dune/grid/config.h"

#include <dune/common/fvector.hh>
#include "dune/grid/uggrid.hh"

#include "fem/advect.hh"
#include <fem/fetransfer.hh>
#include <fem/functionspace.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>

using namespace Kaskade;


// This test case creates a (arbitrary) function on a square and rotates it counterclockwise by pi/2. 
int main(void)
{
  int const order = 2;
  constexpr int m = 1;
  using Grid = Dune::UGGrid<2>;
  Dune::FieldVector<double,2> x0(-5), dx(10);
  GridManager<Grid> gridManager(createRectangle<Grid>(x0,dx,0.5));
  
  typedef Grid::LeafGridView LeafView;
  LeafView leafView = gridManager.grid().leafGridView();

  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView>>;
  H1Space h1Space2(gridManager,leafView,order);
  
  // Create a scalar discontinuos piecewise quadratic function
  H1Space::Element_t<m> f(h1Space2);
  interpolateGloballyWeak(f,makeWeakFunctionView([](auto const& cell, auto const& xi) { auto x = cell.geometry().global(xi); return Dune::FieldVector<double,m>(std::cos(2*x[0]+x[1]*x[1]/10)); }));
//   interpolateGloballyWeak(f,makeWeakFunctionView([](auto const& cell, auto const& xi) { return cell.geometry().global(xi);  }));
  
  writeVTK("original",f,IoOptions().setOrder(2).setDataMode(IoOptions::conforming));
  
  auto velocity = makeWeakFunctionView([](auto const& cell, auto xi) { auto x = cell.geometry().global(xi); auto tmp = x[0]; x[0] = -x[1]; x[1] = tmp; return x; }); // left turning circulation
  
  for (int i=0; i<5; ++i)
  {
    std::cerr << "\nStarting i=" << i << "\n==================\n";
    auto advectedView = makeAdvectedFunctionView(f,velocity,M_PI/2,1+(1<<i));
    
//     Dune::FieldVector<double,2> x(1);
//     auto cell = *findCell(leafView,x);
//     advectedView.value(cell,cell.geometry().local(x));
    
    H1Space::Element_t<m> g(h1Space2);
    interpolateGloballyWeak(g,advectedView);
    writeVTK("advected-"+paddedString(i),g,IoOptions().setOrder(order).setDataMode(IoOptions::conforming),"advected");
    if (m>1)
    {
      g -= f;
      writeVTK("displacement-"+paddedString(i),g,IoOptions().setOrder(order).setDataMode(IoOptions::conforming),"displacement");
    }
  }
  
  return 0;
}

#endif
