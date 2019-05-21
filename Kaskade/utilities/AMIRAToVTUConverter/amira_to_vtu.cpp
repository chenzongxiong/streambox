/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "utilities/kaskopt.hh"
#include "utilities/enums.hh"

#include <iostream>
#include <string>

#include "dune/config.h"
#define HAVE_UG 1

#include "dune/grid/uggrid.hh"

#include "fem/lagrangespace.hh"

#include "io/amirameshreader.hh"
#include "io/vtk.hh"



int main(int argc, char *argv[])
{
  using namespace boost::fusion;

  int verbosity = 1;
  bool dump = false;
  int const dim = 3;
  std::auto_ptr<boost::property_tree::ptree> pt(getKaskadeOptions(argc, argv, verbosity, dump));

  std::string const gridfile = getParameter(*pt, "gridfile", std::string("/datanumerik/bzflubko/Gitter/implant.grid"));
  int const order = getParameter(*pt, "order", 0);
  std::string names[1] = { "name" };

  // define grid
  typedef Dune::UGGrid<dim> Grid;
  typedef Grid::LeafGridView LeafView;

  // define function space
  typedef FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> > Space;
  // define vector holding spaces
  typedef vector<Space const*> Spaces;
  // specify variable details
  typedef vector<VariableDescription<0,1,0> > VariableDescriptions;
  typedef VariableSetDescription<Spaces,VariableDescriptions> VariableSet;

  GridManager<Grid> gridManager( AmiraMeshReader::readGrid<Grid,float>(gridfile) );
  Space space(gridManager,gridManager.grid().leafView(),order);
  Spaces spaces(&space);
  VariableSet variableSet(spaces,names);
  VariableSet::VariableSet x(variableSet);

  IoOptions options;
  options.outputType = IoOptions::ascii;

  std::string const ending_1(".grid"), ending_2(".am"), new_ending("_converted");
  std::string name(gridfile);

  if(name.find(ending_1) != std::string::npos) name.replace(name.find(ending_1), ending_1.length(), new_ending);
  if(name.find(ending_2) != std::string::npos) name.replace(name.find(ending_2), ending_2.length(), new_ending);

  writeVTKFile(gridManager.grid().leafView(),variableSet,x,name,options);

  std::cout << "Convertion finished." << std::endl;

  return 0;
}
