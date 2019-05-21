/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef UTILITIES_ROTATE_HH
#define UTILITIES_ROTATE_HH

#include <array>
#include <vector>
#include <cmath>
#include <dune/common/fvector.hh>
//#include "fem/forEach.hh"
//#include "fem/gridmanager.hh"

namespace Kaskade 
{

/**
   * \brief Rotates a 2D object around it's center (and restores the upper boundary)
   *
   */

void rotate(double alpha, std::array<double,2> *barycenter, std::vector<Dune::FieldVector<double,2> > & nodes);

void find2DBarycenter(std::vector<Dune::FieldVector<double,2> > const & nodes, std::array<double,2>* barycenter);

void translate(double deviation, std::vector<Dune::FieldVector<double,2> > & nodes);

//template<class Grid>
//void rotate(double alpha, std::array<double,2> *barycenter, GridManager<Grid> * gridManager){

//    if(alpha==0)
//        return;

//    const int dim = Grid::dimension;

//    alpha = alpha*M_PI/180;

//    typedef typename Grid::template Codim<Grid::dimension>::LeafIterator VertexIterator;
//    VertexIterator vIt    = gridManager->grid().leafGridView().template begin<dim>();
//    VertexIterator vEndIt = gridManager->grid().leafGridView().template end<dim>();

//    for (; vIt!=vEndIt; ++vIt){
//        Dune::FieldVector<double,dim> foo = Dune::FieldVector<double,dim>();
//        foo[0] = std::get<0>(*barycenter) + cos(alpha)*(vIt->geometry().corner(0)[0] - std::get<0>(*barycenter)) - sin(alpha)*(vIt->geometry().corner(0)[1] - std::get<1>(*barycenter));
//        foo[1] = std::get<1>(*barycenter) + sin(alpha)*(vIt->geometry().corner(0)[0] - std::get<0>(*barycenter)) + cos(alpha)*(vIt->geometry().corner(0)[1] - std::get<1>(*barycenter));

//        typename Grid::template Codim<dim>::EntityPointer ci(vIt);
//        gridManager->grid_non_const().setPosition(ci, foo);
//    }
//}

//template <class Grid>
//void find2DBarycenter(GridManager<Grid> * gridManager, std::array<double,2>* barycenter){

//    typedef typename Grid::template Codim<Grid::dimension>::Entity Vertex;

//    double x_max = std::numeric_limits<double>::lowest();
//    double x_min = std::numeric_limits<double>::max();
//    double y_max = std::numeric_limits<double>::lowest();
//    double y_min = std::numeric_limits<double>::max();

//    forEachVertex(gridManager->grid().leafGridView(), [&x_max,&x_min,&y_max,&y_min](Vertex const& v)
//    {
//        x_max = std::max(x_max, v.geometry().corner(0)[0]);
//        x_min = std::min(x_min, v.geometry().corner(0)[0]);
//        y_max = std::max(y_max, v.geometry().corner(0)[1]);
//        y_min = std::min(y_min, v.geometry().corner(0)[1]);
//    });

//    std::get<0>(*barycenter) = (x_max+x_min)/2;
//    std::get<1>(*barycenter) = (y_max+y_min)/2;
//}


//template <class Grid>
//void translate(double deviation, GridManager<Grid> * gridManager){

//    const int dim = Grid::dimension;

//    typedef typename Grid::template Codim<dim>::LeafIterator VertexIterator;
//    VertexIterator vIt    = gridManager->grid().leafGridView().template begin<dim>();
//    VertexIterator vEndIt = gridManager->grid().leafGridView().template end<dim>();

//    for (; vIt!=vEndIt; ++vIt){
//        Dune::FieldVector<double,dim> foo = Dune::FieldVector<double,dim>();
//        foo[0] = vIt->geometry().corner(0)[0];
//        foo[1] = vIt->geometry().corner(0)[1] - deviation;
//        typename Grid::template Codim<dim>::EntityPointer ci(vIt);
//        gridManager->grid_non_const().setPosition(ci, foo);
//    }
//}

}

#endif

