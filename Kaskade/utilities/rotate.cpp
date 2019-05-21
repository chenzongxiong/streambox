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

#include <array>
#include <cmath>

#include "rotate.hh"
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
namespace Kaskade
{
/**
 * \brief Offers functionalities to rotate a 2D object (i.e. after reading in the nodes, before
 * handing it over to a gridManager) and translate it back to its original upper y-position.
 */
void rotate(double alpha, std::array<double,2> *barycenter, std::vector<Dune::FieldVector<double,2> > & nodes){

    if(alpha==0)
        return;

    alpha *= M_PI/180;

    if(barycenter==nullptr){
        throw std::runtime_error("barycenter in rotate is a nullptr");
    }

    for (int i=0; i<nodes.size(); ++i){
        double foo = std::get<0>(*barycenter) + cos(alpha)*((nodes[i])[0] - std::get<0>(*barycenter)) - sin(alpha)*((nodes[i])[1] - std::get<1>(*barycenter));
        (nodes[i])[1] = std::get<1>(*barycenter) + sin(alpha)*((nodes[i])[0] - std::get<0>(*barycenter)) + cos(alpha)*((nodes[i])[1] - std::get<1>(*barycenter));
        (nodes[i])[0] = foo;
    }
}

void find2DBarycenter(std::vector<Dune::FieldVector<double,2> > const & nodes, std::array<double,2>* barycenter){


    double x_max = std::numeric_limits<double>::lowest();
    double x_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::lowest();
    double y_min = std::numeric_limits<double>::max();

    for (int i=0; i<nodes.size(); ++i)
    {
        x_max = std::max(x_max, (nodes[i])[0]);
        x_min = std::min(x_min, (nodes[i])[0]);
        y_max = std::max(y_max, (nodes[i])[1]);
        y_min = std::min(y_min, (nodes[i])[1]);
    }

    if(barycenter==nullptr){
        throw std::runtime_error("barycenter in find2DBarycenter is a nullptr");
    }

    std::get<0>(*barycenter) = (x_max+x_min)/2;
    std::get<1>(*barycenter) = (y_max+y_min)/2;
}

void translate(double deviation, std::vector<Dune::FieldVector<double,2> > & nodes){
    for (int i = 0; i<nodes.size(); i++){
        (nodes[i])[1] -= deviation;
    }
}

}
