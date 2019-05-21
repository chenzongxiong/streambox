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

#ifndef SPECIAL_QUADRATURE_Hh
#define SPECIAL_QUADRATURE_Hh

#include "dune/grid/common/quadraturerules.hh"

namespace Kaskade
{
  template<class GType, int dim>
  class TrapezoidalRule : public Dune::QuadratureRule<GType, dim>
  {
  public:
    TrapezoidalRule() : Dune::QuadratureRule<GType, dim>()
    {
      Dune::FieldVector<double, dim> pos;
      for(int i=0; i<dim+1; ++i)
      {
        for(int k=0; k<dim; ++k)
        {
          pos[k]=0.0;
          if(i==k) pos[k]=1.0;
        }
        double wgt=1.0/(dim+1);
        this->push_back(Dune::QuadraturePoint<GType, dim>(pos,wgt));
      }
    }
  };


  template<class GType, int dim>
  class NewtonCotesRule3 : public Dune::QuadratureRule<GType, dim>
  {
  public:
    NewtonCotesRule3() : Dune::QuadratureRule<GType, dim>()
    {
      double wgt(0.0);
      //    if(dim==1) wgt=1.0/3.0;
      //    if(dim==2) wgt=1.0/6.0;
      if(dim==3) wgt=-1.0/20.0;

      Dune::FieldVector<double, dim> pos;
      for(int i=0; i<dim+1; ++i)
      {
        for(int k=0; k<dim; ++k)
        {
          pos[k]=0.0;
          if(i==k) pos[k]=1.0;
        }
        this->push_back(Dune::QuadraturePoint<GType, dim>(pos,wgt));
      }

      if(dim==3) wgt=1.0/5.0;
      for(int i=0; i<dim; ++i)
      {
        for(int k=0; k<dim; ++k)
        {
          pos[k]=0.5;
          if(i==k) pos[k]=0.0;
        }
        Dune::QuadratureRule<GType, dim>::push_back(Dune::QuadraturePoint<GType, dim>(pos,wgt));
      }

      for(int i=0; i<dim; ++i)
      {
        for(int k=0; k<dim; ++k)
        {
          pos[k]=0.0;
          if(i==k) pos[k]=0.5;
        }
        Dune::QuadratureRule<GType, dim>::push_back(Dune::QuadraturePoint<GType, dim>(pos,wgt));
      }
    }
  };

  template<class GType, int dim>
  class RefinedTrapezoidal : public Dune::QuadratureRule<GType, dim>
  {
  public:
    RefinedTrapezoidal() : Dune::QuadratureRule<GType, dim>()
    {
      Dune::FieldVector<double, dim> pos;
      assert(dim==2);
      for(int ix=0; ix<2; ++ix)
      {
        for(int iy=0; iy<2-ix; ++iy)
        {
          pos[0]=0.5*ix;
          pos[1]=0.5*iy;
          double wgt=1.0/(dim+1)/4.0;
          if(ix==1 || iy == 1) wgt*=3;
          Dune::QuadratureRule<GType, dim>::push_back(Dune::QuadraturePoint<GType, dim>(pos,wgt));
        }
      }
    }
  };
} // end of namespace Kaskade

#endif
