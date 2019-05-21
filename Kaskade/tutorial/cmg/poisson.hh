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

#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

extern int problemNo;

/// Example simple heat transfer equation
///

template <class RType, class VarSet>
class PoissonFunctional: public FunctionalBase<VariationalFunctional>
{
public:
  typedef RType  Scalar;
  typedef VarSet OriginVars;
  typedef VarSet AnsatzVars;
  typedef VarSet TestVars;

/// \class DomainCache
///

//StartSnippet1
  class DomainCache 
  {
  public:
    DomainCache(PoissonFunctional<RType,AnsatzVars> const&,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      data(vars_) 
    {}

    template <class Entity>
    void moveTo(Entity const &entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>
                          xglob = e->geometry().global(x);

      T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
      dT = at_c<0>(data.data).gradient(at_c<TIdx>(evaluators))[0];
      if (problemNo==5)
        f = -1000.0*exp(-1000.0*((xglob[0]-0.5)*(xglob[0]-0.5)+
                                (xglob[1]-0.5)*(xglob[1]-0.5)));
      else
        f = -1.0;
      if (problemNo==6)
        {
          if ((xglob[0]>0.3)&&(xglob[0]<0.7)&&
              (xglob[1]>0.3)&&(xglob[1]<0.7))
            lap = 1.0e6;
          else
            lap = 1.0;
        }
      else
        lap = 1.0;
    }

    Scalar
    d0() const 
    {
      return dT*dT/2-f*T;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return dT*arg.gradient[0] - f*arg.value; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1,
                               VariationalArg<Scalar,dim> const &arg2) const 
    {
      return lap*arg1.gradient[0]*arg2.gradient[0];
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Scalar T, f, lap;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> dT;
  };
//StopSnippet1

/// \class BoundaryCache
///

//StartSnippet2
  class BoundaryCache 
  {
  public:
    static const bool hasInteriorFaces = false;
    typedef typename AnsatzVars::Grid::template Codim<0>::
                     Entity::LeafIntersectionIterator FaceIterator;

    BoundaryCache(PoissonFunctional<RType,AnsatzVars> const&,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      data(vars_), e(0)
    {}

    void moveTo(FaceIterator const& entity)
    {
      e = &entity;
      penalty = 1.0e30;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                                      AnsatzVars::Grid::dimension-1>
                   const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      
      int const TIdx = result_of::value_at_c<typename AnsatzVars::Variables,
                                             0>::type::spaceIndex;
      Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                        AnsatzVars::Grid::dimension>
                          xglob = (*e)->geometry().global(x);

      switch (problemNo)
        {
      case 1:
		  if ((xglob[0]>0.5000000001)&&(xglob[0]<0.9999999999)&&
			  (xglob[1]>0.49)&&(xglob[1]<0.4999999999))
			{
			  isNeumann = true;
			}
		  else
			{
			  isNeumann = false;
			  T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
			  T0 = 0;
			}
		  break;
      case 2:
      case 4:
		  break;
      default:
		  isNeumann = false;
		  T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
		  T0 = 0;
		  break;
      case 3:
		  if ((xglob[0]>0.9999999999)&&(xglob[1]>0.03))
			{
			  isNeumann = false;
			  T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
			  T0 = 1000.00;
			}
		  else if ((xglob[0]>0.9999999999)&&(xglob[1]<-0.03))
			{
			  isNeumann = false;
			  T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
			  T0 = 0.0;
			}
		  else
			{
			  isNeumann = true;
			}
		  break;
      case 7:
		  isNeumann = false;
		  T = at_c<0>(data.data).value(at_c<TIdx>(evaluators));
		  if ((xglob[0]>0.5)&&(xglob[1]>0.5)&&(xglob[2]>0.5))
		    isNeumann = true;
		  else
		    T0 = 0;
		  break;
        }
    }

    Scalar
    d0() const 
    {
      return isNeumann?0.0:penalty*(T-T0)*(T-T0)/2;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return isNeumann?0.0:penalty*(T-T0)*arg.value[0];
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, 
                       VariationalArg<Scalar,dim> const &arg2) const 
    {
      return isNeumann?0.0:penalty*arg1.value*arg2.value;
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    FaceIterator const* e;
    Scalar penalty, T, T0;
    bool isNeumann;
  };
//StopSnippet1

/// \struct D2
///

//StartSnippet3
  template <int row>
  struct D1 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };
  
public:
  template <int row, int col>
  struct D2 
  {
    static bool const present = true;
    static bool const symmetric = true;
    static bool const lumped = false;
  };

/// \fn integrationOrder
///

  template <class Cell>
  int integrationOrder(Cell const& /* cell */,
                       int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder-1;
  }
  
public:
};

/// \example cmgtest.cpp
/// show the usage of PoissonFunctional
///
#endif
