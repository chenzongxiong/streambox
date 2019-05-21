/*
 * neoHooke.hh
 *
 *  Created on: 20.06.2014
 *      Author: bzflubko
 */

#ifndef NEOHOOKE_HH_
#define NEOHOOKE_HH_

#include "linalg/determinant.hh"
#include "linalg/invariants.hh"
#include "utilities/functionTools.hh"

namespace Dune
{
  template <class,int,int> class FieldMatrix;
}

namespace Kaskade
{
  template <int dim, Invariant i = Invariant::Principal, class Direction = void, class VolumetricPenalty = void> class NeoHooke;

  template <int dim, Invariant i, class Direction, class VolumetricPenalty>
  class NeoHooke
  : public StrainBase<double,dim>, 
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> >,
               MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> > >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det> Penalty;
    typedef Sum<Scaled<SInv1>,Penalty> Base;
    using StrainBase<double,dim>::F;
    using StrainBase<double,dim>::S;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    NeoHooke(double c, double d1, double d2, Direction const& dir)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c,SInv1(Inv1(S,dir),dim)), Penalty(Det(F),VolumetricPenalty(d1,d2)))
    { assert(c>0 && d1>0 && d2>0); }

    NeoHooke(NeoHooke const&) = default;
    NeoHooke& operator=(NeoHooke const&) = default;
  };


  template <int dim, Invariant i, class VolumetricPenalty>
  class NeoHooke<dim,i,void,VolumetricPenalty>
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type> >,
               MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> > >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det> Penalty;
    typedef Sum<Scaled<SInv1>,Penalty> Base;
    using StrainBase<double,dim>::S;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;
    using StrainBase<double,dim>::F;

    NeoHooke(double c, double d1, double d2)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c,SInv1(Inv1(S),dim)), Penalty(Det(F),VolumetricPenalty(d1,d2)))
    { assert(c>0 && d1>0 && d2>0); }

    void updateDisplacementGradient(Argument const& arg)
    {
      StrainBase<double,dim>::updateDisplacementGradient(arg);
    }

    NeoHooke(NeoHooke const&) = default;
    NeoHooke& operator=(NeoHooke const&) = default;
  };


  template <int dim, Invariant i>
  class NeoHooke<dim,i,void,void>
  : public StrainBase<double,dim>, 
    public Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type> >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type Inv;
    typedef ShiftedInvariant<Inv> SInv;
    typedef Scaled<SInv> Base;
    using StrainBase<double,dim>::S;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    explicit NeoHooke(double c) : StrainBase<double,dim>(), Base(c,SInv(Inv(S),dim)) { assert(c>0); }

    NeoHooke(NeoHooke const&) = default;
    NeoHooke& operator=(NeoHooke const&) = default;
  };

  template <int dim, Invariant i, class Direction>
  class NeoHooke<dim,i,Direction,void>
  : public StrainBase<double,dim>,
    public Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv;
    typedef ShiftedInvariant<Inv> SInv;
    typedef Scaled<SInv> Base;
    using StrainBase<double,dim>::S;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    NeoHooke(double c, Direction const& d) : StrainBase<double,dim>(), Base(c,SInv(Inv(S,d),dim)) { assert(c>0); }

    NeoHooke(NeoHooke const&) = default;
    NeoHooke& operator=(NeoHooke const&) = default;
  };

}


#endif /* NEOHOOKE_HH_ */
