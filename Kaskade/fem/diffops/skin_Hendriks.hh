#ifndef SKIN_HENDRIKS_HH
#define SKIN_HENDRIKS_HH

#include "linalg/invariants.hh"
#include "utilities/functionTools.hh"

namespace Kaskade
{
  /**
   * Extended Mooney-Rivlin type model according to Hendriks, Brokken, Oomens eta l.: Influence of hydration and experimental length scale on the mechanical response of human skin in vivo using optical coherence tomography, 2004
   * and Hendriks, Brokken, Oomens et al.: The relative contributions of different skin layers to the mechanical behavior of human skin in vivo using suction experiments, 2006
   */
  template <int dim, class VolumetricPenalty = void> class Skin_Hendriks;

  template <int dim, class VolumetricPenalty>
  class Skin_Hendriks
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type> >,
               Scaled<
                      Product<
                              ShiftedInvariant<typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type>,
                              ShiftedInvariant<typename ChooseInvariant<2,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type>
                              >
                      >,
               MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> > >
  {
    typedef typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef typename ChooseInvariant<2,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type Inv2;
    typedef ShiftedInvariant<Inv2> SInv2;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det>Penalty;
    typedef Sum<Scaled<SInv1>,Scaled<Product<SInv1,SInv2> >,Penalty> Base;
    using StrainBase<double,dim>::F;
    using StrainBase<double,dim>::S;
  public:
    typedef double Scalar;

    Skin_Hendriks(double c0, double c1, double c2, double c3)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S),dim)), Scaled< Product<SInv1,SInv2> >(c1,Product<SInv1,SInv2>( SInv1(Inv1(S),dim), SInv2(Inv2(S),dim) )), Penalty(Det(F),VolumetricPenalty(c2,c3)))
    { assert(c0>0 && c1>0 && c2>0); }

    /**
     * @brief Constructor using parameters from the above mentioned literature for human skin.
     * @param c2 parameter for the volumetric part related to inflation
     * @param c3 parameter for the volumetric part related to compression
     */
    Skin_Hendriks(double c2, double c3)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(9.4,SInv1(Inv1(S),dim)), Scaled< Product<SInv1,SInv2> >(82,Product<SInv1,SInv2>( SInv1(Inv1(S),dim), SInv2(Inv2(S),dim) )), Penalty(Det(F),VolumetricPenalty(c2,c3)))
    { assert(c2>0); }


    Skin_Hendriks(Skin_Hendriks const&) = default;
    Skin_Hendriks& operator=(Skin_Hendriks const&) = default;
  };
  
  template <int dim>
  class Skin_Hendriks<dim,void>
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type> >,
               Scaled<
                      Product<
                              ShiftedInvariant<typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type>,
                              ShiftedInvariant<typename ChooseInvariant<2,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type>
                              >
                      >
              >
  {
    typedef typename ChooseInvariant<1,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef typename ChooseInvariant<2,dim,Invariant::Principal,CauchyGreenTensor<double,dim>,void>::type Inv2;
    typedef ShiftedInvariant<Inv2> SInv2;
    typedef Sum<Scaled<SInv1>,Scaled<Product<SInv1,SInv2> > > Base;
    using StrainBase<double,dim>::F;
    using StrainBase<double,dim>::S;
  public:
    typedef double Scalar;

    Skin_Hendriks(double c0, double c1, double c2, double c3)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S),dim)), Scaled< Product<SInv1,SInv2> >(c1,Product<SInv1,SInv2>( SInv1(Inv1(S),dim), SInv2(Inv2(S),dim) )))
    { assert(c0>0 && c1>0 && c2>0); }

    /**
     * @brief Constructor using parameters from the above mentioned literature for human skin.
     * @param c2 parameter for the volumetric part related to inflation
     * @param c3 parameter for the volumetric part related to compression
     */
    Skin_Hendriks()
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(9.4,SInv1(Inv1(S),dim)), Scaled< Product<SInv1,SInv2> >(82,Product<SInv1,SInv2>( SInv1(Inv1(S),dim), SInv2(Inv2(S),dim) )))
    {}

    Skin_Hendriks(Skin_Hendriks const&) = default;
    Skin_Hendriks& operator=(Skin_Hendriks const&) = default;
  };
}

#endif
