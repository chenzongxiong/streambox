#ifndef MOONEY_RIVLIN_HH
#define MOONEY_RIVLIN_HH

#include "linalg/invariants.hh"
#include "utilities/functionTools.hh"

namespace Kaskade
{
  template <int dim, class VolumetricPenalty = void, Invariant i = Invariant::Principal, class Direction = void> class MooneyRivlin;

  /**
   * \param lambda first Lame constant
   * \param mu second Lame constant
   */
  template <class VolumetricPenalty, int dim, Invariant i = Invariant::Principal, class Direction = void>
  MooneyRivlin<dim,VolumetricPenalty,i,Direction> createMooneyRivlinFromLameConstants(double lambda, double mu)
  {
    static_assert(!std::is_same<VolumetricPenalty,void>::value,"not implemented");
    VolumetricPenalty g;
//    std::cout << "g': " << g.d1(1) << ", g'': " << g.d2(1) << std::endl;
    //double c = (lambda + 2*mu)/(g.d2(1)-g.d1(1));
    //double b = -0.5*mu - 0.5*c*g.d1(1);
    //double a = mu + 0.5*c*g.d1(1);

    double rho = g.d1(1)/(-g.d2(1)+g.d1(1));
    double d = (lambda+2.0*mu)/(g.d2(1)-g.d1(1));
    double c = (0.5*rho-0.25)*mu+0.25*rho*lambda;
    if(c > 0.25*mu) c = (rho-0.75)*mu+0.5*rho*lambda;
    double b = -mu + rho*(lambda+2.*mu)-2.*c;
    double a = b + mu;
    double alpha = 0.5*a - b;
    double beta = 0.5*b;
    //    std::cout << "alpha: " << alpha << ", beta: " << beta << ", c: " << c << ", d: " << d << std::endl;

    if(a<0 || b<0 || c<0)
    {
      std::cout << "computed parameters: " << a << ", " << b << ", " << c << ", " << d << std::endl;
      std::cout << "alpha=" << alpha << ", beta=" << beta << std::endl;
      std::cout << "material law is not polyconvex" << std::endl;
      exit(1);
    }

    return MooneyRivlin<dim,VolumetricPenalty,i,Direction>(alpha,beta,c,d);
  }


    /**
   * \param E Young's modulus
   * \param nu Poisson ratio
   */
  template <class VolumetricPenalty, int dim, Invariant i = Invariant::Principal, class Direction = void>
  MooneyRivlin<dim,VolumetricPenalty,i,Direction> createMooneyRivlinFromMaterialConstants(double E, double nu)
  {
    static_assert(!std::is_same<VolumetricPenalty,void>::value,"not implemented");
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));
    return createMooneyRivlinFromLameConstants<VolumetricPenalty,dim,i,Direction>(lambda,mu);
  }


  template <int dim, class VolumetricPenalty, Invariant i, class Direction>
  class MooneyRivlin
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> >,
               Scaled<ShiftedInvariant<typename ChooseInvariant<2,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> >,
               MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> > >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef typename ChooseInvariant<2,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv2;
    typedef ShiftedInvariant<Inv2> SInv2;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det>Penalty;
    typedef Sum<Scaled<SInv1>,Scaled<SInv2>,Penalty> Base;
    using StrainBase<double,dim>::F;
    using StrainBase<double,dim>::S;
  public:
    typedef double Scalar;

    MooneyRivlin(double lambda, double mu) : MooneyRivlin(createMooneyRivlinFromLameConstants<VolumetricPenalty,dim>(lambda,mu)) {}

//    template <class Scalar, class enable = typename std::enable_if<std::is_same<Scalar,double>::value && std::is_same<Direction,void>::value,void>::type>
    MooneyRivlin(double c0, double c1, double c2, double c3)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S),dim)), Scaled<SInv2>(c1,SInv2(Inv2(S),dim)), Penalty(Det(F),VolumetricPenalty(c2,c3)))
    { assert(c0>0 && c1>0 && c2>0); }

    template <class Dir, class enable = typename std::enable_if< std::is_same<Direction,Dir>::value && !std::is_same<void,Direction>::value,void>::type>
    MooneyRivlin(double c0, double c1, double c2, double c3, Dir d)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S,d),dim)), Scaled<SInv2>(c1,SInv2(Inv2(S,d),dim)), Penalty(Det(F),VolumetricPenalty(c2,c3)))
    { assert(c0>0 && c1>0 && c2>0); }

    MooneyRivlin(MooneyRivlin const&) = default;
    MooneyRivlin& operator=(MooneyRivlin const&) = default;
  };

  
  template <int dim, Invariant i,class Direction>
  class MooneyRivlin<dim,void,i,Direction>
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> >,
               Scaled<ShiftedInvariant<typename ChooseInvariant<2,dim,i,CauchyGreenTensor<double,dim>,Direction>::type> > >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv1;
    typedef typename ChooseInvariant<2,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv2;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef ShiftedInvariant<Inv2> SInv2;
    typedef Sum<Scaled<SInv1>,Scaled<SInv2> > Base;
    using StrainBase<double,dim>::S;
  public:
    typedef double Scalar;

    template <class Scalar, class enable = typename std::enable_if<std::is_same<Scalar,double>::value && std::is_same<Direction,void>::value,void>::type>
    MooneyRivlin(double c0, double c1)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S),dim)),Scaled<SInv2>(c1,SInv2(Inv2(S),dim)))
    { assert(c0>0 && c1>0); }

    template <class Scalar, class enable = typename std::enable_if<std::is_same<Scalar,double>::value && !std::is_same<Direction,void>::value,void>::type>
    MooneyRivlin(double c0, double c1, Direction const& d)
    : StrainBase<double,dim>(), Base(Scaled<SInv1>(c0,SInv1(Inv1(S,d),dim)),Scaled<SInv2>(c1,SInv2(Inv2(S,d),dim)))
    { assert(c0>0 && c1>0); }

    MooneyRivlin(MooneyRivlin const&) = default;
    MooneyRivlin& operator=(MooneyRivlin const&) = default;
  };
}

#endif
