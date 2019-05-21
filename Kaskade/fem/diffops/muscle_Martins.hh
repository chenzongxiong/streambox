#ifndef MUSCLE_MARTINS_HH
#define MUSCLE_MARTINS_HH

#include "neoHooke.hh"
#include "linalg/invariants.hh"
#include "utilities/elementary_functions.hh"
#include "utilities/functionTools.hh"

namespace Dune
{
  template <class,int> class FieldVector;
  template <class,int,int> class FieldMatrix;
}

namespace Kaskade
{
  /**
   * Material model for muscle tissue, from Martins et al.: A numerical model of passive active behavior of skeletal muscles, 1998
   */
//  template <int dim, class Direction, class VolumetricPenalty = void> class MuscleTissue_Martins;

//  template <int dim, class Direction>
//  class MuscleTissue_Martins<dim,Direction,void>
//      : public Sum< NeoHooke<dim>, Fung_OneDistributedFiber<dim,Direction> >
//  {
//    typedef NeoHooke<dim> CellFoam;
//    typedef Fung_OneDistributedFiber<dim,Direction> InterlobularSepta;
//    typedef Sum< CellFoam, InterlobularSepta > Base;
//  public:
//    /**
//     * @brief MuscleTissue_Martins
//     * @param cCells parameter for the neo-hookean contribution from the adipocytes (considered as cell foam)
//     * @param k1 stress like parameter for the interlobular septa
//     * @param k2 dimensionless parameter
//     * @param kappa fiber distribution parameter \f$ 0 \le \kappa \le \frac{1}{3} \f$
//     * @param d lightweight object that contains information on fiber directions in terms of a matrix that can be accessed via d0().
//     */
//    MuscleTissue_Martins(double cCells, double k1, double k2, double kappa, Direction d) : Base( CellFoam(cCells), InterlobularSepta(k1,k2,d,kappa))
//    {
//      assert( cCells > 0 && k1 > 0 && k2 > 0);
//    }

//    /**
//     * @brief MuscleTissue_Martins Mean values form the above mentioned publication, Table 2.
//     * @param d lightweight object that contains information on fiber directions in terms of a matrix that can be accessed via d0().
//     */
//    explicit MuscleTissue_Martins(Direction d) : MuscleTissue_Martins(0.15, 0.8, 47.3, 0.09, d)
//    {}

//    void updateDisplacementGradient(Dune::FieldMatrix<double,dim,dim> const& G)
//    {
//      Base::f.updateDisplacementGradient(G);
//      Base::g.updateDisplacementGradient(G);
//    }

//    void updateDeformationGradient(Dune::FieldMatrix<double,dim,dim> const& F)
//    {
//      Base::f.updateDeformationGradient(F);
//      Base::g.updateDeformationGradient(F);
//    }
//  };

  namespace MuscleTissue_Details
  {
    template <int dim, Invariant i>
    using IsotropicExpArg = Scaled< ShiftedInvariant< typename ChoosePrincipalInvariant<1,dim,i,CauchyGreenTensor<double,dim> >::type > >;

    template <int dim, class Direction, Invariant i>
    using FiberShiftedInvariant = ShiftedInvariant<typename ChooseMixedInvariant<3,dim,i,CauchyGreenTensor<double,dim>,Direction>::type>;

    template <int dim, class Direction, Invariant i>
    using FiberExpArg = Scaled< Squared< FiberShiftedInvariant<dim,Direction,i> > >;

    template <int dim, class Direction, Invariant i>
    using FiberModel = Sum<MatrixToScalarFunction<Exp,FiberExpArg<dim,Direction,i> >, ScalarConstant<Dune::FieldMatrix<double,dim,dim> > >;

    template <int dim, Invariant i>
    using IsotropicModel = Sum<MatrixToScalarFunction<Exp,IsotropicExpArg<dim,i> >, ScalarConstant<Dune::FieldMatrix<double,dim,dim> > >;
  }

  template <int dim, class Direction, class VolumetricPenalty, Invariant i = Invariant::Principal, class FiberDistribution = SimpleFiberDistribution<dim> >
  class MuscleTissue_Martins
  : public StrainBase<double,dim>,
    public Sum<
                Scaled<MuscleTissue_Details::IsotropicModel<dim,i> >,
                Scaled<MuscleTissue_Details::FiberModel<dim,Direction,i> >,
                MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> >
              >
  {
    typedef typename ChoosePrincipalInvariant<1,dim,i,CauchyGreenTensor<double,dim> >::type Inv1;
    typedef typename ChooseMixedInvariant<3,dim,i,CauchyGreenTensor<double,dim>,Direction>::type Inv6;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef ShiftedInvariant<Inv6> SInv6;
//    typedef Sum<MatrixToScalarFunction<Exp,MuscleTissue_Details::IsotropicExpArg<dim,Direction,i> >,  ScalarConstant<Dune::FieldMatrix<double,dim,dim> > > IsotropicModel;
//    typedef Sum<MatrixToScalarFunction<Exp,MuscleTissue_Details::FiberExpArg<dim,Direction,i> >, ScalarConstant<Dune::FieldMatrix<double,dim,dim> > > FiberModel;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det> Penalty;
    typedef Sum<Scaled<MuscleTissue_Details::IsotropicModel<dim,i> >,Scaled<MuscleTissue_Details::FiberModel<dim,Direction,i> >,Penalty> Base;
    using StrainBase<double,dim>::S;
    using StrainBase<double,dim>::F;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    MuscleTissue_Martins(double c0iso, double c1iso, double c0fiber, double c1fiber, double d1, double d2, FiberDistribution fdist, Direction d)
    : StrainBase<double,dim>(),
      Base( Scaled<MuscleTissue_Details::IsotropicModel<dim,i> >(c0iso, MuscleTissue_Details::IsotropicModel<dim,i>( MatrixToScalarFunction<Exp,MuscleTissue_Details::IsotropicExpArg<dim,i> >(
                                                            Exp(),
                                                            MuscleTissue_Details::IsotropicExpArg<dim,i>(c1iso, SInv1(Inv1(S),dim) )
                                                          ),
                                                           ScalarConstant<Dune::FieldMatrix<double,dim,dim> >(-1)
                                                        )
                                   ),
            Scaled<MuscleTissue_Details::FiberModel<dim,Direction,i> >(c0fiber, MuscleTissue_Details::FiberModel<dim,Direction,i>( MatrixToScalarFunction<Exp,MuscleTissue_Details::FiberExpArg<dim,Direction,i> >(
                                                      Exp(),
                                                      MuscleTissue_Details::FiberExpArg<dim,Direction,i>(c1iso,Squared<MuscleTissue_Details::FiberShiftedInvariant<dim,Direction,i> >
                                                                                                         (
                                                                                                            MuscleTissue_Details::FiberShiftedInvariant<dim,Direction,i>(Inv6(S,d),1)
                                                                                                         ) )
                                                    ),
                                                     ScalarConstant<Dune::FieldMatrix<double,dim,dim> >(-1)
                                                  )
                              ),
            Penalty(Det(F),VolumetricPenalty(d1,d2))
          )
    { assert(c>0 && d>0); }

    MuscleTissue_Martins(Direction d, double d1, double d2) :
      MuscleTissue_Martins(0.387, 23.46, 0.584, 12.43, d1, d2, FiberDistribution(0.), d)
    {}


    MuscleTissue_Martins(MuscleTissue_Martins const&) = default;
    MuscleTissue_Martins& operator=(MuscleTissue_Martins const&) = default;
  };


/*  template <int dim, class Direction, class VolumetricPenalty, class FiberDistribution = SimpleFiberDistribution>
  class MuscleTissue_Martins
      : public Sum< NeoHooke<dim,Invariant::Principal,void,VolumetricPenalty>, Fung_OneDistributedFiber<dim,Direction> >
  {
    typedef NeoHooke<dim,Invariant::Principal,void,VolumetricPenalty> CellFoam;
    typedef Fung_OneDistributedFiber<dim,Direction> InterlobularSepta;
    typedef Sum< CellFoam, InterlobularSepta> Base;
  public:*/
    /**
     * @brief MuscleTissue_Martins
     * @param cCells parameter for the neo-hookean contribution from the adipocytes (considered as cell foam)
     * @param k1 stress like parameter for the interlobular septa
     * @param k2 dimensionless parameter
     * @param kappa fiber distribution parameter \f$ 0 \le \kappa \le \frac{1}{3} \f$
     * @param d lightweight object that contains information on fiber directions  in terms of a rank-one matrix that can be accessed via d0().
     * @param d1 parameter for volumetric penalty function corresponding to inflation
     * @param d2 parameter for volumetric penalty function corresponding to compression
     */
//    MuscleTissue_Martins(double cCells, double k1, double k2, double kappa, Direction d, double d1, double d2) : Base( CellFoam(cCells,d1,d2), InterlobularSepta(k1,k2,d,kappa))
//    {
//      assert( cCells > 0 && k1 > 0 && k2 > 0);
//    }

    /**
     * @brief MuscleTissue_Martins Mean values form the above mentioned publication, Table 2.
     * @param d lightweight object that contains information on fiber directions in terms of a rank-one matrix that can be accessed via d0().
     * @param d1 parameter for volumetric penalty function corresponding to inflation
     * @param d2 parameter for volumetric penalty function corresponding to compression
     */
/*    MuscleTissue_Martins(Direction d, double d1, double d2) : MuscleTissue_Martins(0.15, 0.8, 47.3, 0.09, d, d1, d2)
    {}


    void updateDisplacementGradient(Dune::FieldMatrix<double,dim,dim> const& G)
    {
      std::cout << "update displacement gradient: " << G << std::endl;
      Base::f.updateDisplacementGradient(G);
      Base::g.updateDisplacementGradient(G);
      std::cout << "f: " << Base::f.F << std::endl;
    }

    void updateDeformationGradient(Dune::FieldMatrix<double,dim,dim> const& F)
    {
      Base::f.updateDeformationGradient(F);
      Base::g.updateDeformationGradient(F);
    }
  };
*/
}

#endif // ADIPOSE_SOMMERHOLZAPFEL_HH
