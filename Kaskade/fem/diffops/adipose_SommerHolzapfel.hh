#ifndef ADIPOSE_SOMMERHOLZAPFEL_HH
#define ADIPOSE_SOMMERHOLZAPFEL_HH

#include "neoHooke.hh"
#include "fung_OneDistributedFiber.hh"
#include "utilities/functionTools.hh"

namespace Dune
{
  template <class,int> class FieldVector;
  template <class,int,int> class FieldMatrix;
}

namespace Kaskade
{
  /**
   * Material model for adipose tissue, from Sommer et al.: Multiaxial mechanical properties and constitutive modeling of human adipose tissue: A basis for preoperative simulations
   * in plastic and reconstructive surgery, 2013
   */
//  template <int dim, class Direction, class VolumetricPenalty = void> class AdiposeTissue_SommerHolzapfel;

//  template <int dim, class Direction>
//  class AdiposeTissue_SommerHolzapfel<dim,Direction,void>
//      : public Sum< NeoHooke<dim>, Fung_OneDistributedFiber<dim,Direction> >
//  {
//    typedef NeoHooke<dim> CellFoam;
//    typedef Fung_OneDistributedFiber<dim,Direction> InterlobularSepta;
//    typedef Sum< CellFoam, InterlobularSepta > Base;
//  public:
//    /**
//     * @brief AdiposeTissue_SommerHolzapfel
//     * @param cCells parameter for the neo-hookean contribution from the adipocytes (considered as cell foam)
//     * @param k1 stress like parameter for the interlobular septa
//     * @param k2 dimensionless parameter
//     * @param kappa fiber distribution parameter \f$ 0 \le \kappa \le \frac{1}{3} \f$
//     * @param d lightweight object that contains information on fiber directions in terms of a matrix that can be accessed via d0().
//     */
//    AdiposeTissue_SommerHolzapfel(double cCells, double k1, double k2, double kappa, Direction d) : Base( CellFoam(cCells), InterlobularSepta(k1,k2,d,kappa))
//    {
//      assert( cCells > 0 && k1 > 0 && k2 > 0);
//    }

//    /**
//     * @brief AdiposeTissue_SommerHolzapfel Mean values form the above mentioned publication, Table 2.
//     * @param d lightweight object that contains information on fiber directions in terms of a matrix that can be accessed via d0().
//     */
//    explicit AdiposeTissue_SommerHolzapfel(Direction d) : AdiposeTissue_SommerHolzapfel(0.15, 0.8, 47.3, 0.09, d)
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

  template <int dim, class Direction, class VolumetricPenalty, Invariant i = Invariant::Principal, class FiberDistribution = SimpleFiberDistribution<dim> >
  class AdiposeTissue_SommerHolzapfel
  : public StrainBase<double,dim>,
    public Sum<Scaled<ShiftedInvariant<typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type> >,
               MultiplyWithIndicator< FODF::Base<dim,i,Direction>, FODF::Inv4<dim,i,Direction> >,
               MatrixToScalarFunction<VolumetricPenalty,Determinant<dim> > >
  {
    typedef typename ChooseInvariant<1,dim,i,CauchyGreenTensor<double,dim>,void>::type Inv1;
    typedef ShiftedInvariant<Inv1> SInv1;
    typedef Determinant<dim> Det;
    typedef MatrixToScalarFunction<VolumetricPenalty,Det> Penalty;
    typedef Sum< Scaled<SInv1>, MultiplyWithIndicator< FODF::Base<dim,i,Direction>, FODF::Inv4<dim,i,Direction> >, Penalty > Base;
    using StrainBase<double,dim>::S;
    using StrainBase<double,dim>::F;
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;

    AdiposeTissue_SommerHolzapfel(double cCells, double k1, double k2, FiberDistribution fdist, Direction d, double d1, double d2)
    : StrainBase<double,dim>(),
      Base( Scaled<SInv1>(cCells,SInv1(Inv1(S),dim)),
            MultiplyWithIndicator< FODF::Base<dim,i,Direction>, FODF::Inv4<dim,i,Direction> >
            (
              FODF::Base<dim,i,Direction>
              (
                k1/k2 ,
                FODF::ExpFun<dim,i,Direction>
                (
                  FODF::ExpFun1<dim,i,Direction>
                  (
                    FODF::ExpArg<dim,i,Direction>
                    (
                      k2,
                      FODF::initProduct<dim,i>(S, d, fdist)
                    )
                  ),
                  Constant< FODF::M<dim> >(-1)
                )
              ),
              FODF::Inv4<dim,i,Direction>(S,d)
            ),
            Penalty(Det(F),VolumetricPenalty(d1,d2)))
    { assert(cCells>0 && d1>0 && d2>0); }

    AdiposeTissue_SommerHolzapfel(double cCells, double k1, double k2, double kappa, Direction d, double d1, double d2)
    : AdiposeTissue_SommerHolzapfel(cCells, k1, k2, FiberDistribution(kappa), d, d1, d2)
    {}

    AdiposeTissue_SommerHolzapfel(Direction d, double d1, double d2) :
      AdiposeTissue_SommerHolzapfel(0.15, 0.8, 47.3, FiberDistribution(0.09), d, d1, d2)
    {}


    AdiposeTissue_SommerHolzapfel(AdiposeTissue_SommerHolzapfel const&) = default;
    AdiposeTissue_SommerHolzapfel& operator=(AdiposeTissue_SommerHolzapfel const&) = default;
  };



/*  template <int dim, class Direction, class VolumetricPenalty, class FiberDistribution = SimpleFiberDistribution>
  class AdiposeTissue_SommerHolzapfel
      : public Sum< NeoHooke<dim,Invariant::Principal,void,VolumetricPenalty>, Fung_OneDistributedFiber<dim,Direction> >
  {
    typedef NeoHooke<dim,Invariant::Principal,void,VolumetricPenalty> CellFoam;
    typedef Fung_OneDistributedFiber<dim,Direction> InterlobularSepta;
    typedef Sum< CellFoam, InterlobularSepta> Base;
  public:*/
    /**
     * @brief AdiposeTissue_SommerHolzapfel
     * @param cCells parameter for the neo-hookean contribution from the adipocytes (considered as cell foam)
     * @param k1 stress like parameter for the interlobular septa
     * @param k2 dimensionless parameter
     * @param kappa fiber distribution parameter \f$ 0 \le \kappa \le \frac{1}{3} \f$
     * @param d lightweight object that contains information on fiber directions  in terms of a rank-one matrix that can be accessed via d0().
     * @param d1 parameter for volumetric penalty function corresponding to inflation
     * @param d2 parameter for volumetric penalty function corresponding to compression
     */
//    AdiposeTissue_SommerHolzapfel(double cCells, double k1, double k2, double kappa, Direction d, double d1, double d2) : Base( CellFoam(cCells,d1,d2), InterlobularSepta(k1,k2,d,kappa))
//    {
//      assert( cCells > 0 && k1 > 0 && k2 > 0);
//    }

    /**
     * @brief AdiposeTissue_SommerHolzapfel Mean values form the above mentioned publication, Table 2.
     * @param d lightweight object that contains information on fiber directions in terms of a rank-one matrix that can be accessed via d0().
     * @param d1 parameter for volumetric penalty function corresponding to inflation
     * @param d2 parameter for volumetric penalty function corresponding to compression
     */
/*    AdiposeTissue_SommerHolzapfel(Direction d, double d1, double d2) : AdiposeTissue_SommerHolzapfel(0.15, 0.8, 47.3, 0.09, d, d1, d2)
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
