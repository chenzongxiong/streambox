#ifndef MATRIX_INVARIANTS
#define MATRIX_INVARIANTS

#include <cmath>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "linalg/cofactorMatrix.hh"
#include "linalg/determinant.hh"
#include "fem/fixdune.hh"
#include "linalg/tensorProduct.hh"
#include "linalg/wrappedMatrix.hh"
#include "utilities/functionTools.hh"
#include "utilities/geometry/geomtools.hh"

namespace Kaskade
{
  // forward declarations
  template <size_t i, class Scalar, size_t n, class Source = WrappedMatrix<Scalar,n> > class PrincipalInvariant;
  template <size_t i, class Scalar, size_t n, class Source = WrappedMatrix<Scalar,n> > class ModifiedPrincipalInvariant;
  template <size_t i, class Scalar, size_t n, class Source = WrappedMatrix<Scalar,n>, class Direction = WrappedMatrix<Scalar,n> > class MixedInvariant;
  template <size_t i, class Scalar, size_t n, class Source = WrappedMatrix<Scalar,n>, class Direction = WrappedMatrix<Scalar,n> > class ModifiedMixedInvariant;

  enum class Invariant { Principal, Mixed, ModifiedPrincipal, ModifiedMixed, Modified };

  namespace Invariant_Details
  {
    template <class T> struct IsSquareFieldMatrix{ static constexpr bool value = false; };

    template <class Scalar, int n>
    struct IsSquareFieldMatrix<Dune::FieldMatrix<Scalar,n,n> >{ static constexpr bool value = true; };

    template <class T>
    constexpr bool isSquareFieldMatrix() { return IsSquareFieldMatrix<T>::value; }
  }

    struct OverThirdRoot
    {
      OverThirdRoot() = default;
      
      explicit OverThirdRoot(double c0) : c(c0) {}
      
      OverThirdRoot(OverThirdRoot const&) = default;
      OverThirdRoot& operator=(OverThirdRoot const&) = default;
      
      double d0(double t) const { return c/pow(t,1./3.); }
      
      double d1(double t) const { return -c/(3*pow(t,4./3.)); }
      
      double d2(double t) const { return 4*c/(9*pow(t,7./3.)); }
      
      double d3(double t) const { return -28*c/(27*pow(t,10./3.)); }
      
    private:
      double c = 1.;
    };
    
    
    struct OverThirdRootSquared
    {
      OverThirdRootSquared() = default;
      
      explicit OverThirdRootSquared(double c0) : c(c0) {}
      
      OverThirdRootSquared(OverThirdRootSquared const&) = default;
      OverThirdRootSquared& operator=(OverThirdRootSquared const&) = default;
      
      double d0(double t) const { return c/pow(t,2./3.); }
      
      double d1(double t) const { return -2*c/(3*pow(t,5./3.)); }
      
      double d2(double t) const { return 10*c/(9*pow(t,7./3.)); }
      
      double d3(double t) const { return -70*c/(27*pow(t,10./3.)); }
      
    private:
      double c = 1.;
    };

    template <size_t n, class Source = WrappedMatrix<double,n> >
    using DetOverThirdRoot = MatrixToScalarFunction<OverThirdRoot,Determinant<n,Source> >;

    template <size_t n, class Source = WrappedMatrix<double,n> >
    using DetOverThirdRootSquared = MatrixToScalarFunction<OverThirdRootSquared,Determinant<n,Source> >;

    template <int i, int dim,Invariant I,class Source,class Direction=void> struct ChooseInvariant;

    template <int i, int dim, class Source>
    struct ChooseInvariant<i,dim,Invariant::Principal,Source,void> { typedef PrincipalInvariant<i,double,dim,Source> type; };

    template <int i, int dim, class Source>
    struct ChooseInvariant<i,dim,Invariant::ModifiedPrincipal,Source,void> { typedef ModifiedPrincipalInvariant<i,double,dim,Source> type; };
    
    template <int i, int dim, class Source, class Direction>
    struct ChooseInvariant<i,dim,Invariant::Mixed,Source,Direction> { typedef MixedInvariant<i,double,dim,Source,Direction> type; };

    template <int i, int dim, class Source, class Direction>
    struct ChooseInvariant<i,dim,Invariant::ModifiedMixed,Source,Direction> { typedef ModifiedMixedInvariant<i,double,dim,Source,Direction> type; };


    template <size_t i, int dim, Invariant I, class Source> struct ChoosePrincipalInvariant;
    
    template <size_t i, int dim, class Source> 
    struct ChoosePrincipalInvariant<i,dim,Invariant::Principal,Source> { typedef PrincipalInvariant<i,double,dim,Source> type; };

    template <size_t i, int dim, class Source> 
    struct ChoosePrincipalInvariant<i,dim,Invariant::Modified,Source> { typedef ModifiedPrincipalInvariant<i,double,dim,Source> type; };
    
    
    template <size_t i, int dim, Invariant I, class Source, class Direction> struct ChooseMixedInvariant;
    
    template <size_t i, int dim, class Source, class Direction> 
    struct ChooseMixedInvariant<i,dim,Invariant::Principal,Source,Direction> { typedef MixedInvariant<i,double,dim,Source,Direction> type; };

    template <size_t i, int dim, class Source, class Direction> 
    struct ChooseMixedInvariant<i,dim,Invariant::Modified,Source,Direction> { typedef ModifiedMixedInvariant<i,double,dim,Source,Direction> type; };
    
  
   /*
    * The first three principal matrix invariants.
    */
  template <class Scalar, size_t n, class Source_>
  class PrincipalInvariant<1,Scalar,n,Source_>
  {
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;
    
    explicit PrincipalInvariant(Source const& f_) : f(f_) {}

    PrincipalInvariant(PrincipalInvariant const&) = default;
    PrincipalInvariant& operator=(PrincipalInvariant const&) = default;
    
    ReturnType d0() const { return trace(f.d0()); }
    
    ReturnType d1(Argument const& dF1) const { return trace(f.d1(dF1)); }
    
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return trace(f.d2(dF1,dF2)); }
    
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return trace(f.d3(dF1,dF2,dF3)); }
    
  private:
    Source f;
  };
  
  template <class Scalar, size_t n, class Source_>
  class PrincipalInvariant<2,Scalar,n,Source_>
  {
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;
    
    explicit PrincipalInvariant(Source const& F) : cof(F) {}

    PrincipalInvariant(PrincipalInvariant const&) = default;
    PrincipalInvariant& operator=(PrincipalInvariant const&) = default;
    
    ReturnType d0() const { return trace(cof.d0()); }
    
    ReturnType d1(Argument const& dF) const { return trace(cof.d1(dF)); }
    
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return trace(cof.d2(dF1,dF2)); }
    
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return trace(cof.d3(dF1,dF2,dF3)); }
    
  private:
    CofactorMatrix<n,Source> cof;
  };    
  
  template <class Scalar, size_t n, class Source> class PrincipalInvariant<3,Scalar,n,Source> : public Determinant<3,Source>
  {
  public:
    PrincipalInvariant(Source const& f) : Determinant<3,Source>(f) {}

    PrincipalInvariant(PrincipalInvariant const&) = default;
    PrincipalInvariant& operator=(PrincipalInvariant const&) = default;
  };    
  
   /*
    * Modified principal invariants
    */
  template <class Scalar, size_t n, class Source_>
  class ModifiedPrincipalInvariant<1,Scalar,n,Source_> : public Product<PrincipalInvariant<1,Scalar,n,Source_>,DetOverThirdRoot<n,Source_> >
  {
    typedef PrincipalInvariant<1,Scalar,n,Source_> Inv;
    typedef DetOverThirdRoot<n,Source_> Scaling;
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedPrincipalInvariant(Source const& s) : Product<Inv,Scaling>(Inv(s),Scaling(OverThirdRoot(1),Determinant<n,Source>(s))) {}
    ModifiedPrincipalInvariant(ModifiedPrincipalInvariant const&) = default;
    ModifiedPrincipalInvariant& operator=(ModifiedPrincipalInvariant const&) = default;
  };  
  
  template <class Scalar, size_t n, class Source_>
  class ModifiedPrincipalInvariant<2,Scalar,n,Source_> : public Product<PrincipalInvariant<2,Scalar,n,Source_>,DetOverThirdRootSquared<n,Source_> >
  {    
    typedef PrincipalInvariant<2,Scalar,n,Source_> Inv;
    typedef DetOverThirdRootSquared<n,Source_> Scaling;
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedPrincipalInvariant(Source const& s) : Product<PrincipalInvariant<2,Scalar,n,Source>,DetOverThirdRootSquared<n,Source> >(Inv(s),OverThirdRootSquared(1.),Determinant<n,Source>(s)) {}
    ModifiedPrincipalInvariant(ModifiedPrincipalInvariant const&) = default;
    ModifiedPrincipalInvariant& operator=(ModifiedPrincipalInvariant const&) = default;
  };  

  template <class Scalar, size_t n, class Source_>
  class ModifiedPrincipalInvariant<3,Scalar,n,Source_> : public PrincipalInvariant<3,Scalar,n,Source_>
  {
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedPrincipalInvariant(Source const& s) : PrincipalInvariant<3,Scalar,n,Source>(s) {}
    ModifiedPrincipalInvariant(ModifiedPrincipalInvariant const&) = default;
    ModifiedPrincipalInvariant& operator=(ModifiedPrincipalInvariant const&) = default;
  }; 
  
   /*
    * Mixed invariants for anisotropic materials
    */
  template <class Scalar, size_t n, class Source_, class Direction>
  class MixedInvariant<1,Scalar,n,Source_,Direction>
  {
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;
    
    MixedInvariant(Source const& f_, Direction d_) : f(f_), d(d_) {}
    
    /**
      * @param F_ deformation gradient
      * @param v (fibre) direction
      */ 
    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    MixedInvariant(Source const& f_, Vector const& v) : MixedInvariant(f_,tensorProduct(v,v)) {}
    
    MixedInvariant(MixedInvariant const&) = default;
    MixedInvariant& operator=(MixedInvariant const&) = default;
    
    ReturnType d0() const { return trace(f.d0()*d.d0()); }
    
    ReturnType d1(Argument const& dF1) const { return trace(f.d1(dF1)*d.d0()); }
    
    ReturnType d2(Argument const& dF1, Argument const& dF2) const { return trace(f.d2(dF1,dF2)*d.d0()); }
    
    ReturnType d3(Argument const& dF1, Argument const& dF2, Argument const& dF3) const { return trace(f.d3(dF1,dF2,dF3)*d.d0()); }
    
  private:
    Source f;
    Direction d;
  };
  
  template <class Scalar, size_t n, class Source_, class Direction>
  class MixedInvariant<2,Scalar,n,Source_,Direction> : public MixedInvariant<1,Scalar,n,Product<Source_,Source_>,Direction>
  {
  public:
    typedef Source_ Source;
    
    MixedInvariant(Source const& f, Direction const& d) : MixedInvariant<1,Scalar,n,Product<Source,Source>,Direction>(Product<Source,Source>(f,f),d) {}
    
    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    MixedInvariant(Source const& f, Vector const& v) : MixedInvariant(f,tensorProduct(v,v)) {}
    
    MixedInvariant(MixedInvariant const&) = default;
    MixedInvariant& operator=(MixedInvariant const&) = default;
  };

  template <class Scalar, size_t n, class Source_, class Direction>
  class MixedInvariant<3,Scalar,n,Source_,Direction> : public MixedInvariant<1,Scalar,n,Source_,Product<Direction,Direction> >
  {
  public:
    typedef Source_ Source;
    
    MixedInvariant(Source const& f, Direction const& d) : MixedInvariant<1,Scalar,n,Source,Product<Direction,Direction> >(f,Product<Direction,Direction>(d,d)) {}
    
    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    MixedInvariant(Source const& f, Vector const& v) : MixedInvariant(f,tensorProduct(v,v)) {}
    
    MixedInvariant(MixedInvariant const&) = default;
    MixedInvariant& operator=(MixedInvariant const&) = default;
  };

  /**
   * Modified mixed invariants for anisotropic materials
   */
  template <class Scalar, size_t n, class Source_, class Direction>
  class ModifiedMixedInvariant<1,Scalar,n,Source_,Direction> : public Product<MixedInvariant<1,Scalar,n,Source_,Direction>,DetOverThirdRoot<n,Source_> >
  {
    typedef MixedInvariant<1,Scalar,n,Source_,Direction> F;
    typedef DetOverThirdRoot<n,Source_> G;
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedMixedInvariant(Source const& s, Direction const& d) : Product<F,G>(F(s,d),G(s)) {}

    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    ModifiedMixedInvariant(Source const& s, Vector const& v) : ModifiedMixedInvariant(s,tensorProduct(v,v)) {}

    ModifiedMixedInvariant(ModifiedMixedInvariant const&) = default;
    ModifiedMixedInvariant& operator=(ModifiedMixedInvariant const&) = default;
  };  
  
  template <class Scalar, size_t n, class Source_, class Direction>
  class ModifiedMixedInvariant<2,Scalar,n,Source_,Direction> : public Product<MixedInvariant<1,Scalar,n,Squared<Source_>,Direction>,DetOverThirdRootSquared<n,Source_> >
  {
    typedef MixedInvariant<1,Scalar,n,Squared<Source_>,Direction> F;
    typedef DetOverThirdRootSquared<n,Source_> G;
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedMixedInvariant(Source const& s, Direction const& d) : Product<F,G>(F(s,d),G(s)) {}

    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    ModifiedMixedInvariant(Source const& s, Vector const& v) : ModifiedMixedInvariant(s,tensorProduct(v,v)) {}  

    ModifiedMixedInvariant(ModifiedMixedInvariant const&) = default;
    ModifiedMixedInvariant& operator=(ModifiedMixedInvariant const&) = default;
  };  

  template <class Scalar, size_t n, class Source_, class Direction>
  class ModifiedMixedInvariant<3,Scalar,n,Source_,Direction> : public Product<MixedInvariant<1,Scalar,n,Source_,Squared<Direction> >,DetOverThirdRoot<n,Source_> >
  {
    typedef MixedInvariant<1,Scalar,n,Source_,Squared<Direction> > F;
    typedef DetOverThirdRoot<n,Source_> G;
  public:
    typedef Source_ Source;
    typedef Scalar ReturnType;
    typedef typename Source::Argument Argument;

    ModifiedMixedInvariant(Source const& s, Direction const& d) : Product<F,G>(F(s,d),G(s)) {}

    template <class Vector, class enable = typename std::enable_if<std::is_same<Direction,WrappedMatrix<Scalar,Vector::dimension,true> >::value,void>::type>
    ModifiedMixedInvariant(Source const& s, Vector const& v) : ModifiedMixedInvariant(s,tensorProduct(v,v)) {}  

    ModifiedMixedInvariant(ModifiedMixedInvariant const&) = default;
    ModifiedMixedInvariant& operator=(ModifiedMixedInvariant const&) = default;
  };  


  /**
   * Helper functions for convenient access of invariants
   */
  template <int i, class Scalar, size_t n>
  PrincipalInvariant<i,Scalar,n> getPrincipalInvariant(Dune::FieldMatrix<Scalar,n,n> const& F)
  {
    return PrincipalInvariant<i,Scalar,n>(F);
  }
  
  template <int i, class Scalar, size_t n>
  ModifiedPrincipalInvariant<i,Scalar,n> getModifiedPrincipalInvariant(Dune::FieldMatrix<Scalar,n,n> const& F)
  {
    return ModifiedPrincipalInvariant<i,Scalar,n>(F);
  }
  
  template <int i, class Scalar, size_t n>
  MixedInvariant<i,Scalar,n> getMixedInvariant(Dune::FieldMatrix<Scalar,n,n> const& F, Dune::FieldMatrix<Scalar,n,n> const& M)
  {
    return MixedInvariant<i,Scalar,n>(F,M);
  }
  
  template <int i, class Scalar, size_t n>
  MixedInvariant<i,Scalar,n> getMixedInvariant(Dune::FieldMatrix<Scalar,n,n> const& F, Dune::FieldVector<Scalar,n> const& v)
  {
    return MixedInvariant<i,Scalar,n>(F,v);
  }
  

  template <int i, class Scalar, size_t n>
  ModifiedMixedInvariant<i,Scalar,n> getModifiedMixedInvariant(Dune::FieldMatrix<Scalar,n,n> const& F, Dune::FieldMatrix<Scalar,n,n> const& M)
  {
    return ModifiedMixedInvariant<i,Scalar,n>(F,M);
  }
  
  template <int i, class Scalar, size_t n>
  ModifiedMixedInvariant<i,Scalar,n> getModifiedMixedInvariant(Dune::FieldMatrix<Scalar,n,n> const& F, Dune::FieldVector<Scalar,n> const& v)
  {
    return ModifiedMixedInvariant<i,Scalar,n>(F,v);
  }

  /**
   * Subtract integer valued offset from invariant.
   */
  template <class Invariant>
  class ShiftedInvariant : public Sum<Invariant,Constant<typename Invariant::Argument,typename Invariant::ReturnType> >
  {
  public:
    typedef typename Invariant::Source Source;
    typedef typename Invariant::Argument Argument;
    typedef typename Invariant::ReturnType ReturnType;
    typedef Constant<Argument,ReturnType> Offset;

    explicit ShiftedInvariant(Invariant const& i, ReturnType const& offset) : Sum<Invariant,Offset>(i,Offset(-offset)) {}
    explicit ShiftedInvariant(Source const& s, ReturnType const& offset) : Sum<Invariant,Offset>(Invariant(s),Offset(-offset)) {}

    ShiftedInvariant(ShiftedInvariant const&) = default;
    ShiftedInvariant& operator=(ShiftedInvariant const&) = default;
  };


  /**
   * One fiber direction, transforming with \f$ \text{Cof}(\nabla\Phi) \f$.
   */
  template <class Scalar, int dim>
  class DeformingVector
  {
  public:
    typedef Dune::FieldMatrix<double,dim,dim> Argument;
    typedef Argument ReturnType;

    DeformingVector(Dune::FieldVector<Scalar,dim> const& v_, Dune::FieldMatrix<Scalar,dim,dim> const& dphi_) : dphi(dphi_), v(v_)
    {}

    ReturnType d0() const
    {
      auto w = dphi * v;
      GeomTools::normalize(w);
      return tensorProduct(w,w);
    }

    ReturnType d1(Argument const& arg) const
    {
      auto w = arg.gradient * v;
      GeomTools::normalize(w);
      return tensorProduct(w,w);
    }

    ReturnType d2(Argument const&, Argument const&) const
    {
      return ReturnType(0);
    }

    ReturnType d3(Argument const&, Argument const&, Argument const&) const
    {
      return ReturnType(0);
    }

    DeformingVector(DeformingVector const&) = default;
    DeformingVector& operator=(DeformingVector const&) = default;

  private:
    Argument const& dphi;
    Dune::FieldVector<Scalar,dim> v;
  };

  /**
   * One fiber direction, transforming with \f$ \text{Cof}(\nabla\Phi) \f$.
   */
  template <class Scalar, int dim>
  class DeformingSurfaceVector
  {
  public:
    DeformingSurfaceVector(Dune::FieldVector<Scalar,dim> const& v_, Dune::FieldMatrix<Scalar,dim,dim> const& dphi_) : cof(dphi_), v(v_)
    {}

    Dune::FieldVector<Scalar,dim> d0() const { return cof * v; }

    Scalar operator*(Dune::FieldVector<Scalar,dim> const& w) const { return d0()*w; }

    template <class Arg>
    Dune::FieldVector<Scalar,dim> d1(Arg const& arg) const { return cof.d1(if_(arg.gradient,cof.d0())) * v; }

    template <class Arg1, class Arg2>
    Dune::FieldVector<Scalar,dim> d1(Arg1 const& arg1, Arg2 const& arg2) const { return cof.d2(if_(arg1.gradient,cof.d0()),if_(arg2.gradient,cof.d0())) * v; }

    template <class Arg1, class Arg2, class Arg3>
    Dune::FieldVector<Scalar,dim> d1(Arg1 const& arg1, Arg2 const& arg2, Arg3 const& arg3) const { return cof.d3(if_(arg1.gradient,cof.d0()),if_(arg2.gradient,cof.d0()),if_(arg3.gradient,cof.d0())) * v; }

    DeformingSurfaceVector(DeformingSurfaceVector const&) = default;
    DeformingSurfaceVector& operator=(DeformingSurfaceVector const&) = default;

  private:
    CofactorMatrix<dim> cof;
    Dune::FieldVector<Scalar,dim> const& v;
  };

//  template <class Scalar, int dim>
//  Scalar operator*(Dune::FieldVector<Scalar,dim> const& v, DeformingSurfaceVector<Scalar,dim> const& w) { return w*v; }

}

#endif
