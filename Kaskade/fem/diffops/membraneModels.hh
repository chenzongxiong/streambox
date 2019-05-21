/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2013 Zuse Institute Berlin                            */
/*  Copyright (C) 2013-2013 U Milan                                          */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MEMBRANEMODELS_HH
#define MEMBRANEMODELS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>


namespace Kaskade {
  
  /**
   * \ingroup diffops
   * \defgroup membranemodels Membrane models
   * \brief Electrophysiological membrane models for use in cardiac simulations
   * 
   * This module contains electrophysiological cell membrane models defining the transmembrane ionic current depending
   * on the transmembrane voltage and so-called gating variables, which define the internal state of ion channels 
   * in the membrane or ion concentrations inside certain cell compartments. These states and concentrations 
   * evolve themselves, again depending on transmembrane voltage and gating variables. This evolution is defined by
   * a membrane model as well.
   * 
   * Additionally, models for active tensile stress generation in muscle cells are defined here.
   * 
   * A lot of membrane models are available at the CellML site (http://www.cellml.org/).
   * 
   * Membrane models describe ODEs for the transmembrane voltage \f$ u \f$ and gating variables \f$ v \f$:
   * \f[ \begin{aligned} \dot u &= I(u,v) \\ \dot v &= f(u,v) \end{aligned} \f]
   * Here, \f$ I \f$ is the transmembrane current and \f$ f \f$ the gating dynamic. These two functions and their
   * deriviatives are defined by membrane models.
   * 
   * Active stress generation models define ann ODE for the tensile stress generated in muscle cells. 
   * \f[ \begin{aligned} \dot s &= f(s,t) \\ a &= g(s) \end{aligned} \f]
   * Here, \f$ s \f$ is the internal state of the stress generation system, \f$ t \f$ are trigger variables, e.g., 
   * calcium ion concentration, and \f$ a \f$ is the active stress. The functions \f$ f \f$ and \f$ g \f$ and
   * their derivatives are defined by active stress models.
   */
  
  /**
   * \ingroup membranemodels
   * \brief Convenience base class for membrane models providing numerical differentiation.
   */
  template <class Derived, int n>
  struct MembraneModelBase
  {
    /**
     * \brief Number of gating variables.
     */
    static int const nGating = n;
    
    /**
     * \brief Vector type holding gating variables.
     */
    typedef Dune::FieldVector<double,nGating>         Gating;
    
    /**
     * \brief Matrix type for the Jacobian of the gating dynamics.
     */
    typedef Dune::FieldMatrix<double,nGating,nGating> GatingJacobian;
    
    /**
     * \brief Default constructor.
     * 
     * Name and resting state are default initialized, i.e. set to "" and 0.0, respectively.
     */
    MembraneModelBase() {}
    
    /**
     * \brief Constructor specifying the model data.
     * \param name the human-readable model name
     * \param uFix the resting state of the transmembrane voltage
     * \param vFix the resting state of the gating variables
     */
    MembraneModelBase(std::string const& name, double uFix, Gating const& vFix): nm(name), fix(uFix,vFix) {}
    
    /**
     * \brief Human-readable name of the membrane model.
     */
    std::string const& name() const { return nm; }
    
    /**
     * \brief Value of the resting state fixed point.
     */
    std::pair<double,Gating> const& restState() const { return fix; }
    
    /**
     * \brief The transmembrane ion current derivative w.r.t. the transmembrane voltage.
     * 
     * This is a default implementation based on numerical differentiation. Overload this if a 
     * better implementation is available.
     * \param u transmembrane voltage [V]
     * \param v gating variables      
     * \param h step size for numerical differentiation
     */
    double current_du(double u, Gating const& v, double h=1e-5) const
    {
      double fuv = static_cast<Derived const&>(*this).current(u-h,v);
      double fuiv = static_cast<Derived const&>(*this).current(u+h,v);
      
      return (fuiv-fuv)/(2*h);
    }
    
    /**
     * \brief The transmembrane ion current derivative w.r.t. the gating variables.
     * 
     * This convenience method performs numerical differentiation and should be overwritten
     * in derived classes by a better implementation.
     * 
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating current_dv(double u, Gating const& v) const
    {
      double f = static_cast<Derived const&>(*this).current(u,v);
      Gating fv_uv,dir;
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i] = 1e-5;
        fv_uv[i] = (static_cast<Derived const&>(*this).current(u,v+dir)-f)/1e-5;
      }
      return fv_uv;
    }
        
    /**
     * \brief The derivative of the right hand side for the evolution of gating variables w.r.t. the transmembrane voltage.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating gatingRhs_du(double u, Gating const& v) const
    {
      Gating guA = static_cast<Derived const&>(*this).gatingRhs(u-1e-5,v);
      Gating guB = static_cast<Derived const&>(*this).gatingRhs(u+1e-5,v);
      return (guB-guA)/2e-5;
    }
    
    /**
     * \brief The derivative of the right hand side for the evolution of gating variables w.r.t. the gating variables.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    GatingJacobian gatingRhs_dv(double u, Gating const& v) const
    {
      GatingJacobian gv_uv;
      Gating col, dir, guv = static_cast<Derived const&>(*this).gatingRhs(u,v);
      for (int i=0; i<nGating; ++i){
        dir = 0;
        dir[i]=1e-5;
        col = (static_cast<Derived const&>(*this).gatingRhs(u,v+dir)-guv)/1e-5;
        for (int j=0; j<nGating; ++j) gv_uv[j][i] = col[j];
      }
      return gv_uv;
    }
    
  private:
    std::string nm;
    std::pair<double,Gating> fix;
  };
  
  //----------------------------------------------------------------------------------------
  
  
  /**
   * \ingroup membranemodels
   * \brief Phenomenological model by Aliev and Panfilov.
   * 
   * A phenomenological model with just one gating variables published in 
   * Chaos, Solitions, and Fractals 7(3):293-301, 1996. 
   * 
   * The transmembrane voltage covers a physiological range [0,1].
   */
  struct AlievPanfilov: public MembraneModelBase<AlievPanfilov,1>
  {
    /**
     * \brief number of gating variables
     */
    static int const nGating = 1;
    
    typedef Dune::FieldVector<double,nGating> Gating;
    typedef Dune::FieldMatrix<double,nGating,nGating> GatingJacobian;
    
    /**
     * \brief Default constructor.
     */
    AlievPanfilov();
    
    /**
     * \brief The transmembrane ion current.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    double current(double u, Gating const& v) const;
    
    /**
     * \brief The transmembrane ion current derivative w.r.t. the transmembrane voltage.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    double current_du(double u, Gating const& v) const;
    
    /**
     * \brief The transmembrane ion current derivative w.r.t. the gating variables.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating current_dv(double u, Gating const& v) const;
    
    /**
     * \brief The right hand side for the evolution of gating variables.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating gatingRhs(double u, Gating const& v) const;
    
    /**
     * \brief The derivative of the right hand side for the evolution of gating variables w.r.t. the transmembrane voltage.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating gatingRhs_du(double u, Gating const& v) const;
    
    /**
     * \brief The derivative of the right hand side for the evolution of gating variables w.r.t. the gating variables.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    GatingJacobian gatingRhs_dv(double u, Gating const& v) const;
    
  private:
    double a, eps1, ga, gs, mu1, mu2;
  };
  
  //-----------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup membranemodels
   * \brief Physiologcial model by tenTusscher et al.
   * 
   * A detailed physiological model with 18 gating variables.
   */
  struct TenTusscher: public MembraneModelBase<TenTusscher,18>
  {
    /**
     * \brief Default constructor.
     */
    TenTusscher();
    
    /**
     * \brief The transmembrane ion current.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    double current(double u, Gating const& v) const;
    
    /**
     * \brief The right hand side for the evolution of gating variables.
     * \param u transmembrane voltage
     * \param v gating variables
     */
    Gating gatingRhs(double u, Gating const& v) const;
  };

  //-----------------------------------------------------------------------------------------------------

  /**
   * \ingroup membranemodels
   * \brief Convenience base class for active stress generation models providing numerical differentiation.
   * 
   * This models the active stress generation of myocardial tissue by pointwise ODEs of the form
   * \f[ \dot s = f(s,e), \quad a = g(s), \f]
   * where \f$ s \f$ are the internal states, \f$ e \f$ are the trigger variables driving the activation,
   * and \f$ a \f$ is the actual stress generated.
   * 
   * \tparam n the number of internal state variables
   * \tparam m the number of external trigger variables
   */
  template <class Derived, int n, int m>
  struct ActiveStressModelBase
  {
    /**
     * \brief Number of internal state variables.
     */
    static int const nState = n;
    
    /**
     * \brief Number of trigger variables.
     */
    static int const nTrigger = m;
    
    /**
     * \brief Vector type holding the internal variables.
     */
    typedef Dune::FieldVector<double,nState>  State;
    
    /**
     * \brief Vector type holding trigger variables.
     */
    typedef Dune::FieldVector<double,nTrigger> Trigger;
    
    /**
     * \brief Matrix type for the Jacobian of the State dynamics.
     */
    typedef Dune::FieldMatrix<double,nState,nState> StateJacobian;
    
    /**
     * \brief Matrix type for the derivative of the state dynamics wrt the triggers.
     */
    typedef Dune::FieldMatrix<double,nState,nTrigger> StateTriggerJacobian;
    
    /**
     * \brief Default constructor.
     */
    ActiveStressModelBase() {}
    
    /**
     * \brief Constructor specifying the model data.
     * \param name the human-readable model name
     * \param fix the resting state of the state variables
     */
    ActiveStressModelBase(std::string const& name, State const& fix_): nm(name), fix(fix_) {}
    
    /**
     * \brief Human-readable name of the active stress model.
     */
    std::string const& name() const { return nm; }
    
    /**
     * \brief Value of the resting state fixed point.
     */
    State const& restState() const { return fix; }
    
    /**
     * \brief The deriviative of internal state dynamics wrt state variables.
     * \param s state variables
     * \param t trigger variables
     */
    StateJacobian state_ds(State const& s, Trigger const& t) const
    {
      auto const& me = static_cast<Derived const&>(*this);
      StateJacobian ds;
      State r = me.state(s,t), dr, h;
      for (int i=0; i<nState; ++i)
      {
        h = s;
        h[i] += 1e-5;
        dr = (me.state(h,t)-r)/1e-5;
        for (int j=0; j<nState; ++j)
          ds[j][i] = dr[j];
      }
      return ds;
    }
    
    /**
     * \brief The deriviative of internal state dynamics wrt trigger variables.
     * \param s state variables
     * \param t trigger variables
     */
    StateTriggerJacobian state_dt(State const& s, Trigger const& t) const
    {
      auto const& me = static_cast<Derived const&>(*this);
      StateTriggerJacobian ds;
      State r = me.state(s,t), dr;
      Trigger h;
      for (int i=0; i<nTrigger; ++i)
      {
        h = t;
        h[i] += 1e-5;
        dr = (me.state(s,h)-r)/1e-5;
        for (int j=0; j<nState; ++j)
          ds[j][i] = dr[j];
      }
      return ds;
    }
    
    /**
     * \brief The derivative of stress wrt state variables.
     * \param s state variables
     * \param t trigger variables
     */
    State stress_ds(State const& s) const
    {
      auto const& me = static_cast<Derived const&>(*this);
      State da, h;
      double r = me.stress(s);
      for (int i=0; i<nState; ++i)
      {
        h = s;
        h[i] += 1e-5;
        da[i] = (me.stress(h)-r)/1e-5;
      }
      return da;
    }
    
  private:
    std::string nm;
    State fix;
  };
  
  /**
   * \brief Active stress generation from XX
   */
  class XXstress: public ActiveStressModelBase<XXstress,1,1>
  {
  public:
    XXstress():  ActiveStressModelBase<XXstress,1,1>("XXstress",State(0)) {}
    
    double stress(State const& s) const
    {
      return s[0];
    }
    
    State state(State const& s, Trigger const& t) const
    {
      double e0 = .5; 
      double kta = 47.9; 
//       double un = t/120;
      double un = t;
      return e0 * (kta*un-s) / (1+exp(-100*(un-.05))); 
    }
  };
  
  
}

#endif