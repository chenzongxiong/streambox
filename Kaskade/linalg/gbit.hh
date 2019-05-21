/*
        GBIT - Good Broyden - Error based iterative linear solver

 *  Written by        L. Weimann 
 *  Purpose           Iterative solution of large linear systems
 *  Category          ???. - Linear systems
 *  Keywords          large linear system, iterative solver, Finite elements
 *  Version           2.0
 *  Revision          April 2014
 *  Latest Change     May 2014
 *  Library           NewtonLib
 *  Code              C++, Double Precision
 *  Environment       Standard C environment on PC's,
                      workstations and hosts.
 *  Copyright     (c) Konrad-Zuse-Zentrum fuer
                      Informationstechnik Berlin (ZIB)
                      Takustrasse 7, D-14195 Berlin-Dahlem
                      phone : + 49/30/84185-0
                      fax   : + 49/30/84185-125
 *  Contact           Lutz Weimann
                      ZIB, Division Scientific Computing, 
                           Department Numerical Analysis and Modelling
                      phone : + 49/30/84185-185
                      fax   : + 49/30/84185-107
                      e-mail: weimann@zib.de
 
   ---------------------------------------------------------------
 
 * Licence
    KASKADE 7 is distributed under the terms of the ZIB Academic License.
    see $KASKADE/academic.txt 
 
       ------------------------------------------------------------
 
 *    Class description
      =================
 
      The calling interface looks as follows:

      template<class LinearSpace>
      class Gbit
      
      The constructor argument list looks like following:
      Gbit( AssembledGalerkinOperator& a, Precond& preconditioner,
            double tol, int maximumNoIterations, int verbosity, int maxRestart)
      
      AssembledGalerkinOperator a:
      The linear operator obtained after assemblation (see Kaskade7 documentation).
      
      Preconditioner precon:
      The preconditioner, which will be used by Gbit.

      double tol: (dune-istl compatibility parameter, see also method setTolerance)
      The prescribed precision for the final iterate y^i.
      The default value for tol is 1.0e-10

      int maximumNoIterations: (dune-istl compatibility parameter, see also method
                                setMaximumNoIterations)
      The maximum number of allowed iterations. 
      The default value for maximumNoIterations is MAXITER_DEFAULT (=1000).

      int verbosity: (dune-istl compatibility parameter, see also method
                      setMonitorLevel)
      The integer value which corresponds to the monitorLevel PrintLevel.
      The verbosity default value is 0.
      
      int maxRestart:
      The maximum number of Gbit iterations before a restart will be done.
      The amount of internal storage used by Gbit mainly depend on this parameter.
      The default value of maxRestart is 9.

      Public types are:
      
      struct GbitInfo
      {
       double precision, normDx;
       int iterationCount, rcode, noMatVecMult, noPrecondCalls;
      };
      
      typedef enum {none=0, minimum=1, verbose=2, debug=3 } PrintLevel;      
      
      The solver must be called by calling the method apply:
      gbit.apply(CoefficientVectors& y, CoefficientVectors const& b,
                 struct GbitInfo& res)

      The parameters are: 
      
      CoefficientVectors& y :
      .
      y must contain on input an initial guess of the linear system
      solution, which is used as the start-vector of the iteration.
      On output, the pointed array contains an approximate solution vector y*,
      which fits the relative error condition epsilon <= rho*||y^i||*opt->tol,
      where ||y||_A = sqrt( sum_1 to n( (y_i/scale_i)^2 ) / n ) denotes a
      scaled Euclidian norm, epsilon is an estimate of the quantity 
      ||y^solution-y^i||, and rho <= 1.0 denotes a safety factor.
            
      CoefficientVectors const& b :
      b must hold the right hand side of the linear system to solve.
     
      struct GbitInfo& info:
      The Gbit info structure. The fields of the structure are set to output
      info of Gbit.
      
      info.precision is of type double and is set to the relative error
      estimate epsilon/(rho*||y^i||). For detailed info, refer to the 
      description of the parameter y above.
      
      info.iterationCount is set to number of iteration steps done.
      
      info.noMatVecMult is set to the number of done calls to the matrix times
      vector multiplication routine a->apply. 
      
      info.noPrecondCalls is set to the number of done calls to the preconditioner
      routine precon->apply.
      
      info.rcode is set to the return-code of GBIT. A return-code 0
      means that Gbit has terminated sucessfully. For the meaning of a
      nonzero return-code, see the error messages list below.

      The following error conditions may occur: (returned via info.rcode)
      --------------------------------------------------------------------
      
         2 Maximum number of iterations (as set by opt->maxiter) exceeded.
         3 No convergence, error estimate epsilon exceeded the limit 1.0e20. 
        20 Nonpositive input for dimensional parameter n.
        21 Nonpositive value for tol supplied.
        22 Negative scaling value for some component of vector scale supplied.
         
                                 
      Further methods for setting parameters are:

      void setTolerance( double tol )
        tol is of type double and must contain the error threshold
        which the relative error of the final iterate y^i must fit.
      
      void setMaximumNoIterations( int maximumIter )
        maximumIter is of type int and must contain the maximum number of allowed
        iterations. if a zero or negative value is supplied, then the maximum
        iteration count will be set to MAXITER_DEFAULT (=1000).
      
      void setSafetyFactor( double safetyFactor )
        safetyFactor is of type double and must contain a safety factor <= 1.0 by
        which tol will be multiplied. The default value of safetyFactor is 0.25 .
      
      void setRescaleMode ( bool rescale )
        rescale is of type bool.
        If set to true, the scaling values, which are used in the norm 
        computation, will be adjusted at the start and at each restart
        to the current iterate.
        If set to false, the initial scaling values are used unmodified 
        throughout the whole iteration.
        The default of this parameter is true.
      
      void setGiantSimpleTermination( bool giantSimple )
        The flag giantSimple is by default set to false. If it is set to true, a
        special termination criterium for computing the simplified Newton correction
        by Giant will be used.
      
      void setErrorStream( std::ostream& errorStream )
        errorStream is of type pointer to std::ostream, as defined by the <iostream>
        header file. It must be set either to a NULL pointer, std::cout, std::cerr,
        or to another file pointer which has been initialized by a fopen call.
        If it is set to NULL, errorStream will be set to std:cerr. The error 
        messages will be printed to errorStream.
        
      void setMonitorStream( std::ostream& monitorStream )
        monitorStream is of type pointer to std::ostream, as defined by the <iostream>
        header file. It must be set either to a NULL pointer, std::cout, std::cerr,
        or to another file pointer which has been initialized by a fopen call.
        If it is set to NULL, monitorStream will be set to std::cout. The monitor 
        output will be printed to monitorStream.
        
      void setErrorLevel( PrintLevel errorLevel )
        errorLevel is of type PrintLevel. If it is set to level None,
        then no error message will be printed if an error condition occurs.
        If it is set to level Minimum or any higher level, then error messages
        will be printed, if appropriate.
        
      void setMonitorLevel( PrintLevel monitorLevel )
        monitorLevel is of type PrintLevel. If it is set to level None,
        then no monitor output will be printed.
        If it is set to level Minimum, a few infomation will be printed.
        If set to level Verbose, then some infomation about each iteration
        step, fitting into a single line, will be printed. The higher level Debug
        is reserved for future additional information output.
                
      void setScalingVector( std::vector<double>* scale ) 
        scale is a pointer to a std::vector<double> of size n, which must hold
        nonnegative scaling threshold values. If zero values are supplied,
        the zeros are replaced by the (possibly adjusted) value of tol.
        If a zero pointer is supplied, then also all scaling threshold components
        are set to the value of tol. 
  
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      2.0      2014/04/28  Special version for use with Kaskade7.2.
      1.0      2006/11/30  Initial release.
      
*/
#ifndef GBIT_HH
#define GBIT_HH
#include <cstdlib>
#include <fstream>
#include <cmath>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#ifndef MACHCONST_
#define MACHCONST_
#define SMALL  DBL_EPSILON
#define EPMACH DBL_MIN
#endif

namespace Kaskade
{
  /**
   * \class Gbit
   * \brief Iterative solution of a large linear system using an error based 
   *   termination-criterium. Includes some special termination-criterium 
   *   options which are used by the giantGbit code.
   * \tparam LinearSpace The underlying linear space of the unknown vector and
   *   right-hand-side.
   */
  template<class LinearSpace>
  class Gbit
  {
    #define TAU_MIN 1.0e-8
    #define TAU_MAX 1.0e2
    #define EPSILON_LIMIT 1.0e20
    
    #define MAXITER_DEFAULT 1000
    
    public:
    typedef enum {
                   none=0,     /*!<no output */
                   minimum=1,  /*!<summary output */
                   verbose=2,  /*!<detailed output - info about each iterate */
                   debug=3     /*!<additional debug output - reserved for future use */
                 } 
                 PrintLevel ;
    
    /**
     * \struct GbitInfo
     *   A structure to hold information returned of the solvers run.
     * \param precision The achieved precision for the final iterate.
     * \param normDx The scaled Euclidian norm of the last iterate-correction.
     * \param iterationCount The number of iteration steps done.
     * \param rcode The return-code of apply.  0 indicates a successfull
     *   termination, a nonzero code indicates some error condition. 
     * \param noMatVecMult The number of done calls to the matrix times
     *   vector multiplication routine a->apply.
     * \param noPrecondCalls The number of done calls to the preconditioner
     *   routine precon->apply.
     */
    struct GbitInfo
    {
       double precision, normDx;
       int iterationCount, rcode, noMatVecMult, noPrecondCalls;
    };
    
    private:
    typedef enum { initial=1,intermediate=2,solution=3,final=4 } DataMode ;

    struct GbitData
    {
      double tau, t, normDx;
      DataMode mode;
    };
    
    public:
    typedef Dune::LinearOperator<LinearSpace,LinearSpace> AssembledGalerkinOperator;
    typedef Dune::Preconditioner<LinearSpace,LinearSpace> Precond;
    
    /**
     * \brief The Gbit constructor routine
     * \param a_ The Dune-istl compatible linear operator A of the linear system 
     *   \f$ A \dot x = b \f$.
     * \param precon_ The Dune-istl compatible preconditioner.
     * \param iteEps_ The prescribed precision for the solution. This is an optional
     *   argument for Dune-istl compatibility.
     * \param iteSteps_ The maximum number of allowed iteration-steps. This is an
     *   optional argument for Dune-istl compatibility.
     * \param verbosity_ The int equivalent of the monitor-output level. This is an
     *   optional argument for Dune-istl compatibility.
     * \param maxRestart_ The maximum number of iteration-steps before a restart
     *   occurs. Since this parameter mainly determines the amount of storage
     *   used by Gbit, there is no alternative set-method for this parameter available.
     */
    Gbit( AssembledGalerkinOperator& a_, Precond& precon_, double iteEps_=1.0e-10,
          int iteSteps_=MAXITER_DEFAULT, int verbosity_=0, int maxRestart_=9):
          tol(iteEps_), maximumNoIterations(iteSteps_), maxRestart(maxRestart_)
    {
      a=&a_;
      precon = &precon_;
      rho = 0.25;
      rescale = true;
      giantSimplifiedCorrection = false;
      errorLevel = verbose;
      monitorLevel = (PrintLevel) verbosity_;
      errorStream = &(std::cerr);
      monitorStream = &(std::cout);
      scale = NULL;
      delta.resize(maxRestart+2);
      for (int i=0;i<maxRestart+2;i++) delta[i]=NULL;
    }
    
    ~Gbit()
    {
      for (int i=0;i<maxRestart+2;i++) delete delta[i];
     }
    
    /**
     * \brief Sets the error tolerance which the final approximate solution \f$x_k\f$ must fit.
     *   The default value is 1.0e-10.
     * \param tol_  The desired error for the solution. 
     */
    void setTolerance( double tol_ ) { tol=tol_; }
    
    /**
     * \brief Sets the maximum allowed number of Newton-iterations. The default value is the
     *   value of the constructor argument iteSteps_, e.g. MAXITER_DEFAULT, if the iteSteps_
     *   argument is not set.
     * \param maximumIter_  The desired maximum number of iterations.
     */
    void setMaximumNoIterations( int maximumIter_ ) { maximumNoIterations=maximumIter_; }
    
    /**
     * \brief Sets the safety factor for the precision control. Must be a positive number <= 1.0 by
     *  which the tolerance will be multiplied.
     *   The default value is 0.25.
     * \param safetyFactor_  The desired safety factor. 
     */
    void setSafetyFactor( double safetyFactor_ ) { rho=safetyFactor_; }

    /**
     * \brief Sets the rescaling mode. The default value is true.
     * \param rescale_ The desired rescaling mode.
     *   If set to true, the scaling values, which are used in the norm 
     *   computation, will be adjusted at the start and at each restart
     *   to the current iterate. \n
     *   If set to false, the initial scaling values are used unmodified 
     *   throughout the whole iteration.
     *  
     */
    void setRescaleMode ( bool rescale_ ) { rescale=rescale_; }
    
    /**
     * \brief Sets the termination-criterium mode. The default value is false.
     * \param giantSimple_ The desired termination-criterium mode.
     *   If it is set to true, a special termination criterium for computing
     *   the simplified Newton correction by giantGbit will be applied. \n
     *   If set to false, the standard termination-criterium will be applied.
     */
    void setGiantSimpleTermination( bool giantSimple_ ) { giantSimplifiedCorrection=giantSimple_; }
    
    /**
     * \brief Sets the error-output stream.
     *   The default error-output stream is std::cerr.
     * \param errorStream_  The desired error-output stream.
     */
    void setErrorStream( std::ostream& errorStream_ ) { errorStream=&errorStream_; }
    
    /**
     * \brief Sets the monitor-output stream.
     *   The default monitor-output stream is std::cout.
     * \param monitorStream_  The desired monitor-output stream.
     */
    void setMonitorStream( std::ostream& monitorStream_ ) { monitorStream=&monitorStream_; }
    
    /**
     * \brief Sets the error-output level.
     *   The default error-output level is verbose.
     * \param errorLevel_  The desired error-output level.
     */
    void setErrorLevel( PrintLevel errorLevel_ ) { errorLevel=errorLevel_; }
    
    /**
     * \brief Sets the monitor-output level.
     *   The default monitor-output level is the value of the constructur argument
     *   verbosity_ (casted to the corresponding PrintLevel value), e.g. 0, if
     *   the verbosity_ argument is not set.
     * \param monitorLevel_  The desired monitor-output level.
     */
    void setMonitorLevel( PrintLevel monitorLevel_ ) { monitorLevel=monitorLevel_; }
    
    /**
     * \brief Sets a user-defined scaling-threshold vector.
     * \param *scale_  A pointer to the selected scaling-threshold vector.
     */
    void setScalingVector( std::vector<double>* scale_ ) 
    { 
      if ( scale )  delete scale;
      if ( scale_ ) scale=scale_;
    }
    
    /**
     * \brief Calls the giantGbit inexact damped Newton-iteration solver.
     * \param y On input, the start-vector for the iteration. \n
     *   On output, the approximate solution vector \f$ x_k \f$, which fits the error 
     *   condition \f$ \parallel x_k - x\star \parallel <= tol \f$ ,
     *   where \f$ x\star \f$  denotes the exact solution, and \f$  \parallel \ldots \parallel \f$
     *   is a scaled Euclidian norm.
     * \param b The right hand side vector of the linear system to be solved.
     * \param res An instance of type GbitInfo to hold the iteration info on termination
     *   of the apply routine.
     */
    void apply(LinearSpace& y, LinearSpace const& b, struct GbitInfo& res)
    { 
      n=b.dim(); 
      LinearSpace qDomain(y),qRange(y);
      int i, iterationCounter=0, noMatVecMult=0, noPrecondCalls=0;
      struct GbitData* data=new struct GbitData;
      if (!data)
      { 
        fprintf(stderr,"\n could not allocate struct data\n");
        res.rcode=-994; return; 
      };
      data->mode      = initial;
      res.rcode = gbitParameterCheckAndPrint();
      if ( res.rcode !=0 ) 
      { 
        if (data) free(data);
        return;
      };
      for (int i=0;i<maxRestart+2;i++) 
        if (!delta[i]) delta[i]=new std::vector<double>(n);
      std::vector<double> *zeta, *delta0, *deltaOfM, *deltaOfMplusOne, *deltai, *deltaip1;
      std::vector<double> sigma(maxRestart+2), t(maxRestart+2), q(n), z(n), ythresh(n), y0(n),
                          ydata(n), bdata(n);
      double gamma, tau, factor, OneMinusTm, ti, epsilon, normOfyIterPlusOne, normdiff, errtol=tol;
      bool stopIteration;
      if ( !errorStream && errorLevel>0 )     errorStream   = &(std::cerr);
      if ( !monitorStream && monitorLevel>0 ) monitorStream = &(std::cout);
      if ( maximumNoIterations <= 0 ) maximumNoIterations = MAXITER_DEFAULT; 
      if ( monitorLevel > 0 ) *monitorStream << std::endl << " GBIT - Version  2.0" << std::endl;
      zeta = &z;
      y.write(ydata.begin());
      b.write(bdata.begin());
      for (int j=0;j<n;j++) y0[j]=ydata[j];
      if ( monitorLevel > 1 )
        *monitorStream << std::endl << " Iter       sigma    NormIter         tau           t" << 
                          std::endl << std::endl;
      for (int j=0;j<n;j++) ythresh[j] = (*scale)[j];
      delta0 = delta[0];
   
      /* Initialization */
      bool doRestart=true, doRestartNow;
      do
      {
        doRestartNow=false;
        if (rescale) gbitRescale(*scale,ydata,ythresh);
        if ( iterationCounter > 0 && monitorLevel > 1 ) *monitorStream << " > RESTART" << std::endl;
        normOfyIterPlusOne = scaledNorm2(ydata,*scale);
        a->apply(y,qRange);  noMatVecMult++; qRange.write(z.begin());
        for (int j=0;j<n;j++) z[j] = bdata[j]-z[j];
        qRange.read(z.begin());
        precon->apply(qDomain,qRange);  
        noPrecondCalls++; 
        qDomain.write(delta0->begin());
        sigma[0] = scaledScalarProduct(*(delta[0]),*(delta[0]),*scale);
  
        /* Iteration loop (restarts excluded) */
        stopIteration = false;
        for (i=0; i<maxRestart+1 && !stopIteration && iterationCounter < maximumNoIterations+1; i++)
        { 
          qDomain.read(delta[i]->begin());
          a->apply(qDomain,qRange);  noMatVecMult++;;
          precon->apply(qDomain,qRange);  noPrecondCalls++; qDomain.write(zeta->begin());
          for (int m=0;m<i;m++)
          {
            OneMinusTm = 1.0-t[m];  
            factor = scaledScalarProduct(*(delta[m]),*zeta,*scale)/sigma[m];
            deltaOfM = delta[m];  deltaOfMplusOne = delta[m+1];
            for (int j=0;j<n;j++) (*zeta)[j] += factor*((*deltaOfMplusOne)[j]-OneMinusTm*(*deltaOfM)[j]);
          };
          gamma = scaledScalarProduct(*(delta[i]),z,*scale);
          if ( gamma != 0.0 ) tau = sigma[i]/gamma;
          else                tau = TAU_MAX;
          if ( tau < TAU_MIN ) 
          { 
            if ( monitorLevel > 1 ) *monitorStream << " > tau = " << tau << std::endl;
            if ( i>0 ) 
            { 
              doRestartNow=true; 
              break; 
            }
            else
            { 
              t[i] = 1.0;
              if ( monitorLevel > 0 )
                 *monitorStream << " > Restart condition ignored, set t[0] = " << t[i] << std::endl;
            };
          }
          else
          {
            if ( tau <= TAU_MAX ) 
              t[i] = tau;
            else
              t[i] = 1.0;
          };
          ti = t[i];  deltai = delta[i];
          if ( monitorLevel > 1 )
            *monitorStream << " " << std::setw(4) << iterationCounter << " " << std::setw(11)
                           << sigma[i] << " " << std::setw(11) << normOfyIterPlusOne << " " 
                           << std::setw(11) << tau << " " << std::setw(11)
                           << t[i] << std::endl;
          data->normDx = sqrt(sigma[i]);  data->tau = tau;  data->t = t[i];
          //itlin_dataout(iterationCounter,n,y,data);
          data->mode = intermediate;
          for (int j=0;j<n;j++) ydata[j] += ti*(*deltai)[j];
          y.read(ydata.begin());
          deltaip1 = delta[i+1];
          factor = 1.0-ti+tau;
          for (int j=0;j<n;j++) (*deltaip1)[j] = factor*(*deltai)[j]-tau*z[j];
          sigma[i+1] = scaledScalarProduct(*deltaip1,*deltaip1,*scale);
          if ( i>=1 ) epsilon = 0.5*sqrt(sigma[i-1]+2.0*sigma[i]+sigma[i+1]);
          else        epsilon = sqrt(sigma[1]);
          normOfyIterPlusOne = scaledNorm2(ydata,*scale);
          if (giantSimplifiedCorrection)
          { 
            for (int j=0;j<n;j++) q[j]=ydata[j]-y0[j];
            normdiff = scaledNorm2(q,*scale);
            stopIteration = (epsilon/normdiff < rho*errtol);
          }
          else 
            stopIteration = (epsilon < rho*normOfyIterPlusOne*errtol);
          iterationCounter++;
          if ( !stopIteration && epsilon > EPSILON_LIMIT ) 
          { 
            res.rcode = 3; break;
          };
        };
        if ( res.rcode != 0 ) break;
        if ( doRestartNow || ( i>=maxRestart && !stopIteration && iterationCounter<maximumNoIterations) )
          continue;
        else
          doRestart = false;
      }
      while ( doRestart ) ;
      if ( (res.rcode==0) && (iterationCounter >= maximumNoIterations) ) res.rcode = 2;
      
      if ( monitorLevel > 1 )
        *monitorStream << " " << std::setw(4) << iterationCounter << " " << std::setw(11)
                           << sigma[i] << " " << std::setw(11) << normOfyIterPlusOne << " " 
                           << std::setw(11) << tau << " " << std::setw(11)
                           << t[i-1] << std::endl;
      data->mode = ( res.rcode==0 ? solution : final );
      data->normDx = sqrt(sigma[i]);  data->tau = tau;  data->t = t[i-1];
      // itlin_dataout(iterationCounter,n,y,data);
      if ( errorLevel > 0 && res.rcode != 0 )
        {
          switch ( res.rcode )
           {
            case     2:
              *errorStream << "\n Error - Maximum allowed number of iterations exceeded\n" ;
              break;
            case     3:
              *errorStream << "\n Error - no convergence, correction too large\n" ;
              break;
            default   :
              *errorStream << "\n Error, code=" << res.rcode << std::endl;
           };
        };
      delete data;   
      res.precision = epsilon/(rho*normOfyIterPlusOne);
      res.normDx    = data->normDx;
      res.iterationCount = iterationCounter;
      res.noMatVecMult   = noMatVecMult;
      res.noPrecondCalls = noPrecondCalls;
    }
    
    private:

    void gbitRescale(std::vector<double>& scale, std::vector<double> const& y, 
                      std::vector<double> const& ythresh)
    {  
      int const n = scale.size();
      assert ( (y.size()==n) && (ythresh.size()==n) );
      for (int i=0;i<n;i++) 
        scale[i] = std::max(ythresh[i],std::max(fabs(y[i]),SMALL));
    }
    
    int gbitParameterCheckAndPrint()
    #define TOLMIN 1.0e-15
    #define TOLMAX 0.5
    {  
      double large = 1.0/SMALL, default_scale;
      if ( n<=0 ) 
      { 
        if ( errorLevel>0 )
          *errorStream << "\n Error - Number of Equations/unknowns must be >0\n" ;
        return 20;
      };
      if ( tol <= 0 )
      { 
        if ( errorLevel>0 )
          *errorStream << "\n Error - tol must be positive\n" ;
        return 21;
      }
      else
        { 
          if ( tol > TOLMAX ) 
          {  
            tol = TOLMAX;
             if ( errorLevel>1 )
               *errorStream <<
               "\n User prescribed RTOL decreased to reasonable largest value RTOL=" <<
               tol << std::endl;
          }
          else if ( tol < TOLMIN ) 
          { 
            tol = TOLMIN;
            if ( errorLevel>1 )
               *errorStream <<
               "\n User prescribed RTOL increased to reasonable smallest value RTOL=" <<
               tol << std::endl;
          };
         };
      default_scale = 1;
      if ( !scale )
      {
        scale = new std::vector<double>(n);
        for (int i=0;i<n;i++) (*scale)[i]=default_scale;
      }
      else
      { 
        for (int i=0;i<n;i++) 
        { 
          if ( (*scale)[i] < 0.0 ) 
          { 
            if ( errorLevel>0 )
              *errorStream <<
                "\n Error - negative value (*scale)[" << i << "] supplied\n";
            return 22;
          }
          else if ( (*scale)[i] == 0.0 ) 
            (*scale)[i] = default_scale;
          else if ( (*scale)[i] < SMALL )
          { 
            if ( errorLevel>1 )
              *errorStream <<
              "\n Warning: (*scale)[" << i << "] too small - increased to "<< SMALL << "\n";
            (*scale)[i] = SMALL;
          }
          else if ( (*scale)[i] > large )
          { 
            if ( errorLevel>1 )
              *errorStream <<
              "\n Warning: (*scale)[" << i << "] too large - decreased to %" << large << "\n";
            (*scale)[i] = large;
          };
        };
      };

      if ( maxRestart <= 0 ) 
      { 
        maxRestart = std::min(10,n);
        if ( errorLevel > 1 ) 
          *errorStream <<
          " Warning: maxRestart not set; is reset to maxRestart=" << maxRestart << std::endl;
      }
      else if ( maxRestart > n ) 
      { 
        maxRestart = std::min(10,n);
        if ( errorLevel > 1 ) 
          *errorStream <<
          " Warning: maxRestart is greater than n; is reset to maxRestart=" << maxRestart << std::endl;
      };

      if ( rho <= 0.0 )
      { 
        rho = 1.0;
        if ( errorLevel > 1 ) 
          *errorStream <<
          " Warning: rho not set; is reset to rho=" << rho << std::endl;
      };

      if ( monitorLevel==0 ) return 0;
      *monitorStream << "\n Problem dimension: n = " << n << std::endl;
      *monitorStream << "\n Prescribed relative precision: " << tol << std::endl;
      *monitorStream << " The maximum permitted number of iteration steps is: " <<
                      maximumNoIterations << std::endl;
      *monitorStream << " The maximum number of iterations before restart is: " << 
                        maxRestart << std::endl;
      *monitorStream << " The safety factor is rho = " << rho << std::endl;                
      return 0;
    }
    
    double scaledNorm2(std::vector<double> const& v,std::vector<double> const& scale)
    {
      int const n=scale.size();
      assert( v.size()==n );
      double t, rval = 0.0;
      for (int i=0;i<n;i++) {t=v[i]/scale[i]; rval += t*t;};
      return sqrt( rval / (double)n );
    }
    
    double scaledScalarProduct( std::vector<double> const& v1, std::vector<double> const& v2,
                                std::vector<double> const& scale)
    {
      int const n=scale.size();
      assert( (v1.size()==n) && (v2.size()==n) );
      double t1, t2, rval = 0.0;
      for (int i=0;i<n;i++) 
        {t1=v1[i]/scale[i]; t2=v2[i]/scale[i]; rval += t1*t2;};
      return rval / (double)n ;
    }
    
    int n;
    double tol, rho;
    std::vector<std::vector<double>*> delta;
    std::vector<double>* scale;
    int maxRestart, maximumNoIterations;
    bool rescale, giantSimplifiedCorrection;
    PrintLevel errorLevel, monitorLevel, dataLevel;
    std::ostream *errorStream, *monitorStream;
    GbitInfo info;
    AssembledGalerkinOperator* a;
    Precond* precon;
  };
}
#endif
