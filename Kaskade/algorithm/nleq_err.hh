/*
    nleqErr - Error-based Damped Newton algorithm
    
 *  Written by        L. Weimann 
 *  Purpose           solution of systems of nonlinear equations
 *  Method            ERRor-based Damped Newton algorithm
                      (see reference below)
 *  Category          F2a. - Systems of nonlinear equations
 *  Keywords          Nonlinear equations, Newton methods, Finite elements
 *  Version           2.0
 *  Revision          April 2014
 *  Latest Change     May 2014
 *  Library           Kaskade7
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
 
 *    References:
 
      /1/ P. Deuflhard:
          Newton Methods for Nonlinear Problems. -
          Affine Invariance and Adaptive Algorithms.
          Series Computational Mathematics 35, Springer (2004)
 
   ---------------------------------------------------------------
 
 * Licence
    KASKADE 7 is distributed under the terms of the ZIB Academic License.
    see $KASKADE/academic.txt 

      ------------------------------------------------------------
 
 *    Class description
      =================
 
      The calling interface looks as follows:

      template<class Grid, class Equation, class VariableSet, class Spaces>
      class NleqSolver

      Public types are:
      
      typedef enum {none=0, minimum=1, verbose=2, debug=3} PrintLevel

      typedef enum {mildlyNonlinear=2, highlyNonlinear=3, extremelyNonlinear=4}
        NonlinProblemType

      struct NleqInfo
      {
         double precision, normDx;
         std::vector<double> *fx, *dx;
         int noIterations, returnCode, subCode, noFunctionEvaluations, noJacobianEvaluations;
      }

      Public methods are:

      void setTolerance( double tol )
        tol is of type double and must be the error tolerance which
        the final approximate solution x^k must fit.

      void setNonlinProblemType(NonlinProblemType nonlinType )
        Sets the problem type of the nonlinear problem.
        The following classifications may be used:
        mildlyNonlinear: The problem is considered to be mildly nonlinear and
                          NleqSolver starts up with dampingfactor=1.
        highlyNonlinear: The problem is considered to be highly nonlinear and
                          NleqSolver starts up with dampingfactor=1.0e-4.
        extremelyNonlinear: The problem is considered to be extremely nonlinear
                          and NleqSolver starts up with dampingfactor=1.0e-6.
                          Moreover, restricted is set automatically to true.

      void setMaximumNoIterations( int maximumIter )
        maximumIter is of type int and must contain the maximum number of allowed
        Newton iterations. if a zero or negative value is supplied, then the maximum
        iteration count will be set to 50.
      
      void setPreconFillLevel( int preconFillLevel )
        preconFillLevel is of type int and must contain the PrecondType::ILUK preconditioners
        filllevel (fill-in level).
    
      void setRestricted( bool restricted )
        restricted is of type bool.
        If set to true, then the restricted monotonicity test will be applied for
        determination whether the next iterate (and the associate damping factor
        lambda) will be accepted. This means, with
        theta = ||dxBar(k+1)|| / ||dx(k)||,
        (where dx(k) denotes the k-th Newton correction, and dxBar(k+1) denotes
         the sucessive simplified Newton correction)
        the condition theta <= 1.0 - lambda/4 must be fit.
        If set to false, then the standard monotonicity test will be applied, i.e.
        the following condition must be fit:
        theta < 1.0.
      
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
        errorLevel is of type PrintLevel. If it is set to level none,
        then no error message will be printed if an error condition occurs.
        If it is set to level minimum or any higher level, then error messages
        will be printed, if appropriate.
      
      void setMonitorLevel( PrintLevel monitorLevel )
        monitorLevel is of type PrintLevel. If it is set to level none,
        then no monitor output will be printed.
        If it is set to level minimum, a few infomation will be printed.
        If set to level verbose, then some infomation about each Newton iteration
        step, fitting into a single line, will be printed. The higher level debug
        is reserved for future additional information output.
              
      void setDataLevel( PrintLevel dataLevel )
        dataLevel is of type PrintLevel. 
        If it is set to level none, then no data output will be done.
        If it is set to level minimum, the initial values and the final
        values of the iteration vector x will be output to a file named 
        "<outFilePrefix>mmm.vtu" in vtu-format.
        If it is set to level verbose or any higher level, then each Newton iterate
        x will be output to a file named "<outFilePrefix>nnn.vtu" in vtu-format.
       
      void setOutFilePrefix( const std::string& outFilePrefix )
        outFilePrefix is a filename prefix used for the .vtu output file
        for each iterate. The complete filename for the k-th iterate output
        is outFilePrefinnn.vtu, where nnn is the value of k filled with leading 
        zeros for k<100.
      
      void setScalingVector( std::vector<double>* scale )
        scale is of type pointer to a std::vector<double> of size n. 
        *scale must, if scale is not the NULL pointer, contain positive scaling
        values, which are used in computations of scaled norms, as follows:
        || x || = squareroot(sum(1 to n) ( x_i/scale_i )^2)
        The pointer may be initialized with a NULL pointer. In this case, all
        scaling values are internally set to a default value, which depends on
        the value of nonlinType.

      struct NleqInfo getInfo()
        returns some info about the Newton iteration and final Newton correction
        and residuum vector.
        In detail: 
        info.precision is of type double and is set to the achieved scaled norm
        of the error of the final iterate.
        
        info.normDx is of type double and is set to the unscaled norm of the
        last Newton correction.
        
        info.fx is a pointer to a double array of size n, which contains the
        final residuum vector.
        
        info.noIterations is set to number of Newton iteration steps done.

        info.noLiniterIterations ist set to the total number of linear solver iterations
        done.
        
        info.noFunctionEvaluations is set to the number of done evaluations of the right-
        hand side of the linearized system.
        
        info.noJacobianEvaluations is set to the number of computations of the Jacobian
        of the linearized system.
        
        info.returnCode is set to the return-code of nleq_err. A return-code 0
        means that nleq_err has terminated sucessfully. For the meaning of a
        nonzero return-code, see the error messages list below.
        
        info.subCode is set for certain failure conditions to the error code
        which has been returned by a routine called from nleq_err.
  
        The following error conditions may occur: (returned via info.returnCode)
        --------------------------------------------------------------------
        
           1 Singular Jacobian matrix (not yet realized),
             nleqErr cannot proceed the iteration.
           2 Maximum number of Newton iteration (as set by maximumNoIterations) exceeded.
           3 No convergence of Newton iteration, damping factor became too small.
           4 Linear solver Dune::BiCGSTABSolver failed to converge when computing
             the ordinary Newton-correction.
           5 Linear solver Dune::BiCGSTABSolver failed to converge when computing
             the simplified Newton-correction.
          21 Nonpositive value for xtol supplied.
          22 Negative scaling value for some component of vector *scale
             supplied.
  
      void nleqErr(GridManager<Grid>& gridManager,
                   Equation& eq,
                   VariableSet const& variableSet, 
                   typename VariableSet::VariableSet* x,
                   Spaces const& spaces)

      A short description of the parameters of nleqErr follows: 
            
      GridManager<Grid>& gridManager:
      The grid obtained by discretization (see Kaskade7 documentation).
      
      Equation& eq:
      The structure eq defining the Weak-formulation equation and boundary
      conditions mainly by the class-menbers DomainCache::d1, DomainCache::d2,
      BoundaryCache::d1 and BoundaryCache::d2. For details, consult the 
      Kaskade7 documentation.
      
      VariableSet const& variableSet: 
      For details, consult the Kaskade7 documentation.
      
      VariableSet::VariableSet* x:
      x must contain on input an initial guess of the problems solution,
      which is used as the start-vector of the damped Newton iteration.
      On output, the pointed array contains an approximate solution vector x^k,
      which fits the error condition
      || x^k - x* || <= xtol,
      where x* denotes the exact solution, and ||...|| is a scaled 
      Euclidian norm.
      For the class member VariableSet::VariableSet, consult the Kaskade7 
      documentation.
      
      Spaces const& spaces:
      For details, consult the Kaskade7 documentation.
               
      Summary of changes:
      -------------------
      
      Version  Date        Changes
      2.0      2014/04/29  Special version for use with Kaskade7.2
      1.1.1    2006/06/06  Missing int return code in function 
                           nleqErrInitScale, bug fixed.
      1.1      2006/06/02  Added the output data files iterfile, resfile and
                           miscfile, where optional data output is provided,
                           a single line, starting with the iteration number,
                           for each iteration step. 
      1.0      2006/05/30  initial release.
      
*/

#ifndef NLEQERR_HH
#define NLEQERR_HH

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/gridmanager.hh"
#include "linalg/iluprecond.hh"
#include "io/vtk.hh"

#ifndef MACHCONST_
#define MACHCONST_
#define SMALL  DBL_EPSILON
#define EPMACH DBL_MIN
#endif

namespace Kaskade
{
/** \ingroup alg
   * \class NleqSolver
   * \brief Solution of a nonlinear system of equations by a damped Newton-method.
   * \tparam Grid The underlying Grid, eg. UGGrid, AlbertaGrid or some type of ALUGrid.
   * \tparam Equation The name of the equation class, which defines the subclass DomainCache
   *   with d1, d2 and BoundaryCache with d1, d2.
   * \tparam VariableSet The VariableSet of which representation-type is the Newton-iteration
   *   vector.
   * \tparam Spaces The spaces used in the problem formulation.
   */
  template<class Grid, class Equation, class VariableSet, class Spaces>
  class NleqSolver
  {
    public:
    typedef enum {
                   none=0,     /*!<no output */
                   minimum=1,  /*!<summary output */
                   verbose=2,  /*!<detailed output - info about each iterate */
                   debug=3     /*!<additional debug output - reserved for future use */
                 } 
                 PrintLevel ;
    /**
     * \enum PrintLevel
     *   The level of error-, monitor- and data-output.
     */
    typedef enum {
                   mildlyNonlinear=2,     /*!< problem is mildly nonlinear */
                   highlyNonlinear=3,     /*!< problem is highly nonlinear */
                   extremelyNonlinear=4   /*!< problem is extremely nonlinear */
                 }
                 NonlinProblemType ;
    /**
     * \enum NonlinProblemType
     *   The nonlinearity-level of the problem to be solved.
     */
        
    /**
     * \struct NleqInfo
     *   A structure to hold information returned of the solvers run.
     * \param precision The achieved precision for the final iterate.
     * \param normDx The scaled Euclidian norm of the last Newton-correction.
     * \param *fx A pointer to the std::vector<double> which holds the final residuum. 
     * \param *dx A pointer to the std::vector<double> which holds the last Newton-correction.
     * \param noIterations The number of Newton-iterations done.
     * \param noLiniterIterations The total number of linear solver iterations done.
     * \param returnCode The return-code of the nleqErr solver. 0 indicates a successfull
     *   termination, a nonzero code indicates some error condition.
     * \param subCode May hold some addtional integer information when a nonzero return-code
     *   occurs.
     * \param noFunctionEvaluations The number of done function- (eg. d1-) evaluations.
     * \param noJacobianEvaluations The number of done Jacobian- (eg. d2-) evaluations.
     */
    struct NleqInfo
    {
       double precision, normDx;
       std::vector<double> *fx, *dx;
       int noIterations, noLiniterIterations, returnCode, subCode,
           noFunctionEvaluations, noJacobianEvaluations;
    };
    
    private:
    typedef enum {initial=1,intermediate=2,solution=3,final=4} OutputMode;
    struct NleqData
    {
      std::vector<double> *fx, *dx;
      double normF, normDx, normDxBar, lambda, theta;
      OutputMode mode;
    };
    
    typedef typename Grid::LeafGridView LeafView;
    typedef typename VariableSet::VariableSet VariableSetRepresentation;
    
    #define THETA_MAX 0.5
    #define MAX_ITER_DEFAULT 50
    #define LAMBDA_START_DEFAULT 1.0e-2
    #define LAMBDA_START_EXTREMELY_DEFAULT 1.0e-4
    #define LAMBDA_MIN_DEFAULT 1.0e-4
    #define LAMBDA_MIN_EXTREMELY_DEFAULT 1.0e-8
    
    public:

    NleqSolver()
    {
      xtol = 1.0e-5;
      nonlinType = highlyNonlinear;
      maximumNoIterations = 50;
      restricted = false;
      preconFillLevel=0;
      errorStream = &(std::cerr);
      monitorStream = &(std::cout);
      errorLevel = verbose;
      monitorLevel = none;
      dataLevel = none;
      outFilePrefix = "iter_";
      scale = new std::vector<double>(0);
    }

    /**
     * \brief Sets the error tolerance which the final approximate solution \f$x_k\f$ must fit.
     *   The default value is 1.0e-5.
     * \param tol_  The desired error for the solution. 
     */
    void setTolerance( double tol_ ) { xtol=tol_; }
    
    /**
     * \brief Set the classification of the nonlinear problem.
     *   The default value is highlyNonlinear.
     * \param nonlinType_  The following classifications may be used: \n
     *    mildlyNonlinear: The problem is considered to be mildly nonlinear and
     *                     NleqSolver starts up with dampingfactor=1. \n
     *    highlyNonlinear: The problem is considered to be highly nonlinear and
     *                     NleqSolver starts up with dampingfactor=1.0e-4. \n
     *    extremelyNonlinear: The problem is considered to be extremely nonlinear
     *                        and NleqSolver starts up with dampingfactor=1.0e-6.
     *                        Moreover, restricted is set automatically to true.
     */
    void setNonlinProblemType(NonlinProblemType nonlinType_ ) { nonlinType=nonlinType_; }

    /**
     * \brief Sets the maximum allowed number of Newton-iterations. The default value is 50.
     * \param maximumIter_  The desired maximum number of Newton-iterations.
     */
    void setMaximumNoIterations( int maximumIter_ ) { maximumNoIterations=maximumIter_; }
    
    /**
     * \brief Sets the fillin-level of the preconditioner PrecondType::ILUK. The default fillin-level is 0.
     * \param preconFillLevel_  The desired PrecondType::ILUK fillin-level.
     */
    void setPreconFillLevel( int preconFillLevel_ ) { preconFillLevel=preconFillLevel_; }
    
    /**
     * \brief Sets the restricted monotonicity-test flag.
     *   The default value depends on the selected NonlinProblemType.
     * \param restricted_  true means the the restricted damping-strategy is applied. \n
     *   This means, that the restricted monotonicity test will be applied for
     *   determination whether the next iterate (and the associate damping factor
     *   lambda) will be accepted. With
     *   \f$ \theta = \parallel  \overline{\Delta x}_{k+1} \parallel / \parallel \Delta x_{k} \parallel \f$,
     *   (where \f$ \Delta x_{k} \f$ denotes the k-th Newton correction, and \f$ \overline{\Delta x}_{k+1} \f$
     *    denotes the sucessive simplified Newton correction)
     *   the condition \f$ \theta <= 1 - \lambda / 4 \f$ must be fit. \n
     *   If set to false, then the standard monotonicity test will be applied, i.e.
     *   the following condition must be fit:
     *   \f$ \theta < 1 \f$ .
     */
    void setRestricted( bool restricted_ ) { restricted=restricted_; }
    
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
     *   The default monitor-output level is none.
     * \param monitorLevel_  The desired monitor-output level.
     */
    void setMonitorLevel( PrintLevel monitorLevel_ ) { monitorLevel=monitorLevel_; }
    
    /**
     * \brief Sets the plot-data-output level.
     *   The default plot-data-output level is none.
     * \param dataLevel_  The desired plot-data-output level.
     */
    void setDataLevel( PrintLevel dataLevel_ ) { dataLevel=dataLevel_; }
    
    /**
     * \brief Sets the prefix of the output-files named \<prefix\>nnn.vtu, 
     *   where n denotes a decimal digit.
     *   The default-prefix is iter_ .
     * \param outFilePrefix_  The desired filename-prefix.
     */
    void setOutFilePrefix( const std::string& outFilePrefix_ ) { outFilePrefix=outFilePrefix_; }

    /**
     * \brief Sets a user-defined scaling-threshold vector.
     * \param *scale_  A pointer to the selected scaling-threshold vector.
     */
    void setScalingVector( std::vector<double>* scale_ ) 
    { 
      delete scale;
      if (scale_) scale=scale_;
    }
    
    /**
     * \brief After calling nleqErr, getInfo returns the info structure contents.
     */
    struct NleqInfo getInfo() { return info; }

    /**
     * \brief Calls the nleqErr damped Newton-iteration solver.
     * \param gridManager The GridManager where the Grid is associated to.
     * \param eq The weak formulation class instance with the problem- and Jacobian functions.
     * \param variableSet The underlying VariableSet.
     * \param x On input, the start-vector for the damped Newton-iteration. \n
     *   On output, the approximate solution vector \f$ x_k \f$, which fits the error 
     *   condition \f$ \parallel x_k - x\star \parallel <= tol \f$ ,
     *   where \f$ x\star \f$  denotes the exact solution, and \f$  \parallel \ldots \parallel \f$
     *   is a scaled Euclidian norm.
     * \param spaces The underlying problem spaces.
     */
    void nleqErr(GridManager<Grid>& gridManager,
                  Equation& eq,
                  VariableSet const& variableSet, 
                  typename VariableSet::VariableSet* x,
                  Spaces const& spaces)
    {
      double lambda, lambdaNew, mue, normFk, normFkPlusOne, normDxkMinusOne, 
             normDxBar, normDx, theta, s;
      double lambdaMin;
      int const nvars = Equation::AnsatzVars::noOfVariables;
      int const n     = variableSet.degreesOfFreedom(0,nvars);
      int iterationCounter=0, fail=0;
      bool reducted;
      std::vector<double> dxdata(n), dxBar(n), xThresh(n), w(n), xData(n);
      std::vector<double> *fxk = new std::vector<double>(n);                 
      int noRHSComputations=0, noJacobianComputations=0, noLinearSolverIterations=0;
      struct NleqData *data = new struct NleqData;
      if (!data)
      { 
        std::cerr << std::endl << " could not allocate struct data" << std::endl;
        info.returnCode=-994; return;
      };
      data->theta     = 0.0;
      data->mode      = initial;

      if ( monitorLevel > 0 ) 
        *monitorStream << std::endl << " nleqErr - Version 2.0" << std::endl;
      if ( maximumNoIterations <= 0 ) maximumNoIterations = MAX_ITER_DEFAULT;
      if      ( nonlinType==mildlyNonlinear ) 
      { 
        lambda = 1.0; lambdaMin = LAMBDA_MIN_DEFAULT;
      }
      else if ( nonlinType==highlyNonlinear ) 
      { 
        lambda = LAMBDA_START_DEFAULT; lambdaMin = LAMBDA_MIN_DEFAULT;
      }
      else if ( nonlinType==extremelyNonlinear ) 
      { 
        lambda = LAMBDA_START_EXTREMELY_DEFAULT;
        lambdaMin = LAMBDA_MIN_EXTREMELY_DEFAULT;
        restricted = true;
      };

      info.returnCode = nleqParcheckAndPrint(n);
      if ( info.returnCode !=0 ) 
      { 
        if (data) delete data;
        return;
      };

      LeafView leafGridView = gridManager.grid().leafGridView();
      Dune::InverseOperatorResult res;
      typedef VariationalFunctionalAssembler<LinearizationAt<Equation> > Assembler;
      Assembler assembler(gridManager,spaces);
      typename VariableSet::VariableSet dx(variableSet);
      typename VariableSet::VariableSet tmp(variableSet);
      typename VariableSet::VariableSet xkPlusOne(variableSet);
      int const neq = Equation::TestVars::noOfVariables;
      size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
      if ( monitorLevel > 0 )
        *monitorStream << " nvars=" << nvars << ", neq=" << neq <<  ", nnz=" << nnz << "\n";
      typedef typename VariableSet::template CoefficientVectorRepresentation<0,neq>::type CoefficientVectors;
      CoefficientVectors theSolution(VariableSet::template CoefficientVectorRepresentation<0,neq>::init(spaces));
      typedef VariationalFunctionalAssembler<LinearizationAt<Equation> > Assembler;
    
      x->write(xData.begin());
      if ( monitorLevel > 1 ) 
         *monitorStream << std::endl
            << "|            Newton iteration                 |iterative linear solver|" << std::endl
            << " iter      normScl(dx)      norm(fk)    lambda   res  niter      rate  "
            << std::endl << std::endl;
      normDx = 0.0;  normDxBar = 0.0;
      w=xData;   xThresh=*scale;
      info.returnCode = 2;
    
      do 
      { 
        nleqErrRescale(*scale,xData,w,xThresh);
        if ( iterationCounter>0 )
        // recompute norms after rescaling
        { 
          normDxkMinusOne = nleqScaledNorm2(dxdata,*scale);
          normDxBar = nleqScaledNorm2(dxBar,*scale);
        };
        assembler.assemble(linearization(eq,*x));
        noRHSComputations++;  noJacobianComputations++;
        CoefficientVectors fx(assembler.rhs());
        fx.write(fxk->begin());
        normFk = nleqNorm2(*fxk);
        AssembledGalerkinOperator<Assembler,0,neq,0,neq> A(assembler);
    
        int iteSteps = 1000;
        double iteEps = xtol*1.0e-2;
        typedef typename VariableSet::template CoefficientVectorRepresentation<0,neq>::type LinearSpace;
        ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,neq> > iluk(A,preconFillLevel,0);
        Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,0);
        theSolution = 0;
        cg.apply(theSolution,fx,res);
        noLinearSolverIterations += res.iterations;
        if ( !res.converged ) { info.returnCode=4; break; };
        dx.data = theSolution.data;
        dx.write(dxdata.begin());
        normDx = nleqScaledNorm2(dxdata,*scale);
        if ( monitorLevel > 1 ) 
          nleqErrMonitor(iterationCounter,normDx,normFk,lambda,' ',res);
        if ( normDx <= xtol ) /* solution found, if condition is true */
        { 
          data->normF     = normFk;     data->normDx = normDx;
          data->normDxBar = normDxBar;  data->lambda = lambda;
          nleqDataOut(iterationCounter,leafGridView,*x,data);
          *x += dx;  iterationCounter++; 
          info.returnCode=0; break; 
        };  
        if ( iterationCounter>0 )
        { 
          for (int i=0;i<n;i++) w[i] = dxBar[i]-dxdata[i];
          s = nleqScaledNorm2(w,*scale);
          mue = (normDxkMinusOne*normDxBar)/(s*normDx) * lambda; 
          lambda = std::min(1.0,mue);
        };
        reducted = false;

        bool doCorrectorLoop=true;
        do
        {
          if ( lambda < lambdaMin ) { info.returnCode=3; break; };  // stop, convergence failure!
          tmp = dx;
          tmp *= -lambda;
          xkPlusOne=*x; xkPlusOne += tmp; /* new trial iterate */
          assembler.assemble(linearization(eq,xkPlusOne),assembler.RHS); 
          noRHSComputations++;
          CoefficientVectors fxq(assembler.rhs());
          fxq.write(fxk->begin());
          
          normFkPlusOne = nleqNorm2(*fxk);
          dxBar=*fxk;
          theSolution = 0;
          cg.apply(theSolution,fxq,res);
          noLinearSolverIterations += res.iterations;
          if ( !res.converged ) { info.returnCode=5; break; };
          tmp.data = theSolution.data;
          tmp.write(dxBar.begin());
          normDxBar = nleqScaledNorm2(dxBar,*scale);
          if ( monitorLevel > 1 ) 
            nleqErrMonitor(iterationCounter,normDxBar,normFkPlusOne,lambda,'*',res);
          theta = normDxBar/normDx;
          s = 1.0-lambda;
          for (int i=0;i<n;i++)  w[i] = dxBar[i]-s*dxdata[i];
          mue = (0.5*normDx*lambda*lambda) / nleqScaledNorm2(w,*scale);
          if ( ( !restricted && theta >= 1.0 ) || 
               (  restricted && theta > 1.0-lambda/4.0) )
          { 
            lambdaNew = std::min(mue,0.5*lambda);
            if ( lambda <= lambdaMin ) lambda = lambdaNew;
            else                        lambda = std::max(lambdaNew,lambdaMin);
            reducted = true;
            continue;
          };
          lambdaNew = std::min(1.0,mue);
          if ( lambda==1.0 && lambdaNew==1.0 )
          {
            if ( normDxBar <= xtol )
            {
              *x += tmp;
              info.returnCode=0;
              break; 
            };
          }
          else
          { 
            if( lambdaNew >= 4.0*lambda && !reducted ) 
            {
              lambda=lambdaNew; continue;
            };
          };
          doCorrectorLoop=false;
        }
        while ( doCorrectorLoop );
        if ( info.returnCode !=2 ) break;
        data->normF     = normFk;     data->normDx = normDx;
        data->normDxBar = normDxBar;  data->lambda = lambda;
        nleqDataOut(iterationCounter,leafGridView,*x,data);
        data->mode = intermediate;
        // save previous iterate for scaling purposes and accept new iterate
        x->write(w.begin()); *x=xkPlusOne; x->write(xData.begin()); 
        // next step
        iterationCounter++;
        normFk = normFkPlusOne;
      }
      while ( iterationCounter <= maximumNoIterations && info.returnCode == 2 );
    
      data->normF     = normFkPlusOne;   data->normDx = normDx;
      data->normDxBar = normDxBar;  data->lambda = lambda;
      data->mode = ( info.returnCode==0 ? solution : final );
      nleqDataOut(iterationCounter,leafGridView,*x,data);
        
      if ( (monitorLevel > 0) && (info.returnCode==0) )
        *monitorStream << std::endl <<
          " solution of nonlinear system obtained within " << iterationCounter << " Newton iteration steps" <<
          std::endl << " using " << noJacobianComputations << " Jacobian assemblies and " << noRHSComputations <<
          " function assemblies" << std::endl;
      
      if ( (errorLevel > 0) && (info.returnCode != 0) )
      {
        switch ( info.returnCode )
        {
          case     2:
            *errorStream << std::endl << 
              " Error - Maximum allowed number of iterations exceeded" << std::endl;
            break;
          case     3:
            *errorStream << std::endl << 
              " Error - no convergence, damping factor became too small" << std::endl;
            break;
          case     4:
            *errorStream << std::endl << 
              " Error - no convergence of iterative linear solver for dx computation" <<
              std::endl;
            break;
          case     5:
            *errorStream << std::endl << 
              " Error - no convergence of iterative linear solver for dxq computation" << std::endl;
            break;
          default   :
            *errorStream << std::endl << " Error, code=" << info.returnCode << 
              ",  subCode=" << fail << std::endl;
        };
      };
      info.subCode = fail;
      delete data;
    
      info.precision = normDx;
      info.normDx    = normDx;
    
      // store info
      info.noIterations       = iterationCounter; 
      info.noFunctionEvaluations = noRHSComputations;
      info.noJacobianEvaluations = noJacobianComputations;
      info.noLiniterIterations = noLinearSolverIterations;
      info.fx         = fxk;
    }
    
    private:
    
    void nleqErrRescale(std::vector<double>& scale, std::vector<double> const & x, 
                         std::vector<double> const& xa, std::vector<double> const& xThresh)
    {
      int const n=scale.size();
      assert ( (x.size()==n) && (xa.size()==n) && (xThresh.size()==n) );
      for (int i=0;i<n;i++) 
        scale[i] = std::max(xThresh[i],std::max(0.5*(std::abs(x[i])+std::abs(xa[i])),SMALL));
    }
    
    void nleqErrMonitor(int const iterationCounter, double const normDx,
                         double const normF, double const lambda, char const indicator,
                         Dune::InverseOperatorResult res)
    {
      *monitorStream << " " << std::setw(4) << iterationCounter << "  " << std::setw(1) << indicator << "  "
        << std::setw(12) << std::setprecision(5) << std::scientific << normDx << "  "
        << std::setw(12) << std::setprecision(5) << normF << "  ";
      monitorStream->unsetf(std::ios::fixed | std::ios::scientific);
      *monitorStream << std::setw(8) << std::fixed << lambda << "  " << std::setw(4)
         << (res.converged?"conv":"fail") << "  " << std::setw(5) << res.iterations << "  "
         << std::setw(8) << res.conv_rate << std::endl;
      monitorStream->unsetf(std::ios::fixed | std::ios::scientific);
    }
    
    double nleqScaledNorm2(std::vector<double> const& v, std::vector<double> const& scale)
    {
      int const n=v.size();
      assert ( scale.size()==n );
      double t, rval = 0.0;
      for (int i=0;i<n;i++) {t=v[i]/scale[i]; rval += t*t;};
      return sqrt( rval / (double)n );
    }
    
    double nleqNorm2(std::vector<double> const& v)
    {
      int const n=v.size();
      double rval = 0.0;
      for (int i=0;i<n;i++) rval += v[i]*v[i];
      return sqrt( rval / (double)n );
    }
    
    void nleqDataOut(int const iterationCounter, LeafView const leafGridView, VariableSetRepresentation const x,
                     struct NleqData const *data)
    {
      if ( (dataLevel==none) || ( (dataLevel==minimum) && (data->mode==intermediate) ) )
        return;
      // output of solution in VTK format for visualization,
      // the data are written as ascii stream into file temperature.vtu,
      // possible is also binary
      std::ostringstream fname;
      fname << outFilePrefix;
      fname.width(3);
      fname.fill('0');
      fname.setf(std::ios_base::right,std::ios_base::adjustfield);
      fname << iterationCounter;
      fname.flush();
      writeVTKFile(x,fname.str(),IoOptions());
    }
    
    int  nleqParcheckAndPrint(int const n)
    #define TOLMIN EPMACH*1.0e2
    #define TOLMAX 1.0e-1
    {  
      double large = 1.0/SMALL, small = SMALL, defaultScale;
      if ( n<=0 ) 
      { 
        if ( errorLevel>0 )
          *errorStream << std::endl << " Error - Number of Equations/unknowns must be >0" << std::endl;
        return 20;
      };
      if ( xtol <= 0 )
      { 
        if ( errorLevel>0 )
            *errorStream << std::endl << " Error - xtol must be positive" << std::endl;
          return 21;
      }
      else
      { 
        if ( xtol > TOLMAX ) 
        {
          xtol = TOLMAX;
          if ( errorLevel>1 )
            *errorStream << std::endl <<
              " User prescribed RTOL decreased to reasonable largest value RTOL=" << 
              xtol << std::endl;
        }
        else if ( xtol < TOLMIN ) 
        { 
          xtol = TOLMIN;
          if ( errorLevel>1 )
            *errorStream << std::endl <<
              "User prescribed RTOL increased to reasonable smallest value RTOL=" <<
              xtol << std::endl;
        };
      };
      if ( nonlinType >= highlyNonlinear ) 
        defaultScale = xtol;
      else
        defaultScale = 1.0;
    
      if ( scale->empty() ) 
      { 
        scale->resize(n); 
        for (int i=0;i<n;i++) (*scale)[i]=defaultScale;
      }
      else
      {
        for (int i=0;i<n;i++) 
        {
          if ( (*scale)[i] < 0.0 ) 
          { 
            if ( errorLevel>0 )
              *errorStream << std::endl <<
                " Error - negative value scale[" << i << "] supplied" << std::endl;
            return 22;
          }
          else if ( (*scale)[i] == 0.0 ) 
            (*scale)[i] = defaultScale;
          else if ( (*scale)[i] < SMALL )
          {
            if ( errorLevel>1 )
              *errorStream << std::endl <<
              " Warning scale[" << i << "] too small - increased to " << SMALL << std::endl;
            (*scale)[i] = small;
          }
          else if ( (*scale)[i] > large )
          {
            if ( errorLevel>1 )
              *errorStream << std::endl <<
                " Warning scale[" << i << "] too large - decreased to " << large << std::endl;
            (*scale)[i] = large;
          };
        };
      };
      if ( monitorLevel==0 ) return 0;
      *monitorStream << std::endl << " Problem dimension: n = " << n << std::endl;
      *monitorStream << std::endl << " Prescribed relative precision: " << xtol 
        << std::endl;
      *monitorStream << " The problem is specified as being ";
      switch ( nonlinType )
      { 
        case mildlyNonlinear:    *monitorStream << "mildly nonlinear" << std::endl; 
                                  break;
        case highlyNonlinear:    *monitorStream << "highly nonlinear" << std::endl;
                                  break;
        case extremelyNonlinear: *monitorStream << "extremely nonlinear" << std::endl; 
                                  break;
      };
      if ( restricted )
        *monitorStream << " The restricted monotonicity test will be applied" << std::endl;
      else
        *monitorStream << " The standard monotonicity test will be applied" << std::endl;
    
     *monitorStream << " The maximum permitted number of iteration steps is: " <<
       maximumNoIterations << std::endl;
     return 0;
    }
    
    double xtol;
    int maximumNoIterations, preconFillLevel;
    bool restricted;
    std::ostream *errorStream, *monitorStream;
    enum {StandardScale=0,StartValueScale=1} scaleopt;
    PrintLevel errorLevel, monitorLevel, dataLevel;
    NonlinProblemType nonlinType;
    std::vector<double>* scale;
    std::string outFilePrefix;
    struct NleqInfo info;
    
  };
}
#endif
