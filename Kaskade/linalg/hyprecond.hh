/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BOOMERAMG_HH
#define BOOMERAMG_HH

#include <math.h>
#include "globheads.h"
// #include "defs.h"   <--- defines min() and max() macros colliding with std::min() etc. Do NOT include in header files.
#include "protospp.h" 
#ifndef MAX_HBNAME
#include "ios.h"
#endif

#include "HYPRE_sstruct_ls.h"
// delete preprocessor macros min and max which were defined in itsol-1/hypre-2.6.0b headers,
// since these macros may interfere with subsequently included Kaskade7 headers or program code
#undef min
#undef max

#include <vector>

#include "dune/istl/preconditioners.hh"

#include "linalg/triplet.hh"
#include "linalg/umfpack_solve.hh"

namespace Kaskade
{
//---------------------------------------------------------------------
//---------------------------------------------------------------------

/**
 * \ingroup linalg
 * 
 * PrecondType::ILUT Preconditioner from Saads ITSOL library
 */
template <class Op>
class BoomerAMG: public Dune::Preconditioner<typename Op::range_type, typename Op::range_type>
{
  typedef typename Op::range_type range_type;
  typedef range_type               domain_type;
  typedef typename Op::field_type field_type;

public:
  static int const category = Dune::SolverCategory::sequential;

  BoomerAMG(Op& op, int steps_=1, int coarsentype=21, int interpoltype=0, double tol=1.0e-8,
            int cycleType=1, int relaxType=2, double strongThreshold=0.3,
            int variant=0, int overlap=1, int syseqn=1, int verbosity_=2):
            op(op), steps(steps_), verbosity(verbosity_), eqnIndex(0)
            {
    MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
    n = A.nrows();
    int nnz = A.nnz(); /* ierr; */
    std::vector<int> Ap, Ai; 
    std::vector<double> Az;
//    FILE *flog = stdout;

    if ( verbosity>=2 )
    {
      std::cout << "HYPRE_BoomerAMG: n=" << n << ", nnz=" << nnz << 
        " steps=" << steps <<
        ", coarsentype=" << coarsentype << ", interpolationtype=" << 
        interpoltype << "\n";
    };
    Ap.resize(n+1);  Ai.resize(nnz);  Az.resize(nnz);
    umfpack_triplet_to_col(n,n,A.ridx,A.cidx,A.data,
                           Ap,Ai,Az);
    for (int i=0;i<n;i++) Ap[i]=Ap[i+1]-Ap[i];
    ixx   = (int *)Malloc(n*sizeof(int), "main");
    for (int i=0;i<n;i++) ixx[i]=i;

    HYPRE_IJMatrixCreate(0, 0, n-1, 0, n-1, &ij_matrix); 
    HYPRE_IJMatrixSetPrintLevel(ij_matrix,10); /* Debug */
    HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR); 
    HYPRE_IJMatrixInitialize(ij_matrix); 
    /* set matrix coefficients */ 
    HYPRE_IJMatrixSetValues(ij_matrix,n,&Ap[0],ixx,&Ai[0],&Az[0]); 
    HYPRE_IJMatrixAssemble(ij_matrix); 
    HYPRE_IJMatrixGetObject(ij_matrix, (void **) &parcsr_matrix);

    HYPRE_BoomerAMGCreate(&preconditioner);
    HYPRE_BoomerAMGSetPrintLevel(preconditioner,1); /* Debug */
    HYPRE_BoomerAMGSetCoarsenType(preconditioner,coarsentype);
    HYPRE_BoomerAMGSetInterpType(preconditioner,interpoltype);
    HYPRE_BoomerAMGSetMaxIter(preconditioner,steps);
    HYPRE_BoomerAMGSetTol(preconditioner,tol);
    HYPRE_BoomerAMGSetStrongThreshold(preconditioner, strongThreshold); //RR 0.6 3D
    HYPRE_BoomerAMGSetCycleType(preconditioner, cycleType); //RR 2 W-Zyklus
    HYPRE_BoomerAMGSetRelaxType(preconditioner, relaxType); //RR 3 PrecondType::SSOR smoother
    HYPRE_BoomerAMGSetVariant (preconditioner, variant);
    HYPRE_BoomerAMGSetOverlap (preconditioner, overlap) ;
    if (syseqn>1)
      {
         HYPRE_BoomerAMGSetNodal(preconditioner, 1);
         HYPRE_BoomerAMGSetNumFunctions(preconditioner, syseqn);
         int i, k, lng = n/syseqn;
         eqnIndex = (int *)Malloc(n*sizeof(int), "main");
         for (i=0; i<syseqn; i++)
           for (k=0; k<lng; k++)
             {
               eqnIndex[i*lng+k] = i;
             }
         HYPRE_BoomerAMGSetDofFunc(preconditioner, eqnIndex);
      }

//    HYPRE_BoomerAMGSetCycleRelaxType(preconditioner,6,1);
//    HYPRE_BoomerAMGSetCycleRelaxType(preconditioner,6,2);
//    HYPRE_BoomerAMGSetCycleRelaxType(preconditioner,6,3);
    HYPRE_BoomerAMGSetup(preconditioner,parcsr_matrix,0,0);
    
    xxraw = (double *)Malloc(n*sizeof(double), "main");
    yyraw = (double *)Malloc(n*sizeof(double), "main");

    HYPRE_IJVectorCreate(0, 0, n-1, &xxij); 
    HYPRE_IJVectorSetPrintLevel(xxij,10); /* Debug */
    HYPRE_IJVectorSetObjectType(xxij, HYPRE_PARCSR); 
    HYPRE_IJVectorInitialize(xxij); 
    HYPRE_IJVectorAssemble(xxij); 
    HYPRE_IJVectorGetObject(xxij, (void **) &xx); 

    HYPRE_IJVectorCreate(0, 0, n-1, &yyij); 
    HYPRE_IJVectorSetPrintLevel(yyij,10); /* Debug */
    HYPRE_IJVectorSetObjectType(yyij, HYPRE_PARCSR); 
    HYPRE_IJVectorInitialize(yyij); 
    HYPRE_IJVectorAssemble(yyij); 
    HYPRE_IJVectorGetObject(yyij, (void **) &yy); 
    }
  ~BoomerAMG()
    {
      HYPRE_IJMatrixDestroy (ij_matrix);
      HYPRE_BoomerAMGDestroy(preconditioner);
      HYPRE_IJVectorDestroy (xxij);
      HYPRE_IJVectorDestroy (yyij);
      free(xxraw);
      free(yyraw);
      free(ixx);
      if (eqnIndex) free(eqnIndex);
    }

  virtual void pre (domain_type&, range_type&) {}
  virtual void post (domain_type&) {}
  
  virtual void apply (domain_type& x, range_type const& y) {
/* set vector values */ 
    y.write(yyraw);
    for (int i=0; i<n; i++)
      xxraw[i] = 0.0;

    HYPRE_IJVectorSetValues(yyij, n, ixx, yyraw); 
    HYPRE_IJVectorSetValues(xxij, n, ixx, xxraw); 
    HYPRE_BoomerAMGSolve(preconditioner,parcsr_matrix,yy,xx);
    if (steps>1)
      {
        int num_iterations;
        double rel_resid_norm;
        HYPRE_BoomerAMGGetNumIterations (preconditioner,&num_iterations); 
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm (preconditioner, &rel_resid_norm);
        if ( verbosity>=2 )
        {
          std::cout << "HYPRE_BoomerAMG1: " << num_iterations << " steps, rel_resid_norm="
                    << rel_resid_norm << std::endl;
        }
      }
    HYPRE_IJVectorGetValues(xxij,n,ixx,xxraw);
    x.read(xxraw);
  }

private:
    Op& op;
    HYPRE_Solver preconditioner;
    HYPRE_IJMatrix ij_matrix; 
    HYPRE_ParCSRMatrix parcsr_matrix;
    HYPRE_IJVector xxij, yyij;
    HYPRE_ParVector xx, yy;
    int n, steps, verbosity;
    double *xxraw, *yyraw;
    int *ixx, *eqnIndex;
};

//---------------------------------------------------------------------

/**
 * \ingroup linalg
 * 
 * PrecondType::ILUT Preconidoner from Saads ITSOL library
 */
template <class Op>
class Euclid: public Dune::Preconditioner<typename Op::range_type, typename Op::range_type>
{
  typedef typename Op::range_type range_type;
  typedef range_type               domain_type;
  typedef typename Op::field_type field_type;

public:
  static int const category = Dune::SolverCategory::sequential;

  Euclid(Op& op, int level=1, double droptol=0.0, int printlevel=0, int bj=0, 
         int verbosity=2): op(op)
  {
    MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
    n = A.nrows();
    int nnz = A.nnz(); /* ierr; */
    std::vector<int> Ap, Ai; 
    std::vector<double> Az;
//    FILE *flog = stdout;

    if ( verbosity>=2 )
    {
      std::cout << "PrecondType::EUCLID: n=" << n << ", nnz=" << nnz << ", level=" << 
          level <<  ", droptol="  << droptol << ", printlevel=" << printlevel <<
          ", bj=" << bj << "\n";
    };
    Ap.resize(n+1);  Ai.resize(nnz);  Az.resize(nnz);
    umfpack_triplet_to_col(n,n,A.ridx,A.cidx,A.data,
                           Ap,Ai,Az);
    for (int i=0;i<n;i++) Ap[i]=Ap[i+1]-Ap[i];
    ixx   = (int *)Malloc(n*sizeof(int), "main");
    for (int i=0;i<n;i++) ixx[i]=i;

    HYPRE_IJMatrixCreate(0, 0, n-1, 0, n-1, &ij_matrix); 
    HYPRE_IJMatrixSetPrintLevel(ij_matrix,10); /* Debug */
    HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR); 
    HYPRE_IJMatrixInitialize(ij_matrix); 
    /* set matrix coefficients */ 
    HYPRE_IJMatrixSetValues(ij_matrix,n,&Ap[0],ixx,&Ai[0],&Az[0]); 
    HYPRE_IJMatrixAssemble(ij_matrix); 
    HYPRE_IJMatrixGetObject(ij_matrix, (void **) &parcsr_matrix);

    HYPRE_EuclidCreate(0,&preconditioner);
    HYPRE_EuclidSetStats(preconditioner,printlevel);
    HYPRE_EuclidSetMem(preconditioner,printlevel);
    HYPRE_EuclidSetLevel(preconditioner,level);
//    HYPRE_EuclidSetSparseA(preconditioner,droptol);
    HYPRE_EuclidSetRowScale(preconditioner,0);
    if (bj)
      HYPRE_EuclidSetBJ(preconditioner,bj);
    if (droptol!=0.0)
      HYPRE_EuclidSetILUT(preconditioner,droptol);
    HYPRE_EuclidSetup(preconditioner,parcsr_matrix,0,0); 

    xxraw   = (double *)Malloc(n*sizeof(double), "main");
    yyraw   = (double *)Malloc(n*sizeof(double), "main");

    HYPRE_IJVectorCreate(0, 0, n-1, &xxij); 
    HYPRE_IJVectorSetPrintLevel(xxij,10); /* Debug */
    HYPRE_IJVectorSetObjectType(xxij, HYPRE_PARCSR); 
    HYPRE_IJVectorInitialize(xxij); 
    HYPRE_IJVectorAssemble(xxij); 
    HYPRE_IJVectorGetObject(xxij, (void **) &xx); 

    HYPRE_IJVectorCreate(0, 0, n-1, &yyij); 
    HYPRE_IJVectorSetPrintLevel(yyij,10); /* Debug */
    HYPRE_IJVectorSetObjectType(yyij, HYPRE_PARCSR); 
    HYPRE_IJVectorInitialize(yyij); 
    }
    
  ~Euclid()
    { free(xxraw);
      free(yyraw);
      free(ixx);
      HYPRE_EuclidDestroy(preconditioner);
      HYPRE_IJMatrixDestroy(ij_matrix);
      HYPRE_IJVectorDestroy(xxij);
      HYPRE_IJVectorDestroy(yyij);
    }

  virtual void pre (domain_type&, range_type&) {}
  virtual void post (domain_type&) {}
  
  virtual void apply (domain_type& x, range_type const& y) {
/* set vector values */ 
    y.write(yyraw);
    HYPRE_IJVectorSetValues(yyij, n, ixx, yyraw); 
    HYPRE_IJVectorAssemble(yyij); 
    HYPRE_IJVectorGetObject(yyij, (void **) &yy); 
    
    HYPRE_EuclidSolve(preconditioner,parcsr_matrix,yy,xx);
    
    HYPRE_IJVectorGetValues(xxij,n,ixx,xxraw);
    x.read(xxraw);
  }

private:
    Op& op;
    HYPRE_Solver preconditioner;
    HYPRE_IJMatrix ij_matrix; 
    HYPRE_ParCSRMatrix parcsr_matrix;
    HYPRE_IJVector xxij, yyij;
    HYPRE_ParVector xx, yy;
    int n;
    double *xxraw, *yyraw;
    int *ixx;
};

}  // namespace Kaskade
#endif
