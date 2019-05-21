/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ILUPRECOND_HH
#define ILUPRECOND_HH

#include "globheads.h"
#include "protospp.h" 
#include "ios.h"
extern "C" void set_arms_pars(io_t* io,  int Dscale, int *ipar, double *tolcoef, int *lfil);

#include <boost/timer/timer.hpp>
#include "dune/istl/preconditioners.hh"
#include "linalg/triplet.hh"

//---------------------------------------------------------------------
//---------------------------------------------------------------------

namespace Kaskade
{
/**
 * \ingroup linalg
 * 
 * PrecondType::ILUT Preconditioner from Saads ITSOL library
 */
template <class Op>
class ILUTPreconditioner: public Dune::Preconditioner<typename Op::Range, typename Op::Range>
{
  typedef typename Op::Range Range;
  typedef Range               Domain;
  typedef typename Op::field_type field_type;

public:
  static int const category = Dune::SolverCategory::sequential;

  ILUTPreconditioner(Op& op, int lfil=240, double tol=1.0e-2, int verbosity=2) {
    boost::timer::cpu_timer iluTimer;
    std::cout << "ilut constructor: ";
    MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
    int n = A.nrows(), nnz = A.nnz(), ierr;
    FILE *flog = stdout;

    csmat = (csptr)Malloc( sizeof(SparMat), "main" );
    ierr = COOcs(n, nnz,  &A.data[0], &A.cidx[0],  &A.ridx[0], csmat);
    if( ierr != 0 )
      {
        printf(" *** ILU error - COOcs code %d csmat=%p\n", ierr, csmat);
        throw ierr;
      }

    lu = (iluptr)Malloc(sizeof(ILUSpar), "main");
    ierr = ilut(csmat, lu, lfil, tol, flog);
    if( ierr != 0 )
      {
        printf(" *** PrecondType::ILUT error - ilut code %d lu=%p\n", ierr, lu);
        throw ierr;
      }
    
    xx   = (double *)Malloc(n*sizeof(double), "main");
    yy   = (double *)Malloc(n*sizeof(double), "main");
    
    if ( verbosity>=2 )
    {
	std::cout << "PrecondType::ILUT: n=" << n << ", nnz=" << nnz << " dropTol=" << tol << ", fillfac=" 
	          << 1.0*nnz_ilu(lu)/(nnz+1.0) << ", lfil=" << lfil << ", time=" << (double)(iluTimer.elapsed().user)/1e9 << "s\n";
	}
  }
  ~ILUTPreconditioner()
    {
      cleanILU(lu);
      cleanCS(csmat);
      free(xx);
      free(yy);
    }

  virtual void pre (Domain&, Range&) {}
  virtual void post (Domain&) {}
  
  virtual void apply (Domain& x, Range const& y) {
    y.write(yy);
    lusolC(yy, xx, lu);
    x.read(xx);
  }

private:
  csptr csmat;
  iluptr lu;
  double *xx, *yy;
};

//---------------------------------------------------------------------

/**
 * \ingroup linalg
 * 
 * PrecondType::ILUT Preconditioner from Saads ITSOL library
 */
template <class Op>
class ILUKPreconditioner: public Dune::Preconditioner<typename Op::Range, typename Op::Range>
{
  typedef typename Op::Range Range;
  typedef Range               Domain;
  typedef typename Op::field_type field_type;

public:
  static int const category = Dune::SolverCategory::sequential;

  ILUKPreconditioner(Op& op, int fill_lev=3, int verbosity=2) {
    boost::timer::cpu_timer iluTimer;
    MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
    int n = A.nrows(), nnz = A.nnz(), ierr;
    FILE *flog = stdout;

    csmat = (csptr)Malloc( sizeof(SparMat), "main" );
    ierr = COOcs(n, nnz,  &A.data[0], &A.cidx[0],  &A.ridx[0], csmat);
    if( ierr != 0 )
      {
        printf(" *** ILU error - COOcs code %d csmat=%p\n", ierr, csmat);
        throw ierr;
      }

    lu = (iluptr)Malloc(sizeof(ILUSpar), "main");
    ierr = ilukC(fill_lev, csmat, lu, flog);
    if( ierr != 0 )
      {
        printf(" *** PrecondType::ILUK error - ilut code %d lu=%p\n", ierr, lu);
        throw ierr;
      }

    xx   = (double *)Malloc(n*sizeof(double), "main");
    yy   = (double *)Malloc(n*sizeof(double), "main");

    if ( verbosity>=2 )
    {
	std::cout << "PrecondType::ILUK: n=" << n << ", nnz=" << nnz << ", fillfac="
	          << 1.0*nnz_ilu(lu)/(nnz+1.0) << ", fill_lev=" 
	          << fill_lev << ", time=" << (double)(iluTimer.elapsed().user)/1e9 << "s\n";
    }
  }
  ~ILUKPreconditioner()
    {
      cleanILU(lu);
      cleanCS(csmat);
      free(xx);
      free(yy);
    }

  virtual void pre (Domain&, Range&) {}
  virtual void post (Domain&) {}
  
  virtual void apply (Domain& x, Range const& y) {
    y.write(yy);
    lusolC(yy, xx, lu);
    x.read(xx);
  }

private:
  csptr csmat;
  iluptr lu;
  double *xx, *yy;
};
//---------------------------------------------------------------------

/**
 * \ingroup linalg
 * 
 * PrecondType::ILUT Preconidoner from Saads ITSOL library
 */
template <class Op>
class ARMSPreconditioner: public Dune::Preconditioner<typename Op::Range, typename Op::Range>
{
  typedef typename Op::Range Range;
  typedef Range               Domain;
  typedef typename Op::field_type field_type;

public:
  static int const category = Dune::SolverCategory::sequential;

  ARMSPreconditioner(Op& op, int lfil=4, double tol=1.0e-2, int lev_reord=1, double tolind = 0.2, int verbosity=2) {
    boost::timer::cpu_timer iluTimer;
    MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
    int n = A.nrows(), nnz = A.nnz(), ierr, diagscal = 1;

    int lfil_arr[7]; 
    double droptol[7], dropcoef[7];
    int ipar[18];
    io_t io;
    FILE *flog = stdout;

    csmat = (csptr)Malloc( sizeof(SparMat), "main" );
    ierr = COOcs(n, nnz,  &A.data[0], &A.cidx[0],  &A.ridx[0], csmat);
    if( ierr != 0 )
      {
        printf(" *** PrecondType::ARMS error - COOcs code %d csmat=%p\n", ierr, csmat);
        throw ierr;
      }

    memset(&io, 0, sizeof(io) );
    io.perm_type = 0;
    io.Bsize = 400;
    set_arms_pars(&io, diagscal, ipar, dropcoef, lfil_arr);
    for (int j=0; j<7; j++)
      {
         lfil_arr[j] = lfil*((int) nnz/n); 
         droptol[j] =  tol*dropcoef[j];
      }
    ipar[1] = lev_reord; 
    ArmsSt = (arms) Malloc(sizeof(armsMat),"main:ArmsSt");
    setup_arms(ArmsSt);

    ierr = arms2(csmat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);
    if( ierr != 0 )
      {
        printf(" *** PrecondType::ARMS error - arms2 code %d ArmsSt=%p\n", ierr, ArmsSt);
        throw ierr;
      }

    yy   = (double *)Malloc(n*sizeof(double), "main");

    if ( verbosity>=2 )
    {
	std::cout << "PrecondType::ARMS: n=" << n << ", nnz=" << nnz << ", dropTol=" << tol
	          << ", lev_reord=" << lev_reord
	          << ", lfil=" << lfil << ", time=" << (double)(iluTimer.elapsed().user)/1e9 << "s\n";
    }
  }
  ~ARMSPreconditioner()
    {
      cleanARMS(ArmsSt); 
      cleanCS(csmat);
      free(yy);
    }

  virtual void pre (Domain&, Range&) {}
  virtual void post (Domain&) {}
  
  virtual void apply (Domain& x, Range const& y) {
    y.write(yy);
    armsol2(yy, ArmsSt);
    x.read(yy);
  }

private:
  csptr csmat;
  arms ArmsSt;
  double *yy;
};

}  // namespace Kaskade
#endif
