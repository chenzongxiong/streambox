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

#ifndef KASKADE_ENUMS_HH_
#define KASKADE_ENUMS_HH_

#include <iostream>

enum class PrecondType { NONE, JACOBI, PARTIAL, ILUT, ILUK, ARMS, INVERSE, ADDITIVESCHWARZ,
                   BOOMERAMG, EUCLID, TRILINOSML, SSOR, ICC0, ICC, ILUKS, HB, DIRECT };

std::ostream& operator<<(std::ostream& out, PrecondType direct);

enum class IterateType { CG, BICGSTAB, GMRES, PCG, APCG, SGS };

std::ostream& operator<<(std::ostream& out, IterateType direct);

/**
 * \ingroup direct
 * \brief Available direct solvers for linear equation systems.
 */
enum class DirectType  { ANY=-1, UMFPACK, PARDISO, MUMPS, SUPERLU, UMFPACK3264, UMFPACK64 };

std::ostream& operator<<(std::ostream& out, DirectType direct);

/**
 * \ingroup linalg
 * \brief Characterizations of sparse matrix properties.
 */
enum class MatrixProperties { GENERAL, SYMMETRICSTRUCTURE, SYMMETRIC, POSITIVEDEFINITE };

std::ostream& operator<<(std::ostream& out, MatrixProperties direct);

#endif /* KASKADE_ENUMS_HH_ */
