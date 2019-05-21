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

#ifndef UMFPACK_TOOLS_HH_
#define UMFPACK_TOOLS_HH_

#include <vector>

/**
 * \ingroup linalgbasic
 * \brief Triplet to compressed column storage conversion .
 * This is just a frontent to DirectType::UMFPACK utility routines.
 * \param rows number of rows
 * \param cols number of columns
 * \param[in] ridx row indices of matrix entries
 * \param[in] cidx column indices of matrix entries
 * \param[in] values matrix entries
 * \param[out] Ap column start array
 * \param[out] Ai row indices
 * \param[out] Az matrix entries
 * 
 * The output arrays \arg Ap, \arg Ai, and \arg Az are resized as needed.
 */
void umfpack_triplet_to_col(
    int rows,
    int cols,
    std::vector<long> const& ridx,
    std::vector<long> const& cidx,
    std::vector<double> const& values,
    std::vector<long>& Ap,
    std::vector<long>& Ai,
    std::vector<double>& Az
);

/**
 * \ingroup linalgbasic
 * \brief Triplet to compressed column storage conversion .
 * This is just a frontent to DirectType::UMFPACK utility routines.
 */
void umfpack_triplet_to_col(
    int rows,
    int cols,
    std::vector<int> const& ridx,
    std::vector<int> const& cidx,
    std::vector<double> const& values,
    std::vector<int>& Ap,
    std::vector<int>& Ai,
    std::vector<double>& Az
);

/**
 * \ingroup linalgbasic
 * \brief Compressed column storage to triplet storage conversion.
 * This is just a frontent to DirectType::UMFPACK utility routines.
 */
void umfpack_col_to_triplet(
    std::vector<long> const& Ap,
    std::vector<long>& Ti
);

/**
 * \ingroup linalgbasic
 * \brief Compressed column storage to triplet storage conversion.
 * This is just a frontent to DirectType::UMFPACK utility routines.
 */
void umfpack_col_to_triplet(
    std::vector<int> const& Ap,
    std::vector<int>& Ti
);

#endif /* UMFPACK_TOOLS_HH_ */
