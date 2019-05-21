/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * \file
 * \author Martin Weiser, modified by Felix Lehmann
 */

#ifndef CONJUGATION_HH
#define CONJUGATION_HH

#include <cassert>
#include <memory>

#include <boost/timer/timer.hpp>

#include "dune/istl/bcrsmatrix.hh"
#include "dune/istl/matrixindexset.hh"

#include "linalg/localMatrices.hh"
#include "linalg/threadedMatrix.hh"

#include "utilities/timing.hh"

namespace Kaskade
{
  /**
   * \brief Creates the sparsity pattern of \f$ P^T A P\f$.
   *
   * Out now: onlyLowerTriangle. If parameter is set true, we will only touch the lower triangle of A
   * and only create the lower triangle of the resulting matrix P^T A P
   */
  template <class Scalar, class Entry>
  std::unique_ptr<Dune::MatrixIndexSet> conjugationPattern(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
							   Dune::BCRSMatrix<Entry> const& A ,
							   bool onlyLowerTriangle = false)
      {
    assert(A.N()==A.M());
    assert(A.N()==P.N());

    typedef Dune::BCRSMatrix<Entry> MatA;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > MatP;

    // If C = P^T A P, we have that C_{ij} = \sum_{k,l} P_{ki} P_{lj} A_{kl}. Hence, the entry A_{kl} contributes to
    // all C_{ij} for which there are nonzero entries P_{ki} and P_{lj} in the rows k and l of P. Thus we can simply
    // run through all nonzeros A_{kl} of A, look up the column indices i,j of rows k and l of P, and flag C_{ij}
    // as nonzero.

    std::unique_ptr<Dune::MatrixIndexSet> nzC(new Dune::MatrixIndexSet(P.M(),P.M()));


    if(onlyLowerTriangle == false)
    {
      // Step through all entries of A
      for (int k=0; k<A.N(); ++k)
        for (typename MatA::ConstColIterator ca=A[k].begin(); ca!=A[k].end(); ++ca)
        {
          int const l = ca.index();
          // Step through all entries of rows k and l of P and add entry
          for (typename MatP::ConstColIterator cpk=P[k].begin(); cpk!=P[k].end(); ++cpk)
            for (typename MatP::ConstColIterator cpl=P[l].begin(); cpl!=P[l].end(); ++cpl)
              nzC->add(cpk.index(),cpl.index());
        }
    }
    else
    {
      // Step through all entries of A
      for (int k=0; k<A.N(); ++k)
        for (typename MatA::ConstColIterator ca=A[k].begin(); ca!=A[k].end(); ++ca)
        {
          int const l = ca.index();
          // Step through all entries of rows k and l of P and add entry
          for (typename MatP::ConstColIterator cpk=P[k].begin(); cpk!=P[k].end(); ++cpk)
            for (typename MatP::ConstColIterator cpl=P[l].begin(); cpl!=P[l].end(); ++cpl)
            {
              if( cpk.index() >= cpl.index() )
                nzC->add(cpk.index(),cpl.index());
              else
                nzC->add(cpl.index(),cpk.index());
            }
        }
    }

    return nzC;

    // An alternative way of computing the sparsity pattern would be to use that the nonzero entries j in column i
    // of C are exactly those for which there is k with (nonzero P_{jk} and there is l with (nonzero P_{il} and A_{lk})).
    // Hence we can obtain the column index set J directly by the following steps:
    // (i) find all l with P_{li} nonzero -> L  [requires to access columns of P - compute the transpose patterns once]
    // (ii) find all k with A_{lk} nonzero for some l in L -> K  [probably sorting K and removing doubled entries would be a good idea here]
    // (iii) find all j with P_{kj} nonzero for some k in K -> J
    // Compared to the above implementation this would have the following (dis)advantages
    // + easy to do in parallel (since write operations are separated)
    // + fewer scattered write accesses to memory
    // - more complex implementation
    // - requires the transpose pattern of P
      }

  /**
   * \brief Creates the conjugation product \f$ P^T A P\f$.
   *
   * Note that for typical multigrid Galerkin projections the memory access patterns of the conjugation are quite
   * bad for performance. Consider assembling the projected matrix directly on the coarser discretization.
   *
   * \param onlyLowerTriangle. If true, we will only touch the lower triangle of A
   * and only create the lower triangle of the resulting matrix P^T A P
   */
  template <class IndexP, class EntryP, class IndexA, class EntryA>
  NumaBCRSMatrix<EntryA,IndexA>
  conjugation(NumaBCRSMatrix<EntryP,IndexP> const& P, NumaBCRSMatrix<EntryA,IndexA> const& A, bool onlyLowerTriangle = false)
  {
    assert(A.N()==A.M());
    assert(A.N()==P.N());

    Timings& timer = Timings::instance();


    // First create the sparsity pattern
    timer.start("conjugation pattern");
    NumaCRSPatternCreator<IndexA> creator(P.M(),P.M(),onlyLowerTriangle);

    // If C = P^T A P, we have that C_{ij} = \sum_{k,l} P_{ki} P_{lj} A_{kl}. Hence, the entry A_{kl} contributes to
    // all C_{ij} for which there are nonzero entries P_{ki} and P_{lj} in the rows k and l of P. Thus we can simply
    // run through all nonzeros A_{kl} of A, look up the column indices i,j of rows k and l of P, and flag C_{ij}
    // as nonzero.

    { // just a new scope
      // helper routine for extracting all column indices of a row
      auto getColumnIndices = [] (auto const& row, std::vector<IndexP>& ci)
      {
        ci.clear();
        for (auto i=row.begin(); i!=row.end(); ++i)       // indices i for which Pki != 0
          ci.push_back(i.index());
      };
      std::vector<IndexP> is, js;

      // Step through all entries of A
      for (IndexA k=0; k<A.N(); ++k)
      {
        getColumnIndices(P[k],is);                          // indices i for which Pki != 0

        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          IndexA const l = ca.index();
          getColumnIndices(P[l],js);                        // indices j for which Plj != 0

          // add all combinations i,j
          creator.addElements(std::begin(is),std::end(is),std::begin(js),std::end(js));
          if (onlyLowerTriangle && k>l)                     // subdiagonal entry (k,l) of A -> entry (l,k) must be treated implicitly:
            // add all combinations (j,i)
            creator.addElements(std::begin(js),std::end(js),std::begin(is),std::end(is));
        }
      }
    }
    timer.stop("conjugation pattern");

    // An alternative way of computing the sparsity pattern would be to use that the nonzero entries j in column i
    // of C are exactly those for which there is k with (nonzero P_{jk} and there is l with (nonzero P_{il} and A_{lk})).
    // Hence we can obtain the column index set J directly by the following steps:
    // (i) find all l with P_{li} nonzero -> L  [requires to access columns of P - compute the transpose patterns once]
    // (ii) find all k with A_{lk} nonzero for some l in L -> K  [probably sorting K and removing doubled entries would be a good idea here]
    // (iii) find all j with P_{kj} nonzero for some k in K -> J
    // Compared to the above implementation this would have the following (dis)advantages
    // + easy to do in parallel (since write operations are separated)
    // + fewer scattered write accesses to memory
    // - more complex implementation
    // - requires the transpose pattern of P

    // Create the sparse matrix.
    timer.start("matrix creation");
    NumaBCRSMatrix<EntryA,IndexA> pap(creator);
    timer.stop("matrix creation");

    // Fill the sparse matrix PAP. This is done as before by stepping through all Akl entries and scatter
    // Pki*Plj*Akl into PAPij.
    //
    // An alternative way of computing P^TAP would be a gather operation with inverted loop order:
    // Cij = sum_kl Pki*Plj*Akl. This requires P^T for efficient determination of required kl indices.
    // While the transpose construction is efficient, the gather implementation ist not (tested 2016-01-17),
    // presumably because A is larger than PAP and the scattered accesses have a worse locality.
    // Sequential performance was more than 10-fold slower than the scatter implementation below, so
    // we stick to the scatter.
    //
    // A second alternative is to create a triplet matrix first. This appears to be a factor 3 slower in
    // sequential implementation, and incurs a high memory footprint as several entries are duplicate.

    auto getEntryValues = [] (auto const& row, std::vector<EntryP>& vi)
    {
      vi.clear();
      for (auto i=row.begin(); i!=row.end(); ++i)       // indices i for which Pki != 0
        vi.push_back(*i);
    };

    auto getColumnIndices = [] (auto const& row, auto& ci)     // computes (global,local) column index pairs of row k of P
    {
      ci.clear();
      int idx = 0;
      for (auto i=row.begin(); i!=row.end(); ++i, ++idx)       // indices i for which Pki != 0
        ci.push_back(std::make_pair(i.index(),idx));
    };

    timer.start("conjugation scatter");
    parallelFor([&](size_t block, size_t nBlocks)
    {
      size_t rowStart = uniformWeightRangeStart(block,nBlocks,A.N());
      size_t rowEnd   = uniformWeightRangeStart(block+1,nBlocks,A.N());
      std::vector<EntryP> pki, plj;
      using SortedIndices = std::vector<std::pair<IndexP,int>>;
      SortedIndices is, js;
      LocalMatrices<EntryA,false,SortedIndices,SortedIndices> localMatrices(pap);
      for (IndexA k=rowStart; k<rowEnd; ++k)
      {
        auto Pk = P[k];
        getColumnIndices(Pk,is);                          // indices i for which Pki != 0
        getEntryValues(Pk,pki);

        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          IndexA const l = ca.index();                    // column index of Akl
          auto Pl = P[l];
          getColumnIndices(Pl,js);                        // indices j for which Plj != 0
          getEntryValues(Pl,plj);

          // For each entry Akl of A we have to create one local matrix.
          localMatrices.push_back(is,js);

          auto Akl = *ca;
          for (int i=0; i<is.size(); ++i)                 // just scatter Akl to all affected PAP entries
          {
            auto Pki_Akl = transpose(pki[i]) * Akl;
            for (int j=0; j<plj.size(); ++j)
              localMatrices.back()(i,j) = Pki_Akl * plj[j];
          }

          if (onlyLowerTriangle && k>l)
          {
            // treat Alk entry
            abort();
          }
        }
      }
    },8*NumaThreadPool::instance().cpus());

    timer.stop("conjugation scatter");

    return pap;
  }

  /**
   * \brief Computes the triple sparse matrix product \f$ C = C + P^T A P \f$.
   *
   * \param C has to have a sparsity pattern that is a superset of the sparsity pattern of \f$ P^T A P \f$.
   *
   */
  template <class Scalar, class Entry>
  void conjugation(Dune::BCRSMatrix<Entry>& C,
		   Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& P,
		   Dune::BCRSMatrix<Entry> const& A,
		   bool onlyLowerTriangle = false )
  {
    assert(A.N()==A.M());
    assert(A.N()==P.N());
    assert(C.N()==P.M());
    assert(C.M()==P.M());

    typedef Dune::BCRSMatrix<Entry> MatA;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > MatP;

    if(onlyLowerTriangle == false )
    {
      // Step through all entries of A
      for (int k=0; k<A.N(); ++k)
        for (typename MatA::ConstColIterator ca=A[k].begin(); ca!=A[k].end(); ++ca)
        {
          int const l = ca.index();
          // Step through all entries of rows k and l of P and add entry
          for (typename MatP::ConstColIterator cpk=P[k].begin(); cpk!=P[k].end(); ++cpk)
            for (typename MatP::ConstColIterator cpl=P[l].begin(); cpl!=P[l].end(); ++cpl)
              C[cpk.index()][cpl.index()].axpy((*cpl) * (*cpk),(*ca));
        }
    }
    else
    {
      // Step through all entries of A
      for (int k=0; k<A.N(); ++k)
        for (typename MatA::ConstColIterator ca=A[k].begin(); ca!=A[k].end(); ++ca)
        {
          int const l = ca.index();
          for (typename MatP::ConstColIterator cpk=P[k].begin(); cpk!=P[k].end(); ++cpk)
            for (typename MatP::ConstColIterator cpl=P[l].begin(); cpl!=P[l].end(); ++cpl)
            {
              if( cpk.index() >= cpl.index() )
                C[cpk.index()][cpl.index()].axpy((*cpl) * (*cpk), (*ca) );
            }
          if(k>l)
          {
            for (typename MatP::ConstColIterator cpl=P[l].begin(); cpl!=P[l].end(); ++cpl)
              for (typename MatP::ConstColIterator cpk=P[k].begin(); cpk!=P[k].end(); ++cpk)
                if( cpl.index() >= cpk.index() )
                  C[cpl.index()][cpk.index()].axpy((*cpl) * (*cpk) , (*ca) );
          }
        }
    }
  }
}
#endif


