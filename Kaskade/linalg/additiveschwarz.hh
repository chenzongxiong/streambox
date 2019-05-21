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

#ifndef ADDSCHWARZ_HH
#define ADDSCHWARZ_HH

#include <iostream>
#include <vector>
#include <boost/timer/timer.hpp>
#include <dune/istl/preconditioners.hh>

#include "linalg/triplet.hh"

namespace Kaskade
{
  enum schwarzType { DIAGONAL, COMPLETE };

  extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
  extern "C" void dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv,
      double *b, int *ldb, int *info);

  //---------------------------------------------------------------------

  // TODO: document this
  template <class Op>
  class AdditiveSchwarzPreconditioner: public Dune::Preconditioner<typename Op::domain_type, typename Op::range_type>
  {
    typedef typename Op::domain_type domain_type;
    typedef typename Op::domain_type range_type;
    typedef typename domain_type::field_type field_type;

  public:
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \arg op the assembled operator
     */
    AdditiveSchwarzPreconditioner(Op const& op, schwarzType select, int verbosity=2)
    {
      MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
      size_t k, n = A.nrows();
      switch (select)
      {
      case DIAGONAL:
        iBlocks.resize(n+1);
        iBlocksColIndices.resize(n);
        for (k=0; k<n; k++)
        {
          iBlocks[k] = k;
          iBlocksColIndices[k] = k;
        }
        iBlocks[n] = n;
        break;
      case COMPLETE:
        iBlocks.resize(2);
        iBlocksColIndices.resize(n);
        iBlocks[0] = 0;
        iBlocks[1] = n;
        for (k=0; k<n; k++)
        {
          iBlocksColIndices[k] = k;
        }
        break;
      }
      init(A,verbosity);
    }

    AdditiveSchwarzPreconditioner(Op const& op, size_t first, size_t last, int verbosity=2)
    {
      MatrixAsTriplet<field_type> A = op.template get<MatrixAsTriplet<field_type> >();
      size_t i, k, nBlocks = last-first, cnt, nnz = A.nnz();
      iBlocks.resize(nBlocks+1);
      // printf("AdditiveSchwarzPreconditioner: n=%ld, nBlocks=%ld\n", n, nBlocks);
      std::vector<size_t> nCounts(nBlocks);
      for (k=0; k<nBlocks; k++)
        nCounts[k] = 1;
      cnt = nBlocks;

      for (i=0; i<nnz; i++)
      {
        size_t ri = A.ridx[i];
        size_t ci = A.cidx[i];
        if ((ri>=first)&&(ri<last))
        {
          if ((ci<first)||(ci>=last))
          {
            nCounts[ri-first]++;
            cnt++;
          }
        }
        //             else
          //              {
          //                 if ((ci>=first)||(ci<last))
        //                   {
        //                     nCounts[ci-first]++;
        //                     cnt++;
        //                   }
        //               }
      }
      // printf("nCounts:  cnt=%ld ", cnt);
      // for(k=0; k<nBlocks; k++) printf("%3ld ", nCounts[k]);
      // printf("\n");

      iBlocksColIndices.resize(cnt);
      iBlocks[0] = 0;
      for (k=1; k<=nBlocks; k++)
      {
        iBlocks[k] = iBlocks[k-1]+nCounts[k-1];
        iBlocksColIndices[iBlocks[k-1]] = first+k-1;
        nCounts[k-1] = 1;
      }
      assert(iBlocks[nBlocks]==cnt);
      // printf("iBlocks:  ");
      // for(k=0; k<=nBlocks; k++) printf("%3ld ", iBlocks[k]);
      // printf("\n");
      for (i=0; i<nnz; i++)
      {
        size_t ri = A.ridx[i];
        size_t ci = A.cidx[i];
        if ((ri>=first)&&(ri<last))
        {
          if ((ci<first)||(ci>=last))
          {
            iBlocksColIndices[iBlocks[ri-first]+nCounts[ri-first]] = ci;
            nCounts[ri-first]++;
          }
        }
        //             else
          //              {
          //                 if ((ci>=first)||(ci<last))
        //                   {
        //                     iBlocksColIndices[iBlocks[ci-first]+nCounts[ci-first]] = ri;
        //                     nCounts[ci-first]++;
        //                   }
        //               }
      }

      init(A,verbosity);
    }

    void init(MatrixAsTriplet<field_type> &A, int verbosity=2)
    {
      size_t i, k, n = A.nrows(), nnz = A.nnz();
      // computing the amount of storage needed vo all sub-matrices

      size_t nBlocks = iBlocks.size()-1;
      // printf("AdditiveSchwarzPreconditioner: n=%ld, nBlocks=%ld\n", n, nBlocks);
      // for (k=0; k<nBlocks; k++) {
      // printf("%3ld: ", k);
      // for (i=iBlocks[k]; i<iBlocks[k+1]; i++) printf("%3ld ", iBlocksColIndices[i]);
      // printf("\n");
      // }
      blocks.resize(nBlocks);
      ipiv.resize(nBlocks);
      size_t maxBlockSize = -1;
      size_t memSize = 0, inverseIndexSize = 0;
      for (k=0; k<nBlocks; k++)
      {
        int blockSize = iBlocks[k+1]-iBlocks[k];
        if (blockSize>maxBlockSize)
          maxBlockSize = blockSize;
        memSize += blockSize*blockSize;
        inverseIndexSize += blockSize;
      }
      // printf("    memSize=%ld, inverseIndexSize=%ld\n", memSize, inverseIndexSize);

      // allocating the memory

      blocks[0] = (field_type*)malloc(memSize*sizeof(field_type));
      ipiv[0] = (int*)malloc(inverseIndexSize*sizeof(int));
      for (k=1; k<nBlocks; k++)
      {
        int blockSize = iBlocks[k]-iBlocks[k-1];
        blocks[k] = blocks[k-1]+blockSize*blockSize;
        ipiv[k] = ipiv[k-1]+blockSize;
      }
      // printf("blocks allocated: \n");
      // for (k=0; k<nBlocks; k++) {
      // int blockSize = iBlocks[k+1]-iBlocks[k];
      // printf("%8p %d ", blocks[k], blockSize*blockSize);
      // }
      // printf("\n");

      // which index is used in which submatrix? (inverse)

      std::vector<size_t> usedCnts(n);
      std::vector<size_t> usedIn(n+1);
      std::vector<size_t> usedInIndices(inverseIndexSize);

      for (k=0; k<n; k++)
        usedCnts[k] = 0;
      for (k=0; k<nBlocks; k++)
      {
        for (i=iBlocks[k]; i<iBlocks[k+1]; i++)
        {
          usedCnts[iBlocksColIndices[i]]++;
        }
      }
      // printf("usedCnts: ");
      // for (k=0; k<n; k++) printf(" %3ld", usedCnts[k]);
      // printf("\n");
      size_t cnt = 0;
      for (k=0; k<n; k++)
      {
        usedIn[k] = cnt;
        cnt += usedCnts[k];
      }
      if (cnt!=inverseIndexSize) printf("Assert cnt=%ld!=%ld=inverseIndexSize\n", cnt, inverseIndexSize);
      usedIn[n] = inverseIndexSize;

      for (k=0; k<n; k++)
        usedCnts[k] = 0;
      for (k=0; k<nBlocks; k++)
      {
        for (i=iBlocks[k]; i<iBlocks[k+1]; i++)
        {
          size_t ii = iBlocksColIndices[i];
          usedInIndices[usedIn[ii]+usedCnts[ii]] = k;
          usedCnts[ii]++;
        }
      }
      // printf("Inverse index: n=%ld, nBlocks=%ld\n", n, nBlocks);
      // for (k=0; k<n; k++) {
      // printf("%3ld: ", k);
      // for (i=usedIn[k]; i<usedIn[k+1]; i++) printf("%3ld ", usedInIndices[i]);
      // printf("\n");
      // }

      // filling the submatrices

      for (i=0; i<nnz; i++)
      {
        size_t ri = A.ridx[i];
        size_t ci = A.cidx[i];
        field_type vali = A.data[i];
        for (k=usedIn[ri]; k<usedIn[ri+1]; k++)
        {
          size_t kk = usedInIndices[k];
          // printf("A(%ld,%ld)=%14.5e in Block %ld\n", ri, ci, vali, kk);
          size_t subi, subk, bsize = iBlocks[kk+1]-iBlocks[kk];
          // now we have found the blocknummer kk
          for (subi=iBlocks[kk]; subi<iBlocks[kk+1]; subi++)
            if (ri==iBlocksColIndices[subi]) break;
          if (subi==iBlocks[kk+1]) continue;
          subi -= iBlocks[kk];
          for (subk=iBlocks[kk]; subk<iBlocks[kk+1]; subk++)
            if (ci==iBlocksColIndices[subk]) break;
          if (subk==iBlocks[kk+1]) continue;
          subk -= iBlocks[kk];
          // now we have the indices for the submatrix
          // printf("subA(%ld,%ld) %8p set to%14.5e\n", subi, subk, &blocks[kk][subi*bsize+subk], vali);
          blocks[kk][subi*bsize+subk] = vali;    // or adding? preset blocks with zero!
        }
      }

      // calling lapacks to compute the inverses inverseIndexSize

      for (k=0; k<nBlocks; k++)
      {
        size_t bsize = iBlocks[k+1]-iBlocks[k];
        //             printf("block=%ld %8p, bsize=%ld, ipiv=%8p\n", k, blocks[k], bsize, ipiv[k]);
        // 			for (int j=0; j<bsize; j++)
        // 			  {
        // 				printf("%8p:", &blocks[k][j*bsize]);
        // 				for (int l=0; l<bsize; l++)
        // 				  printf("%14.5e ", blocks[k][j*bsize+l]);
        // 				printf("\n");
        // 			  }
        if (bsize==1)
        {
          if (fabs(blocks[k][0])<1.0e-20)
          {
            printf("error diagonal element to small %e\n", blocks[k][0]);
            exit(1);
          }
          blocks[k][0] = 1.0/blocks[k][0];
        }
        else
        {
          int info, f_n = bsize;
          dgetrf_(&f_n, &f_n, blocks[k], &f_n, ipiv[k], &info);
          if (info!=0)
          {
            printf("error from dgetrf_: info = %d, block=%ld %8p, bsize=%ld\n", info,
                k, blocks[k], bsize);
            exit(1);
          }
        }
      }

      if ( verbosity>=2 )
      {
        std::cout << "PrecondType::ADDITIVESCHWARZ: n=" << n << ", nnz=" << nnz << ", time=" << boost::timer::format(timer.elapsed()) << "\n";
      }
    }

    ~AdditiveSchwarzPreconditioner()
    {
      free(blocks[0]);
      free(ipiv[0]);
    }

    virtual void pre(domain_type&, range_type&) {}
    virtual void post (domain_type&) {}

    virtual void apply (domain_type& x, range_type const& y)
    {
      size_t i, k;
      int n = y.dim();
      std::vector<field_type> b(n), bb(n), bbb(n);
      y.write(b.begin());
      int nBlocks = iBlocks.size()-1;
      // printf("apply: n=%d, nBlock=%d\n", n, nBlocks);
      for (k=0; k<nBlocks; k++)
      {
        int nn, nrhs = 1, lda, ldb, info, cnt;
        char trans[1] = {'N'};
        cnt = 0;
        for (i=iBlocks[k]; i<iBlocks[k+1]; i++)
        {
          size_t ii = iBlocksColIndices[i];
          // printf("bb[%d] %14.5e = b[%ld] %14.5e\n", cnt, bb[cnt], ii, b[ii]);
          bb[cnt] = b[ii];
          cnt++;
        }
        nn = lda = ldb = iBlocks[k+1]-iBlocks[k];
        // printf("apply%3ld: cnt=%d %8p, lng=%d\n", k, cnt, ipiv[k], lda);
        if (nn==1)
        {
          bb[0] = bb[0]*blocks[k][0];
        }
        else
        {
          dgetrs_(trans, &nn, &nrhs, blocks[k], &lda, ipiv[k], &bb[0], &ldb, &info);
          if (info!=0)
          {
            printf("error from dgetrs_: info =%d\n", info);
            exit(1);
          }
        }
        cnt = 0;
        for (i=iBlocks[k]; i<iBlocks[k+1]; i++)
        {
          size_t ii = iBlocksColIndices[i];
          // printf("b[%ld] %14.5e = bb[%d] %14.5e\n", ii, b[ii], cnt, b[cnt]);
          bbb[ii] += bb[cnt];
          cnt++;
        }
      }
      x.read(bbb.begin());
    }

  private:
    std::vector<size_t> iBlocks, iBlocksColIndices;
    std::vector<field_type*> blocks;
    std::vector<int*> ipiv;
    boost::timer::cpu_timer timer;
  };
} // namespace Kaskade
#endif
