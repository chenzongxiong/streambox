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

#ifndef LINALG_PATCHDOMAINDECOMPOSITIONPPRECONDITIONER_HH
#define LINALG_PATCHDOMAINDECOMPOSITIONPPRECONDITIONER_HH

#include <vector>

#include <boost/timer/timer.hpp>

#include <dune/common/dynvector.hh>
#include "dune/istl/preconditioner.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/symmetricOperators.hh"
#include "fem/supports.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  /**
   * \ingroup iterative
   * \brief An additive overlapping domain decomposition type preconditioner for higher order finite elements applied to elliptic equations.
   * 
   * The subdomains are the stars around grid vertices. The ansatz functions with a support contained in a subdomain
   * are treated as a block, and the corresponding submatrix is inverted using dense linear algebra. The subdomain
   * corrections are added up.
   * 
   * \tparam Space a finite element space (FEFunctionSpace)
   * \tparam m the number of components in the equation
   */
  template <class Space, int m>
  class PatchDomainDecompositionPreconditioner: public SymmetricPreconditioner<typename Space::template Element<m>::type::StorageType,
                                                                               typename Space::template Element<m>::type::StorageType>
  {
    using Scalar = typename Space::Scalar;
    using Entry = Dune::FieldMatrix<Scalar,m,m>;
    using VEntry = Dune::FieldVector<Scalar,m>;
    using PEntry = float;
  public:
    using domain_type = typename Space::template Element<m>::type::StorageType;
    using range_type = domain_type;
    
    template <class Matrix>
    PatchDomainDecompositionPreconditioner(Space const& space, Matrix const& A)
    : patchMatrices(space.gridView().size(space.gridView().dimension)),
      patchDofs(patchMatrices.size())
    {
      assert(A.N()>0);
      assert(A.N()==A.M());
      Timings& timer = Timings::instance();
      timer.start("make Schwarz preconditioner");
      
      auto patches = computePatches(space.gridView());      
      assert(patches.size()==patchMatrices.size());
      
      // We can create more tasks than CPU cores available - these will just wait in the thread pool queue.
      // Having more tasks can be beneficial if the loads are unbalanced (e.g. in unstructured grids
      // or when in a multi-user environement different processes occupy some cores).
      // On a single-CPU machine, however, a larger number of tasks can be slower, probably because
      // more threads than cores compete for cache and memory bandwith.
      int const taskFactor = 1;
      
      parallelFor([&patches,this,&A,&space](size_t j, size_t s)
      {
        DynamicMatrix<Entry> aloc;
        for (size_t i=uniformWeightRangeStart(j,s,patches.size()); i<uniformWeightRangeStart(j+1,s,patches.size()); ++i)
        {
          this->patchDofs[i] = algebraicPatch(space,patches[i]);
          auto const& idx = this->patchDofs[i];
          int n = idx.size();
          auto& aloc = this->patchMatrices[i];
          aloc.resize(n*m,n*m);
//          aloc = 0;     
          for (int i=0; i<n*m;i++)                     // workaround for a bug in
            for (int j=0; j<n*m;j++) aloc[i][j]=0;     // Macintosh clang++
          for (int k=0; k<n; ++k)         // extract the submatrix of current patch
          {
            auto Ak = A[idx[k]];
            for (int l=0; l<n; ++l)       // TODO: scanning through the row, ticking off the desired entries, is 
            {                             //       probably more efficient (idx is sorted anyway)
              auto const& akl = Ak[idx[l]];
              for (int r=0; r<m; ++r)
                for (int s=0; s<m; ++s)
                  aloc[k*m+r][l*m+s] = akl[r][s];    
            }
          }
//          line below temporarely outcommented, because it depends on further changes which
//          cannot yet be checked in (L.W.)
//          invertSpd(aloc); // Jacobi style: invert the local block
        }
      },taskFactor*NumaThreadPool::instance().cpus());

      timer.stop("make Schwarz preconditioner");
    }
    
    virtual Scalar applyDp(domain_type& x, range_type const& y)
    {
      apply(x,y);
      return x*y;
    }
    
    virtual void apply(domain_type& x, range_type const& y)
    {
      // TODO: consider parallelization, e.g. maintain groups of
      //       patches where the patches in each group are algebraically 
      //       disjoint and can be processed in parallel. Of course, this
      //       requires a higher setup cost.
      std::vector<PEntry> xi, yi;
      for (size_t i=0; i<patchMatrices.size(); ++i)
      {
        auto const& idx = patchDofs[i];
        xi.resize(m*patchMatrices[i].N());
        yi.resize(m*patchMatrices[i].N());
        for (int j=0; j<idx.size(); ++j)   // extract the local residual values
          for (int k=0; k<m; ++k)
            yi[j*m+k] = y[idx[j]][k];
        patchMatrices[i].mv(yi,xi);        // solve the local system (using BLAS ssymv does not improve performance)
        for (int j=0; j<idx.size(); ++j)   // copy the solution values back, adding up
          for (int k=0; k<m; ++k)
            x[idx[j]][k] += xi[j*m+k];     // all the corrections
      }
    }
    
    virtual bool requiresInitializedInput() const
    {
      return true;
    }
    
  private:
    std::vector<DynamicMatrix<PEntry>> patchMatrices;
    std::vector<std::vector<size_t>>   patchDofs;
  };
}

#endif
