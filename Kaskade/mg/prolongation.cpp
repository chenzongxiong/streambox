/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <limits>
#include <tuple>

#include "dune/grid/config.h"

#include "fem/firstless.hh"
#include "mg/prolongation.hh"

namespace Kaskade
{
  MGProlongation::MGProlongation(std::vector<ProlongationDetail::Node> const& parents, std::vector<size_t> const& indexInCoarse,
                                 size_t nc_, int fineLevel)
  : entries(parents.size()), nc(nc_)
  {
    for (size_t row=0; row<parents.size(); ++row)
      if (parents[row].level == fineLevel)                                      // this is a fine grid node appearing on level l
      {
        entries[row][0] = indexInCoarse[parents[row].p];                        // define its value
        entries[row][1] = indexInCoarse[parents[row].q];                        // as the average of its parent nodes
        if (entries[row][0]>entries[row][1])                                    // sort the two entries ascendingly
          std::swap(entries[row][0],entries[row][1]);
      }
      else                                                                      // this already existed on coarse level, just copy value
        entries[row][0] = entries[row][1] = indexInCoarse[row];                 // (we use copy as twice the half of the value for uniform access)
  }

  // ----------------------------------------------------------------------------------------------

  namespace ProlongationDetail
  {
    template <class IndexInCoarse>
    void mapToCoarseIndices(Node& node, IndexInCoarse const& indexInCoarse, size_t range)
    {
      using namespace std;
      if (node.p < range)
        node.p = indexInCoarse[node.p];
      if (node.q < range)
        node.q = indexInCoarse[node.q];
    }


    std::vector<MGProlongation> makeProlongationStack(std::vector<Node> parents, int maxLevel, size_t minNodes)
    {
      std::vector<MGProlongation> prolongations;

      // For moving between levels, create one prolongation matrix.
      for (int l=maxLevel; l>0 && parents.size()>minNodes; --l)
      {
        // Find all nodes on level l-1. Assign consecutive coarse level indices to them.
        std::vector<Node> coarseNodes;
        std::vector<size_t> indexInCoarse(parents.size(),std::numeric_limits<size_t>::max()); // fill with sentinel
        for (size_t i = 0; i<parents.size(); ++i)
        {
          if (parents[i].level < l)
          {
            indexInCoarse[i] = coarseNodes.size();
            coarseNodes.push_back(parents[i]);
          }
        }

        // Create prolongation and put it to the lowest level.
        prolongations.insert(begin(prolongations),MGProlongation(parents,indexInCoarse,coarseNodes.size(),l));


        // Now map all parent node indices to the coarse level numbering.
        for (auto& node: coarseNodes)
          mapToCoarseIndices(node,indexInCoarse,parents.size());

        // Move to next lower level.
        parents = std::move(coarseNodes);
      }

      return prolongations;
    }
  }

  std::ostream& operator<<(std::ostream& out, MGProlongation const& p)
  {
    for (int i=0; i<p.N(); ++i)
    {
      out << i << ' ' << p.parents(i)[0] << ' ' << 0.5 << "\n";
      out << i << ' ' << p.parents(i)[1] << ' ' << 0.5 << "\n";
    }
    return out;
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  // Algebraic multigrid

  namespace
  {
    template <class Entry, class Index>
    MGProlongation getAlgebraicProlongation(NumaBCRSMatrix<Entry,Index> const& A)
    {
      std::vector<ProlongationDetail::Node> parents(A.N(),ProlongationDetail::Node{0,0,-1});

      // For each node, decide whether this is a coarse grid or a fine grid node. If a
      // node has a coarse grid neighbor, it becomes a fine grid node (level 1), otherwise a coarse
      // grid node (level 0).
      for (Index i=0; i<A.N(); ++i)
      {
        auto const& row = A[i];
        parents[i].level = 0;                 // assume we'll be a coarse grid node
        for (auto ci=row.begin(); ci!=row.end() && ci.index()<i; ++ci)
          if (parents[ci.index()].level==0)   // found a neighboring coarse grid node
          {
            parents[i].level = 1;             // -> turn this into a fine grid node
            break;                            // ...decision made
          }
      }

      // For each fine grid node, find (at most) two neighboring coarse grid nodes, from which the
      // node's value can be interpolated. The number is limited to two in order to keep the resulting
      // projected Galerkin matrices very sparse.
      std::vector<std::pair<double,Index>> tmp;
      for (Index i=0; i<A.N(); ++i)
        if (parents[i].level > 0)
        {
          int count = 0;
          Index ps[2];
          auto const& row = A[i];
          // A naive strategy is implemented here: Take the (up to) two neighboring coarse
          // grid nodes with "strongest" connection to current fine grid node.
          tmp.clear();
          for (auto ci=row.begin(); ci!=row.end(); ++ci)                        // look at each neighbor
            if (parents[ci.index()].level==0)                                   // consider only coarse grid nodes
              tmp.push_back(std::make_pair(ci->frobenius_norm2(),ci.index()));  // store their connection strength and index
          std::sort(begin(tmp),end(tmp),FirstGreater());                        // sort by connection strength

          if (tmp.size()==1)                                                    // there is at least one (by construction)
            parents[i].p = parents[i].q = tmp[0].second;                        // if only one, count this one doubly
          else
          {
            parents[i].p = tmp[0].second;                                       // otherwise, extract the two "strongest"
            parents[i].q = tmp[1].second;                                       // coarse grid neighbors
          }
        }

      // For all nodes flagged as coarse, find their index in the coarse grid.
      std::vector<Index> indexInCoarse(A.N());
      Index coarseIndex = 0;
      for (Index i=0; i<A.N(); ++i)
        if (parents[i].level==0)
          indexInCoarse[i] = coarseIndex++;

      // Create the prolongation.
      return MGProlongation(parents,indexInCoarse,coarseIndex,1);
    }

  } // end of anonymous namespace


  template <class Entry, class Index>
  MultiGridStack<MGProlongation,Entry,Index> makeAlgebraicMultigridStack(NumaBCRSMatrix<Entry,Index>&& A, Index n, bool onlyLowerTriangle)
  {
    n = std::max((Index)5,n);

    Timings& timer = Timings::instance();
    timer.start("algebraic multigrid stack");

    std::vector<MGProlongation> prolongations;
    std::vector<NumaBCRSMatrix<Entry,Index>> galerkinMatrices;

    galerkinMatrices.push_back(std::move(A));

    std::cerr << "galerkinMatrices[0].N()=" << galerkinMatrices[0].N() << " vs n=" << n << "\n";
    while (galerkinMatrices[0].N() > n)
    {
      // Create a prolongation from the matrix entries
      timer.start("prolongation creation");
      prolongations.insert(begin(prolongations),getAlgebraicProlongation(galerkinMatrices[0]));
      timer.stop("prolongation creation");
      timer.start("matrix projection");
      galerkinMatrices.insert(begin(galerkinMatrices),prolongations[0].galerkinProjection(galerkinMatrices[0],onlyLowerTriangle));
      timer.stop("matrix projection");
    }
    timer.stop("algebraic multigrid stack");

    return MultiGridStack<MGProlongation,Entry,Index>(std::move(prolongations),std::move(galerkinMatrices));
  }

  // explicit instantiation for scalar problems
  template MultiGridStack<MGProlongation,Dune::FieldMatrix<double,1,1>,size_t>
           makeAlgebraicMultigridStack(NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>,size_t>&& A, size_t n, bool onlyLowerTriangle);
  template MultiGridStack<MGProlongation,Dune::FieldMatrix<double,2,2>,size_t>
           makeAlgebraicMultigridStack(NumaBCRSMatrix<Dune::FieldMatrix<double,2,2>,size_t>&& A, size_t n, bool onlyLowerTriangle);
  template MultiGridStack<MGProlongation,Dune::FieldMatrix<double,3,3>,size_t>
           makeAlgebraicMultigridStack(NumaBCRSMatrix<Dune::FieldMatrix<double,3,3>,size_t>&& A, size_t n, bool onlyLowerTriangle);
}

// ----------------------------------------------------------------------------------------------



#ifdef UNITTEST

#include <iostream>

#include "dune/grid/uggrid.hh"

#include "fem/fetransfer.hh"
#include "fem/spaces.hh"
#include "io/vtk.hh"
#include "utilities/gridGeneration.hh"


using namespace Kaskade;

int main(void)
{
  using Grid = Dune::UGGrid<2>;
  Dune::FieldVector<double,2> x0(0.0), dx(1.0);

  GridManager<Grid> gridman(createRectangle<Grid>(x0,dx,2.0));
  H1Space<Grid> h1Space(gridman,gridman.grid().leafGridView(),1);
  H1Space<Grid> h1Space2(gridman,gridman.grid().leafGridView(),2);

  gridman.globalRefine(2);
  H1Space<Grid>::Element_t<1> f(h1Space);
  H1Space<Grid>::Element_t<1> f2(h1Space2);

  // Check h multigrid prolongations
  auto pStack = prolongationStack(gridman.grid());

  Dune::BlockVector<Dune::FieldVector<double,1>> coeff(pStack[0].M());
  for (int i=0; i<coeff.N(); ++i)
    coeff[i] = (double)(i);

  for (int i=0; i<pStack.size(); ++i)
  {
    auto const& p = pStack[i];
    auto tmp = coeff;
    tmp.resize(p.N());
    tmp = 0;
    p.umv(coeff,tmp);
    swap(tmp,coeff);
  }
  // Now f should be piecewise linear with integer values
  // at the "coarse" grid nodes.
  f.coefficients() = coeff;
  writeVTK(f,"f",IoOptions(),"f");

  // Check p multigrid prolongation
  auto pMat = prolongation(h1Space,h1Space2);
  pMat.mv(f.coefficients(),f2.coefficients());
  writeVTK(f,"fquad",IoOptions().setOrder(2),"f");


  return 0;
}

#endif
