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

#ifndef GRIDCOMBINATORICS_HH
#define GRIDCOMBINATORICS_HH

#include <array>
#include <cassert>
#include <numeric>
#include <utility>

#include "dune/geometry/type.hh"

#include "fem/firstless.hh"

namespace Kaskade
{
  /**
   * \ingroup fem
   * \brief For each codimension (i.e., type of subentity) computes the number
   * of global ansatz functions as well as an accumulated index into an
   * array of all ansatz functions.
   *
   * \param dim the dimension of the grid
   * \param map an associative container with key type Dune::GeometryType and int values.
   * \param is  the index set defining the entities
   * \param localDof a callable object returning the number of local dofs
   *        associated to an entity of given geometry type
   */
  template <class Map, class IndexSet, class LocalDofs>
  size_t computeDofStartIndices(int dim, Map& map, IndexSet const& is, LocalDofs& localDof)
  {
    map.clear();
    size_t n = 0;

    for (int codim=0; codim<=dim; ++codim) {
      for (size_t i=0; i<is.geomTypes(codim).size(); i++) {
        Dune::GeometryType gt = is.geomTypes(codim)[i];
        map[gt] = n;

        n += is.size(gt) * localDof(gt);
      }
    }

    return n;
  }

  //---------------------------------------------------------------------

  /// \cond internals
  namespace GridCombinatoricsDetail {

    /**
     * \brief Computes the vertex indices of the vertices incident to the k-th
     * subentity of codimension c in the reference simplex of dimension d.
     * 
     * \tparam d the spatial dimension
     * \tparam Idx an output iterator type
     * 
     * \param c the codimension
     * \param k the subentity number
     * \param[out] idx an output iterator (storage for number of vertices of the subentity is required)
     */
    template <int d, class Idx>
    void vIds(int c, int k, Idx idx)
    {
      //   Dune::GeometryType gt(Dune::GeometryType::simplex);
      //   Dune::ReferenceSimplex<double,d> const& s = Dune::GenericReferenceElements<double,d>::simplices(gt);
      Dune::ReferenceElement<double,d> const &s = Dune::ReferenceElements<double,d>::simplex();
      assert(0<=k && k<s.size(c));
      assert(s.size(k,c,d)==d-c+1);
      for (int i=0; i<s.size(k,c,d); ++i, ++idx)
        *idx = s.subEntity(k,c,i,d);
    }

    // Computes the local number of the vertex that lies in coordinate
    // direction coord. The (barycentric) coordinate direction dim+1
    // points towards the origin. See the Dune simplex numbering scheme.
    inline int localVertexNumber(int dim, int coord)
    {
      return (coord+1) % (dim+1);
    }
    
    /**
     * \ingroup grid
     * \brief Computes the barycentric coordinate direction in which the vertex with given number lies.
     * This is the inverse of \ref localVertexNumber.
     */
    inline int barycentricDirection(int dim, int vertexNumber)
    {
      return (vertexNumber+dim) % (dim+1);
    }

  } // End of namespace GridCombinatoricsDetail
  /// \endcond

  /**
   * \brief Computes the d-c+1 vertex indices of the vertices incident to the
   * k-th subentity of codimension c in the reference simplex of
   * dimension d.
   */
  inline void vertexids(int d, int c, int k, int* idx)
  {
    switch(d) {
    case 0:
      idx[0] = 0;
      break;
    case 1:
      GridCombinatoricsDetail::vIds<1>(c,k,idx);
      break;
    case 2:
      GridCombinatoricsDetail::vIds<2>(c,k,idx);
      break;
    case 3:
      GridCombinatoricsDetail::vIds<3>(c,k,idx);
      break;
    default:
      assert("Not implemented for d>3!"==0);
      abort();
    }
  }
  
  /**
   * \ingroup grid
   * \brief A class for computing permutations of local vertex numbers of simplex subentities to a globally unique ordering.
   * 
   * Use this to compute globally unique numberings of degrees of freedom, as necessary for global continuity of finite element
   * functions.
   */
  template <int dimension>
  class GlobalBarycentricPermutation
  {
  public:
    template <class IndexSet, class Cell>
    GlobalBarycentricPermutation(IndexSet const& is, Cell const& cell)
    {
      for (int i=0; i<=dimension; ++i)              // obtain global indices of the 
        gvid[i] = is.subIndex(cell,i,dimension);    // cell's vertices
    }
    
    /**
     * \brief Computes a permutation of barycentric coordinates to globally unique ordering.
     * \tparam codim the codimension of the considered subentity
     * \param e the subentity number in the cell (among those of given codimension)
     * \returns the index permutation vector \f$ \pi \f$
     * 
     * If \f$ i \f$ is a barycentric coordinate direction according to the globally unique ordering,
     * \f$ \pi(i) \f$ is its cell-local barycentric coordinate direction in the reference element.
     */
    template <int codim>
    std::array<int,dimension+1-codim> barycentricSubsetPermutation(int e) const
    {
      using namespace GridCombinatoricsDetail;
      
      std::array<int,dimension+1-codim> vs;
      vIds<dimension>(codim,e,begin(vs));                                 // get the cell-local vertex ids of our subentity's vertices
      
      std::array<std::pair<size_t,int>,dimension+1-codim> vids;           // the global and subentity-local vertex indices, sorted by global
      for (int i=0; i<dimension+1-codim; ++i)                             // for each vertex of the subentity...
        vids[i] = std::make_pair(gvid[vs[i]],                             // ... get the global index...
                                 barycentricDirection(dimension,vs[i]));  // ... and the barycentric coordinate direction
      std::sort(begin(vids),end(vids),FirstLess());                       // sort according to global index
      
      std::array<int,dimension+1-codim> pi;                               // array for the permutation
      std::iota(begin(pi),end(pi),0);                                     // initialized by the identity
      for (int i=0; i<dimension+1-codim; ++i)                             // enter the computed permutation
        pi[i] = vids[i].second;
      
      return pi;
    }
    
    /**
     * \brief A dynamic interface to \ref barycentricSubsetPermutation
     * \param e the subentity number in the cell (among those of given codimension)
     * \param codim the codimension of the subentity
     * \param out points to at least dimension+1-codim entries to be filled
     */
    template <class OutIter>
    void barycentricSubsetPermutation(int e, int codim, OutIter out) const
    {
      switch(codim)
      {
        case 0: { auto pi0 = barycentricSubsetPermutation<0>(e); std::copy(begin(pi0),end(pi0),out); break; }
        case 1: { auto pi1 = barycentricSubsetPermutation<1>(e); std::copy(begin(pi1),end(pi1),out); break; }
        case 2: { auto pi2 = barycentricSubsetPermutation<2>(e); std::copy(begin(pi2),end(pi2),out); break; }
// compilation error on line below in experiments/moreExamples/atp-1D
//        case 3: { auto pi3 = barycentricSubsetPermutation<3>(e); std::copy(begin(pi3),end(pi3),out); break; }
        default: assert("Not implemented for codim>3!"==0); abort();
      }
    }
    
  private:
    std::array<int,dimension+1> gvid; // the global vertex ids of the cell's vertices
  };
  
  
} /* end of namespace Kaskade */
#endif
