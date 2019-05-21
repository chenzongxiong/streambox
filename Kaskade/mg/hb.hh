/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
* \file
* \author Martin Weiser
* This file contains an implementation of an additive hierarchical basis preconditioner.
*/

#ifndef HB_HH
#define HB_HH

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <boost/utility/result_of.hpp>

#include "dune/grid/common/grid.hh"
#include "dune/istl/preconditioners.hh"

#include "fem/barycentric.hh"


namespace Kaskade {
  
  /**
  * \ingroup linalgsolution
  * \brief A hierarchical basis preconditioner.
  *
  * This preconditioner realizes an approximative inverse of a Laplace operator with constant diffusion coefficient
  * discretized by linear finite elements in nodal basis over the leaf view of the provided grid.
  *
  * This implementation is geometry- and coefficient-agnostic, i.e. it assumes that the diagonal of the stiffness matrix is
  * just the identity. The advantage of this approach are its extreme simplicity, both in terms of interface and implementation,
  * a minimal memory footprint, and very cheap set-up phase and application. The drawback is, of course, a worse resulting
  * condition number in case non-uniform coarse grids, anisotropic coarse grids, or spatially varying or anisotropic
  * diffusion coefficients are to be treated.
  *
  * Extending the implementation to take more structure (stiffness matrices, non-uniform grids) into account is probably
  * not worth the while, as the iteration count grows anyway in proportion to the number of refinement levels in 2D
  * or even faster in 3D. If the simplicity of the PrecondType::HB preconditioner does'nt pay off, have a look at \ref MultigridSolver
  * or \ref BPXYPreconditioner.
  *
  * We refer to Deuflhard/Weiser Chapter 7.3.
  *
  * \tparam Grid the grid on which the Laplace operator is discretized. This implementation assumes that the grid 
  *              is purely simplicial and that on refinement new nodes are created only at edge midpoints.
  * \tparam Domain the type of domain space vectors
  * \tparam Range the type of range space vectors
  */



  template <class Grid, class Domain, class Range>
  class HierarchicalBasisPreconditioner: public Dune::Preconditioner<Domain,Range> {
    private:
      typedef typename Grid::LevelGridView::IndexSet::IndexType Index;
      typedef std::vector<Index>                                NodeSet;
	    
      static int const dim = Grid::dimension;
      
    public:
      enum { category = Dune::SolverCategory::sequential };
	    
      /**
      * \brief Constructor.
      * \param grid the grid 
      * \param mat_ a matrix containing (at least) the diagonal of the stiffness matrix
      */
      explicit HierarchicalBasisPreconditioner(Grid const& grid) 
      {
	// make sure the grid satisfies our assumptions: a purely simplicial grid
	assert(grid.leafIndexSet().geomTypes(0).size()==1 && grid.leafIndexSet().geomTypes(0)[0].isSimplex()
	      && "HierarchicalBasisPreconditioner currently requires purely simplicial grids");

	      
	// extract and store relevant sizes
	j = grid.maxLevel();
	n = grid.leafIndexSet().size(dim);
	
	
	typedef typename Grid::template Codim<0>::LevelIterator  CellLevelIterator;
	typedef typename Grid::template Codim<dim>::LevelIterator LeafVertexIterator;
	
	// extract the coarse grid nodes
	nodesOnLevel.resize(j+1); 
	for (LeafVertexIterator i=grid.template lbegin<dim>(0); i!=grid.template lend<dim>(0); ++i)
	  nodesOnLevel[0].push_back(grid.leafIndexSet().index(*i));
	
	// extract the coarse grid elements' vertices
	coarseElement.reserve(grid.levelIndexSet(0).size(0));
	for (CellLevelIterator ci=grid.template lbegin<0>(0); ci!=grid.template lend<0>(0); ++ci)
	{
	  Dune::FieldVector<Index,dim+1> idx;
	  for (int i=0; i<=dim; ++i)
	    idx[i] = grid.leafIndexSet().index(*(ci->template subEntity<dim>(i)));
	  coarseElement.push_back(idx);
	}
	
	// compute the parents of each node on level >0
	parents.resize(n,std::make_pair(n,n)); // invalid sentinel (n,n) denotes not yet processed nodes.
	for (int k=1; k<=j; ++k) 
	{
	  // On each level, look out where nodes of this level are located in their father cell. Their barycentric coordinates
	  // should have exactly two nonzero entries of size 0.5. The corresponding corners of the father cell are the parents.
	  for (CellLevelIterator ci=grid.template lbegin<0>(k); ci!=grid.template lend<0>(k); ++ci)
	  {
	    typedef typename Grid::template Codim<0>::EntityPointer CellPointer;
	    CellPointer father = ci->father();
	    
	    // Extract indices of father's vertices
	    // @TODO Variable is never read! Remove or keep for future extension of the preconditioner??
//            Index fatherCornerIndex[dim+1];
//            for (int i=0; i<=dim; ++i)
//              fatherCornerIndex[i] = grid.leafIndexSet().index(*(father->template subEntity<dim>(i)));
	    
	    // consider each vertex of current cell in turn
	    for (int i=0; i<=dim; ++i)
	    {
	      typename Grid::template Codim<dim>::EntityPointer vertexPointer = ci->template subEntity<dim>(i);
	      Index idx = grid.leafIndexSet().index(*vertexPointer);
	      
	      
	      // Edge midpoints can be reached from multiple cells depending on the spatial dimension. The parent nodes 
	      // are independent of from which cell we look. Thus we check wether we've already done the work.
	      if (parents[idx].first==n && vertexPointer->level()>0) // not yet processed (but has parents)
	      {
		// Obtain barycentric coordinates of child corner in father. For child nodes located on an edge,
		// there will be exactly two nonzero entries associated to the vertices. Note that the barycentric
		// coordinates as implemented in fem/barycentric.hh are shifted by index 1 compared to the 
		// reference element vertex numbering.
		Dune::FieldVector<typename Grid::ctype,dim+1> relativeCoord = barycentric(ci->geometryInFather().corner(i));
		
		// Extract the two coordinates with value 0.5
		int coords[dim+1];
		int pos = 0;
		for (int m=0; m<=dim; ++m)
		  if (relativeCoord[m] > 0.4)
		    coords[pos++] = (m+1)%(dim+1);
		  
		if (pos==2)
		{
		  // that's an edge midpoint: obtain the leaf indices of those father corners
		  parents[idx].first  = grid.leafIndexSet().index(*(father->template subEntity<dim>(coords[0])));
		  parents[idx].second = grid.leafIndexSet().index(*(father->template subEntity<dim>(coords[1])));
		  
		  // write down that this node appeared first on level k
		  nodesOnLevel[k].push_back(idx);
		}
	      }
	    }
	  }
	}
	
	// perform sanity checks
	#ifndef NDEBUG
	// check that each level >0 node has parents
	int count = nodesOnLevel[0].size();
	for (int k=1; k<=j; ++k) 
	{
	  for (typename NodeSet::const_iterator i=nodesOnLevel[k].begin(); i!=nodesOnLevel[k].end(); ++i)
	    assert(parents[*i].first<n && parents[*i].second<n);
	  count += nodesOnLevel[k].size();
	}
	assert(count==n);
	#endif
      }
      
      /** 
      * \brief Has to be called before the first call to apply
      */
      virtual void pre (Domain& x, Range& b) 
      {
	// actually, does nothing :)
      }

      /**
      * \brief applies the preconditioner to the residual \arg d, which results in the correction \arg v
      */
      virtual void apply (Domain& v, const Range& d) 
      {
	using namespace boost::fusion;
      // We modify the residual (in-place), hence we have to copy it here.
	Range r = d;
	
	// apply smoother and recursively restrict the residual downwards through the mesh hierarchy
	for (int k=j; k>0; --k) 
	{
	  // correct scaling for diagonal coefficients
	  typename Domain::field_type a = std::pow(2.0,(Grid::dimension-2.0)*(k-j));
	  
	  //  for all nodes on level k
	  for (typename NodeSet::const_iterator i=nodesOnLevel[k].begin(); i!=nodesOnLevel[k].end(); ++i)
	  {
	    // apply one Jacobi step, care for correct scaling of lower level hierarchical basis functions
	    // Beispiel at_c<0>(estSol.data)[*j];
	    
	    at_c<0>(v.data)[*i] = at_c<0>(r.data)[*i];
	    at_c<0>(v.data)[*i] /= a; 
	    
	    // restrict the residual
	    at_c<0>(r.data)[parents[*i].first]  += 0.5 * at_c<0>(r.data)[*i];
	    at_c<0>(r.data)[parents[*i].second] += 0.5 * at_c<0>(r.data)[*i];
	  }
	}
	
	// "solve" on coarse grid. 
	coarseGridSolution(v,r);
	
	// prolongate the correction recursively upwards through the mesh hierarchy, adding up
	// all the corrections from different levels
	for (int k=1; k<=j; ++k) 
	  // the nodes on level k get contributions of 0.5 from each of their parent nodes - that's all
	  for (typename NodeSet::const_iterator i=nodesOnLevel[k].begin(); i!=nodesOnLevel[k].end(); ++i)
	    at_c<0>(v.data)[*i] += 0.5*(at_c<0>(v.data)[parents[*i].first]+at_c<0>(v.data)[parents[*i].second]);
      }

      /**
      * \brief Has to be called after the last call to apply
      */
      virtual void post (Domain& x)
      {
	// actually, does nothing :)
      }

    private:
      std::vector<NodeSet>                       nodesOnLevel;   // hierarchical basis node sets
      std::vector<std::pair<Index,Index> >       parents;        // parent nodes of child nodes on level >0
      int                                        j;              // maximum grid level
      Index                                      n;              // total number of nodes
      std::vector<Dune::FieldVector<Index,dim+1> > coarseElement;
      
      // This coarse grid solution employs a simple Jacobi iteration, where the 
      // matrix is simply patched together from identical elemental stiffness matrices.
      // TODO: maybe it would be better to explicitly form the coarse grid matrix and use a direct solver
      void coarseGridSolution(Domain& v, Range& r) const {
	using namespace boost::fusion;
	// get the correct scaling for the coarse grid
	typename Domain::field_type a = std::pow(2.0,-(Grid::dimension-2.0));
	
	// Initialize correction to zero
	for (typename NodeSet::const_iterator i=nodesOnLevel[0].begin(); i!=nodesOnLevel[0].end(); ++i)
	  at_c<0>(v.data)[*i] = 0;
	Domain dv = v; // TODO: that's a looooong vector (overkill)
	
	// perform a couple of Jacobi iterations
	for (int k=0; k<20; ++k)
	{
	  // add correction: solve with diagonal and set residual to zero
	  for (typename NodeSet::const_iterator i=nodesOnLevel[0].begin(); i!=nodesOnLevel[0].end(); ++i)
	  {
	    at_c<0>(dv.data)[*i] = at_c<0>(r.data)[*i]; 
	    at_c<0>(v.data)[*i] += at_c<0>(dv.data)[*i] / a; 
	  }
	  
	  // update the residual. We are geometry agnostic here and use a rough approximation of the Laplacian: On each
	  // coarse grid element we assume we have the following element matrices in 1D, 2D, 3D, respectively:
	  // [ 1 -1 ]       [ 2 -1 -1 ]        [  3 -1 -1 -1 ]
	  // [ -1 1 ], 1/12 [-1  2 -1 ], 1/60  [ -1  3 -1 -1 ]
	  //                [-1 -1  2 ]        [ -1 -1  3 -1 ]
	  //                                   [ -1 -1 -1  3 ]
	  // Obviously, the matrixes are just  s*(-e*e^T + (d+1)*I), where e is the vector with all entries zero. In order
	  // to have an invertible matrix, we add a little bit more of the identity, i.e. we end up with element matrices
	  // s*(-e*e^T + (d+1.1)*I) with s=1, 1/12, 1/60 for d=1,2,3, respectively. In this form, multiplication with the
	  // total stiffness matrix is just summing up the contributions of the elemental matrices that can efficiently
	  // been computed.
	  typename Dune::Preconditioner<Domain,Range>::field_type s = dim==1? 1.0: dim==2? 1/12.0 : 1/60.0;
	  for (int i=0; i<coarseElement.size(); ++i) 
	  {
	    //typename Domain::block_type eTdv = 0;
	    double eTdv = 0;
	    for (int m=0; m<=dim; ++m)
	      eTdv += at_c<0>(dv.data)[coarseElement[i][m]];
	    for (int m=0; m<=dim; ++m)
	      at_c<0>(r.data)[coarseElement[i][m]] -= a*s*(-eTdv+(dim+1.1)*at_c<0>(dv.data)[coarseElement[i][m]]);
	  }
	}
      }
  };
  
}

#endif
