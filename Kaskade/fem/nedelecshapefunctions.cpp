/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <numeric>  // needed for clang++

#include "dune/grid/config.h"

#include "fem/nedelecshapefunctions.hh"


namespace Kaskade 
{
  
  template <class ctype, int dimension, class T>
  HdivSimplexShapeFunctionSet<ctype,dimension,T>::HdivSimplexShapeFunctionSet(int order)
  : ShapeFunctionSet<ctype,dim,T,dim>(Dune::GeometryType(Dune::GeometryType::simplex,dim))
  {
    assert(order >= 1);
     auto const& simplex = Dune::ReferenceElements<ctype,dim>::simplex();
    
    // We define the shape functions as linear combinations of the Lagrangian shape functions and formulate
    // the Kronecker condition of the normal components as an (underdefined) interpolation problem.
    auto const& lsfs = lagrangeShapeFunctionSet<ctype,dim,T>(Dune::GeometryType(Dune::GeometryType::simplex,dim),order);
    
    int const n = dim * SimplexLagrangeDetail::size(dim,order);               // total number of vectorial shape functions
    int const nFaces = simplex.size(1);                                       // number of faces 
    int const nFaceNodes = SimplexLagrangeDetail::size(dim-1,order);          // number of Lagrangian nodes on each face
    assert(n == dim*lsfs.size());
    
    // We partition our vectorial shape functions into first face functions and second cell-interior functions.
    // Each we represent as a column vector of coefficients for linear combination of Lagrangian shape functions.
    // Let these column vectors form the matrix C = [ Cf Ci ] of first face and second interior shape functions.
    // As the shape functions form a basis, all columns of C are linearly independent, i.e. C is invertible.
    // Evaluating the normal components of the shape functions at the Lagrangian nodes of the faces (as matrix K) 
    // gives the Kronecker condition K C = [ I 0 ]. Note that K has full rank. Choosing Cf = K^+ (pseudoinverse) and
    // Ci as a basis of the nullspace of K satisfies the Kronecker condition. This can be obtained by computing the 
    // SVD K = U [S 0] [V1 V2]^T and setting Cf = V1 S^{-1} U^T and Ci = V2.
    
    // This approach is quite implicit, and we cannot easily deduce how the individual shape functions look like 
    // (except for their normal component at the faces), but it's valid and easy to implement for arbitrary order.
    
    // Create the matrix K of size nFaces*nFaceNodes x n.
    DynamicMatrix<Dune::FieldMatrix<T,1,1>> K(nFaces*nFaceNodes,n);
    for (int f=0; f<nFaces; ++f)
    {
      auto normal = simplex.integrationOuterNormal(f); // obtain the unit outer
      normal /= normal.two_norm();                // normal of the corresponding face
      
      for (int j=0; j<nFaceNodes; ++j)
      {
        auto idx = SimplexLagrangeDetail::tupleIndex<dim-1>(order,j);                                   // interpolation grid position on face
        auto x = simplex.template geometry<1>(f).global(SimplexLagrangeDetail::nodalPosition<dim-1>(idx,order));    
        for (int k=0; k<lsfs.size(); ++k)
        {
          double phi = lsfs[k].evaluateFunction(x);
          for (int d=0; d<dim; ++d)
            K[f*nFaceNodes+j][k*dim+d] = normal[d] * phi;
        }
      }
    }
    
    // Compute SVD.
    DynamicMatrix<Dune::FieldMatrix<double,1,1>> U, V;
    std::vector<double> sigma;
    std::tie(U,sigma,V) = svd(K);
    assert(sigma.back()>0); // check that K has full rank.
    
    // Form C = [Cf Ci].
    DynamicMatrix<Dune::FieldMatrix<double,1,1>> C(n,n);
    
    // Cf = V_1 * sigma^{-1} * U^T
    for (int j=0; j<nFaces*nFaceNodes; ++j)
      for (int i=0; i<V.rows(); ++i)
      {
        C[i][j] = 0;
        for (int k=0; k<nFaces*nFaceNodes; ++k)
          C[i][j] += V[i][k] / sigma[k] * U[j][k]; // <- that's U^T[k][j]
      }
      
    // Ci = V2, do partial copy.
    for (int j=nFaces*nFaceNodes; j<n; ++j)
      for (int i=0; i<V.rows(); ++i)
        C[i][j] = V[i][j];
      
    // Now create the vectorial shape functions corresponding to the columns of C.
    sf.reserve(n);
    std::vector<Dune::FieldVector<double,dim>> coeff(lsfs.size());
    for (int j=0; j<n; ++j)
    {
      for (int i=0; i<coeff.size(); ++i) // extract j-th column of C
        for (int d=0; d<dim; ++d)
          coeff[i][d] = C[i*dim+d][j];
        
      if (j<nFaces*nFaceNodes) // face functions
        sf.push_back(value_type(coeff,std::make_tuple(order,1,j/nFaceNodes,j%nFaceNodes)));
      else                     // interior cell functions
        sf.push_back(value_type(coeff,std::make_tuple(order,0,0,j-nFaces*nFaceNodes)));
    }
    
    this->order_ = order;
    this->size_ = n;
  }
  
  // explicit instantiaion of method
  template HdivSimplexShapeFunctionSet<double,2,double>::HdivSimplexShapeFunctionSet(int order);
  template HdivSimplexShapeFunctionSet<double,3,double>::HdivSimplexShapeFunctionSet(int order);
}


#ifdef UNITTEST

#include <iostream>


using namespace Kaskade;

// output to standard out in gnuplot syntax
// inspect the shape functions with
// plot "test.gnu" index 1 using 1:2:3:4 with vectors


int main(void) 
{
  int const dim = 2;
  int const n = 30;
  for (int order=1; order<3; ++order)
  {
    auto const& hdivsfs = hdivShapeFunctionSet<double,dim,double>(Dune::GeometryType(Dune::GeometryType::simplex,dim),order);
    std::cerr << "number of shape functions: " << hdivsfs.size() << "\n"; 
    for (int ix=0; ix<=n; ++ix)
    {
      double x = (double)ix/n;
      for (int iy=0; iy<=n; ++iy)
      {
        double y = (double)iy/n;
        std::cout << x << ' ' << y << ' ';
        for (int k=0; k<hdivsfs.size(); ++k)
        {
          Dune::FieldVector<double,dim> xi; xi[0] = x; xi[1] = y;
          auto phi = hdivsfs[k].evaluateFunction(xi) / 35.0;
          if (x+y>1) // outside the unit simplex...
            phi = 0;
          std::cout << phi[0] << ' ' << phi[1] << ' ';
        }
        std::cout << "\n";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

#endif
