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

#include "dune/grid/config.h"

#include "fem/hierarchicshapefunctions.hh"

namespace Kaskade
{
  namespace HierarchicDetail
  {

    /**
     * Returns the number of subentites of given codimension.
     */
    template <int d>
    int nSubent(int c)
    {
      Dune::ReferenceElement<double,d> const& s = Dune::ReferenceElements<double,d>::simplex();
      return s.size(c);
    }

    int nSubentities(int d, int c)
    {
      assert(0<=d && d<=3);

      switch(d) {
      case 0:
        return 1;
      case 1:
        return nSubent<1>(c);
      case 2:
        return nSubent<2>(c);
      case 3:
        return nSubent<3>(c);
      default:
        abort();
      }
      return 0;
    }

    //---------------------------------------------------------------------

    /**
     * Given a sequence vGlobal of d+1 different integers, substitute
     * 0,..,d such that i -> vGlobal[vNormal[i]] is monotonely increasing,
     * i.e. 0 is in vNormal located where the minimal value in vGlobal is,
     * etc.
     *
     * This is essentially a monotone, i.e. ordering coserving,
     * renumbering or compression to the index range 0,..,d.
     */
    void normalizeIndices(int d, int const vGlobal[], int vNormal[])
    {
      std::copy(vGlobal,vGlobal+d+1,vNormal);
      bool done[d+1];
      for (int i=0; i<=d; ++i) done[i] = false;

      for (int i=0; i<=d; ++i) {
        int min = std::numeric_limits<int>::max();
        int idx = -1;
        for (int j=0; j<=d; ++j)
          if (!done[j] && vNormal[j]<min) {
            idx = j;
            min = vNormal[j];
          }
        vNormal[idx] = i;
        done[idx] = true;
      }
    }

    //---------------------------------------------------------------------

    int cdimSize(int d, int p, int c)
    {
      if (p<0) return 0;
      if (p==0)
      {
        if (c==d) return d+1;
        else      return 0;
      }
      if (p==1) return 0;

      // On edges (which are essentially 1D, there is just one polynomial
      // of given degree (>=2). Thus the number of edges in the simplex is
      // the number of edge functions.
      if (c==d-1) return nSubentities(d,d-1);

      // Recursive application. On each subentity of given codimension,
      // the number of shape functions is determined by the recursive
      // structure.
      return nSubentities(d,c)*size(d-c,p-(d-c+1));
    }

    //---------------------------------------------------------------------

    int size(int d, int p)
    {
      // Just sum over all subentity codimensions.
      int k=0;
      for (int c=0; c<=d; ++c)
        k += cdimSize(d,p,c);
      return k;
    }

    //---------------------------------------------------------------------

    /* Octave/Matlab program for computing coefficients for three-term recurrence.
function [as,bs,cs] = op(K)

% first order polynomial
x = poly([0]);

pm1 = [];
p = [1];

as = [];
bs = [];
cs = 1/sprod(p,p);

spm1 = 0;

for k=1:K
  sp = sprod(p,p)
  spx = sprod(p,conv(p,x));

  pm2 = pm1;
  pm1 = p;
  spm2 = spm1;
  spm1 = sp;

  a = -spx/spm1;
  b = -spm1/spm2;

  p = polyadd(conv(pm1,poly(-a)),b*pm2);
  sp = sprod(p,p);

  as = [as a];
  bs = [bs b];
  cs = [cs 1/sp];
end

% 
function y = sprod(p1,p2)

% weight function defining the scalar product
w = poly([-1,1]);

p1 = polyderiv(conv(p1,w));
p2 = polyderiv(conv(p2,w));
p = conv(p1,p2);
n = length(p);

% integration of polynomials on [-1,1] by scalar product with
int = 2*mod(n:-1:1,2)./(n:-1:1);
% Remember that polynomial coefficients in octave are reverse!

y = int*p';

function p = polyadd(p1,p2)
n1 = length(p1);
n2 = length(p2);

if n1>n2
  p2 = [zeros(1,n1-n2) p2];
end

if n2>n1
  p1 = [zeros(1,n2-n1) p1];
end

p = p1+p2;

     */

    // And the resulting  coefficients.
    int const nCoefficients = 7;
    double kernel_a[nCoefficients] = { 0, 0, 0, 0, 0, 0, 0 };
    double kernel_b[nCoefficients] = { 0, -0.6,  -0.723809523809523,  -1.357894736842107,  -2.180478008000302,
        -2.939943419665647,  -3.986343856928839 };

    double kernel_c[nCoefficients+1] = { 0.375,
        0.62499999999999978,   0.86348684210526350,   0.63590116279069686,   0.29163383462595738,
        0.09919709089473709,   0.02488422836939123, 0.00499952847651311};

    /**
     * One-dimensional hierarchic polynomials on the interval [-1,1].
     *
     * These are orthogonal polynomials computed by three term
     * recurrece. They are orthogonal w.r.t. the scalar product \f[ (p,q)
     * = \int_{-1}^1 ((1+x)(1-x)p(x))' ((1+x)(1-x)q(x))' \, dx. \f] This
     * improves the condition number of the local stiffness matrix when
     * used as kernel function in hierarchic shape functions.
     */
    double kernel(int p, double x)
    {
      assert(p<=nCoefficients);

      if (p<0) return 0;
      if (p==0) return 1;

      double fa = 0, fb = 1;
      for (int pp=0; pp<p; ++pp) {
        double t = (x+kernel_a[pp])*fb + kernel_b[pp]*fa;
        fa = fb;
        fb = t;
      }
      return fb*kernel_c[p];
    }

    double dKernel(int p, double x)
    {
      assert(p<=nCoefficients);

      if (p<=0) return 0;

      double fa = 0, fb = 1, dfa = 0, dfb = 0;
      for (int pp=0; pp<p; ++pp) {
        double t = (x+kernel_a[pp])*fb + kernel_b[pp]*fa;
        double dt = fb+(x+kernel_a[pp])*dfb + kernel_b[pp]*dfa;
        fa = fb; dfa = dfb;
        fb = t; dfb = dt;
      }
      return dfb*kernel_c[p];
    }

    //---------------------------------------------------------------------

    std::pair<int,int> codim(int d, int p, int k)
            {
      assert(0<=k);

      // codim==dim is equivalent to p==0 here.
      if (p==0) {
        assert(k<=d);
        return std::make_pair(d,k);
      }

      assert(p!=1);

      // Check for higher dimensional subentities.
      for (int c=d-1; c>=0; --c) {
        if (k<cdimSize(d,p,c))
          return std::make_pair(c,k);
        k -= cdimSize(d,p,c);
      }

      assert("k too large!"==0);
      return std::make_pair(-1,-1);
            }

    //---------------------------------------------------------------------

    std::tuple<int,int,int> tupleIndex(int d, int p, int k)
            {
      int c=codim(d,p,k).first, m=codim(d,p,k).second;
      //   int c, m;
      //   std::tie(c,m) = codim(d,p,k);

      if (c>=d-1) // at most one shape function on subentity
        return std::make_tuple(c,m,0);
      else {
        // recursive call
        int n = size(d-c,p-(d-c+1));
        return std::make_tuple(c,m/n,m%n);
      }
            }

    //---------------------------------------------------------------------

    /**
     * Computes sequential index from tuple index (see \ref tupleIndex).
     */
    int sequentialIndex(int d, int p, int c, int e, int  k)
    {
      assert(d>=0 && p>=0);
      assert(0<=c && c<=d);
      assert(e<nSubentities(d,c));
      assert(k<cdimSize(d,p,c)/nSubentities(d,c));

      int idx = 0;

      for (int co=d; co>c; --co)
        idx += cdimSize(d,p,co);
      assert(idx<size(d,p));

      int m = cdimSize(d,p,c) / nSubentities(d,c);
      idx += e*m+k;
      assert(idx<size(d,p));
      assert(std::make_tuple(c,e,k)==tupleIndex(d,p,idx));

      return idx;
    }

    //---------------------------------------------------------------------

    double shapeFunction(int d, int p, int k, double const b[])
    {
      assert(p>=0);

      if (p==0) {
        assert(0<=k && k<=d);
        return b[(k+d)%(d+1)];
      }

      assert(p!=1);

      // Get codimension and number in codimension.
      std::pair<int,int> cdim = codim(d,p,k);
      int c = cdim.first;
      k = cdim.second;
      assert(c<d); // Check that for p>0 no vertex functions are present.

      if (c==d-1) { // edge functions are treated separately
        // Get vertex ids of the incident vertices.
        int idx[2];
        vertexids(d,c,k,idx);
        // Get associated barycentric coordinates
        idx[0] = (idx[0]+d)%(d+1);
        idx[1] = (idx[1]+d)%(d+1);

        // compute and return value
        return 4*b[idx[0]]*b[idx[1]]*kernel(p-2,b[idx[0]]-b[idx[1]]);
      } else { // Higher dimensional subentities
        // Get number of shape functions on subentity, subentity id and shape function inside subentity.
        int n = size(d-c,p-(d-c+1));
        int e = k/n;
        int ksf = k%n;
        // Get vertex ids of the incident vertices.
        int idx[d-c+1];
        vertexids(d,c,e,idx);
        // Get associated barycentric coordinates and compute value.
        double lambda[d-c+1];
        for (int i=0; i<=d-c; ++i)
          lambda[(i+d-c)%(d-c+1)] = b[(idx[i]+d)%(d+1)];

        double phi = pow(d-c+1,d-c+1); // compensate for small value of barycentric product such that max phi = 1
        for (int i=0; i<=d-c; ++i)
          phi *= lambda[i];

        double val = phi*shapeFunction(d-c,p-(d-c+1),ksf,lambda);

        return val;
      }
    }

    //---------------------------------------------------------------------

    void dShapeFunction(int d, int p, int k, double const b[], double* grad)
    {
      //  int D = d, P = p, K = k;

      assert(p>=0);

      if (p==0) {
        assert(0<=k && k<=d);
        std::fill(grad,grad+d+1,0);
        grad[(k+d)%(d+1)] = 1;
        return;
      }

      assert(p!=1);


      // Get codimension and number in codimension.
      std::pair<int,int> cdim = codim(d,p,k);
      int c = cdim.first;
      k = cdim.second;

      if (c==d-1) { // edge functions are treated separately
        // Get vertex ids of the incident vertices.
        int idx[2];
        vertexids(d,c,k,idx);
        // Get associated barycentric coordinates
        idx[0] = (idx[0]+d)%(d+1);
        idx[1] = (idx[1]+d)%(d+1);
        // compute and return derivative value
        for (int i=0; i<=d; ++i) {
          double db0 = idx[0]==i? 1: 0;
          double db1 = idx[1]==i? 1: 0;

          grad[i] = (db0*b[idx[1]]+b[idx[0]]*db1) * kernel(p-2,b[idx[0]]-b[idx[1]]);
          grad[i] += b[idx[0]]*b[idx[1]] * dKernel(p-2,b[idx[0]]-b[idx[1]])*(db0-db1);
          grad[i] *= 4;
        }
      } else { // Higher dimensional subentities
        // Get number of shape functions on subentity, subentity id and shape function inside subentity.
        int n = size(d-c,p-(d-c+1));
        int e = k/n;
        int ksf = k%n;
        // Get vertex ids of the incident vertices.
        int idx[d+1];
        vertexids(d,c,e,idx);
        // Get associated barycentric coordinates
        double lambda[d-c+1];
        for (int i=0; i<=d-c; ++i)
          lambda[(i+d-c)%(d-c+1)] = b[(idx[i]+d)%(d+1)];

        // Compute product of subentity barycentric coordinates (lambda).
        double phi = pow(d-c+1,d-c+1); // compensate for small value of barycentric product such that max phi = 1

        double dphi_dlambda[d-c+1];
        std::fill(dphi_dlambda,dphi_dlambda+d-c+1,0.0);

        for (int i=0; i<=d-c; ++i) {
          for (int j=0; j<=d-c; ++j)
            dphi_dlambda[j] *= lambda[i];
          dphi_dlambda[i] += phi;

          phi *= lambda[i];
        }

        double sf = shapeFunction(d-c,p-(d-c+1),ksf,lambda);

        // Compute recursive shapefunction derivative w.r.t. lambda.
        double dSF_dlambda[d-c+1];
        dShapeFunction(d-c,p-(d-c+1),ksf,lambda,dSF_dlambda);

        // Compute derivative of phi*sf w.r.t.  lambda.
        double gr[d-c+1];
        for (int i=0; i<=d-c; ++i)
          gr[i] = dphi_dlambda[i]*sf + phi*dSF_dlambda[i];

        // Compute derivative of phi*sf w.r.t. complete barycentric coordinates.
        for (int i=0; i<=d; ++i)
          grad[i] = 0;
        for (int i=0; i<=d-c; ++i)
          grad[(idx[i]+d)%(d+1)] = gr[(i+d-c)%(d-c+1)];

        // // Implementation-error free by numerical differentiation
        // for (int i=0; i<=D; ++i) {
        //   double h = 1e-5;
        //   double bl[d+1], br[d+1];
        //   std::copy(b,b+D+1,bl);
        //   std::copy(b,b+D+1,br);
        //   bl[i] -= h;
        //   br[i] += h;
        //   double diff = (shapeFunction(D,P,K,br)-shapeFunction(D,P,K,bl))/(2*h);
        //   if (std::abs(grad[i]-diff)>1e-3)
        //     std::cout << "i=" << i << " alg=" << grad[i] << " num=" << diff << "\n";
        //   grad[i] = diff;
        // }
        return;
      }
    }

    //---------------------------------------------------------------------

    int actualOrder(int d, int p, int k)
    {
      assert(p>=0 && d>=0 && k>=0);

      if (p==0) return 1;
      assert(p!=1);

      std::pair<int,int> cdim = codim(d,p,k);
      int c = cdim.first;
      k = cdim.second;

      if (c==d-1) return p; // edge function always have nominal order
      else {                // recursive call
        int n = size(d-c,p-(d-c+1));
        int ksf = k%n;
        return (d-c+1)+actualOrder(d-c,p-(d-c+1),ksf);
      }
    }

    //---------------------------------------------------------------------

    int entityPermutation(int d, int c, int e, int const vGlobal[])
    {
      // normalize global vertex numbers to range 0,..,d
      int vGl[d+1];
      normalizeIndices(d,vGlobal,vGl);
      int dsub = d-c;


      // extract normalized global vertex numbers of given subentity
      int vIdx[dsub+1];
      vertexids(d,c,e,vIdx);
      for (int i=0; i<dsub+1; ++i)
        vIdx[i] = vGl[vIdx[i]];

      // Find a codim c subentity with local vertex numbers as given in
      // vIdx. To simplify the lookup, sort our vertex indices first.
      std::sort(vIdx,vIdx+dsub+1);

      for (int i=0; i<nSubentities(d,c); ++i) {
        int idx[dsub+1];
        vertexids(d,c,i,idx);
        std::sort(idx,idx+dsub+1);
        if (std::equal(vIdx,vIdx+dsub+1,idx)) {
          return i;
        }
      }

      assert("Never get here!"==0);
      return -1;
    }

    //---------------------------------------------------------------------

    std::tuple<int,int> sfPermutation(int d, int p, int c, int e, int k, int const vGlobal[])
    {
      if (c==d || p<2) // vertex functions match here: there is only one on each vertex, hence no
        return std::make_tuple(0,1); // permutation is necessary


      if (c==d-1) {
        // edge functions, only one per order on edge, hence no index permutation
        // But sign changes for odd order if vertices are permuted.
        int idx[2];
        vertexids(d,c,e,idx);
        bool permuted = vGlobal[idx[1]] < vGlobal[idx[0]];
        return std::make_tuple(0,(permuted && p%2==1)? -1: 1);
      } else {
        int dsub = d-c;         // dimension of recursive call shape function
        int psub = p-(d-c+1);   // order of recursive call shape function

        // Get global ids of subentity vertices
        int idx[dsub+1];
        vertexids(d,c,e,idx);
        for (int i=0; i<=dsub; ++i)
          idx[i] = vGlobal[idx[i]];

        // We are the shapefunction with local number k on the dsub-dimensional simplex
        // with order psub. Now let's look up the codimension (csub), subentity number (esub),
        // and index in subentity (ksub) with respect to this dsub-dimensional simplex.
        int csub, esub, ksub;
        std::tie(csub,esub,ksub) = tupleIndex(dsub,psub,k);

        // The global indices of the dsub-dimensional simplex are contained in idx. Get the
        // globally unique number based on these global vertex indices.
        int sign;
        std::tie(ksub,sign) = sfPermutation(dsub,psub,csub,esub,ksub,idx);

        // Now the subentity number esub will be used to compute the total index k in our
        // dsub-dimensional subsimplex, however, it is still local and hence non-unique. We
        // correct this by selecting a globally unique (but otherwise arbitrary) subentity index.
        esub = entityPermutation(dsub,csub,esub,idx);

        // Now both esub and ksub are globally unique. Compute the total index in the
        // dsub-dimensional subsimplex.
        k = sequentialIndex(dsub,psub,csub,esub,ksub);
        assert(k<size(dsub,psub));

        return std::make_tuple(k,sign);
      }
    }

    //---------------------------------------------------------------------

  } // End of namespace HierarchicDetail
} // End of namespace Kaskade

