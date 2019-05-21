/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <math.h>
#include <float.h>

typedef double coord[3];


int main(int argc, char *argv[])
  {
    const char *fName = "data.out", *outName = "v.dat";
    if (argc>1) fName = argv[1];
    if (argc>2) outName = argv[2];
    FILE *f = fopen(fName, "r");
    if (f==0)
      {
        printf("Cann't open  %s\n", fName);
        return -1;
      }

    int i, k, n = -1, rc;
    rc = fscanf(f, "%d", &n);
    if ((rc!=1)&&(n>0)&&(n<100000))
      {
        printf("Data file should start with a valid integer\n");
        return -1;
      }

    double *freq = new double[n];
    coord *u1 = new coord[n], *v1 = new coord[n];
    coord *u2 = new coord[n], *v2 = new coord[n];

    for (k=0; k<n; k++)
      {
        rc = fscanf(f, "%le", &freq[k]);
        if (rc!=1) break;
        for (i=0; i<3; i++)
          {
            rc = fscanf(f, "%le", &u1[k][i]);
            if (rc!=1) break;
          }
        if (rc!=1) break;
        for (i=0; i<3; i++)
          {
            rc = fscanf(f, "%le", &v1[k][i]);
            if (rc!=1) break;
          }
        if (rc!=1) break;
        for (i=0; i<3; i++)
          {
            rc = fscanf(f, "%le", &u2[k][i]);
            if (rc!=1) break;
          }
        if (rc!=1) break;
        for (i=0; i<3; i++)
          {
            rc = fscanf(f, "%le", &v2[k][i]);
            if (rc!=1) break;
          }
        if (rc!=1) break;
      }
    if (k!=n)
      {
        printf("Not enough data\n");
        n = k;
      }
    printf("%d lines read from %s\n",  n, fName);
    fclose(f);

    fopen(outName, "w");
    double t, dmax, r, r1, r2;
    const double PI2 = 2*3.1415926535;
    for (k=0; k<n; k++)
      {
        fprintf(f,"%10.02f", freq[k]);
        dmax = 0.0;
        for (i=0; i<1000; i++)
          {
            t = i*PI2/n;
            r1 = u1[k][0]*u1[k][0] + u1[k][1]*u1[k][1] + u1[k][2]*u1[k][2];
            r2 = v1[k][0]*v1[k][0] + v1[k][1]*v1[k][1] + v1[k][2]*v1[k][2];
            r = cos(t)*sqrt(r1)+sin(t)*sqrt(r2);
            if (r>dmax) dmax = r;
          }
        fprintf(f," %16.8e", dmax);
        dmax = 0.0;
        for (i=0; i<1000; i++)
          {
            t = i*PI2/n;
            r1 = u2[k][0]*u2[k][0] + u2[k][1]*u2[k][1] + u2[k][2]*u2[k][2];
            r2 = v2[k][0]*v2[k][0] + v2[k][1]*v2[k][1] + v2[k][2]*v2[k][2];
            r = cos(t)*sqrt(r1)+sin(t)*sqrt(r2);
            if (r>dmax) dmax = r;
          }
        fprintf(f," %16.8e\n", dmax);
      }
    printf("%d lines written to %s\n",  n, outName);
    fclose(f);
  }
