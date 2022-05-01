// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#include "mex.h"
#include "matrix.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include "Header.h"
using namespace std;

double SmatrixLP(double *eVector, double *rVector, double *nVector, int v, double k0, int nopN, int nopR, double *S)
{
	int kx, k1,a,nIndex,rIndex;
	double uw, n;
	double Dmrm[4],invDmrm[4],Dm1rm[4],invDm1rm[4],cmp1[2],tt[2];

    for(kx = 0;kx < nopN;kx++)
	{
		n = nVector[kx];

        // m
        nIndex = 0;
        rIndex = 1;
		uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
        DmatrixLP(eVector[nIndex],rVector[rIndex],n,v,k0,uw,Dmrm,invDmrm);
        
        // m + 1
        nIndex = nIndex+1;
        rIndex = rIndex;
		uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
		DmatrixLP(eVector[nIndex],rVector[rIndex],n,v,k0,uw,Dm1rm,invDm1rm);
    
        // C
		cmp1[1-1]   = invDm1rm[1-1] * Dmrm[1-1] + invDm1rm[2-1] * Dmrm[3-1];
		cmp1[2-1]   = invDm1rm[3-1] * Dmrm[1-1] + invDm1rm[4-1] * Dmrm[3-1];
		for(k1 = 2;k1 < nopR;k1++)
		{
            // m
            nIndex = nIndex;
            rIndex = rIndex+1;
			uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
			DmatrixLP(eVector[nIndex],rVector[rIndex],n,v,k0,uw,Dmrm,invDmrm);
        
            // m + 1
            nIndex = nIndex+1;
            rIndex = rIndex;
			uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
			DmatrixLP(eVector[rIndex],rVector[rIndex],n,v,k0,uw,Dm1rm,invDm1rm);
            
            // C
            tt[0] = cmp1[0]; tt[1] = cmp1[1];        
			cmp1[0] = invDm1rm[0] * ( Dmrm[0]*tt[0] + Dmrm[1]*tt[1] ) + invDm1rm[1] * ( Dmrm[2]*tt[0] + Dmrm[3]*tt[1] );
			cmp1[1] = invDm1rm[2] * ( Dmrm[0]*tt[0] + Dmrm[1]*tt[1] ) + invDm1rm[3] * ( Dmrm[2]*tt[0] + Dmrm[3]*tt[1] );
        }
		
        S[kx] = cmp1[0];
	}
	return S[0];
}