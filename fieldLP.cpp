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

void fieldLP(int profType, pPar profPar, double dr, double lambda, int V, int U, double pc, double rho[maxLrho], double nXX[maxLrho], double Ex[maxLrho], int &lRho) {
    int k1, k2, k3, kk0, noMb, a, noM;
    double aux1, aux2, S[100], absS[100], x0[2], modesPc[100], e0, u0, c0, wArb;
    double pi = 3.141592653589793, k0;
    double nCo, nCl;
    double n,w0, eVector[maxLrho], xGe[maxLrho], xF[maxLrho], nX[maxLrho][maxLwArb], rVector[maxLrho], deltan[maxLrho];
    int kx,nIndex,rIndex,nopR;
	double uw;
	double DmrmX[4],invDmrmX[4],Dmrm[4],invDmrm[4],Dm1rm[4],invDm1rm[4],cmp1[2],tt[2],r1,r2,r3;
    
    // Physical constans
    e0      = 8.854187817e-12;
    u0      = 1.25663706e-6;
    c0      = 1/sqrt(e0*u0);
    
    //
    w0 = profPar.w0;
    wArb = 2*pi*c0/lambda;
    k0   = (2*pi)/lambda;
    
    // Generate Profile
    r1 = profPar.w1;
    r2 = profPar.w1 + profPar.w2;
    r3 = profPar.w1 + profPar.w2 + profPar.w3;

    gen_rho(profType,profPar,dr,rho,lRho);
    gen_deltan_prof(w0,profType,profPar,rho,deltan,lRho);

    refIndProfile(w0, &wArb, 1, rho, deltan, lRho, nX, xGe, xF);
    for(k1 = 0;k1 < lRho;k1++){
        rVector[k1] = rho[k1]*profPar.radius;
        eVector[k1] = nX[k1][0]*nX[k1][0];
        nXX[k1] = nX[k1][0];
    }
       
    nopR = lRho;
    n    = pc;
    
    // m
    nIndex = 0;
    rIndex = 1;
    uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
    DmatrixLP(eVector[nIndex],rVector[rIndex],n,V,k0,uw,Dmrm,invDmrm);
    
    //Field r = 0
    DmatrixLP(eVector[nIndex],rVector[rIndex-1],n,V,k0,uw,DmrmX,invDmrmX);
    Ex[0] = DmrmX[3-1];
    
    // m + 1
    nIndex = nIndex+1;
    rIndex = rIndex;
    uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
    DmatrixLP(eVector[nIndex],rVector[rIndex],n,V,k0,uw,Dm1rm,invDm1rm);
    
    // C
    cmp1[1-1]   = invDm1rm[1-1] * Dmrm[1-1] + invDm1rm[2-1] * Dmrm[3-1];
    cmp1[2-1]   = invDm1rm[3-1] * Dmrm[1-1] + invDm1rm[4-1] * Dmrm[3-1];
    
    // Field r = 0 + dr
    Ex[1] = Dm1rm[3-1] * cmp1[1-1] + Dm1rm[4-1] * cmp1[2-1];

    for(k1 = 2;k1 < nopR;k1++)
    {   
        // m
        nIndex = nIndex;
        rIndex = rIndex+1;
        uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
        DmatrixLP(eVector[nIndex],rVector[rIndex],n,V,k0,uw,Dmrm,invDmrm);
        
        // m + 1
        nIndex = nIndex+1;
        rIndex = rIndex;
        uw = sqrt(abs(k0*k0*eVector[nIndex] - k0*k0*n*n));
        DmatrixLP(eVector[rIndex],rVector[rIndex],n,V,k0,uw,Dm1rm,invDm1rm);
        
        // C
        tt[0] = cmp1[0]; tt[1] = cmp1[1];
        cmp1[0] = invDm1rm[0] * ( Dmrm[0]*tt[0] + Dmrm[1]*tt[1] ) + invDm1rm[1] * ( Dmrm[2]*tt[0] + Dmrm[3]*tt[1] );
        cmp1[1] = invDm1rm[2] * ( Dmrm[0]*tt[0] + Dmrm[1]*tt[1] ) + invDm1rm[3] * ( Dmrm[2]*tt[0] + Dmrm[3]*tt[1] );

        // Field
        Ex[rIndex] = Dm1rm[3-1] * cmp1[1-1] + Dm1rm[4-1] * cmp1[2-1];
    }
}