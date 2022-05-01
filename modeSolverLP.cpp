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

void modeSolverLP(int profType,pPar profPar,double dr,double lambda,modes *modesLP,int &tnom)
{
    const int nopN = 40;
	int k1,k2,k3,kk0,v,noMb,a,noM,lRho;
	double aux1, aux2,S[100],absS[100],x0[2],modesPc[100],e0,u0,c0,wArb;
	double pi = 3.141592653589793,k0;
    double V,nCo, nCl,nVector[nopN+1];
    double w0,rho[maxLrho],eVector[maxLrho],xGe[maxLrho],xF[maxLrho],deltan[maxLrho],n[maxLrho][maxLwArb],rVector[maxLrho];
    
     // Physical constans
    e0      = 8.854187817e-12;
    u0      = 1.25663706e-6;
    c0      = 1/sqrt(e0*u0);
    
    //
    w0 = profPar.w0;
    wArb = 2*pi*c0/lambda;
    
    // Generate Profile
    gen_rho(profType,profPar,dr,rho,lRho);
    gen_deltan_prof(w0,profType,profPar,rho,deltan,lRho);
    refIndProfile(w0,&wArb,1,rho,deltan,lRho,n,xGe,xF);
    for(k1 = 0;k1 < lRho;k1++){
        rVector[k1] = rho[k1]*profPar.radius;
        eVector[k1] = n[k1][0]*n[k1][0];
    }
    
    
    nCo = sqrt(eVector[0]);
    nCl = sqrt(eVector[lRho-1]);
    V = 2*pi/lambda*profPar.w1*sqrt(nCo*nCo-nCl*nCl);
    
    for(k1 = 0;k1 < nopN+1;k1++)
        nVector[k1] = nCl + (nCo-nCl)/double(nopN)*double(k1);
    nVector[0]    = nVector[0] + 1e-9;
    nVector[nopN] = nVector[nopN] - 1e-9;

	tnom = 0;
	k0       = (2*pi)/lambda;
	x0[0] = nVector[0]; x0[1] = nVector[nopN];

	for(k1 = 0;k1 <= 8;k1++)
	{
		noM = 0;
		noMb = noM;
		v = k1;
		aux1 = SmatrixLP(eVector,rVector,&x0[0],v,k0,1,lRho,S);
		aux2 = SmatrixLP(eVector,rVector,&x0[1],v,k0,1,lRho,S);

		if((aux1/abs(aux1) != aux2/abs(aux2)) && V < 8)
		{
			x0[0] = nVector[0]; x0[1] = nVector[nopN];
			modesPc[noM] = zeroFinderMode(eVector, rVector, v, k0,lRho,S,x0,100);
			noM++;tnom++;
		}else{
			SmatrixLP(eVector,rVector,nVector,v,k0,nopN+1,lRho,S);
			for(k2 = 0;k2<=nopN;k2++){
				absS[k2] = S[k2]/abs(S[k2]);
			}
			for(k3 = 1;k3 <= nopN;k3++)
			{
				if(absS[k3] != absS[k3-1])
				{
					x0[0] = nVector[k3-1];x0[1] = nVector[k3];
					modesPc[noM] = zeroFinderMode(eVector, rVector, v, k0,lRho,S,x0,100);
					noM++;tnom++;
				}
			}
		}
		for(k3=noM-1;k3>=0;k3--){
			modesLP[tnom-k3-1].vi = v;
			modesLP[tnom-k3-1].ui = noM-k3;
			modesLP[tnom-k3-1].pc = modesPc[k3];
		}
		if(noMb == noM){
			break;
		}
	}
}