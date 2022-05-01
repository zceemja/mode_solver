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

void gen_deltan_prof(double w0,int profType,pPar profPar,double *rho,double *deltan,int lRho)
{
    if(profType == 1){
        double delta1,delta2,delta3,delta1E,alpha1,alpha2,alpha3,n_rho[maxLrho],r1,r2,r3,nCl,nCo,nCo1,nCo2,nClD,*aux1;
        int k1;
        
        delta1 = profPar.sI1;
        delta1E = profPar.eI1;
        delta2 = profPar.sI2;
        delta3 = profPar.sI3;

        alpha1 = profPar.a1;
        alpha2 = profPar.a2;
        alpha3 = profPar.a3;

        r1 = profPar.w1;
        r2 = profPar.w1 + profPar.w2;
        r3 = profPar.w1 + profPar.w2 + profPar.w3;

        sellmeierEquation(&w0,1,0,0,&nCl);
        nCo1        = nCl / sqrt(1-2*delta1);
        nCo2        = nCl / sqrt(1-2*delta2);
        nClD        = nCl / sqrt(1-2*delta3);

        n_rho[0] = nCo1;
        deltan[0] = (n_rho[0]*n_rho[0]-nCl*nCl)/(2*n_rho[0]*n_rho[0]);
        for(k1=1;k1<lRho;k1++)
        {
            if(rho[k1] <= r1/r1){
                n_rho[k1] = nCo1*(1 - (delta1-delta1E)*(pow((rho[k1])/(r1/r1),alpha1)));
            }else if(rho[k1] > r1/r1 && rho[k1] <= r2/r1)
                n_rho[k1] = nCo2*(1 - delta2*(pow(rho[k1]/(r2/r1),alpha2)));
            else if(rho[k1] > r2/r1 && rho[k1] <= r3/r1)
                n_rho[k1] = nClD;
            else if(rho[k1] > r3/r1)
                n_rho[k1] = nCl;
            deltan[k1] = (n_rho[k1]*n_rho[k1]-nCl*nCl)/(2*n_rho[k1]*n_rho[k1]);//(n_rho[k1]*n_rho[k1]-nCl*nCl)/(2*n_rho[k1]*n_rho[k1]);
        }
    }
}