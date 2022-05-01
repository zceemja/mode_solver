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

void gen_rho(int profType, pPar profPar, double incMax,double *rho,int &lRho)
{
    if(profType == 1){
        double r1, r2, r3, *aux1;
        int k1;       
        r1 = profPar.w1;
        r2 = profPar.w1 + profPar.w2;
        r3 = profPar.w1 + profPar.w2 + profPar.w3;
        
        rho[0] = 0;
        k1 = 1;
        while(1){
            rho[k1] = rho[k1-1] + incMax;

            if(rho[k1] > r3/r1 + 2){
                lRho = k1 + 1;
                break;
            }
            k1   = k1 + 1;
            lRho = k1;
        }
    }
}