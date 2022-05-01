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

void refIndProfile(double w0, double *wArb, int lw, double *rho, double *deltan, int lRho, double n[maxLrho][maxLwArb], double *xGe, double *xF) {
    double e0, u0, c0, Bi0[3], li0[3], dBiGe[3], dliGe[3], dBiF[3], dliF[3], nCl;
    double x0[2], Bi[maxLrho][3], li[maxLrho][3], wi[maxLrho][3], er[maxLrho][maxLwArb];
    double r1, r2, r3, delta, nCo;
    int k0, k1, k2, a;
    double pi = 3.141592653589793, aux;
    
    // Physical constans
    e0      = 8.854187817e-12;
    u0      = 1.25663706e-6;
    c0      = 1/sqrt(e0*u0);
    
    // The Sellmeir parameters for Pure Silica
    Bi0[0] = 0.697668; Bi0[1] = 0.407339; Bi0[2] = 0.889883;//     Bi0   = [ 0.697668  0.407339  0.889883];
    li0[0] = 0.070861; li0[1] = 0.113600; li0[2] = 9.784231;//     li0   = [ 0.070861  0.113600  9.784231];
    
    // The Sellmeir parameters for Germanium dopants
    dBiGe[0] = 0.031510; dBiGe[1] = 0.267302; dBiGe[2] = -0.012950;//     dBiGe   = [ 0.031510  0.267302  -0.012950];
    dliGe[0] = 0.001677; dliGe[1] = 0.032138; dliGe[2] = 0.318034;//     dliGe   = [ 0.001677  0.032138   0.318034];
    
    // The Sellmeir parameters for Germanium dopants
    dBiF[0] = -3.234366; dBiF[1] = 0.164911; dBiF[2] = 1.369490;//     dBiF   = [ -3.234366  0.164911  1.369490];
    dliF[0] = -1.108703; dliF[1] = 0.752919; dliF[2] = 2.906858;//     dliF   = [ -1.108703  0.752919  2.906858];
    
    // Dopants Concentration Profile Determination
    sellmeierEquation(&w0, 1, 0, 0, &nCl);
    
    for(k1=0;k1<lRho;k1++)//:length(fN.rho(fN.rho < r3/r1))
    {        
        delta = deltan[k1];
        nCo = nCl / sqrt(1-2*delta);
        
        if(delta == 0){
            xF[k1]  = 0;
            xGe[k1] = 0;
        }else if((float)delta > 0){
            x0[0] = 0;x0[1] = 0.2;
            xGe[k1] = zeroFinderGE(w0, nCo, x0, 100);
            xF[k1]  = 0;
        }else if((float)delta < 0){
            x0[0] = 0;x0[1] = 0.2;
            xF[k1]   = zeroFinderF(w0, nCo, x0, 100);
            xGe[k1]  = 0;
        }
    }

    // Sellmeier Equation Parameters
    for(k1=0;k1<lRho;k1++) {
        for(k2=0;k2<3;k2++) {
            Bi[k1][k2] = Bi0[k2] + xGe[k1]*dBiGe[k2] + xF[k1]*dBiF[k2];
            li[k1][k2] = li0[k2] + xGe[k1]*dliGe[k2] + xF[k1]*dliF[k2];
            wi[k1][k2] = 2*pi*c0/(li[k1][k2]*1e-6);
        }
    }

    // Equation evaluation
    for(k0=0;k0<lRho;k0++) {
        for(k1=0;k1<lw;k1++)//length(wArb)
        {
            er[k0][k1] = 1;
            for(k2=0;k2<3;k2++)// = 1:3
            {
                er[k0][k1] = er[k0][k1] + Bi[k0][k2]*wi[k0][k2]*wi[k0][k2]/(wi[k0][k2]*wi[k0][k2] - wArb[k1]*wArb[k1]);
            }
            n[k0][k1] = sqrt(er[k0][k1]);
        }
    }
}
