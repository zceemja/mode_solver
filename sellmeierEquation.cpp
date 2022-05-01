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

void sellmeierEquation(double *w,int wl,double xGe,double xF,double *n){
    double e0,u0,c0,Bi0[3],li0[3],dBiGe[3],dliGe[3],dBiF[3],dliF[3],Bi[3],li[3],wi[3],er[maxLwArb];
    int k1,k2;
    double pi = 3.141592653589793;
    
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

    // Sellmeier Equation Parameters
    for(k1=0;k1<3;k1++)
    {
        Bi[k1] = Bi0[k1] + xGe*dBiGe[k1] + xF*dBiF[k1];
        li[k1] = li0[k1] + xGe*dliGe[k1] + xF*dliF[k1];
    }
    for(k1=0;k1<3;k1++){
        wi[k1] = 2*pi*c0/(li[k1]*1e-6);
    }

    // Equation evaluation
    for(k1=0;k1<wl;k1++)
    {
        er[k1] = 1;
        for(k2=0;k2<3;k2++){
            er[k1] = er[k1] + Bi[k2]*wi[k2]*wi[k2]/(wi[k2]*wi[k2] - w[k1]*w[k1]);
        }
        // Output
        n[k1] = sqrt(er[k1]);
    }
}

    
