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

void dispersionLP(int profType,pPar profPar,double dr,double lambda,double dl,double *ng,double *vg,double *dmd,int &tnom,modes *modesLP1)
{
    const int nopN = 40;
	int a,k1,k2,tnom1,tnom2,tnom3;
	double e0,u0,c0,dw,wArb[3],lArb[3],w0,der1neff0_w[20],delay[20],delayLP01=0;
	double pi = 3.141592653589793, nVector[nopN+1];
    modes modesLP2[20],modesLP3[20];//modesLP1[100],

	// Physical constans
	e0      = 8.854187817e-12;
	u0      = 1.25663706e-6;
	c0      = 1/sqrt(e0*u0);

	// Frequency Range
	w0       = 2*pi*c0/lambda;
	dw       = 2*pi*c0/(lambda*lambda)*dl;
	wArb[0]  = w0-dw;//wArb[0]  = 2*pi*c0/1.625e-6;
	wArb[1]  = w0;//wArb[1]  = 2*pi*c0/1.625e-6 + 2*dw;
	wArb[2]  = w0+dw;//wArb[2]  = 2*pi*c0/1.530e-6;
	lArb[0]  = 2*pi*c0/(wArb[0]);
	lArb[1]  = 2*pi*c0/(wArb[1]);
	lArb[2]  = 2*pi*c0/(wArb[2]);

	modeSolverLP(profType,profPar,dr,lArb[0],modesLP1,tnom1);
	modeSolverLP(profType,profPar,dr,lArb[1],modesLP2,tnom2);
	modeSolverLP(profType,profPar,dr,lArb[2],modesLP3,tnom3);

	if(tnom1 == tnom2 && tnom2 == tnom3)
	{
		for(k2 = 0;k2 < tnom1;k2++)
		{
			der1neff0_w[k2] = (modesLP3[k2].pc - modesLP1[k2].pc)/(2*dw);
			ng[k2]          =  modesLP2[k2].pc + w0*der1neff0_w[k2];
            vg[k2]          = c0 / ng[k2];
			delay[k2]		= 1/vg[k2];

            if(modesLP2[k2].vi == 0 && modesLP2[k2].ui == 1)
				delayLP01 = delay[k2];
		}
        
		for(k2 = 0;k2 < tnom1;k2++){
			dmd[k2] = delay[k2]-delayLP01;// s/m
        }
		tnom = tnom1;
	}else{
		for(k2 = 0;k2 <= tnom1-1;k2++)
		{
			dmd[k2] = 1e6;// s/m
		}
		tnom = -1;
	}
}