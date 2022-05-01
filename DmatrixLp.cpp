// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#include "mex.h"
#include "matrix.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include "Header.h"
#include <gsl/gsl_sf_bessel.h>
using namespace std;

void DmatrixLP(double epsi,double r,double n,int v,double k0,double u, double* D, double* invD)
{
	double a1, a2, a3, a4;
	int a;

	if(n >= sqrt(epsi))
	{
		a1 = +u/2*(gsl_sf_bessel_In(v-1,u*r) + gsl_sf_bessel_In(v+1,u*r));
		a2 = -u/2*(gsl_sf_bessel_Kn(v-1,u*r) + gsl_sf_bessel_Kn(v+1,u*r));
		a3 = gsl_sf_bessel_In(v,u*r);
		a4 = gsl_sf_bessel_Kn(v,u*r);
	}else{
		a1 = +u/2*(jn(v-1,u*r) - jn(v+1,u*r));
		a2 = +u/2*(yn(v-1,u*r) - yn(v+1,u*r));
		a3 = jn(v,u*r);
		a4 = yn(v,u*r);
	}

	D[1-1] = a1;
	D[2-1] = a2;
	D[3-1] = a3;
	D[4-1] = a4;
   
 	invD[1-1] = a4/(a1*a4 - a2*a3);
	invD[2-1] = -a2/(a1*a4 - a2*a3);
	invD[3-1] = -a3/(a1*a4 - a2*a3);
	invD[4-1] = a1/(a1*a4 - a2*a3);
}
