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

double zeroFinderMode(double *eVector, double *rVector, int v, double k0, int nopR, double *S,double *x0,int n)
{
	int i;
	double x_left,x_right,x_middle,f_left,f_right,f_middle,error,new_error;
	// abscissa interval
	x_left  = x0[1-1];
	x_right = x0[2-1];

	// ordinate interval
	f_left   = SmatrixLP(eVector,rVector,&x_left,v,k0,1,nopR,S);
	f_right  = SmatrixLP(eVector,rVector,&x_right,v,k0,1,nopR,S);

	// lets split the input interval
	error = abs(abs(f_left) - abs(f_right));
	for(i = 1-1;i<=n-1;i++)
	{
		x_middle    = (x_left+x_right)/2;
		f_middle    = SmatrixLP(eVector,rVector,&x_middle,v,k0,1,nopR,S);//f( x_middle );
    
		if(f_middle/abs(f_middle) == f_right/abs(f_right))
		{
			x_right = x_middle;
			f_right = f_middle;//f(x_right);
		}else{
			x_left = x_middle;
			f_left = f_middle;//f(x_left);
		} 
		new_error = abs(abs(f_left) - abs(f_right));
		if(abs(new_error)/abs((f_left + f_right))*100 < 0.1 || abs(new_error-error) == 0)
		{
			return (x_left + x_right)/2;
		}else{
			error = new_error;
		}
	}

	return (x_left + x_right)/2;
}
