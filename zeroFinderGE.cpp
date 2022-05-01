// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#include <math.h>
#include <cmath>
#include "Header.h"
using namespace std;

double zeroFinderGE(double w0,double nCo,double *x0,int n)
{
	int i;
	double x_left,x_right,x_middle,f_left,f_right,f_middle,error,new_error;
    double aux1,aux2;
    const int wl=1;
    
	// abscissa interval
	x_left  = x0[1-1];
	x_right = x0[2-1];

	// ordinate interval
    sellmeierEquation(&w0,wl,x_left,0,&aux1);
    sellmeierEquation(&w0,wl,x_right,0,&aux2);
	f_left   = (nCo - aux1);
	f_right  = (nCo - aux2);

	// lets split the input interval
	error = abs(abs(f_left) - abs(f_right));
	for(i = 1-1;i<=n-1;i++)
	{
		x_middle    = (x_left+x_right)/2;
        sellmeierEquation(&w0,wl,x_middle,0,&aux1);
        f_middle   = (nCo - aux1);
    
		if(f_middle/abs(f_middle) == f_right/abs(f_right))
		{
			x_right = x_middle;
			f_right = f_middle;
		}else{
			x_left = x_middle;
			f_left = f_middle;
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
