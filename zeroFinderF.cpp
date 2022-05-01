// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#include <math.h>
#include <cmath>
#include "Header.h"
using namespace std;

double zeroFinderF(double w0,double nCl,double *x0,int n)
{
	int i;
	double x_left,x_right,x_middle,f_left,f_right,f_middle,error,new_error;
    double aux1,aux2;
    
	// abscissa interval
	x_left  = x0[1-1];
	x_right = x0[2-1];

	// ordinate interval
    sellmeierEquation(&w0,1,0,x_left,&aux1);
    sellmeierEquation(&w0,1,0,x_right,&aux2);
	f_left   = (nCl - aux1);
	f_right  = (nCl - aux2);

	// lets split the input interval
	error = abs(abs(f_left) - abs(f_right));
	for(i = 1-1;i<=n-1;i++)
	{
		x_middle    = (x_left+x_right)/2;
        sellmeierEquation(&w0,1,0,x_middle,&aux1);
        f_middle   = (nCl - aux1);
    
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
