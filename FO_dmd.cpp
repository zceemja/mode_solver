// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#include "matrix.h"
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

void FO_dmd(int profType,pPar profPar,double dr,double *lArb,int lArbL,double dl,double &FO,int &tnom,modes *modesLP,double dmd[20][10])
{
	int k1,k2,tnom0=0;
	double ng[20],vg[20],dmdAux[20],dmdMax[20],FO_dmd_var=0,meanDMD=0,sFO = 0;
	
	for(k1=0;k1<lArbL;k1++){
		dmdMax[k1] = 0;
		dispersionLP(profType,profPar,dr,lArb[k1],dl,ng,vg,dmdAux,tnom,modesLP);
        printf("tnom: %d, dmdAux[1]: %f\n",tnom,dmdAux[1]*1e15);
        if(k1==0)
            tnom0 = tnom;
        else if(tnom != tnom0){
            printf("ERROR - tnom varies with lambda!");
            FO = 1e6;
            tnom = -1;
        }
        
		for(k2 = 0;k2 < tnom;k2++)
		{
            dmd[k2][k1] = dmdAux[k2];
			if(abs(dmdAux[k2])>dmdMax[k1])
				dmdMax[k1] = (dmdAux[k2]);
		}
	}
	for(k1=0;k1<lArbL;k1++){
        printf("dmdMax[%d]: %f\n",k1,dmdMax[k1]*1e15);
        sFO  = sFO + (dmdMax[k1]); 
		FO_dmd_var  = FO_dmd_var + abs(dmdMax[k1]); 
	}
	FO = FO_dmd_var/lArbL*(sFO/abs(sFO));
    printf("FO: %f\n",FO*1e15);
	return;
}
