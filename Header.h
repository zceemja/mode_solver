// AUTHOR:  Filipe M. Ferreira (filipef@ieee.org)
#define maxLwArb 10
#define maxLrho 20000

typedef struct{
    int vi;
    int ui;
    double pc;
}modes;

typedef struct{
    // graded core, double trench
    double w1;
    double w2;
    double w3;
    double sI1;
    double sI2;
    double sI3;
    double eI1;
    double a1;
    double a2;
    double a3;
    double w0;
    double radius;
}pPar;

void DmatrixLP(double epsi,double r,double n,int v,double k0,double u, double* D, double* invD);
double SmatrixLP(double *eVector, double *rVector, double *nVector, int v, double k0, int nopN, int nopR, double *S);
double zeroFinderF(double w0,double nCl,double *x0,int n);
double zeroFinderGE(double w0,double nCo,double *x0,int n);
void gen_deltan_prof(double w0,int profType,pPar profPar,double *rho,double *deltan,int lRho);
void refIndProfile(double w0,double *wArb,int lw,double *rho,double *deltan,int lRho,double n[maxLrho][maxLwArb],double *xGe,double *xF);
void sellmeierEquation(double *w,int wl,double xGe,double xF,double *n);
void gen_rho(int profType, pPar profPar, double incMax,double *rho,int &lRho);
void modeSolverLP(int profType,pPar profPar,double dr,double lambda,modes *modesLP,int &tnom);
double zeroFinderMode(double *eVector, double *rVector, int v, double k0, int nopR, double *S,double *x0,int n);
void FO_dmd(int profType,pPar profPar,double dr,double *lArb,int lArbL,double dl,double &FO,int &tnom,modes *modesLP,double dmd[20][10]);
void dispersionLP(int profType,pPar profPar,double dr,double lambda,double dl,double *ng,double *vg,double *dmd,int &tnom,modes *modesLP1);
void fieldLP(int profType, pPar profPar, double dr, double lambda, int V, int U, double pc, double rho[maxLrho], double nRho[maxLrho], double Ex[maxLrho], int &lRho);

