#include <stdio.h>
#include <math.h>
#include "mlf_cintf.h"

#define NG 32
#define NR 64
#define NPAR 16

double params[NG][NR][2], weights[NG][NR];

double W[NPAR];

int fBasis_invExp(const double *x, const double *rpar, double *y, int nX, int nPar, int sPar, void *ptr){
  for(int i=0; i<nPar; i++) {
    double a = rpar[sPar*i], b = rpar[sPar*i+1];
    for(int j=0; j<nX; j++) {
      double z = 2.0/(1.0+exp(-exp(b)*x[j]-a))-1;
      y[i*nX+j] = 1-z*z;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  mlf_init();
  MLF_OBJ *fobj = mlf_basisfunction(fBasis_invExp, NULL);
  for(int i=0; i<NG; i++) {
    double vGi = 4.0*((double) i)/(NG-1.0);
    double x = ((double) i) / ((double) NG-1);
    double w1 = 1.0*(1.0-x)+0.2*x;
    double w2 = 0.2*(1.0-x)+1.0*x;
    for(int j=0; j<NR; j++) {
      double vRj = -6.0 + 6.0*((double) i) / ((double) NR-1);
      params[i][j][0]=vGi;
      params[i][j][1]=vRj;
      double y = ((double) j) / ((double) NR-1);
      weights[i][j] = w1*(1.0-y)+w2*y;
    }
  }
  double alpha = 0.0;
  MLF_OBJ *fbasis = mlf_funBasisInit(fobj, 2, alpha, 0.0, 3.0, (double *) params, NR*NG, NPAR, 8192, (double *) weights);
  MLF_OBJ *f2dgrid = mlf_2dgridModelInit(fbasis, 0.0, 4.0, -6.0, 0.0, NPAR, 1024, 1024);
  //MLF_OBJ *f2dgrid = mlf_2dgridModelInit(fbasis, 0.0, 4.0, 0.0, 4.0, NPAR, 2048, 2048);
  double Y[2];
  for(Y[0] = 0.0; Y[0] < 4.0; Y[0] += 0.3)
    for(Y[1] = -6; Y[1] < 0; Y[1] += 0.4) {
      printf("Function %f %f:\n", Y[0], Y[1]);
      //mlf_getProj(fbasis, (double*) Y, (double*) W, 1, 2, NPAR);
      mlf_getProj(f2dgrid, (double*) Y, (double*) W, 1, 2, NPAR);
      for(double x = 0.1; x<3.0; x+=0.3) {
        double vx = mlf_getValue(fbasis, (double*) W, x, NPAR), vy;
        fBasis_invExp(&x, (double*) Y, &vy, 1, 1, 2, NULL);
        printf("x: %f vx: %f vy: %f error: %f\n", x, vx, vy, 2.0*(vx-vy)/(fabs(vx)+fabs(vy)));
      }
    }
  mlf_dealloc(f2dgrid);
  mlf_dealloc(fobj);
  mlf_dealloc(fbasis);
  mlf_quit();
  return 0;
}

