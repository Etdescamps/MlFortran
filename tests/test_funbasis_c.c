#include <stdio.h>
#include <math.h>
#include "mlf_model.h"
#include "mlf_funintf.h"

#define nG 16
#define nR 32

double params[nG][nR][2], weights[nG][nR];

double Y[2] = {0.5, 0.5}, W[16];

int fBasis_invExp(const double *x, const double *rpar, double *y, int nX, int nPar, int sPar, void *ptr){
  for(int i=0; i<nPar; i++) {
    double a = rpar[sPar*i], b = rpar[sPar*i+1];
    for(int j=0; j<nX; j++) {
      y[i*nX+j] = 1.0/(1.0+a*exp(-b*x[j]))-1.0;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  mlf_init();
  MLF_OBJ *fobj = mlf_basisfunction(fBasis_invExp, NULL);
  for(int i=0; i<nG; i++) {
    double vGi = exp(-((double) i)/(nG-1.0));
    double x = ((double) i) / ((double) nG-1);
    double w1 = 1.0*(1.0-x)+0.01*x;
    double w2 = 0.2*(1.0-x)+1.*x;
    for(int j=0; j<nR; j++) {
      double vRj = 5.0*exp(-((double) j)/(nR-1.0)*log(2.0));
      params[i][j][0]=vGi;
      params[i][j][1]=vRj;
      double y = ((double) j) / ((double) nR-1);
      weights[i][j] = w1*(1.0-y)+w2*y;
    }
  }
  MLF_OBJ *fbasis = mlf_funBasisInit(fobj, 2, 1.0, 0.0, HUGE_VAL, (double *) params, nR*nG, 16, 2000000, (double *) weights);
  mlf_getProj(fbasis, (double*) Y, (double*) W, 1);
  double x = 2.0;
  double vx = mlf_getValue(fbasis, (double*) W, x), vy;
  fBasis_invExp(&x, (double*) Y, &vy, 1, 1, 2, NULL);
  printf("vx: %f vy: %f error: %f\n", vx, vy, 2.0*(vx-vy)/(fabs(vx)+fabs(vy)));
  mlf_dealloc(fobj);
  mlf_dealloc(fbasis);
  mlf_quit();
  return 0;
}

