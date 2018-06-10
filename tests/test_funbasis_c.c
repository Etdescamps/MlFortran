#include <stdio.h>
#include <math.h>
#include "mlf_model.h"
#include "mlf_funintf.h"

int fBasis_invExp(const double *x, const double *rpar, double *y, int nX, int nPar, int nY, void *ptr){
  for(int i=0; i<nPar; i++) {
    double a = rpar[nPar*i], b = rpar[nPar*i+1];
    for(int j=0; j<nX; j++) {
      y[i*nX*nY+j*nY] = 1.0/(1.0+a*exp(-b*x[j]))-1.0;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  mlf_init();
  MLF_OBJ *fobj = mlf_basisfunction(fBasis_invExp, NULL);
  mlf_dealloc(fobj);
  mlf_quit();
  return 0;
}

