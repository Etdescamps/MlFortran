#include <math.h>

int fun_objfun(const double *x, double *y, int nD, int nY, int lambda, void *ptr) {
  for(int i=0; i<lambda; ++i) {
    double y0 = 0;
    const double *x0 = &x[i*nD];
    for(int j=0; j<nD-1; j++) {
      double v = 10.0*(x0[j+1]-x0[j])*x0[j], w = x0[j]-1;
      y0 += v*v+w*w;
    }
    y[i*nY] = y0;
  }
  return 0;
}


