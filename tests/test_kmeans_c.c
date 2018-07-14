/*-
 * Copyright (c) 2017-2018 Etienne Descamps
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation and/or
 *    other materials provided with the distribution.

 * 3. Neither the name of the copyright holder nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "mlf_cintf.h"
#include <stdio.h>
#define nX 1024
#define nY 16
#define nC 6

double X[nC][nX][nY], Mu[nC][nY];

int main(int argc, char **argv) {
  double t;

  mlf_init();

  mlf_randN((double*) Mu, nC*nY);
  for(int i=0; i<nC; i++)
    for(int j=0; j<nY; j++)
      Mu[i][j] *= 10.0;
  mlf_printMat((double *) Mu, nY, nC);
  mlf_randN((double*) X, nX*nY*nC);
  for(int i=0; i<nC; i++)
    for(int j=0; j<nX; j++)
      for(int k=0; k<nY; k++)
        X[i][j][k] += Mu[i][k];
  MLF_OBJ *obj=mlf_kmeans_c((double *) X, nX*nC, nY, nC, (double *) NULL);
  MLF_DT dt;
  int64_t nstep = 32;
  int info = mlf_step(obj, &t, &nstep);
  if(info<0)
    return info;
  int rank, dim0[18];
  double *V = (double *) mlf_getrsc(obj, 3, &dt, &rank, dim0, NULL);
  printf("%f\n", t);
  mlf_printMat(V, nY, nC);
  mlf_dealloc(obj);

  mlf_quit();
}

