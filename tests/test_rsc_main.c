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

#include <stdio.h>
#include "test_rsc.h"

int main(int argc, char **argv) {
  //mlf_init();

  MLF_OBJ *o = c_alloctest();
  char *desc = mlf_getinfo(o, 1, mlf_NAME);
  if(desc)
    printf("%s\n", desc);
  MLF_DT dt;
  int rank, dim0[5];
  void *data = mlf_getrsc(o, 1, &dt, &rank, (int*) dim0, NULL);
  printf("%d %d %d %d\n", dt.dt, dt.access, rank, dim0[0]);
  if(data && dt.dt == mlf_DOUBLE) {
    double *v = (double *) data;
    for(int i=0; i<dim0[0]; i++)
      printf("%f ", v[i]);
    printf("\n");
  }
  int ret = mlf_dealloc(o);

  //mlf_quit();

  return ret;
}

