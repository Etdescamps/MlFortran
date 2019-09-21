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

#ifdef __linux
#define _GNU_SOURCE
#endif

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "mlf_cfuns.h"

typedef union {
  int64_t i;
  double d;
} _dbl;

inline static int _vcompare(const double *v, int id1, int id2, uint32_t ND, uint32_t L) {
  const double *v1=&v[((long)id1)*((long)L)];
  const double *v2=&v[((long)id2)*((long)L)];
  for(int i=0; i<ND; i++) {
    if(v1[i] != v2[i]) {
      if(v1[i]<v2[i])
        return -1;
      else
        return 1;
    }
  }
  return 0;
}

inline static uint32_t _genMask(int N) {
  uint32_t k = mlf_countZ((uint32_t) N);
  return 0xFFFFFFFF >> k;
}

typedef struct {
  int64_t mask1;
  const double *v;
  uint32_t ND, N, L;
} _fqsortS;

#ifdef _GNU_SOURCE // GNU interface for qsort_r different from those of BSD and Windows
static int _fcompare(const void *a, const void *b, void *s) {
#else
static int _fcompare(void *s, const void *a, const void *b) {
#endif
  _dbl *va = (_dbl *) a;
  _dbl *vb = (_dbl *) b;
  _fqsortS *u = (_fqsortS*) s;
  int64_t mask1 = u->mask1;
  _dbl ua = {.i= va->i&mask1}, ub = {.i= vb->i&mask1};
  double d = ua.d-ub.d;
  // Shortcut that avoid cache miss
  if(d)
    return d<0 ? -1 : 1;
  // Shall use the array to determine the comparison
  int ki = va->i & (~mask1);
  int kj = vb->i & (~mask1);
  return _vcompare(u->v, ki, kj, u->ND, u->L);
}

static int _icompare(const void *a, const void *b) {
  int *va = (int *) a;
  int *vb = (int *) b;
  if(*va<*vb)
    return -1;
  if(*va>*vb)
    return 1;
  return 0;
}

static _dbl *_c_mlf_qsort(const double *v, int64_t mask0, int N, int ND, int L, int mu) {
  // TODO: avoid complete sort when mu << N
  
  // Use a mask trick in order to operate the sort on a single array of 64-bits words
  // Best strategy to avoid cache miss on big array
  int64_t mask1 = ~mask0;
  _fqsortS s = {mask1, v, ND, N, L};
  _dbl *u = (_dbl*) malloc(sizeof(_dbl)*N);
  if(ND==1) {
    memcpy(u, v, sizeof(_dbl)*N);
    for(int i=0; i<N; i++)
      u[i].i = (u[i].i & mask1) | i;
  }
  else {
    for(int i=0; i<N; i++) {
      u[i].d=v[i*L];
      u[i].i = (u[i].i & mask1) | i;
    }
  }
// Determine three "flavors" of qsort with external pointer
#if _WIN32 // Windows use a mix between GNU and BSD
  qsort_s((void *) u,  N, sizeof(_dbl), _fcompare, (void*) &s);
#elif defined(_GNU_SOURCE) // GNU libc interface ("GNU"/Linux)
  qsort_r((void *) u,  N, sizeof(_dbl), _fcompare, (void*) &s);
#else // For *BSD, macOS/iOS and other Unixes
  qsort_r((void *) u,  N, sizeof(_dbl), (void*) &s, _fcompare);
#endif
  return u;
}

/* void mlf_qsort(const double *v, int *idSorted, int N, int ND, int L, int mu) {
  int64_t mask0 = (int64_t) _genMask(N-1);
  _dbl *u = _c_mlf_qsort(v, mask0, N, ND, L, mu);
  for(int i=0; i<mu; i++)
    idSorted[i] = ((int) (u[i].i & mask0))+1;
  free(u);
} */

// Sort elements of a matrix using lexicographical order, and when equals, uses the index
void mlf_qsort(const double *v, int *idSorted, int N, int ND, int L, int mu) {
  int64_t mask0 = (int64_t) _genMask(N-1);
  int64_t mask1 = ~mask0;
  _dbl *u = _c_mlf_qsort(v, mask0, N, ND, L, mu);
  _dbl k1 = {.i = u[0].i & mask1};
  int j=0, l;
  idSorted[0] = ((int) (u[0].i & mask0))+1;
  for(int i=1; i<N; i++) {
    l = (int) (u[i].i & mask0);
    idSorted[i] = l+1;
    _dbl k2 = {.i = u[i].i & mask1};
    if(k1.d == k2.d && _vcompare(v, idSorted[j]-1, l, ND, L) == 0)
      continue;
    // else: sort the equal elements using the index
    if(i-1>j)
      qsort((void *) &idSorted[j],  i-j, sizeof(int), _icompare);
    k1 = k2;
    j = i;
    if(j>=mu-1)
      break;
  }
  if(j<mu-1)
    qsort((void *) &idSorted[j],  mu-j, sizeof(int), _icompare);
  free(u);
}

// Sort elements of a matrix using lexicographical order, and when equals, put negative id to the other elements
void mlf_qsort_neg(const double *v, int *idSorted, int N, int ND, int L, int mu) {
  int64_t mask0 = (int64_t) _genMask(N-1);
  int64_t mask1 = ~mask0;
  _dbl *u = _c_mlf_qsort(v, mask0, N, ND, L, mu);
  _dbl k1 = {.i = u[0].i & mask1};
  int j=0, l;
  idSorted[0] = ((int) (u[0].i & mask0))+1;
  for(int i=1; i<N; i++) {
    l = (int) (u[i].i & mask0);
    idSorted[i] = l+1;
    _dbl k2 = {.i = u[i].i & mask1};
    if(k1.d == k2.d && _vcompare(v, idSorted[j]-1, l, ND, L) == 0)
      continue;
    // else: sort the equal elements using the index
    if(i-1>j) {
      qsort((void *) &idSorted[j],  i-j, sizeof(int), _icompare);
      for(int k=j+1; k<i; k++)
        idSorted[k] = -idSorted[k];
    }
    k1 = k2;
    j = i;
    if(j>=mu-1)
      break;
  }
  if(j<mu-1) {
    qsort((void *) &idSorted[j],  mu-j, sizeof(int), _icompare);
    for(int k=j+1; k<mu; k++)
      idSorted[k] = -idSorted[k];
  }
  free(u);
}

// Remove double entries in the index array (useful for Pareto front)
int mlf_qsort_unify(const double *v, int *idSorted, int N, int ND, int L, int mu) {
  int64_t mask0 = (int64_t) _genMask(N-1);
  int64_t mask1 = ~mask0;
  _dbl *u = _c_mlf_qsort(v, mask0, N, ND, L, mu);
  _dbl k1 = {.i = u[0].i & mask1};
  int j = 1, l;
  idSorted[0] = ((int) (u[0].i & mask0))+1;
  for(int i=1; i<N; i++) {
    if(j >= mu)
      break;
    l = (int) (u[i].i & mask0);
    _dbl k2 = {.i = u[i].i & mask1};
    if(k1.d == k2.d && _vcompare(v, idSorted[j-1]-1, l, ND, L) == 0)
      continue;
    // else
    idSorted[j] = l+1;
    k1 = k2;
    j++;
  }
  free(u);
  return j;
}

