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

#pragma once
#include <stdlib.h>
// Function that counts the number of zero bits on top of an uint32
#if __LZCNT__
// Use the x64 instruction (works on gcc with the good march flags)
#include <x86intrin.h>
#define mlf_countZ __lzcnt32
#else
// Fallback function for other architectures
inline static uint32_t mlf_countZ(uint32_t k) {
  for(uint32_t i =0; i<32; i++) {
    if(k&0x80000000)
      return i;
    else
      k = (k<<1);
  }
  return 32;
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Implement a qsort routine for double format that output identifier (instead of sorting an array).
// This methods as the goal to minimize cache default using the low bits from the double for
// storing identifier.
void mlf_qsort(const double *v, int *idx, int N, int ND, int L, int mu);

// Revert index function that is used after qsort for rearranging arrays
int mlf_revert_idx(const double *v, const int *idPos, int *revId, int *newPos, int N, int ND, int L, int mu);
#ifdef __cplusplus
}
#endif

