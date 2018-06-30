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

#include "mlf_cintf.h"

#ifdef __cplusplus
extern "C" {
#endif

// Objective function type used for multi/mono objective optimisation
typedef int (*mlf_objective_fun)(const double *x, double *y, int nD, int nY, int lambda, void *ptr);

MLF_OBJ *mlf_objfunction(mlf_objective_fun f, void *ptr, mlf_objective_fun cstr);

// Basis function type used for dimension reduction
typedef int (*mlf_basis_fun)(const double *x, const double *rpar, double *y, int nX, int nPar, int sPar, void *ptr);

MLF_OBJ *mlf_basisfunction(mlf_basis_fun f, void *ptr);

// Init/free function for .so binding using a data filename
typedef void *(*mlf_init_fun)(const char *filename);

typedef int (*mlf_free_fun)(void *data);

typedef char *(*mlf_getinfo_fun)(int type);

#ifdef __cplusplus
}
#endif

