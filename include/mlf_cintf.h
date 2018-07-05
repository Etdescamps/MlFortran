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

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MLF_MAXRANK 15 // Maximum rank in Fortran 2008

typedef void MLF_OBJ;

typedef enum {mlf_UNK = 0, mlf_BOOL = 1, mlf_INT = 2, mlf_INT64 = 3, mlf_SIZE = 4, mlf_FLOAT = 5, mlf_DOUBLE = 6, mlf_SIZEPARAM = 7, mlf_RAW = 8} MLF_DATATYPE;
typedef enum {mlf_DIRECT = 1, mlf_INDIRECT = 2, mlf_READONLY = 3, mlf_COPYONLY = 4, mlf_WRITEONLY = 5} MLF_ACCESSTYPE;
typedef enum {mlf_NAME = 1, mlf_DESC = 2, mlf_FIELDS = 3, mlf_FUNINFO} MLF_INFOTYPE;
typedef enum {mlf_OK = 0, mlf_UNINIT = -1, mlf_FUNERROR = -2, mlf_WRONGTYPE = -3, mlf_WRONGRANK = -4, mlf_OTHERERROR = -5} MLF_ERRORTYPE;

typedef struct {
  short dt; // MLF_DATATYPE
  short access; // MLF_ACCESSTYPE
} MLF_DT;

// Initialisation and exit function
// Shall be run respectivly before and after the use of the library
int mlf_init();
int mlf_quit();

// Universal deallocator for C-FORTRAN interface
int mlf_dealloc(MLF_OBJ *obj);

// Function for getting ressources
// Return pointer if dt == mlf_DIRECT, mlf_INDIRECT or mlf_READONLY
void *mlf_getrsc(MLF_OBJ *obj, int id, MLF_DT *dt, int *rank, int *dim0, void *data);

char *mlf_getinfo(MLF_OBJ *obj, int id, int type);

int mlf_updatersc(MLF_OBJ *obj, int id, void *data);

int64_t mlf_step(MLF_OBJ *obj, double *dt, int nstep);

// Utility functions
void mlf_randN(double *v, int N);

// Display a matrix using FORTRAN representation
void mlf_printMat(double *M, int nL, int nC);

#ifdef __cplusplus
}
#endif

