! Copyright (c) 2017-2018 Etienne Descamps
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation and/or
!    other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Module mlf_cfuns
  use ieee_arithmetic
  use iso_c_binding
  IMPLICIT NONE

  PRIVATE
  Public :: c_qsort, c_qsort_neg, c_qsort_unify, c_memcpy, c_strncmp, c_strlen
  Public :: mlf_stringFromC, c_cbrt, c_erf, c_erf_inv

  Interface
    Subroutine c_qsort(V, idSorted, N, ND, L, mu) Bind(C, NAME="mlf_qsort")
      Use iso_c_binding
      real(c_double), intent(in),dimension(*) :: V
      integer(c_int), intent(out),dimension(*) :: idSorted
      integer(c_int), value :: N, ND, L, mu
    End Subroutine c_qsort

    Subroutine c_qsort_neg(V, idSorted, N, ND, L, mu) Bind(C, NAME="mlf_qsort_neg")
      Use iso_c_binding
      real(c_double), intent(in),dimension(*) :: V
      integer(c_int), intent(out),dimension(*) :: idSorted
      integer(c_int), value :: N, ND, L, mu
    End Subroutine c_qsort_neg

    Function c_qsort_unify(V, idSorted, N, ND, L, mu) Bind(C, NAME="mlf_qsort_unify")
      Use iso_c_binding
      real(c_double), intent(in),dimension(*) :: V
      integer(c_int), intent(out),dimension(*) :: idSorted
      integer(c_int), value :: N, ND, L, mu
      integer(c_int) :: c_qsort_unify
    End Function c_qsort_unify

    Function c_memcpy(dest, src, n) Bind(C, NAME="memcpy")
      Use iso_c_binding
      type(c_ptr), value :: src, dest
      integer(c_size_t), value :: n
      type(c_ptr) :: c_memcpy
    End Function c_memcpy

    Function c_strncmp(s1, s2, n) Bind(C, NAME="strncmp")
      Use iso_c_binding
      type(c_ptr), value :: s1, s2
      integer(c_size_t), value :: n
      integer(c_int) ::  c_strncmp
    End Function c_strncmp

    Function c_strlen(s1) Bind(C, NAME="strlen")
      Use iso_c_binding
      type(c_ptr), value :: s1
      integer(c_size_t) ::  c_strlen
    End Function c_strlen

    Pure Function c_cbrt(r) Result(x) Bind(C, NAME="cbrt")
      Use iso_c_binding
      real(c_double), value :: r
      real(c_double) :: x
    End Function c_cbrt

    Function c_erf(r) Result(x) Bind(C, NAME="cpp_erf")
      Use iso_c_binding
      real(c_double), value :: r
      real(c_double) :: x
    End Function c_erf

    Function c_erf_inv(r) Result(x) Bind(C, NAME="cpp_erf_inv")
      Use iso_c_binding
      real(c_double), value :: r
      real(c_double) :: x
    End Function c_erf_inv
  End Interface
Contains
  Subroutine mlf_stringFromC(cptr, string)
    character(len=:, kind=c_char), allocatable, target, intent(out) :: string
    type(c_ptr), intent(in) :: cptr
    type(c_ptr) :: cout
    integer(c_size_t) :: sl
    sl = c_strlen(cptr)
    ALLOCATE(character(len=sl) :: string)
    cout = c_memcpy(c_loc(string), cptr, int(sl, kind=c_size_t))
  End Subroutine mlf_stringFromC
End Module mlf_cfuns

