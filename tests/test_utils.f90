! Copyright (c) 2019 Etienne Descamps
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

Program test_utils
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  !CALL test_sel(128, 65)
  !CALL test_sel(129, 65)
  !CALL test_sel(127, 64)
  CALL test_median([5d0, 3d0, 2d0, 8d0, 1d0, 0d0, 6d0]) 
  CALL test_median([5d0, 4d0, 2d0, 8d0, 1d0, 0d0, 6d0, 7d0]) 
Contains
  Subroutine test_median(V)
    real(c_double) :: V(:)
    integer, allocatable :: K(:)
    ALLOCATE(K(SIZE(V)))
    CALL QSortIdx(V, K)
    PRINT *, Median(V, K)
  End Subroutine test_median
  !Subroutine test_sel(N, i)
  !  integer, intent(in) :: N, i
  !  real(c_double), allocatable :: V(:)
  !  real(c_double) :: x
  !  integer, allocatable :: K(:)
  !  ALLOCATE(V(N), K(N))
  !  CALL RANDOM_NUMBER(V)
  !  x = mlf_QuickSelect(V, i)
  !  CALL QSortIdx(V, K)
  !  PRINT *, x, V(K(i)), x == V(K(i))
  !End Subroutine test_sel
End
