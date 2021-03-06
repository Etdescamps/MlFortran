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

Program test_kmeans
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_utils
  Use mlf_rand
  Use mlf_kmeans_naive
  Use mlf_kmeans
  IMPLICIT NONE
  CALL testkmeans()
CONTAINS
  Subroutine testkmeans()
    integer, parameter :: nX = 1024, nY = 16, nC = 6
    real(c_double), allocatable, target :: X(:,:), Mu(:,:)
    integer, allocatable :: idx(:)
    real(c_double) :: dt
    integer(kind=8) :: nstep = 32
    type(mlf_algo_kmeans_naive) :: km_algo
    integer :: i, info

    info = mlf_init()

    ALLOCATE(X(nY,nX*nC), Mu(nY, nC), idx(NC))
    CALL RandN(Mu, 10d0)
    CALL PrintMatrix(Mu)
    Do i=1,nC
      CALL RandN(X(:,(1+(i-1)*nX):(i*nX)), 1d0, Mu(:,i))
    End Do
    info = km_algo%init(X, nC)
    info = km_algo%step(dt, nstep)
    PRINT *, dt
    CALL mlf_match_points(Mu, km_algo%Mu, idx)
    CALL PrintMatrix(km_algo%Mu(:,idx))
    CALL km_algo%finalize()
    info = mlf_quit()
  End Subroutine testkmeans
End program test_kmeans

