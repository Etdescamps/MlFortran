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

Program test_beta
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_gamma_dist
  Use mlf_beta_dist
  Use mlf_histogram
  Use mlf_hdf5
  Use mlf_intf
  integer :: info
  type(mlf_hdf5_file) :: h5f
  real(8) :: x
  info = mlf_init()
  x = InverseIncompleteBeta(0.5d0, 0.7d0, 0.98d0)
  info = h5f%createFile("test_beta.h5")
  CALL test_incompleteBeta(h5f, 'a2b2', 2d0, 2d0)
  CALL test_incompleteBeta(h5f, "a0.2b1.2", 0.2d0, 1.2d0)
  CALL test_incompleteBeta(h5f, "a0.5b0.7", 0.5d0, 0.7d0)
  !CALL test_random(h5f, "a2b2", 2d0, 2d0, 1024, 1024)
  !CALL test_random(h5f, "a0.2b1.2", 0.2d0, 1.2d0, 1024, 1024)
  !CALL test_random(h5f, "a0.5b0.7", 0.5d0, 0.7d0, 1024, 1024)
  CALL h5f%finalize()
  !CALL test_likelihood(10000000, 2d0, 2d0)
  !CALL test_likelihood(10000000, 5d0, 1.5d0)
  !CALL test_likelihood(10000000, 2d0, 0.5d0)
  !CALL test_likelihood(10000000, 0.3d0, 0.6d0)
  !CALL test_likelihood(10000000, 60d0, 50d0)
  info = mlf_quit()
Contains
  Subroutine test_incompleteBeta(h5f, rname, alpha, beta)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: alpha, beta
    real(c_double), allocatable :: points(:,:)
    integer, parameter :: N = 1000
    integer :: i
    ALLOCATE(points(3, N))
    FORALL(i=1:N) points(1, i) = REAL(i-1,8)/REAL(N-1,8)
    points(2, :) = IncompleteBeta(alpha, beta, points(1, :))
    Do i = 1,N
      points(3, i) = InverseIncompleteBeta(alpha, beta, points(2, i))
    End Do
    info = h5f%pushData(points, rname)
  End Subroutine test_incompleteBeta

  Subroutine test_likelihood(N, alpha, beta)
    real(c_double), intent(in) :: alpha, beta
    integer, intent(in) :: N
    real(c_double), allocatable :: X(:), W(:)
    real(c_double) :: a0, b0
    integer :: i, info
    ALLOCATE(X(N), W(N))
    Do i = 1, N
      Do
        X(i) = RandomBeta(alpha,beta)
        If(X(i) /= 0d0 .AND. X(i) /= 1d0) EXIT
      End Do
    End Do
    CALL RANDOM_NUMBER(W)
    W = 1-W
    info = MaxLikelihoodBeta(X, a0, b0, W)
    PRINT *, alpha, a0, beta, b0
  End Subroutine test_likelihood

  Subroutine test_random(h5f, rname, a, b, k, N)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: a, b
    integer, intent(in) :: k, N
    real(c_double), allocatable :: X(:), points(:,:), betaV(:,:)
    integer, allocatable :: idx(:)
    type(mlf_histogram_int) :: hist
    integer :: info, i, j
    integer, parameter :: sizeX = 16384
    info = hist%init(0d0, 1d0, k)
    ALLOCATE(X(sizeX), idx(sizeX), points(2,k+1), betaV(2, 64*k))
    Do j = 1, N
      Do i = 1, sizeX
        X(i) = RandomBeta(a, b)
      End Do
      CALL QSortIdx(X, idx)
      X = X(idx)
      CALL hist%addVectSorted(X)
    End Do
    info = hist%get(points(2,:), points(1,:))
    points(2,:) = REAL(k-1, 8)*points(2,:)
    info = h5f%pushData(points, rname)
    FORALL(i=1:64*k) betaV(1,i) = REAL(i-1,8)/REAL(64*k-1,8) 
    betaV(2,:) = BetaDensity(a, b, betaV(1,:))
    info = h5f%pushData(betaV, rname//"_real")
  End Subroutine test_random
End
