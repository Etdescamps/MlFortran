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

Program test_gamma
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_gamma_dist
  Use mlf_histogram
  Use mlf_hdf5
  Use mlf_intf
  integer :: info
  type(mlf_hdf5_file) :: h5f
  info = mlf_init()
  info = h5f%createFile("test_gamma.h5")
  CALL test_random(h5f, "g1", 1d0, 1024, 1024)
  CALL test_random(h5f, "g0.5", 0.5d0, 1024, 1024)
  CALL test_random(h5f, "g2", 2d0, 1024, 1024)
  CALL h5f%finalize()
  info = mlf_quit()
Contains
  Subroutine test_random(h5f, rname, a, k, N)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: a
    integer, intent(in) :: k, N
    real(c_double), allocatable :: X(:), points(:,:), gammaV(:,:)
    integer, allocatable :: idx(:)
    type(mlf_histogram_int) :: hist
    integer :: info, i, j
    integer, parameter :: sizeX = 16384
    real(c_double) :: vMax = 10d0
    info = hist%init(0d0, vMax, k)
    ALLOCATE(X(sizeX), idx(sizeX), points(2,k+1), gammaV(2, 64*k))
    Do j = 1, N
      Do i = 1, sizeX
        Do
          X(i) = RandomGamma(a)
          If(X(i) < vMax) EXIT
        End DO
      End Do
      CALL QSortIdx(X, idx)
      X = X(idx)
      CALL hist%addVectSorted(X)
    End Do
    info = hist%get(points(2,:), points(1,:))
    points(2,:) = REAL(k-1, 8)*points(2,:)/vMax
    info = h5f%pushData(points, rname)
    FORALL(i=1:64*k) gammaV(1,i) = vMax*REAL(i-1,8)/REAL(64*k-1,8) 
    gammaV(2,:) = GammaDensity(a, gammaV(1,:))
    info = h5f%pushData(gammaV, rname//"_real")
  End Subroutine test_random
End
