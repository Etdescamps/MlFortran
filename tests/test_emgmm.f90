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

Program test_emgmm
  use ieee_arithmetic
  use iso_c_binding
  use mlf_utils
  use mlf_kmeans
  use mlf_intf
  use mlf_emgmm
  use mlf_matrix
  use mlf_rand
  use mlf_gaussian
  use mlf_hdf5
  use test_common
  implicit none

  integer :: ND, NC, NX, info
  integer(kind=8) :: Nstep
  real(c_double) :: sigmaC, sigmaX

  info = mlf_init()

  if(COMMAND_ARGUMENT_COUNT()<6) then
    print *, "Error missing arguments: ND NC NX Nstep sigmaC sigmaX"
    stop
  endif

  ND = GetIntParameter(1)
  NC = GetIntParameter(2)
  NX = GetIntParameter(3)
  Nstep = GetIntParameter(4)
  sigmaC = GetRealParameter(5)
  sigmaX = GetRealParameter(6)
  print *,"ND: ", ND, " NC: ", NC, " NX: ", NX, " Nstep: ", Nstep, " sigmaC: ", sigmaC, " sigmaX: ", sigmaX
  call testOrtho()
  call testMN()
  call testGMM()

  info = mlf_quit()
contains
  subroutine testOrtho()
    real(c_double) :: C(ND,ND)
    call randOrthogonal(C)
    print *, "Random Ortho:"
    call PrintMatrix(C)
    print *, "Shall be Id:"
    call PrintMatrix(matmul(transpose(C), C))
  end subroutine testOrtho
  subroutine testMN()
    real(c_double), allocatable :: invC(:,:), C12(:,:), C(:,:), X(:,:)
    real(c_double), allocatable :: weight(:), P(:), Mu(:)
    real(c_double) ::  lnd
    integer :: i
    ALLOCATE(weight(ND), P(NX), Mu(ND), X(ND,NX), C12(ND,ND), C(ND,ND))
    P = 1
    weight = (/ (1d0/real(i), i=1,ND) /)
    weight = sqrt(sigmaX*weight/Mean(weight))
    call randOrthogonal(C12, weight)
    print *, "Random Org:"
    call PrintMatrix(C12)
    print *, "Covariance Matrix:"
    C = matmul(C12, transpose(C12))
    call PrintMatrix(C)
    allocate(invC(ND,ND))
    i = InverseSymMatrix(C, invC, lnd, .TRUE.)
    print *, "Inverse covariance matrix: lnd: ", lnd, " i: ", i
    call PrintMatrix(invC, matmul(C,invC))
    call randN(X, C12 = C12)
    lnd =  mlf_MaxGaussian(X, P, Mu, C)
    print *, "Evaluated Mu"
    print *, Mu
    print *, "Evaluated Covariance Matrix:"
    call PrintMatrix(C)
    i = mlf_EvalGaussian(X, P, Mu, C, 1d0, lnd)
    print *, "Evaluated Probabilities: lnd: ", lnd
    print *, P
  end subroutine testMN

  subroutine testGMM()
    real(c_double), allocatable :: XC(:,:), X(:,:), C12(:,:,:), weight(:)
    real(c_double) :: dt
    type(mlf_algo_emgmm) :: em_algo
    integer :: i, info, r
    integer, allocatable :: idx(:), idc(:)
    type(mlf_hdf5_file) :: mh5
    ALLOCATE(XC(ND,NC), X(ND,NX*NC), C12(ND,ND, NC), weight(ND), idx(NX*NC), idc(NC))
    ! Generates random centres
    call RandN(XC, sigmaC)
    ! Weights used as diagonal matrix (multiplied by orthogonal matrix)
    weight = (/ (1d0/real(i), i=1,ND) /)
    weight = sqrt(sigmaX*weight/Mean(weight))
    do i = 1,NC
      call randOrthogonal(C12(:,:,i), weight)
      call randN(X(:,1+(i-1)*NX:i*NX), X0 = XC(:,i), C12 = C12(:,:,i))
    end do
    idx = (/ (i, i=1,NX*NC) /)
    call randPerm(idx)
    X = X(:,idx)
    info = em_algo%init(X, nC)
    if(info<0) then
      print *,"Error initialization"
      RETURN
    endif
    r = em_algo%step(dt, Nstep)
    print *,"Time: ", dt
    print *, "Matrix Mu: (determined by GMM_EMAlgo)/ XC"
    call mlf_match_points(em_algo%Mu, XC, idc)
    call PrintMatrix(em_algo%Mu, XC(:, mlf_reverseid(idc)))
    do i=1,NC
      print *, "Matrix Covar:", i
      call PrintMatrix(em_algo%Cov(:,:,i), matmul(C12(:,:,idc(i)), transpose(C12(:,:,idc(i)))))
    end do
    info = mh5%createFile("emgmm.h5")
    info = mh5%pushState(em_algo)
    call mh5%finalize()
  end subroutine testGMM
End Program test_emgmm
