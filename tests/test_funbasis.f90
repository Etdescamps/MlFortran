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

Program test_funbasis
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_utils
  Use mlf_matrix
  Use mlf_fun_intf
  Use mlf_funbasis
  Use test_functions
  IMPLICIT NONE

  integer :: infoX
  infoX = mlf_init()
  infoX = RunTest()

  infoX = mlf_quit()

Contains
  integer function RunTest() result(info)
    type(mlf_basis_test) :: f_test
    integer, parameter :: nR = 32, nG = 16
    real(c_double) :: YB(2,nG,nR), a0, a1, b0, b1, YT(2,2), WB(nG,nR), W1, W0, X, Y
    integer :: i, j
    real(c_double), parameter :: Ri(nR) = [(5d-2*exp(i/real(nR-1)*log(2d1)), i=0,(NR-1))]
    real(c_double), parameter :: Gi(nG) = [(exp(-i/real(nG-1)), i=0,(nG-1))]

    info = mlf_init()
 
    a0 = -1.5; a1 = 2.0
    b0 = 0.1; b1 = 10.0
    forall(i=1:nG) YB(2,i,:) = Ri
    forall(i=1:nR) YB(1,:,i) = Gi
    do i=1,nG
      X = (i-1d0)/real(nG-1)
      W0 = 1d0*(1d0-X)+0.01*X
      W1 = 0.2d0*(1d0-X)+ 1d0*X
      do j=1,nR
        Y = (j-1d0)/real(nR-1)
        WB(i,j) = W0*(1d0-Y)+W1*Y
      End do
    End do
    ! print *, YB
    YT(:,1) = YB(:,1,2)
    YT(:,2) = YB(:,nG-5,nR-2)
    f_test%evalT => FExpInv  
    info = testfb(f_test, reshape(YB, [2,nR*nG]), YT, reshape(sqrt(WB), [nR*nG]))

  End function RunTest

  integer function testfb(f_test, Ybasis, Ytests, WP) result(info)
    class(mlf_basis_test), intent(inout) :: f_test
    type(mlf_algo_funbasis) :: fb
    real(c_double), intent(in) :: Ybasis(:,:), Ytests(:,:), WP(:)
    real(c_double) :: infty
    real(c_double), parameter :: ferr = 1d-8, alpha = 1d0
    real(c_double), allocatable :: W(:,:), AError(:,:), U1(:,:), U2(:,:)
    integer, parameter :: nparam = 8, nx=16, nxp=65536
    integer :: nt, i, j
    real(c_double) :: X0(nx) =  [(-1d0/alpha*log(1d0-i/real(nx)),i=0,(nx-1))]
    infty = ieee_value(infty, ieee_positive_inf)
    info = fb%initF(f_test, alpha, 0d0, infty, Ybasis, nparam, nxp, WP)
    if(info /= 0) then
      print *, "fb%init error ", info
      RETURN
    endif
    nt = size(Ytests,1)
    allocate(W(nparam,nt), Aerror(nt,2), U1(nx,nt), U2(nx,nt))
    print *,ferr
    info = fb%getProj(Ytests, W, Aerror)
    print *, Aerror(:,1)
    print *, Aerror(:,2)
    call FExpInv(X0, Ytests, U1(:,:))
    forall(i=1:nt, j=1:nx) U2(j,i) = (fb%getValue(W(:,i), X0(j))-U1(j,i))*exp(-alpha*X0(j))
    print *, U2(:,:)
  end function testfb
End Program test_funbasis

