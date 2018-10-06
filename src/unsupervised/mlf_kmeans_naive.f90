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

Module mlf_kmeans_naive
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_models
  Use mlf_utils
  Use mlf_rand
  Use mlf_errors
  Use mlf_kmeans
  IMPLICIT NONE
  
  PRIVATE

  Type, Public, extends(mlf_algo_kmeans) :: mlf_algo_kmeans_naive
  Contains
    procedure :: stepF => mlf_algo_kmeans_stepF
  End Type mlf_algo_kmeans_naive
  
  Public :: mlf_kmeans_c

Contains
  ! Simple K-means algorithm implementation
  ! Algorithm step function
  integer Function mlf_algo_kmeans_stepF(this, niter) result(info)
    class(mlf_algo_kmeans_naive), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, N
    N = 1
    if(present(niter)) N = niter
    do i=1,N
      if(.NOT. this%initialized) then
        call mlf_InitKMeans(this%X, this%Mu)
        this%initialized = .TRUE.
      else
        call mlf_EvaluateClass(this%X, this%Mu, this%cl, this%minDist)
        this%meanDist = mean(this%minDist)
        call EvaluateCentres(this%X, this%Mu, this%cl)
      endif
    end do
    info = 0
  End Function mlf_algo_kmeans_stepF

  Subroutine EvaluateCentres(X, Mu, cl)
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Mu(:,:)
    integer, intent(in) :: cl(:)
    integer :: NX, NC, i, k, Ncl(size(cl,1))
    NX = size(X,2)
    NC = size(Mu,2)
    Mu = 0
    Ncl = 0
    Do i=1, NX
      k = cl(i)
      Ncl(k) = Ncl(k)+1
      Mu(:,k) = Mu(:,k)+X(:,i)
    End Do
    Do i=1,NC
      if(Ncl(i) == 0) then
        ! Set to the farest position
        Mu(:,i) = mlf_GetFarPoint(X,Mu)
      else
        Mu(:,i) = Mu(:,i)/real(Ncl(i), kind=c_double)
      endif
    End Do
  End Subroutine EvaluateCentres

  ! C interface to k-means algorithm
  type(c_ptr) Function mlf_kmeans_c(pX, nX, nY, nC, pMu) result(cptr) bind(C, name="mlf_kmeans_c")
    type(c_ptr), value :: pX, pMu
    real(c_double), pointer :: X(:,:), Mu(:,:)
    integer(c_int), value :: nX, nY, nC
    integer :: info
    type(mlf_algo_kmeans_naive), pointer :: this
    class (*), pointer :: obj
    ALLOCATE(this)
    call C_F_POINTER(pX, X, [nY, nX])
    if(C_ASSOCIATED(pMu)) then
      call C_F_POINTER(pMu, Mu, [nC, nX])
      info = this%init(X, Mu = Mu)
    else
      info = this%init(X, nC = nC)
    endif
    obj => this
    cptr = c_allocate(obj)
  End Function mlf_kmeans_c


End Module mlf_kmeans_naive

