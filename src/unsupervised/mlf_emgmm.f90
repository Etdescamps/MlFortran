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

Module mlf_emgmm
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_models
  Use mlf_utils
  Use mlf_rand
  Use mlf_matrix
  Use mlf_gaussian
  Use mlf_kmeans
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE
  
  ! Algorithm structure of Expectation-Maximization on Gaussian Mixture Model
  Type, public, extends(mlf_class_proba_model) :: mlf_algo_emgmm
    real(c_double), pointer :: X(:,:), Mu(:,:), Cov(:,:,:), lambda(:), sumLL
    real(c_double), allocatable :: ProbaC(:,:)
    integer(c_int32_t), pointer :: cl(:)
    logical :: initialized
  Contains
    procedure :: init => mlf_algo_emgmm_init
    procedure :: reinit => mlf_algo_emgmm_reinit
    procedure :: stepF => mlf_algo_emgmm_stepF
    procedure :: getProba => mlf_algo_emgmm_getLikelihood
    procedure :: getNumClasses => mlf_algo_emgmm_getNumClasses
  End Type mlf_algo_emgmm

Contains
  integer Function mlf_algo_emgmm_init(this, X, nC, data_handler) result(info)
    class(mlf_algo_emgmm), intent(out), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), target :: X(:,:)
    integer, intent(inout) :: nC
    integer :: nrsc
    integer(c_int64_t) :: nipar, nrpar, nY, nX, ndMu(2), ndC(3), nC2
    nipar = 0; nrpar = 1; nrsc = 3
    info = mlf_step_obj_init(this, nipar, nrpar, nrsc, C_CHAR_"", C_CHAR_"sumLL;", data_handler = data_handler)
    this%X => X
    if(info < 0) RETURN
    nY = size(X,1); nX = size(X,2); nC2 = nC
    ndMu = [nY, nC2]; ndC = [nY, nY, nC2]
    info = this%add_rmatrix(nrsc, ndMu, this%Mu, C_CHAR_"Mu", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    if(CheckF(info, "emgmm: error creating Mu")) RETURN
    ndC(3) = ndMu(2)
    info = this%add_rmatrix3d(nrsc+1, ndC, this%Cov, C_CHAR_"Cov", &
      data_handler = data_handler, fixed_dims = [.TRUE., .TRUE., .TRUE.])
    if(CheckF(info, "emgmm: error creating Cov")) RETURN
    nC2 = ndMu(2)
    info = this%add_rarray(nrsc+2, nC2, this%lambda, C_CHAR_"lambda", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(CheckF(info, "emgmm: error creating lambda")) RETURN
    nC = int(nC2, kind=4)
    ALLOCATE(this%ProbaC(nX, nC))
    this%sumLL => this%rpar(nrpar)
    if(present(data_handler)) then
      this%initialized = .TRUE.
    else
      info = this%reinit()
    endif
  End Function mlf_algo_emgmm_init

  ! Reinit algorithm parameters
  integer Function mlf_algo_emgmm_reinit(this) result(info)
    class(mlf_algo_emgmm), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    if(associated(this%Mu)) this%Mu = 0
    if(associated(this%Cov)) this%Cov = 0
    if(allocated(this%ProbaC)) this%ProbaC = 0
    this%initialized = .FALSE.
  End Function mlf_algo_emgmm_reinit

  ! Step function
  integer Function mlf_algo_emgmm_stepF(this, niter) result(info)
    class(mlf_algo_emgmm), intent(inout), target :: this
    type(mlf_algo_kmeans) :: km_algo
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, N, nkm = 16
    integer :: j,nC
    N = 1
    nC = size(this%Mu,2)
    if(present(niter)) N = niter
    do i=1,N
      if(.NOT. this%initialized) then
        info = km_algo%init(this%X, nC = nC)
        if(info<0) RETURN
        info = km_algo%stepF(nkm)
        if(info<0) RETURN
        this%Mu=km_algo%Mu
        this%initialized = .TRUE.
        do j=1,nC
          call IdentityMatrix(this%Cov(:,:,j))
        end do
        this%lambda = 1d0/real(nC, kind=c_double)
      else
        info = ExpStep(this%X, this%ProbaC, this%Mu, this%Cov, this%Lambda, this%sumLL)
        if(info /= 0) RETURN
        call MaxStep(this%X, this%ProbaC, this%Mu, this%Cov, this%Lambda)
      endif
    end do
    info = 0
  End Function mlf_algo_emgmm_stepF

  integer Function mlf_algo_emgmm_getNumClasses(this) result(nC)
    class(mlf_algo_emgmm), intent(in), target :: this
    nC = size(this%Mu,2)
  End Function mlf_algo_emgmm_getNumClasses

  ! Function that determines likelihood of input datapoints
  integer Function mlf_algo_emgmm_getLikelihood(this, X, Proba) result(info)
    class(mlf_algo_emgmm), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Proba(:,:)
    real(c_double) :: sumLL
    if(.NOT. this%initialized) then
      info = mlf_UNINIT
      RETURN
    endif
    info = ExpStep(X, Proba, this%Mu, this%Cov, this%Lambda, sumLL)
  End Function mlf_algo_emgmm_getLikelihood

  ! Maximum likelihood step of the EM algorithm
  subroutine MaxStep(X, Proba, Mu, Covar, Lambda)
    real(c_double), intent(in) :: X(:,:), Proba(:,:)
    real(c_double), intent(out) :: Mu(:,:), Covar(:,:,:), Lambda(:)
    integer :: Nmix, i
    Nmix = size(Proba,2)
    !$OMP PARALLEL DO
    do i = 1,Nmix
      Lambda(i) = mlf_MaxGaussian(X, Proba(:,i), Mu(:,i), Covar(:,:,i))
    end do
    !$OMP END PARALLEL DO
    Lambda(:) = Lambda(:)/sum(Lambda)
  end subroutine MaxStep

  ! Expectation step: compute responsabilities
  integer(c_int) function ExpStep(X, Proba, Mu, Covar, Lambda, sumLL) result(info)
    real(c_double), intent(in) :: X(:,:), Mu(:,:), Covar(:,:,:), Lambda(:)
    real(c_double), intent(out) :: Proba(:,:), sumLL
    integer :: ND, NX, Nmix, i
    real(c_double) :: LL
    ND = size(X,1); NX = size(X,2); Nmix = size(Mu,2)
    sumLL = 0
    info = 0 ! Remove warning
    do i = 1,Nmix
      info = mlf_EvalGaussian(X, Proba(:,i), Mu(:,i), Covar(:,:,i), Lambda(i), LL)
      if(info /= 0) RETURN
      sumLL = sumLL + LL
    end do
    !Normalize probabilities
    forall(i=1:NX) Proba(i,:) = Proba(i,:)/sum(Proba(i,:))
  end function ExpStep


End Module mlf_emgmm
