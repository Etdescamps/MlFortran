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
  Use mlf_kmeans_naive
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE
  
  ! Algorithm structure of Expectation-Maximization on Gaussian Mixture Model
  Type, public, extends(mlf_step_obj) :: mlf_algo_emgmm
    real(c_double), pointer :: X(:,:), Mu(:,:), Cov(:,:,:), lambda(:)
    real(c_double), pointer :: sumLL, minProba
    real(c_double), allocatable :: ProbaC(:,:)
    integer(c_int32_t), pointer :: cl(:)
    logical :: initialized
  Contains
    procedure :: init => mlf_algo_emgmm_init
    procedure :: reinit => mlf_algo_emgmm_reinit
    procedure :: stepF => mlf_algo_emgmm_stepF
  End Type mlf_algo_emgmm

  Type, Public, extends(mlf_class_proba_model) :: mlf_model_emgmm
    class(mlf_algo_emgmm), pointer :: top
  Contains
    procedure :: getProba => mlf_algo_emgmm_getLikelihood
    procedure :: getNumClasses => mlf_algo_emgmm_getNumClasses
  End Type mlf_model_emgmm
Contains
  ! Model initilisator
  Subroutine mlf_model_emgmm_init(this, top)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_algo_emgmm), intent(in), target :: top
    ALLOCATE(mlf_model_emgmm :: this)
    Select Type(this)
    Class is (mlf_model_emgmm)
      this%top => top
    End Select
  End Subroutine mlf_model_emgmm_init

  Integer Function mlf_algo_emgmm_init(this, X, nC, minProba, data_handler) Result(info)
    class(mlf_algo_emgmm), intent(out), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), optional :: minProba
    real(c_double), target :: X(:,:)
    integer, intent(inout) :: nC
    type(mlf_step_numFields) :: numFields
    integer(c_int64_t) :: nY, nX, ndMu(2), ndC(3), nC2
    CALL numFields%initFields(nRVar = 1, nRsc = 3, nIVar = 1)
    info = mlf_step_obj_init(this, numFields, data_handler = data_handler)
    this%X => X
    If(info < 0) RETURN
    nY = SIZE(X,1); nX = SIZE(X,2); nC2 = nC
    ndMu = [nY, nC2]; ndC = [nY, nY, nC2]
    info = this%add_rmatrix(numFields, ndMu, this%Mu, C_CHAR_"Mu", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    If(CheckF(info, "emgmm: error creating Mu")) RETURN
    ndC(3) = ndMu(2)
    info = this%add_rmatrix3d(numFields, ndC, this%Cov, C_CHAR_"Cov", &
      data_handler = data_handler, fixed_dims = [.TRUE., .TRUE., .TRUE.])
    If(CheckF(info, "emgmm: error creating Cov")) RETURN
    nC2 = ndMu(2)
    info = this%add_rarray(numFields, nC2, this%lambda, C_CHAR_"lambda", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    If(CheckF(info, "emgmm: error creating lambda")) RETURN
    nC = INT(nC2, kind=4)
    ALLOCATE(this%ProbaC(nX, nC))
    CALL this%addRVar(numFields, this%sumLL, "sumLL")
    CALL this%addRPar(numFields, this%minProba, "minProba")
    CALL InitOrDefault(this%minProba, 0.0d0, minProba)
    If(PRESENT(data_handler)) Then
      this%initialized = .TRUE.
    Else
      info = this%reinit()
    Endif
    CALL mlf_model_emgmm_init(this%model, this)
  End Function mlf_algo_emgmm_init

  ! Reinit algorithm parameters
  integer Function mlf_algo_emgmm_reinit(this) Result(info)
    class(mlf_algo_emgmm), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    If(associated(this%Mu)) this%Mu = 0
    If(associated(this%Cov)) this%Cov = 0
    If(allocated(this%ProbaC)) this%ProbaC = 0
    this%initialized = .FALSE.
  End Function mlf_algo_emgmm_reinit

  ! Step function
  Integer Function mlf_algo_emgmm_stepF(this, niter) Result(info)
    class(mlf_algo_emgmm), intent(inout), target :: this
    type(mlf_algo_kmeans_naive) :: km_algo
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, N, nkm
    integer :: j,nC
    N = 1; nkm = 16
    nC = SIZE(this%Mu,2)
    If(PRESENT(niter)) N = niter
    Do i=1,N
      If(.NOT. this%initialized) Then
        info = km_algo%init(this%X, nC = nC)
        If(info<0) RETURN
        info = km_algo%stepF(nkm)
        If(info<0) RETURN
        this%Mu=km_algo%Mu
        this%initialized = .TRUE.
        Do j=1,nC
          CALL IdentityMatrix(this%Cov(:,:,j))
        End Do
        this%lambda = 1d0/REAL(nC, KIND=c_double)
      Else
        info = ExpStep(this%X, this%ProbaC, this%Mu, this%Cov, this%Lambda, this%sumLL)
        If(info /= 0) RETURN
        CALL MaxStep(this%X, this%ProbaC, this%Mu, this%Cov, this%Lambda, this%minProba)
      Endif
    End Do
    info = 0
  End Function mlf_algo_emgmm_stepF

  Integer Function mlf_algo_emgmm_getNumClasses(this) Result(nC)
    class(mlf_model_emgmm), intent(in), target :: this
    nC = SIZE(this%top%Mu,2)
  End Function mlf_algo_emgmm_getNumClasses

  ! Function that determines likelihood of input datapoints
  Integer Function mlf_algo_emgmm_getLikelihood(this, X, Proba) Result(info)
    class(mlf_model_emgmm), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Proba(:,:)
    real(c_double) :: sumLL
    If(.NOT. this%top%initialized) Then
      info = mlf_UNINIT
      RETURN
    Endif
    info = ExpStep(X, Proba, this%top%Mu, this%top%Cov, this%top%Lambda, sumLL)
  End Function mlf_algo_emgmm_getLikelihood

  ! Maximum likelihood step of the EM algorithm
  Subroutine MaxStep(X, Proba, Mu, Covar, Lambda, minProba)
    real(c_double), intent(in) :: X(:,:), Proba(:,:), minProba
    real(c_double), intent(out) :: Mu(:,:), Covar(:,:,:), Lambda(:)
    integer :: Nmix, i
    Nmix = size(Proba,2)
    !$OMP PARALLEL DO default(shared) PRIVATE(i)
    Do i = 1,Nmix
      Lambda(i) = mlf_MaxGaussian(X, Proba(:,i), Mu(:,i), Covar(:,:,i))
    End Do
    !$OMP END PARALLEL DO
    If(minProba > 0) Then
      Lambda(:) = MAX(Lambda(:), minProba*MAXVAL(Lambda))
    Endif
    Lambda(:) = Lambda(:)/sum(Lambda)
  End Subroutine MaxStep

  ! Expectation step: compute responsabilities
  Integer(c_int) Function ExpStep(X, Proba, Mu, Covar, Lambda, sumLL) Result(info)
    real(c_double), intent(in) :: X(:,:), Mu(:,:), Covar(:,:,:), Lambda(:)
    real(c_double), intent(out) :: Proba(:,:), sumLL
    integer :: ND, NX, Nmix, i, infoI
    real(c_double) :: LL
    ND = size(X,1); NX = size(X,2); Nmix = size(Mu,2)
    sumLL = 0
    info = 0 ! Remove warning
    !$OMP PARALLEL DO default(shared) PRIVATE(i, LL, infoI) REDUCTION(+:sumLL)
    Do i = 1,Nmix
      infoI = mlf_EvalGaussian(X, Proba(:,i), Mu(:,i), Covar(:,:,i), Lambda(i), LL)
      sumLL = sumLL + LL
    End Do
    !$OMP END PARALLEL DO
    !Normalize probabilities
    FORALL(i=1:NX) Proba(i,:) = Proba(i,:)/SUM(Proba(i,:))
  End Function ExpStep


End Module mlf_emgmm
