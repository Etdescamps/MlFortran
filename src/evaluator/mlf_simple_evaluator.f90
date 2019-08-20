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

Module mlf_simple_evaluator
  !$ Use OMP_LIB
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_errors
  Use mlf_rand
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_supervised_model
  IMPLICIT NONE
  PRIVATE
 
  Type, Public, Abstract, Extends(mlf_objective_fun) :: mlf_simple_stochastic_evaluator
    type(mlf_data_point), allocatable :: dataSet(:)
    real(c_double), allocatable :: ds_target(:,:)
    real(c_double), allocatable :: weights(:,:)
    type(mlf_local_experience_runner), allocatable :: runners(:)
    integer(8) :: nCstr, nOut, nIn
    integer :: nEval
  Contains
    procedure :: setup => mlf_simple_stochastic_setup
    procedure :: eval => mlf_simple_stochastic_eval
  End Type mlf_simple_stochastic_evaluator
Contains
  Integer Function EvalFunMedian(this, dS, X, Y, targ, Cstr, Cstr2, nE) Result(info)
    class(mlf_local_experience_runner), intent(inout) :: this
    type(mlf_data_point), intent(in), target :: dS
    real(c_double), intent(in) :: X(:), targ(:)
    real(c_double), intent(out) :: Y(:), Cstr(:), Cstr2(:)
    integer, intent(in) :: nE
    integer :: i, idX(nE)
    If(ALLOCATED(this%Z)) Then
      If(ANY(SIZE(this%Z) /= [SIZE(Y), nE])) DEALLOCATE(this%Z)
    Endif
    If(.NOT. ALLOCATED(this%Z)) ALLOCATE(this%Z(SIZE(Y), nE))
    Cstr = 0
    Do i = 1, nE
      info = this%runExp(dS%experiment, X, this%Z(:,i), Cstr2)
      Cstr = MAX(Cstr, Cstr2)
      If(info /= 0) Then
        Y = this%Z(:,i)
        RETURN
      Endif
    End Do
    Do i = 1, SIZE(Y)
      CALL QSortIdx(this%Z(i,:), idX)
      Y(i) = dS%weight*(Median(this%Z(i,:), idX) - targ(i))**2
    End Do
    info = 0
  End Function EvalFunMedian

  Integer Function mlf_simple_stochastic_setup(this, experiences, results, &
      initr, nR, nE, weights, ds_weights) Result(info)
    class(mlf_simple_stochastic_evaluator), intent(inout), target :: this
    class(mlf_model_experiment), intent(in) :: experiences(:)
    class(mlf_result_vectReal), intent(in) :: results(:)
    class(mlf_local_initializer), intent(inout) :: initr
    integer, intent(in) :: nR, nE
    real(c_double), intent(in), optional :: weights(:,:), ds_weights(:)
    integer :: i, N
    integer(8) :: M
    info = -1
    N = SIZE(experiences)
    If(SIZE(results) /= N .OR. nR <= 0) RETURN
    this%weights = weights
    M = results(1)%getNOutput()
    ALLOCATE(this%dataSet(N), this%ds_target(M,N), this%runners(nR))
    Do i = 1, nR
      info = initr%init_runner(this%runners(i))
      If(info < 0) RETURN
    End Do
    this%nEval = nE
    CALL this%runners(1)%params%getNParameters(this%nIn, this%nCstr)
    this%nCstr = this%nCstr + this%runners(1)%model%getNCstr()
    If(PRESENT(weights)) this%weights = weights
    Do i = 1, N
      this%dataSet(i)%experiment = experiences(i)
      this%dataSet(i)%results = results(i)
      info = results(i)%vectResults(this%ds_target(:,i))
      If(PRESENT(ds_weights)) Then
        this%dataSet(i)%weight = ds_weights(i)
      Else
        this%dataSet(i)%weight = 1d0
      Endif
    End Do
    info = 0
  End Function mlf_simple_stochastic_setup

  Integer Function mlf_simple_stochastic_eval(this, X, Y) Result(info)
    class(mlf_simple_stochastic_evaluator), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    class(mlf_local_experience_runner), pointer :: runner
    integer(8) :: i, j, Nx, Ny, Nd, Nr, id, nCstr
    real(c_double), allocatable :: Z(:), Cstr(:), Cstr2(:)
    info = -1
    Nx = SIZE(X,2)
    Ny = SIZE(Y,1)
    Nd = SIZE(this%dataSet)
    Nr = SIZE(this%runners)
    nCstr = this%nCstr
    Y = 0d0
    !$OMP PARALLEL num_threads(Nr) private(id, runner, info, Z)
    id = OMP_GET_THREAD_NUM()
    runner => this%runners(id)
    ALLOCATE(Z(this%nOut), Cstr(this%nCstr), Cstr2(this%nCstr))
    !$OMP DO collapse(2)
    Do i = 1, Nx
      Do j = 1, Nd
        info = EvalFunMedian(runner, this%dataSet(j), X(:,i), Z, this%ds_target(:,i), Cstr, Cstr2, this%nEval)
        !$OMP CRITICAL
        If(ALLOCATED(this%weights)) Then
          Y(nCstr+1:, i) = Y(nCstr+1:, i) + MATMUL(this%weights, Z)
        Else
          Y(nCstr+1:, i) = Y(nCstr+1:, i) + Z
        Endif
        Y(1:nCstr, i) = MAX(Y(1:nCstr, i), Cstr)
        !$OMP END CRITICAL
      End Do
    End Do
    !$OMP END DO
    DEALLOCATE(Z, Cstr, Cstr2)
    !$OMP END PARALLEL
    info = 0
  End Function mlf_simple_stochastic_eval
End Module mlf_simple_evaluator

