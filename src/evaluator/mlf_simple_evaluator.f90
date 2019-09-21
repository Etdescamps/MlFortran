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
 
  Type, Public, Extends(mlf_local_experience_runner) :: mlf_simple_stochastic_runner
    real(c_double), allocatable :: Z(:,:), Y(:), Cstr(:), Cstr2(:)
    integer, allocatable :: idX(:)
  Contains
    procedure :: init_r => mlf_simple_stochastic_runner_init
  End Type mlf_simple_stochastic_runner

  Type, Public, Extends(mlf_objective_fun) :: mlf_simple_stochastic_evaluator
    type(mlf_data_point), allocatable :: dataSet(:)
    real(c_double), allocatable :: ds_target(:,:)
    real(c_double), allocatable :: weights(:,:)
    type(mlf_simple_stochastic_runner), allocatable :: runners(:)
    integer :: nEval
  Contains
    procedure :: setup => mlf_simple_stochastic_setup
    procedure :: eval => mlf_simple_stochastic_eval
  End Type mlf_simple_stochastic_evaluator
Contains
  Subroutine mlf_simple_stochastic_runner_init(this, nY, nCstr, nE)
    class(mlf_simple_stochastic_runner), intent(inout) :: this
    integer, intent(in) :: nY, nCstr, nE
    ALLOCATE(this%Z(nY, nE), this%Y(nY), this%idX(nE))
    ALLOCATE(this%Cstr(nCstr), this%Cstr2(nCstr))
  End Subroutine mlf_simple_stochastic_runner_init

  Integer Function EvalFunMedian(this, dS, X, targ) Result(info)
    class(mlf_simple_stochastic_runner), intent(inout) :: this
    type(mlf_data_point), intent(in), target :: dS
    real(c_double), intent(in) :: X(:), targ(:)
    integer :: i, nE
    this%Cstr = 0
    nE = SIZE(this%idX)
    Do i = 1, nE
      info = this%runExp(dS%experiment, X, this%Z(:,i), this%Cstr2)
      this%Cstr = MAX(this%Cstr, this%Cstr2)
      If(info /= 0) Then
        this%Y = HUGE(1d0)
        RETURN
      Endif
    End Do
    Do i = 1, SIZE(this%Y)
      CALL QSortIdx(this%Z(i,:), this%idX)
      this%Y(i) = dS%weight*(Median(this%Z(i,:), this%idX) - targ(i))**2
    End Do
    info = 0
  End Function EvalFunMedian

  Integer Function mlf_simple_stochastic_setup(this, experiences, results, &
      initr, nRunners, nEval, weights, ds_weights, statData) Result(info)
    class(mlf_simple_stochastic_evaluator), intent(inout), target :: this
    class(mlf_model_experiment), intent(in) :: experiences(:)
    class(mlf_result_vectReal), intent(in) :: results(:)
    class(mlf_local_initializer), intent(inout) :: initr
    type(mlf_local_experience_stats), optional, target :: statData(:)
    integer, intent(in) :: nRunners, nEval
    real(c_double), intent(in), optional :: weights(:,:), ds_weights(:)
    integer :: i, N, nOut
    info = -1
    N = SIZE(experiences)
    If(SIZE(results) /= N .OR. nRunners <= 0) RETURN
    this%weights = weights
    nOut = results(1)%getNOutput()
    ALLOCATE(this%dataSet(N), this%ds_target(nOut,N), this%runners(nRunners))
    Do i = 1, nRunners
      info = initr%init_runner(this%runners(i))
      If(info < 0) RETURN
      If(PRESENT(statData)) Then
        this%runners(i)%stats => statData(i)
      Else
        this%runners(i)%stats => NULL()
      Endif
    End Do
    this%nEval = nEval
    CALL this%runners(1)%params%getNParameters(this%nD, this%nC)
    this%nC = this%nC + this%runners(1)%model%getNCstr()
    If(PRESENT(weights)) Then
      this%weights = weights
      this%nY = SIZE(weights, 1)
    Else
      this%nY = nOut
    Endif
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
    Do i = 1, nRunners
      CALL this%runners(i)%init_r(nOut, this%nC, nEval)
    End Do
    info = 0
  End Function mlf_simple_stochastic_setup

  Integer Function mlf_simple_stochastic_eval(this, X, Y) Result(info)
    class(mlf_simple_stochastic_evaluator), intent(inout), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer :: i, j, Nx, Ny, Nd, Nr, id, nCstr
    info = -1
    Nx = SIZE(X,2)
    Ny = SIZE(Y,1)
    Nd = SIZE(this%dataSet)
    Nr = SIZE(this%runners)
    nCstr = this%nC
    Y = 0d0
    !$OMP PARALLEL num_threads(Nr) private(id, info) default(shared)
    id = OMP_GET_THREAD_NUM() + 1
    !$OMP DO collapse(2) schedule(dynamic)
    Do i = 1, Nx
      Do j = 1, Nd
        !PRINT *, "Start id:", id, "dataset:", j, "point:", i
        CALL this%runners(id)%statSetPoint(dataSet = j, paramSet = i)
        info = EvalFunMedian(this%runners(id), this%dataSet(j), X(:,i), this%ds_target(:,i))
        !$OMP CRITICAL
        If(ALLOCATED(this%weights)) Then
          Y(nCstr+1:, i) = Y(nCstr+1:, i) + MATMUL(this%weights, this%runners(id)%Y)
        Else
          Y(nCstr+1:, i) = Y(nCstr+1:, i) + this%runners(id)%Y
        Endif
        Y(1:nCstr, i) = MAX(Y(1:nCstr, i), this%runners(id)%Cstr)
        !PRINT *, "Finished id:", INT(id,2), "dataset:", INT(j,2), "point:", INT(i,2), "dist:", this%runners(id)%Y
        !$OMP END CRITICAL
      End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    info = 0
  End Function mlf_simple_stochastic_eval
End Module mlf_simple_evaluator

