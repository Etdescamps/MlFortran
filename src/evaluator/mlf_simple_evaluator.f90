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
    integer :: nEval, nCstr, nOut
  Contains
    procedure :: eval => mlf_simple_stochastic_eval
  End Type mlf_simple_stochastic_evaluator
Contains
  Integer Function EvalFunMedian(this, dS, X, Y, targ, Cstr, nE) Result(info)
    class(mlf_local_experience_runner), intent(inout) :: this
    type(mlf_data_point), intent(in), target :: dS
    real(c_double), intent(in) :: X(:), targ(:)
    real(c_double), intent(out) :: Y(:), Cstr(:)
    integer, intent(in) :: nE
    real(c_double) :: Cstr2(SIZE(Cstr))
    integer :: i, idX(nE), nCstr
    If(info < 0) RETURN
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
  End Function EvalFunMedian

  Integer Function mlf_simple_stochastic_eval(this, X, Y) Result(info)
    class(mlf_simple_stochastic_evaluator), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    class(mlf_local_experience_runner), pointer :: runner
    integer :: i, j, k, Nx, Ny, Nd, Nr, id, nOut, nCstr
    real(c_double), allocatable :: Z(:), Cstr(:)
    Nx = SIZE(X,2)
    Ny = SIZE(Y,1)
    Nd = SIZE(this%dataSet)
    Nr = SIZE(this%runners)
    nCstr = this%nCstr
    Y = 0d0
    !$OMP PARALLEL num_threads(Nr) private(id, runner, info, Z)
    id = OMP_GET_THREAD_NUM()
    runner => this%runners(id)
    ALLOCATE(Z(this%nOut), Cstr(this%nCstr))
    !$OMP DO collapse(2)
    Do i = 1, Nx
      Do j = 1, Nd
        info = EvalFunMedian(runner, this%dataSet(j), X(:,i), Z, this%ds_target(:,i), Cstr, this%nEval)
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
    DEALLOCATE(Z)
    !$OMP END PARALLEL
  End Function mlf_simple_stochastic_eval
End Module mlf_simple_evaluator

