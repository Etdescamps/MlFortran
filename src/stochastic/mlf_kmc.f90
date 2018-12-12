! Copyright (c) 2018 Etienne Descamps
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

Module mlf_kmc
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_supervised_model
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_kmc_reinit, mlf_kmc_init, mlf_kmc_reinitT
  Public :: mlf_kmc_stepFun
  
  ! Simple Kinetic Monte-Carlo class, used for comparison with hybrid methods.
  ! It is not the most efficient way to simulate this kind of models.
  Type, Public, Abstract, Extends(mlf_experience_model) :: mlf_kmc_model
    real(c_double), pointer :: t, Rates(:)
  Contains
    procedure (mlf_kmc_apply_action), deferred :: applyAction
    procedure (mlf_kmc_transition_rates), deferred :: funTransitionRates
    procedure :: stepF => mlf_kmc_stepFun
    procedure :: reinit => mlf_kmc_reinit
  End Type mlf_kmc_model

  Abstract Interface
    Integer Function mlf_kmc_transition_rates(this, Rates)
      Use iso_c_binding
      import :: mlf_kmc_model
      class(mlf_kmc_model), intent(inout), target :: this
      real(c_double), intent(out), target :: Rates(:)
    End Function mlf_kmc_transition_rates

    Integer Function mlf_kmc_apply_action(this, id, Rate)
      Use iso_c_binding
      import :: mlf_kmc_model
      class(mlf_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: Rate
      integer, intent(in) :: id
    End Function mlf_kmc_apply_action
  End Interface

Contains

  Integer Function mlf_kmc_init(this, numFields, t0, NActions, data_handler) Result(info)
    class(mlf_kmc_model), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: t0
    integer, intent(in), optional :: NActions
    integer(8) :: N
    N = -1
    CALL numFields%addFields(nRsc = 1, nRVar = 1)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info < 0) RETURN
    If(PRESENT(NActions)) N = NActions
    info = this%add_rarray(numFields, N, this%Rates, C_CHAR_"Rates", &
      data_handler = data_handler)
    If(info < 0) RETURN
    CALL this%addRVar(numFields, this%t, "t")
    If(.NOT. PRESENT(data_handler)) info = mlf_kmc_reinitT(this, t0)
  End Function mlf_kmc_init

  Integer Function mlf_kmc_reinit(this) Result(info)
    class(mlf_kmc_model), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    If(info < 0) RETURN
    this%t = 0d0
  End Function mlf_kmc_reinit

  Integer Function mlf_kmc_reinitT(this, t0) Result(info)
    class(mlf_kmc_model), intent(inout), target :: this
    real(c_double), intent(in), optional :: t0
    info = this%reinit()
    If(PRESENT(t0)) this%t = t0
  End Function mlf_kmc_reinitT

  Integer Function mlf_kmc_stepFun(this, niter) Result(info)
    class(mlf_kmc_model), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, niter0
    integer :: N, idAction
    real(c_double) :: U, r
    niter0=1; info = 0
    If(PRESENT(niter)) niter0 = niter
    Do i=1,niter0
      N = this%funTransitionRates(this%Rates)
      If(N <= 0) Then
        info = N
        If(N == 0) info = mlf_ODE_HardCstr
        RETURN
      Endif
      U = SUM(this%Rates(1:N))
      If(U <= 0) Then
        info = mlf_STEP_STOP
        If(U < 0) info = -1
        RETURN
      Endif
      CALL RANDOM_NUMBER(r)
      ! If r == 0 -> we consider the rare case where dt = 0
      If(r > 0) this%t = this%t-log(r)/U
      CALL RANDOM_NUMBER(r)
      r = r*U
      Do idAction = 1,N-1
        r = r - this%Rates(idAction)
        If(r <= 0) Then
          info = this%applyAction(idAction, this%Rates(idAction))
          EXIT
        Endif
      End Do
      If(r > 0) info = this%applyAction(N,this%Rates(N))
    End Do
  End Function mlf_kmc_stepFun

End Module mlf_kmc

