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

Module mlf_ode_class
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_supervised_model
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_ode_algo_init, mlf_ode_algo_reinit, mlf_init_ode_model
  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_ode_algo
    real(c_double), pointer :: X0(:), X(:), K(:,:)
    real(c_double), pointer :: t0, t, tMax, hMax
    real(c_double), allocatable :: Cont(:,:)
    real(c_double) :: lastT
    integer(c_int64_t), pointer :: nFun
    class(mlf_ode_fun), pointer :: fun
    integer :: nK
  Contains
    procedure :: setFun => mlf_ode_algo_setFun
    procedure :: reinit => mlf_ode_algo_reinit
    procedure :: initODE => mlf_ode_algo_initODE
    procedure :: denseEvaluation => mlf_ode_algo_denseEvaluation
    procedure(mlf_ode_algo_initHandler), deferred :: initHandler
  End Type mlf_ode_algo

  Type, Public, Abstract, Extends(mlf_experience_model) :: mlf_ode_model
    class(mlf_ode_algo), pointer :: ode
  End Type mlf_ode_model

  Abstract Interface
    Integer Function mlf_ode_algo_initHandler(this, data_handler)
      import :: mlf_ode_algo, mlf_data_handler
      class(mlf_ode_algo), intent(inout), target :: this
      class(mlf_data_handler), intent(inout) :: data_handler
    End Function mlf_ode_algo_initHandler
  End Interface

Contains
  Integer Function mlf_ode_algo_initODE(this, X0, t0, tMax, hMax) Result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    real(c_double), intent(in), optional, target :: X0(:)
    real(c_double), intent(in), optional :: t0, tMax, hMax
    this%lastT = -HUGE(this%lastT)
    If(PRESENT(X0)) this%X0 = X0
    this%X = this%X0
    If(PRESENT(t0)) Then
      this%t0 = t0
    Else
      this%t0 = 0d0
    Endif
    this%t = this%t0
    If(PRESENT(tMax)) Then
      this%tMax = tMax
    Else
      this%tMax = HUGE(this%tMax)
    Endif
    If(PRESENT(hMax)) Then
      this%hMax = hMax
    Else
      this%hMax = HUGE(this%hMax)
    Endif
    info = 0
  End Function mlf_ode_algo_initODE

  ! Generic function for dense evaluation that use only X, X0 and their derivatives
  ! More advanced methods (such as ODE45) may rewrite this function more appropriatly
  Subroutine DefaultUpdateDense(this)
    class(mlf_ode_algo), intent(inout), target :: this
    real(c_double) :: dt
    dt = this%t-this%t0
    If(.NOT. ALLOCATED(this%Cont)) ALLOCATE(this%Cont(SIZE(this%X), 3))
    ASSOCIATE(A => this%Cont, KX0 => this%K(:,1), KX => this%K(:,this%nK))
      A(:,1) = this%X-this%X0
      A(:,2) = dt*KX0(:)-A(:,1)
      A(:,3) = -dt*KX(:)+A(:,1)-A(:,2)
    END ASSOCIATE
    this%lastT = this%t
  End Subroutine DefaultUpdateDense

  Subroutine mlf_ode_algo_denseEvaluation(this, t, Y, K)
    class(mlf_ode_algo), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(out), target :: Y(:)
    real(c_double), intent(out), optional, target :: K(:)
    real(c_double) :: th, th1
    If(this%lastT < this%t) CALL DefaultUpdateDense(this)
    th = (t-this%t0)/(this%lastT-this%t0)
    th1 = th*(1d0-th)
    ASSOCIATE(X0 => this%X0, A => this%Cont)
      Y = X0+th*A(:,1)+th1*A(:,2)+th*th1*A(:,3)
      If(PRESENT(K)) Then
        K = A(:,1)+(1-2*th)*A(:,2)+th*(2-3*th)*A(:,3)
      Endif
    END ASSOCIATE
  End Subroutine mlf_ode_algo_denseEvaluation

  Integer Function mlf_ode_algo_reinit(this) result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%nFun = 0
    this%hMax = HUGE(this%hMax); this%tMax = HUGE(this%tMax)
    this%X0 = 0; this%X = 0; this%t0 = 0; this%t = 0
    this%lastT = -HUGE(this%lastT)
  End Function mlf_ode_algo_reinit

  Subroutine mlf_ode_algo_setFun(this, fun)
    class(mlf_ode_algo), intent(inout), target :: this
    class(mlf_ode_fun), intent(inout), target :: fun
    this%fun => fun
  End Subroutine mlf_ode_algo_setFun

  Integer Function mlf_ode_algo_init(this, numFields, nK, nX, data_handler) &
      Result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer, intent(in), optional :: nK
    integer(c_int64_t), intent(in), optional :: nX
    integer(c_int64_t) :: N, nK2(2)
    CALL numFields%addFields(nRPar = 2, nRsc = 3, nIVar = 1, nRVar = 2)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info < 0) RETURN
    N = -1; nK2 = -1
    If(PRESENT(nX)) N = nX
    info = this%add_rarray(numFields, N, this%X0, C_CHAR_"X0", data_handler = data_handler)
    If(info < 0) RETURN
    info = this%add_rarray(numFields, N, this%X, C_CHAR_"X", data_handler = data_handler)
    If(info < 0) RETURN
    If(PRESENT(nK)) nK2 = [N, INT(nK, KIND = 8)]
    info = this%add_rmatrix(numFields, nK2, this%K, C_CHAR_"K", data_handler = data_handler)
    If(info < 0) RETURN
    this%nK = INT(nK2(2), KIND = 4)
    ! Integer variables
    CALL this%addIVar(numFields, this%nFun, "nFun")
    ! Real parameters
    CALL this%addRPar(numFields, this%hMax, "hMax")
    CALL this%addRPar(numFields, this%tMax, "tMax")
    ! Real variables
    CALL this%addRVar(numFields, this%t, "t")
    CALL this%addRVar(numFields, this%t0, "t0")
  End Function mlf_ode_algo_init

  Integer Function mlf_init_ode_model(this, numFields, ode, data_handler) Result(info)
    class(mlf_ode_model), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_ode_algo), pointer, intent(inout) :: ode
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_obj), pointer :: obj
    info = -1
    If(.NOT. ASSOCIATED(ode)) RETURN
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info<0) RETURN
    obj => ode
    If(PRESENT(data_handler)) Then
      info = ode%initHandler(data_handler%getSubObject(C_CHAR_"ode"))
      If(info < 0) RETURN
    Endif
    CALL this%add_subobject(C_CHAR_"ode", obj)
    this%ode => ode
    ode => NULL()
  End Function mlf_init_ode_model

End Module mlf_ode_class

