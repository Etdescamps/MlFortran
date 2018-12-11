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

Module mlf_ode_model
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_ode_model_reinitT, mlf_ode_model_init, mlf_ode_model_reinit
  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_ode_algo
    real(c_double), pointer :: X0(:), X(:), K(:,:)
    real(c_double), pointer :: t0, t, tMax, hMax
    real(c_double), allocatable :: Cont(:,:)
    real(c_double) :: lastT
    integer(c_int64_t), pointer :: nFun
    class(mlf_ode_fun), pointer :: fun
    integer :: nK
  Contains
    procedure :: reinit => mlf_ode_model_reinit
    procedure :: denseEvaluation => mlf_ode_model_denseEvaluation
  End Type mlf_ode_algo

Contains
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

  Subroutine mlf_ode_model_denseEvaluation(this, t, Y, K)
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
  End Subroutine mlf_ode_model_denseEvaluation

  Integer Function mlf_ode_model_reinit(this) result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%nFun = 0
    this%hMax = HUGE(this%hMax); this%tMax = HUGE(this%tMax)
    this%X0 = 0; this%X = 0; this%t0 = 0; this%t = 0
    this%lastT = -HUGE(this%lastT)
  End Function mlf_ode_model_reinit

  Integer Function mlf_ode_model_reinitT(this, X0, t0, tMax, hMax) Result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    real(c_double), intent(in), optional :: X0(:), t0, hMax, tMax
    If(PRESENT(X0)) Then
      this%X0 = X0
      this%X = X0
    Endif
    If(PRESENT(t0)) Then
      this%t0 = t0
      this%t = t0
    Endif
    If(PRESENT(hMax)) this%hMax = hMax
    If(PRESENT(tMax)) this%tMax = tMax
    info = 0
  End Function mlf_ode_model_reinitT

  Integer Function mlf_ode_model_init(this, numFields, nK, fun, X0, data_handler) &
      Result(info)
    class(mlf_ode_algo), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_ode_fun), intent(inout), target :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:)
    integer, intent(in), optional :: nK
    integer(c_int64_t) :: N, nK2(2)
    CALL numFields%addFields(nRPar = 2, nRsc = 3, nIVar = 1, nRVar = 2)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info < 0) RETURN
    this%fun => fun; N = -1; nK2 = -1
    If(PRESENT(X0)) N = size(X0)
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
  End Function mlf_ode_model_init


End Module mlf_ode_model

