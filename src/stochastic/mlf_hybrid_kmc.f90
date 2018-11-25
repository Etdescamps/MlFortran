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

Module mlf_hybrid_kmc
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_ode45
  Use mlf_errors
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_hybrid_kmc_init

  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_hybrid_kmc_model
    class(mlf_ode45_obj), pointer :: ode
    real(c_double), pointer :: Rates(:)
  Contains
    procedure(mlf_hybrid_kmc_apply_action), deferred :: applyAction
    procedure(mlf_hybrid_kmc_transition_rates), deferred :: funTransitionRates
    procedure(mlf_kmc_evalOde), deferred :: evalOde
    procedure :: stepF => mlf_hybrid_kmc_stepFun
  End Type mlf_hybrid_kmc_model

  Type, Public, Extends(mlf_ode_funCstr) :: mlf_kmc_odeModel
    class(mlf_hybrid_kmc_model), pointer :: kmc_model
    real(c_double) :: kmc_alpha
  Contains
    procedure :: eval => mlf_kmc_eval
    procedure :: getHMax => mlf_kmc_getHMax
    procedure :: updateCstr => mlf_kmc_update
    procedure :: reachCstr => mlf_kmc_reach
    procedure :: getDerivatives => mlf_kmc_getDerivatives
  End Type mlf_kmc_odeModel

!  Type, Public, Abstract, Extends(mlf_kmc_odeModel) :: mlf_kmc_odeModelCstr
!    integer :: idCstr
!  Contains
!    procedure :: modelCstr => mlf_kmc_modelCstr
!    procedure(mlf_kmc_cstrDerivatives), deferred :: cstrDerivatives
!    procedure :: updateCstr => mlf_kmc_cstr_update
!    procedure :: reachCstr => mlf_kmc_cstr_reach
!    procedure :: getDerivatives => mlf_kmc_cstr_getDerivatives
!  End Type mlf_kmc_odeModelCstr

  Abstract Interface
    Integer Function mlf_kmc_evalOde(this, t, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Function mlf_kmc_evalOde

!    Subroutine mlf_kmc_cstrDerivatives(this, ids, K, C0, C, Q)
!      class(mlf_kmc_odeModelCstr), intent(inout), target :: this
!      real(c_double), intent(in), target :: K(:,:)
!      real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
!      integer, intent(inout), target :: ids(:)
!    End Subroutine mlf_kmc_cstrDerivatives

    Integer Function mlf_hybrid_kmc_transition_rates(this, t, X, F, Rates)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      real(c_double), intent(out), target :: Rates(:)
    End Function mlf_hybrid_kmc_transition_rates

    Integer Function mlf_hybrid_kmc_apply_action(this, id, t, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      integer, intent(in) :: id
    End Function mlf_hybrid_kmc_apply_action
  End Interface
Contains
  Integer Function mlf_hybrid_kmc_init(this, numFields, fun, nActions, X0, t0, tMax, &
      atoli, rtoli, fac, facMin, facMax, hMax, nStiff, data_handler) result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_kmc_odeModel), intent(inout), target, optional :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:), t0, tMax, atoli, rtoli
    real(c_double), intent(in), optional :: fac, facMin, facMax, hMax
    integer(c_int64_t), intent(in), optional :: nStiff
    integer, intent(in), optional :: NActions
    type(mlf_ode45_obj), pointer :: ode
    class(mlf_kmc_odeModel), pointer :: funSelected
    class(mlf_obj), pointer :: obj
    real(c_double), allocatable :: X(:)
    real(c_double) :: r
    integer(8) :: N
    CALL numFields%addFields(nRsc = 1)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info<0) RETURN
    ALLOCATE(ode)
    If(PRESENT(fun)) Then
      funSelected => fun
    Else
      ALLOCATE(mlf_kmc_odeModel :: funSelected)
      obj => funSelected
      CALL this%add_subobject(C_CHAR_"odeFun", obj)
    Endif
    If(PRESENT(X0)) Then
      ALLOCATE(X(SIZE(X0)+1))
      X(2:) = X0
      CALL RANDOM_NUMBER(r)
      X(1) = -log(1d0-r)
      If(PRESENT(data_handler)) Then
        info = ode%init(funSelected, X, t0, tMax, atoli, rtoli, fac, facMin, &
          facMax, hMax, nStiff, data_handler%getSubObject(C_CHAR_"ode"))
      Else
        info = ode%init(funSelected, X, t0, tMax, atoli, rtoli, fac, facMin, &
          facMax, hMax, nStiff)
      Endif
    Else
      If(PRESENT(data_handler)) Then
        info = ode%init(funSelected, X0, t0, tMax, atoli, rtoli, fac, facMin, &
          facMax, hMax, nStiff, data_handler%getSubObject(C_CHAR_"ode"))
      Else
        info = -1
        WRITE (error_unit, *) "mlf_hybrid_kmc_init: No X0 nor handler provided"
        RETURN
      Endif
    Endif
    If(PRESENT(NActions)) N = NActions
    info = this%add_rarray(numFields, N, this%Rates, C_CHAR_"Rates", &
      data_handler = data_handler)
    If(info < 0) RETURN

    If(info < 0) Then
      DEALLOCATE(ode)
      RETURN
    Endif
    obj => ode
    CALL this%add_subobject(C_CHAR_"ode", obj)
    this%ode => ode
    funSelected%kmc_model => this
  End Function mlf_hybrid_kmc_init

  Integer Function mlf_hybrid_kmc_stepFun(this, nIter) result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: nIter
    info = this%ode%stepF(nIter)
  End Function mlf_hybrid_kmc_stepFun

  Integer Function mlf_kmc_eval(this, t, X, F) Result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    ASSOCIATE(Model => this%kmc_model, Rates => this%kmc_model%Rates)
      info = Model%evalOde(t, X(2:), F(2:))
      If(info < 0) RETURN
      info = Model%funTransitionRates(t, X(2:), F(2:), Rates)
      If(info < 0) RETURN
      F(1) = -SUM(Rates)
    END ASSOCIATE
  End Function mlf_kmc_eval

  Real(c_double) Function mlf_kmc_getHMax(this, t, X, F) Result(hMax)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    hMax = -this%kmc_alpha*X(1)/F(1)
  End Function mlf_kmc_getHMax


  Subroutine mlf_kmc_getDerivatives(this, ids, X0, X, K, C0, C, Q)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(inout), target :: ids(:)
    C0(1) = X0(1)
    C(1) = X(1)
    Q(1,:) = K(1,:)
  End Subroutine mlf_kmc_getDerivatives

  Integer Function mlf_kmc_update(this, t, X0, X, F0, F, ids, hMax) &
      result(N)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    integer :: k
    If(X(1) < 0 .OR. X(1) == 0) Then
      N = 1
      ids(1) = 1
      RETURN
    Endif
    N = 0
    hMax = MIN(hMax, mlf_kmc_getHMax(this, t, X, F))
  End Function  mlf_kmc_update

  Integer Function mlf_kmc_reach(this, t, id, X, F) result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(inout), target :: t, X(:), F(:)
    integer, intent(in) :: id
    integer :: idAction, N
    real(c_double) :: r
    ASSOCIATE(Model => this%kmc_model, Rates => this%kmc_model%Rates)
    ! id should be equal to 1
      N = Model%funTransitionRates(t, X(2:), F(2:), Rates)
      If(N <= 0) Then
        info = N
        If(N == 0) info = mlf_ODE_StopTime
        RETURN
      Endif
      CALL RANDOM_NUMBER(r)
      r = r*SUM(this%kmc_model%Rates(1:N))
      Do idAction = 1,N-1
        r = r - Rates(idAction)
        If(r <= 0) Then
          info = Model%applyAction(idAction, t, X, F)
          EXIT
        Endif
      End Do
      If(r > 0) info = Model%applyAction(N, t, X, F)
      If(info < 0 .OR. info == mlf_ODE_HardCstr .OR. info == mlf_ODE_StopTime) RETURN
      info = this%eval(t, X, F)
      CALL RANDOM_NUMBER(r)
      X(1) = -log(1d0-r)
    END ASSOCIATE
  End Function mlf_kmc_reach

!  Subroutine mlf_kmc_cstr_getDerivatives(this, ids, K, C0, C, Q)
!    class(mlf_kmc_odeModelCstr), intent(inout), target :: this
!    real(c_double), intent(in), target :: K(:,:)
!    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
!    integer, intent(inout), target :: ids(:)
!    integer :: k
!    k = size(ids)
!    If(ids(k) /= NCstr+1) Then
!      CALL this%cstrDerivatives(this, ids, K, C0, C, Q)
!      RETURN
!    Endif
!    C0(k) = this%kmc_model%U
!    C(k) = this%kmc_model%Utemp
!    Q(k,:) = K(k,:)
!    If(k > 1) Then
!      CALL this%cstrDerivatives(this, ids(:k-1), K(:k-1,:), &
!        C0(:k-1), C(:k-1), Q(:k-1))
!    Endif
!  End Subroutine mlf_kmc_cstr_getDerivatives
!
!  Real(c_double) Function mlf_kmc_cstr_update(this, t, X, F, ids, N) &
!     result(hMax)
!   class(mlf_kmc_odeModelCstr), intent(inout), target :: this
!    real(c_double), intent(in) :: t
!    real(c_double), intent(in), target :: X(:), F(:)
!    integer, intent(out), optional, target :: ids(:)
!    integer, intent(out), optional :: N
!  End Function  mlf_kmc_cstr_update
!  
!  Integer Function mlf_kmc_cstr_reach(this, t, id, X, F, hMax) result(info)
!    class(mlf_kmc_odeModelCstr), intent(inout), target :: this
!    real(c_double), intent(inout), target :: t, X(:), F(:), hMax
!    integer, intent(in) :: id
!  End Function mlf_kmc_cstr_reach


End Module mlf_hybrid_kmc

