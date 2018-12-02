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

  Public :: mlf_hybrid_kmc_init, mlf_hybrid_kmc_h_init

  Type, Public, Abstract, Extends(mlf_ode_funCstr) :: mlf_hybrid_odeFun
  Contains
    procedure(mlf_hybrid_odeFun_setModel), deferred :: setModel 
  End Type mlf_hybrid_odeFun

  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_hybrid_kmc_model
    class(mlf_ode45_obj), pointer :: ode
    real(c_double), pointer :: Rates(:)
    real(c_double) :: kmc_alpha, lastTNext
  Contains
    procedure(mlf_hybrid_kmc_apply_action), deferred :: applyAction
    procedure(mlf_hybrid_kmc_transition_rates), deferred :: funTransitionRates
    procedure(mlf_kmc_evalOde), deferred :: evalOde
    procedure :: stepF => mlf_hybrid_kmc_stepFun
  End Type mlf_hybrid_kmc_model

  Type, Public, Abstract, Extends(mlf_hybrid_kmc_model) :: mlf_hybrid_kmc_cstrModel
  Contains
    procedure(mlf_hybrid_kmc_getHMax), deferred :: m_getHMax
    procedure(mlf_hybrid_kmc_updateCstr), deferred :: m_updateCstr
    procedure(mlf_hybrid_kmc_reachCstr), deferred :: m_reachCstr
    procedure(mlf_hybrid_kmc_getDerivatives), deferred :: m_getDerivatives
  End Type mlf_hybrid_kmc_cstrModel

  Type, Public, Extends(mlf_hybrid_odeFun) :: mlf_kmc_odeModel
    class(mlf_hybrid_kmc_model), pointer :: kmc_model
  Contains
    procedure :: setModel => mlf_kmc_ode_setModel
    procedure :: eval => mlf_kmc_eval
    procedure :: getHMax => mlf_kmc_getHMax
    procedure :: updateCstr => mlf_kmc_update
    procedure :: reachCstr => mlf_kmc_reach
    procedure :: getDerivatives => mlf_kmc_getDerivatives
  End Type mlf_kmc_odeModel

  Type, Public, Extends(mlf_hybrid_odeFun) :: mlf_kmc_constrModel
    class(mlf_hybrid_kmc_cstrModel), pointer :: kmc_model
  Contains
    procedure :: setModel => mlf_kmc_h_ode_setModel
    procedure :: eval => mlf_kmc_h_eval
    procedure :: getHMax => mlf_kmc_h_getHMax
    procedure :: updateCstr => mlf_kmc_h_update
    procedure :: reachCstr => mlf_kmc_h_reach
    procedure :: getDerivatives => mlf_kmc_h_getDerivatives
  End Type mlf_kmc_constrModel

  Abstract Interface
    Integer Function mlf_kmc_evalOde(this, t, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Function mlf_kmc_evalOde

    Integer Function mlf_hybrid_kmc_transition_rates(this, t, X, F, Rates)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      real(c_double), intent(out), target :: Rates(:)
    End Function mlf_hybrid_kmc_transition_rates

    Integer Function mlf_hybrid_kmc_apply_action(this, id, t, X, F, Rate)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_kmc_model), intent(inout), target :: this
      real(c_double), intent(inout) :: t
      real(c_double), intent(in) :: Rate
      real(c_double), intent(inout), target :: X(:), F(:)
      integer, intent(in) :: id
    End Function mlf_hybrid_kmc_apply_action

    Integer Function mlf_kmc_ode_init(this, kmc_model) Result(info)
      import :: mlf_hybrid_odeFun
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_odeFun), intent(inout), target :: this
      class(mlf_hybrid_kmc_model), intent(in), target :: kmc_model
    End Function mlf_kmc_ode_init

    Function mlf_hybrid_kmc_getHMax(this, t, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_cstrModel
      class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      real(c_double) :: mlf_hybrid_kmc_getHMax
    End Function mlf_hybrid_kmc_getHMax

    Integer Function mlf_hybrid_kmc_reachCstr(this, t, tMin, tMax, id, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_cstrModel
      class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
      real(c_double), intent(inout) :: t
      real(c_double), intent(in) :: tMin, tMax
      integer, intent(in) :: id
      real(c_double), intent(inout), target :: X(:), F(:)
    End Function mlf_hybrid_kmc_reachCstr

    Integer Function mlf_hybrid_kmc_updateCstr(this, t, X0, X, F0, F, ids, hMax)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_cstrModel
      class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
      integer, intent(out), target :: ids(:)
      real(c_double), intent(inout) :: hMax
    End Function mlf_hybrid_kmc_updateCstr

    Subroutine mlf_hybrid_kmc_getDerivatives(this, ids, X0, X, K, C0, C, Q)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_cstrModel
      class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
      real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
      real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
      integer, intent(in), target :: ids(:)
    End Subroutine mlf_hybrid_kmc_getDerivatives

    Integer Function mlf_hybrid_odeFun_setModel(this, kmc_model)
      import :: mlf_hybrid_odeFun
      import :: mlf_hybrid_kmc_model
      class(mlf_hybrid_odeFun), intent(inout), target :: this
      class(mlf_hybrid_kmc_model), intent(in), target :: kmc_model
    End Function mlf_hybrid_odeFun_setModel
  End Interface
Contains
  Integer Function mlf_kmc_ode_setModel(this, kmc_model) Result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    class(mlf_hybrid_kmc_model), intent(in), target :: kmc_model
    this%NCstr = 1
    this%kmc_model => kmc_model
    info = 0
  End Function mlf_kmc_ode_setModel

  Integer Function mlf_kmc_h_ode_setModel(this, kmc_model) Result(info)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    class(mlf_hybrid_kmc_model), intent(in), target :: kmc_model
    Select Type(kmc_model)
    Class is (mlf_hybrid_kmc_cstrModel)
      this%kmc_model => kmc_model
      info = 0
    Class default
      info = -1
    End Select
  End Function mlf_kmc_h_ode_setModel

  Integer Function mlf_hybrid_kmc_h_init(this, numFields, fun, numCstr, nActions, X0, &
      t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax, nStiff, data_handler) Result(info)
    class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_kmc_constrModel), intent(inout), target, optional :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:), t0, tMax, atoli, rtoli
    real(c_double), intent(in), optional :: fac, facMin, facMax, hMax
    integer(c_int64_t), intent(in), optional :: nStiff
    integer, intent(in), optional :: NActions, numCstr
    class(mlf_kmc_constrModel), pointer :: funSelected
    class(mlf_obj), pointer :: obj
    If(PRESENT(fun)) Then
      funSelected => fun
    Else
      ALLOCATE(mlf_kmc_constrModel :: funSelected)
      obj => funSelected
      funSelected%NCstr = 1 + numCstr
      CALL this%add_subobject(C_CHAR_"odeFun", obj)
    Endif
    info = mlf_hybrid_kmc_init(this, numFields, funSelected, nActions, X0, &
      t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax, nStiff, data_handler)
  End Function mlf_hybrid_kmc_h_init

  Integer Function mlf_hybrid_kmc_init(this, numFields, fun, nActions, X0, t0, tMax, &
      atoli, rtoli, fac, facMin, facMax, hMax, nStiff, data_handler) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_hybrid_odeFun), intent(inout), target, optional :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:), t0, tMax, atoli, rtoli
    real(c_double), intent(in), optional :: fac, facMin, facMax, hMax
    integer(c_int64_t), intent(in), optional :: nStiff
    integer, intent(in), optional :: NActions
    type(mlf_ode45_obj), pointer :: ode
    class(mlf_hybrid_odeFun), pointer :: funSelected
    class(mlf_obj), pointer :: obj
    real(c_double), allocatable :: X(:)
    real(c_double) :: r
    integer(8) :: N
    CALL numFields%addFields(nRsc = 1)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info<0) RETURN
    ! Allocate the ODE solver and link it to the object
    ALLOCATE(ode)
    obj => ode
    CALL this%add_subobject(C_CHAR_"ode", obj)
    this%ode => ode
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
      X(1) = -LOG(1d0-r)
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
    If(info < 0) GOTO 10
    If(PRESENT(NActions)) N = NActions
    info = this%add_rarray(numFields, N, this%Rates, C_CHAR_"Rates", &
      data_handler = data_handler)
    If(info < 0) GOTO 10
    If(.NOT. PRESENT(data_handler)) Then
      info = this%reinit()
      If(info < 0) GOTO 10
    Endif
    info = funSelected%setModel(this)
    If(info < 0) GOTO 10
    this%kmc_alpha = 1.5
    this%lastTNext = ieee_value(this%lastTNext, ieee_quiet_nan)
    RETURN
    ! Error handling: finalize object
 10 call this%finalize()
  End Function mlf_hybrid_kmc_init

  Integer Function mlf_hybrid_kmc_stepFun(this, nIter) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: nIter
    info = this%ode%stepF(nIter)
  End Function mlf_hybrid_kmc_stepFun

  Integer Function EvalOdeModel(model, t, X, F) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: model
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    integer ::N
    ASSOCIATE(Rates => model%Rates)
      info = model%evalOde(t, X(2:), F(2:))
      If(info < 0) RETURN
      N = model%funTransitionRates(t, X(2:), F(2:), Rates)
      If(N <= 0) Then
        info = N
        If(info == 0) info = mlf_ODE_StopTime
        RETURN
      Endif
      F(1) = -SUM(Rates(1:N))
      info = mlf_ODE_Continue
    END ASSOCIATE
  End Function 

  Integer Function mlf_kmc_eval(this, t, X, F) Result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    info = EvalOdeModel(this%kmc_model, t, X, F)
  End Function mlf_kmc_eval

  Integer Function mlf_kmc_h_eval(this, t, X, F) Result(info)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    info = EvalOdeModel(this%kmc_model, t, X, F)
  End Function mlf_kmc_h_eval

  Real(c_double) Function ModelGetHMax(this, t, X, F) Result(hMax)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    real(c_double) :: Z
    hMax = HUGE(hMax)
    If(F(1) >= 0) RETURN
    If(.NOT. ieee_is_nan(this%lastTNext) .AND. X(1) > 0) Then
      If(this%lastTNext == t) this%kmc_alpha = this%kmc_alpha * 1.5
    Endif
    Z = -this%kmc_alpha*X(1)/F(1)
    If(Z <= 0) Then
      this%lastTNext = ieee_value(this%lastTNext, ieee_quiet_nan)
    Else
      hMax = Z
      this%lastTNext = t + Z
    Endif
  End Function ModelGetHMax

  Real(c_double) Function mlf_kmc_getHMax(this, t, X, F) Result(hMax)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    hMax = ModelGetHMax(this%kmc_model, t, X, F)
  End Function mlf_kmc_getHMax

  Real(c_double) Function mlf_kmc_h_getHMax(this, t, X, F) Result(hMax)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    hMax = ModelGetHMax(this%kmc_model, t, X, F)
    hMax = MIN(this%kmc_model%m_getHMax(t, X(2:), F(2:)), hMax)
  End Function mlf_kmc_h_getHMax

  Subroutine mlf_kmc_getDerivatives(this, ids, X0, X, K, C0, C, Q)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(in), target :: ids(:)
    C0(1) = X0(1)
    C(1) = X(1)
    Q(1,:) = K(1,:)
  End Subroutine mlf_kmc_getDerivatives

  Subroutine mlf_kmc_h_getDerivatives(this, ids, X0, X, K, C0, C, Q)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(in), target :: ids(:)
    If(ids(1) == 1) Then
      C0(1) = X0(1)
      C(1) = X(1)
      Q(1,:) = K(1,:)
      If(SIZE(ids) > 1) CALL this%kmc_model%m_getDerivatives(ids(2:)-1, &
        X0(2:), X(2:), K(2:,:), C0(2:), C(2:), Q(2:,:))
    Else
      CALL this%kmc_model%m_getDerivatives(ids-1, X0(2:), X(2:), K(2:,:), &
        C0, C, Q)
    Endif
  End Subroutine mlf_kmc_h_getDerivatives

  Integer Function KmcUpdate(X, ids) Result(N)
    real(c_double), intent(in) :: X
    integer, intent(inout) :: ids(:)
    If(X < 0 .OR. X == 0) Then
      N = 1
      ids(1) = 1
      RETURN
    Endif
    N = 0
  End Function KmcUpdate

  Integer Function mlf_kmc_update(this, t, X0, X, F0, F, ids, hMax) &
      Result(N)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    integer :: k
    N = KmcUpdate(X(1), ids)
    hMax = MIN(hMax, ModelGetHMax(this%kmc_model, t, X, F))
  End Function  mlf_kmc_update

  Integer Function mlf_kmc_h_update(this, t, X0, X, F0, F, ids, hMax) &
      Result(N)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    integer :: K
    K = KmcUpdate(X(1), ids)
    N = K+this%kmc_model%m_updateCstr(t, X0(2:), X(2:), F0(2:), F(2:), ids(K+1:), hMax)
    If(N > K) ids(K+1:N) = ids(K+1:N)+1
    If(N == 0) hMax = MIN(hMax, ModelGetHMax(this%kmc_model, t, X, F))
  End Function  mlf_kmc_h_update

  Integer Function KMCReachAction(model, t, tMin, tMax, X, F, Rates) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: model
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    real(c_double), intent(out), target :: Rates(:)
    integer :: idAction, N
    real(c_double) :: r, dt
    info = -1
    dt = -X(1)/F(1)
    dt = MIN(MAX(dt, tMin-t), tMax-t)
    t = t + dt
    X = X + dt*F
    model%lastTNext = ieee_value(model%lastTNext, ieee_quiet_nan)
    N = Model%funTransitionRates(t, X(2:), F(2:), Rates)
    If(N <= 0) RETURN ! Shall not happen
    CALL RANDOM_NUMBER(r)
    r = r*SUM(Rates(1:N))
    Do idAction = 1,N-1
      r = r - Rates(idAction)
      If(r <= 0) Then
        info = Model%applyAction(idAction, t, X(2:), F(2:), Rates(idAction))
        EXIT
      Endif
    End Do
    If(r > 0) info = Model%applyAction(N, t, X(2:), F(2:), Rates(N))
    If(info < 0 .OR. info == mlf_ODE_HardCstr .OR. info == mlf_ODE_StopTime) RETURN
    CALL RANDOM_NUMBER(r)
    X(1) = -LOG(1d0-r)
    info = EvalOdeModel(model, t, X, F)
  End Function KMCReachAction

  Integer Function mlf_kmc_reach(this, t, tMin, tMax, id, X, F) Result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: id
    info = -1
    ! id should be equal to 1
    If(id /= 1) RETURN
    info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
  End Function mlf_kmc_reach

  Integer Function mlf_kmc_h_reach(this, t, tMin, tMax, id, X, F) Result(info)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: id
    real(c_double) :: t0
    If(id == 1) Then
      info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
    Else
      t0 = t
      info = this%kmc_model%m_reachCstr(t, tMin, tMax, id-1, X(2:), F(2:))
      If(t0 /= t) X(1) = X(1) + (t0-t)*F(1)
    Endif
  End Function mlf_kmc_h_reach
End Module mlf_hybrid_kmc

