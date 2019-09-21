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
  Use mlf_errors
  Use mlf_supervised_model
  Use mlf_ode_class
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_hybrid_kmc_init, mlf_hybrid_kmc_h_init

  Integer, Parameter, Public :: mlf_h_FunOdeComplete = 5, mlf_h_StopODE = 6

  Type, Public, Abstract, Extends(mlf_ode_funCstr) :: mlf_hybrid_odeFun
  Contains
    procedure(mlf_hybrid_odeFun_setModel), deferred :: setModel 
  End Type mlf_hybrid_odeFun

  Type, Public, Abstract, Extends(mlf_ode_model) :: mlf_hybrid_kmc_model
    real(c_double), pointer :: Rates(:)
    real(c_double) :: kmc_alpha, lastTNext
    logical :: without_ode
  Contains
    procedure(mlf_hybrid_kmc_apply_action), deferred :: applyAction
    procedure(mlf_hybrid_kmc_transition_rates), deferred :: funTransitionRates
    procedure(mlf_kmc_evalOde), deferred :: evalOde
    procedure :: stepF => mlf_hybrid_kmc_stepFun
    procedure :: initModel => mlf_hybrid_kmc_initModel
  End Type mlf_hybrid_kmc_model

  Type, Public, Abstract, Extends(mlf_hybrid_kmc_model) :: mlf_hybrid_kmc_cstrModel
  Contains
    procedure(mlf_hybrid_kmc_getHMax), deferred :: m_getHMax
    procedure(mlf_hybrid_kmc_updateCstr), deferred :: m_updateCstr
    procedure(mlf_hybrid_kmc_reachCstr), deferred :: m_reachCstr
    procedure(mlf_hybrid_kmc_getDerivatives), deferred :: m_getDerivatives
    procedure :: m_checkCstr => mlf_hybrid_kmc_checkCstr
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
    procedure :: checkCstr => mlf_kmc_h_checkCstr
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

    Integer Function mlf_hybrid_kmc_reachCstr(this, t, tMin, tMax, ids, X, F)
      Use iso_c_binding
      import :: mlf_hybrid_kmc_cstrModel
      class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
      real(c_double), intent(inout) :: t
      real(c_double), intent(in) :: tMin, tMax
      integer, intent(in) :: ids(:)
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

  Integer Function mlf_hybrid_kmc_h_init(this, numFields, ode, fun, numCstr, &
      nActions, data_handler) Result(info)
    class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_ode_algo), pointer, intent(inout) :: ode
    class(mlf_kmc_constrModel), intent(inout), target, optional :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer, intent(in), optional :: nActions, numCstr
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
    this%without_ode = .FALSE.
    info = mlf_hybrid_kmc_init(this, numFields, ode, funSelected, nActions, &
      data_handler)
  End Function mlf_hybrid_kmc_h_init

  Integer Function mlf_hybrid_kmc_initModel(this, X0, t0, tMax, hMax) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    real(c_double), intent(in), target :: X0(:)
    real(c_double), intent(in), optional :: t0, tMax, hMax
    real(c_double) :: r
 20 CALL RANDOM_NUMBER(r)
    If(r == 0) GOTO 20
    this%ode%X0(1) = -LOG(r)
    this%ode%X0(2:) = X0
    info = this%ode%initODE(t0 = t0, tMax = tMax, hMax = hMax)
  End Function mlf_hybrid_kmc_initModel

  Integer Function mlf_hybrid_kmc_init(this, numFields, ode, fun, nActions, &
      data_handler) Result(info)
    class(mlf_hybrid_kmc_model), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_ode_algo), pointer, intent(inout) :: ode
    class(mlf_hybrid_odeFun), intent(inout), target, optional :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer, intent(in), optional :: nActions
    class(mlf_hybrid_odeFun), pointer :: funSelected
    class(mlf_obj), pointer :: obj
    integer(8) :: N
    CALL numFields%addFields(nRsc = 1)
    info = mlf_init_ode_model(this, numFields, ode, data_handler)
    If(info<0) RETURN
    If(PRESENT(fun)) Then
      funSelected => fun
    Else
      ALLOCATE(mlf_kmc_odeModel :: funSelected)
      obj => funSelected
      CALL this%add_subobject(C_CHAR_"odeFun", obj)
    Endif
    CALL this%ode%setFun(funSelected)
    If(PRESENT(NActions)) N = NActions
    info = this%add_rarray(numFields, N, this%Rates, C_CHAR_"Rates", &
      data_handler = data_handler)
    If(info < 0) GOTO 10
    If(.NOT. PRESENT(data_handler)) Then
      info = this%reinit()
      If(info < 0) GOTO 10
    Endif
    Select Type(funSelected)
    Class is (mlf_hybrid_odeFun)
      info = funSelected%setModel(this)
    Class Default
      info = -1
    End Select
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
      If(info < 0) Then
        WRITE (error_unit, *) "EvalOde error: ", info
        RETURN
      Endif
      N = model%funTransitionRates(t, X(2:), F(2:), Rates)
      If(N <= 0) Then
        info = N
        If(info == 0) info = mlf_ODE_StopTime
        RETURN
      Endif
      F(1) = -SUM(Rates(:N))
      If(info == mlf_h_FunOdeComplete .AND. F(1) == 0) Then
        info = mlf_ODE_StopTime
      Else
        info = mlf_ODE_Continue
      Endif
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
    model%lastTNext = ieee_value(model%lastTNext, ieee_quiet_nan)
    N = model%funTransitionRates(t, X(2:), F(2:), Rates)
    If(N <= 0) Then
      WRITE (error_unit, *) "Error model%funTransitionRates:", N
      RETURN ! Shall not happen
    Endif
    CALL RANDOM_NUMBER(r)
    r = r*SUM(Rates(:N))
    Do idAction = 1,N-1
      r = r - Rates(idAction)
      If(r <= 0) Then
        info = model%applyAction(idAction, t, X(2:), F(2:), Rates(idAction))
        GOTO 20
      Endif
    End Do
    info = model%applyAction(N, t, X(2:), F(2:), Rates(N))
 20 If(info < 0 .OR. info == mlf_ODE_StopTime) Then
      If(info < 0) WRITE (error_unit, *) "Error model%applyAction:", info
      RETURN
    Endif
 10 CALL RANDOM_NUMBER(r)
    If(r == 0) GOTO 10 ! Remove the case r == 0
    X(1) = -LOG(r)
    If(info /= mlf_ODE_HardCstr) info = EvalOdeModel(model, t, X, F)
  End Function KMCReachAction

  Integer Function mlf_kmc_reach(this, t, tMin, tMax, ids, X, F) Result(info)
    class(mlf_kmc_odeModel), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: ids(:)
    info = -1
    ! id should be equal to 1
    If(ids(1)/= 1) RETURN
    info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
  End Function mlf_kmc_reach

  Integer Function mlf_hybrid_kmc_checkCstr(this, t, X0, Xd, K0, K1, ids, K) Result(N)
    class(mlf_hybrid_kmc_cstrModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), Xd(:), K0(:), K1(:)
    integer, intent(inout), target :: ids(:)
    integer, intent(in) :: K
    N = K
  End Function mlf_hybrid_kmc_checkCstr

  Integer Function mlf_kmc_h_checkCstr(this, t, X0, Xd, K0, K1, ids, K) Result(N)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), Xd(:), K0(:), K1(:)
    integer, intent(inout), target :: ids(:)
    integer, intent(in) :: K
    If(ids(1) == 1) Then
      N = 1 + this%kmc_model%m_checkCstr(t, X0(2:), Xd(2:), K0(2:), K1(2:), ids(2:), K-1)
    Else If(Xd(1) <= X0(1)*EPSILON(1d0)) Then
      ! Add KMC event (id = 1) as constraint
      ids(2:K+1) = ids(1:K)
      ids(1) = 1
      N = 1 + this%kmc_model%m_checkCstr(t, X0(2:), Xd(2:), K0(2:), K1(2:), ids(2:), K)
    Else
      ! No KMC event
      N = this%kmc_model%m_checkCstr(t, X0(2:), Xd(2:), K0(2:), K1(2:), ids, K)
    Endif
  End Function mlf_kmc_h_checkCstr

  Integer Function mlf_kmc_h_reach(this, t, tMin, tMax, ids, X, F) Result(info)
    class(mlf_kmc_constrModel), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: ids(:)
    real(c_double) :: t0, U0, r
    If(SIZE(ids) == 1 .AND. ids(1) == 1) Then
      info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
      If(info < 0) Then
        WRITE (error_unit, *) "KMCReachAction error", info
        WRITE (error_unit, *) "ids:", ids, "X:", X, "F:", F
        RETURN
      Endif
    Else If(ids(1) == 1) Then
      ! The global events of the system are applied before applying random actions.
      t0 = t
      U0 = F(1)
      info = this%kmc_model%m_reachCstr(t, tMin, tMax, ids(2:)-1, X(2:), F(2:))
      If(info < 0) Then
        WRITE (error_unit, *) "KMC model reachCstr error", info
        WRITE (error_unit, *) "ids:", ids, "X:", X, "F:", F
        RETURN
      Endif
      If(info == mlf_ODE_ReevaluateDer) Then
        info = EvalOdeModel(this%kmc_model, t, X, F)
      Endif
      ! Look if the probability of an event has been inpacted by the constraint reach
      ! Introduce a null event if the probability decreased
      If(F(1) < U0) Then
        CALL RANDOM_NUMBER(r)
        If(r*U0 > F(1)) Then
          ! Case where a null event is reached: reinit X(1)
          10 CALL RANDOM_NUMBER(r)
          If(r == 0) GOTO 10 ! Remove the case r == 0
          X(1) = -LOG(r)
          info = EvalOdeModel(this%kmc_model, t, X, F)
          RETURN
        Endif
      Endif
      info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
    Else
      t0 = t
      info = this%kmc_model%m_reachCstr(t, tMin, tMax, ids-1, X(2:), F(2:))
      If(info < 0) Then
        WRITE (error_unit, *) "KMC model reachCstr error", info
        WRITE (error_unit, *) "ids:", ids, "X:", X, "F:", F
        RETURN
      Endif
      If((info == mlf_ODE_Continue .OR. info == mlf_ODE_SoftCstr) .AND. t0 /= t) Then
        X(1) = X(1) + (t-t0)*F(1)
      Endif
      If(X(1) <= 0) Then
        info = KMCReachAction(this%kmc_model, t, tMin, tMax, X, F, this%kmc_model%Rates)
      Endif
      If(info == mlf_ODE_ReevaluateDer) Then
        info = EvalOdeModel(this%kmc_model, t, X, F)
      Endif
    Endif
  End Function mlf_kmc_h_reach
End Module mlf_hybrid_kmc

