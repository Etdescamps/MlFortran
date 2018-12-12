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

Module test_cstr_kmc_model
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_hybrid_kmc
  Use mlf_supervised_model
  Use mlf_ode45
  Use mlf_ode_class
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_model_real_parameters) :: test_cstr_kmc_parameters
    real(c_double) :: Alpha, Beta, Delta, Kappa, Zeta, Volume
  Contains
    procedure :: set => test_cstr_kmc_parameters_set
    procedure :: getNParameters => test_cstr_kmc_parameters_getNParameters
  End Type test_cstr_kmc_parameters

  Type, Public, Extends(mlf_model_experiment) :: test_cstr_kmc_experiment
    real(c_double) :: a, b, c, d ! Element concentration
    real(c_double) :: T
  End Type test_cstr_kmc_experiment


  ! Values: X -> [c, d, T]
  ! T (with constaint T > 0) -> when T == 0 create A + increase T (beta/Volume)
  ! c, d
  ! c catalyse A -> B (alpha)
  ! reaction c -> d (kappa) decrease T (-delta)
  ! A catalyse d -> c (zeta) increase T (delta)
  Type, Public, Extends(mlf_hybrid_kmc_cstrModel) :: model_cstr_kmc
    integer(c_int64_t), pointer :: NIndiv(:) ! [A, B]
    real(c_double), pointer :: Alpha, Beta, Delta, Kappa, Zeta, Volume
  Contains
    procedure :: init => test_hybrid_init
    procedure :: setExperiment => test_setExperiment
    procedure :: setParameters => test_setParameters
    procedure :: applyAction => test_applyAction
    procedure :: funTransitionRates => test_funTransitionRates
    procedure :: evalOde => test_evalOde
    procedure :: m_getHMax => test_getHMax
    procedure :: m_updateCstr => test_updateCstr
    procedure :: m_reachCstr => test_reachCstr
    procedure :: m_getDerivatives => test_getDerivatives
  End Type model_cstr_kmc
Contains
  Integer(8) Function test_cstr_kmc_parameters_getNParameters(this) Result(N)
    class(test_cstr_kmc_parameters), intent(inout), target :: this
    N = 5
  End Function test_cstr_kmc_parameters_getNParameters

  Integer Function test_cstr_kmc_parameters_set(this, X) Result(info)
    class(test_cstr_kmc_parameters), intent(inout), target :: this
    real(c_double), intent(in) :: X(:)
    this%Alpha = X(1)
    this%Beta  = X(2)
    this%Delta = X(3)
    this%Kappa = X(4)
    this%Zeta  = X(5)
    info = 0
  End Function test_cstr_kmc_parameters_set

  Integer Function test_setExperiment(this, experiment) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    class(mlf_model_experiment), intent(in) :: experiment
    real(c_double) :: X0(3)
    Select Type(experiment)
    Class is (test_cstr_kmc_experiment)
      this%NIndiv(1) = INT(experiment%a*this%Volume, KIND=8)
      this%NIndiv(2) = INT(experiment%b*this%Volume, KIND=8)
      X0(1) = experiment%c
      X0(2) = experiment%d
      X0(3) = experiment%T
      info = this%initModel(X0)
    Class Default
      info = -1
    End Select
  End Function test_setExperiment

  Integer Function test_setParameters(this, param) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    class(mlf_model_parameters), intent(in) :: param
    Select Type(param)
    Class is (test_cstr_kmc_parameters)
      this%Volume = param%Volume
      this%Alpha  = param%Alpha
      this%Beta   = param%Beta
      this%Delta  = param%Delta
      this%Kappa  = param%Kappa
      this%Zeta   = param%Zeta
      info = 0
    Class Default
      info = -1
    End Select
  End Function test_setParameters


  Integer Function test_hybrid_init(this, data_handler) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    type(mlf_step_numFields) :: numFields
    class(mlf_ode_algo), pointer :: ode
    integer(c_int64_t) :: NCat
    CALL numFields%initFields(nRPar = 6, nRsc = 1)
    If(info < 0) RETURN
    ALLOCATE(mlf_ode45_obj :: ode)
    If(.NOT. PRESENT(data_handler)) Then
      Select Type(ode)
      Class is (mlf_ode45_obj)
        info = ode%init(4_8, atoli = 1d-6, rtoli = 1d-6)
        If(info < 0) RETURN
      End Select
    Endif
    info = mlf_hybrid_kmc_h_init(this, numFields, ode, &
      numCstr = 1, NActions = 1, data_handler = data_handler)
    If(info < 0) RETURN
    NCat = 2
    info = this%add_i64array(numFields, NCat, this%NIndiv, C_CHAR_"NIndiv", &
      data_handler = data_handler)
    If(info < 0) RETURN
    CALL this%addRPar(numFields, this%Alpha, "Alpha")
    CALL this%addRPar(numFields, this%Beta, "Beta")
    CALL this%addRPar(numFields, this%Delta, "Delta")
    CALL this%addRPar(numFields, this%Kappa, "Kappa")
    CALL this%addRPar(numFields, this%Zeta, "Zeta")
    CALL this%addRPar(numFields, this%Volume, "Volume")
  End Function test_hybrid_init

  Integer Function test_funTransitionRates(this, t, X, F, Rates) Result(N)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    real(c_double), intent(out), target :: Rates(:) ! [A->B, B->A]
    ASSOCIATE(A => REAL(this%NIndiv(1)), d => X(2), zeta => this%Zeta)
      Rates(1) = zeta*A*d
    END ASSOCIATE
    N = 1
  End Function test_funTransitionRates

  Integer Function test_evalOde(this, t, X, F) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:) ! [c, d]
    real(c_double), intent(out), target :: F(:)
    real(c_double) :: Val
    ASSOCIATE(A => REAL(this%NIndiv(1)), B => REAL(this%NIndiv(2)), c => X(1), &
        d => X(2), kappa => this%Kappa, zeta => this%Zeta, V => this%Volume, delta => this%Delta)
      Val = zeta*d*B/V - kappa*c
      F(1) = Val
      F(2) = -Val
      F(3) = Val*delta
    END ASSOCIATE
    info = 0
  End Function test_evalOde

  Integer Function test_applyAction(this, id, t, X, F, Rate) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    integer, intent(in) :: id
    real(c_double), intent(in) :: Rate
    real(c_double), intent(inout) :: t
    real(c_double), intent(inout), target :: X(:), F(:)
    this%NIndiv = this%NIndiv + [-1, +1]
    info = 0
  End Function test_applyAction

  Real(c_double) Function test_getHMax(this, t, X, F) Result(hMax)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    hMax = HUGE(hMax)
    If(X(3) > 0 .AND. F(3) > 0) hMax = -this%kmc_alpha*X(3)/F(3)
  End Function test_getHMax

  Integer Function test_reachCstr(this, t, tMin, tMax, ids, X, F) Result(info)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    integer, intent(in) :: ids(:)
    real(c_double), intent(inout), target :: X(:), F(:)
    real(c_double) :: dt
    dt = -X(3)/F(3)
    dt = MIN(MAX(dt, tMin-t), tMax-t)
    t = t + dt
    X = X + dt*F
    X(3) = this%Beta/this%Volume
    this%NIndiv(1) = this%NIndiv(1)+1
    info = mlf_ODE_Continue
  End Function test_reachCstr

  Integer Function test_updateCstr(this, t, X0, X, F0, F, ids, hMax) Result(N)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    If(X(3) < 0) Then
      ids(1) = 1
      N = 1
    Else
      N = 0
    Endif
  End Function test_updateCstr

  Subroutine test_getDerivatives(this, ids, X0, X, K, C0, C, Q)
    class(model_cstr_kmc), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(in), target :: ids(:)
    C0(1) = X0(3)
    C(1) = X(3)
    Q(1,:) = K(3,:)
  End Subroutine test_getDerivatives
End Module test_cstr_kmc_model

