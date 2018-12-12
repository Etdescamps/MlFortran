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

Module test_hybrid_kmc_model
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

  Type, Public, Extends(mlf_model_real_parameters) :: test_hybrid_kmc_parameters
    real(c_double) :: Alpha, Beta, Kappa, Zeta, Volume
  Contains
    procedure :: set => test_hybrid_kmc_parameters_set
    procedure :: getNParameters => test_hybrid_kmc_parameters_getNParameters
  End Type test_hybrid_kmc_parameters

  Type, Public, Extends(mlf_model_experiment) :: test_hybrid_kmc_experiment
    real(c_double) :: a, b, c, d ! Element concentration
  End Type test_hybrid_kmc_experiment

  ! Reaction -> rare elements A and B (represented by discrete variables)
  ! more common element c and d (represented by contious variables)
  ! A -> B (alpha)
  ! c catalyze B -> A (beta)
  ! A catalyse c -> d (kappa)
  ! B catalyse d -> c (zeta)
  Type, Public, Extends(mlf_hybrid_kmc_model) :: model_hybrid_kmc
    integer(c_int64_t), pointer :: NIndiv(:) ! [A, B]
    real(c_double), pointer :: Alpha, Beta, Kappa, Zeta, Volume
  Contains
    procedure :: init => test_hybrid_init
    procedure :: setupModel => test_setupModel
    procedure :: applyAction => test_applyAction
    procedure :: funTransitionRates => test_funTransitionRates
    procedure :: evalOde => test_evalOde
  End Type model_hybrid_kmc
Contains
  Integer(8) Function test_hybrid_kmc_parameters_getNParameters(this) Result(N)
    class(test_hybrid_kmc_parameters), intent(inout), target :: this
    N = 4
  End Function test_hybrid_kmc_parameters_getNParameters

  Integer Function test_hybrid_kmc_parameters_set(this, X) Result(info)
    class(test_hybrid_kmc_parameters), intent(inout), target :: this
    real(c_double), intent(in) :: X(:)
    this%Alpha = X(1)
    this%Beta  = X(2)
    this%Kappa = X(3)
    this%Zeta  = X(4)
    info = 0
  End Function test_hybrid_kmc_parameters_set

  Integer Function test_setupModel(this, param, experiment) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    class(mlf_model_parameters), intent(in) :: param
    class(mlf_model_experiment), intent(in) :: experiment
    real(c_double) :: X0(2)
    info = -1
    Select Type(param)
    Class is (test_hybrid_kmc_parameters)
      this%Volume = param%Volume
      this%Alpha  = param%Alpha
      this%Beta   = param%Beta
      this%Kappa  = param%Kappa
      this%Zeta   = param%Zeta
    Class Default
      RETURN
    End Select
    Select Type(experiment)
    Class is (test_hybrid_kmc_experiment)
      this%NIndiv(1) = INT(experiment%a*this%Volume, KIND=8)
      this%NIndiv(2) = INT(experiment%b*this%Volume, KIND=8)
      X0(1) = experiment%c
      X0(2) = experiment%d
      info = this%initModel(X0)
    Class Default
      RETURN
    End Select
  End Function test_setupModel

  Integer Function test_setParameters(this) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this

  End Function test_setParameters

  Integer Function test_hybrid_init(this, data_handler) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    type(mlf_step_numFields) :: numFields
    class(mlf_ode_algo), pointer :: ode
    integer(c_int64_t) :: NCat
    CALL numFields%initFields(nRPar = 5, nRsc = 1)
    If(info < 0) RETURN
    ALLOCATE(mlf_ode45_obj :: ode)
    If(.NOT. PRESENT(data_handler)) Then
      Select Type(ode)
      Class is (mlf_ode45_obj)
        info = ode%init(3_8, atoli = 1d-6, rtoli = 1d-6)
        If(info < 0) RETURN
      End Select
    Endif
    info = mlf_hybrid_kmc_init(this, numFields, ode, NActions = 2, &
      data_handler = data_handler)
    NCat = 2
    info = this%add_i64array(numFields, NCat, this%NIndiv, C_CHAR_"NIndiv", &
      data_handler = data_handler)
    If(info < 0) RETURN
    CALL this%addRPar(numFields, this%Alpha, "Alpha")
    CALL this%addRPar(numFields, this%Beta, "Beta")
    CALL this%addRPar(numFields, this%Kappa, "Kappa")
    CALL this%addRPar(numFields, this%Zeta, "Zeta")
    CALL this%addRPar(numFields, this%Volume, "Volume")
  End Function test_hybrid_init

  Integer Function test_funTransitionRates(this, t, X, F, Rates) Result(N)
    class(model_hybrid_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    real(c_double), intent(out), target :: Rates(:) ! [A->B, B->A]
    ASSOCIATE(A => REAL(this%NIndiv(1)), B => REAL(this%NIndiv(2)), c => X(1), &
        d => X(2), alpha => this%Alpha, beta => this%Beta)
      Rates(1) = alpha*A*d
      Rates(2) = beta*B*c
    END ASSOCIATE
    N = 2
  End Function test_funTransitionRates

  Integer Function test_evalOde(this, t, X, F) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:) ! [c, d]
    real(c_double), intent(out), target :: F(:)
    real(c_double) :: Val
    ASSOCIATE(A => REAL(this%NIndiv(1)), B => REAL(this%NIndiv(2)), c => X(1), &
        d => X(2), kappa => this%Kappa, zeta => this%Zeta, V => this%Volume)
      Val = (zeta*d*B - kappa*c*A)/V
      F(1) = Val
      F(2) = -Val
    END ASSOCIATE
    info = 0
  End Function test_evalOde

  Integer Function test_applyAction(this, id, t, X, F, Rate) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    integer, intent(in) :: id
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: Rate
    real(c_double), intent(inout), target :: X(:), F(:)
    SELECT CASE(id)
    Case(1)
      this%NIndiv = this%NIndiv + [-1, +1]
    Case(2)
      this%NIndiv = this%NIndiv + [+1, -1]
    END SELECT
    info = 0
  End Function test_applyAction
End Module test_hybrid_kmc_model

