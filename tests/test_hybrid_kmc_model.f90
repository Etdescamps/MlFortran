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
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

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
    procedure :: applyAction => test_applyAction
    procedure :: funTransitionRates => test_funTransitionRates
    procedure :: evalOde => test_evalOde
  End Type model_hybrid_kmc
Contains
  Integer Function test_hybrid_init(this, X0, NIndiv, Volume, &
      Alpha, Beta, Kappa, Zeta, data_handler) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    real(c_double), intent(in), optional :: X0(:), Volume, Alpha, Beta, Kappa, Zeta
    integer(c_int64_t), intent(in), optional :: NIndiv(:)
    class(mlf_data_handler), intent(inout), optional :: data_handler
    type(mlf_step_numFields) :: numFields
    integer(c_int64_t) :: NCat
    CALL numFields%initFields(nRPar = 5, nRsc = 1)
    If(info < 0) RETURN
    info = mlf_hybrid_kmc_init(this, numFields, NActions = 2, X0 = X0, &
      atoli = 1d-6, rtoli = 1d-6, data_handler = data_handler)
    NCat = 2
    info = this%add_i64array(numFields, NCat, this%NIndiv, C_CHAR_"NIndiv", &
      data_handler = data_handler)
    If(info < 0) RETURN
    If(PRESENT(NIndiv)) this%NIndiv = NIndiv
    CALL this%addRPar(numFields, this%Alpha, "Alpha")
    CALL this%addRPar(numFields, this%Beta, "Beta")
    CALL this%addRPar(numFields, this%Kappa, "Kappa")
    CALL this%addRPar(numFields, this%Zeta, "Zeta")
    CALL this%addRPar(numFields, this%Volume, "Volume")
    If(PRESENT(Alpha)) this%Alpha = Alpha
    If(PRESENT(Beta)) this%Beta = Beta
    If(PRESENT(Kappa)) this%Kappa = Kappa
    If(PRESENT(Zeta)) this%Zeta = Zeta
    If(PRESENT(Volume)) this%Volume = Volume
  End Function test_hybrid_init

  Integer Function test_funTransitionRates(this, t, X, F, Rates) Result(N)
    class(model_hybrid_kmc), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    real(c_double), intent(out), target :: Rates(:) ! [A->B, B->A]
    ASSOCIATE(A => REAL(this%NIndiv(1)), B => REAL(this%NIndiv(2)), c => X(1), &
        alpha => this%Alpha, beta => this%Beta)
      Rates(1) = alpha*A
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

  Integer Function test_applyAction(this, id, t, X, F) Result(info)
    class(model_hybrid_kmc), intent(inout), target :: this
    integer, intent(in) :: id
    real(c_double), intent(inout) :: t
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

