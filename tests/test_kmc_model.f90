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

Module test_kmc_model
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_kmc
  Use iso_fortran_env
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_ode_fun) :: kmc_ode
    real(c_double) :: Alpha, Beta
  Contains
    procedure :: eval => kmc_ode_fun
  End Type kmc_ode

  Type, Public, Extends(mlf_kmc_model) :: kmc_model
    integer(c_int64_t), pointer :: NIndiv(:)
    real(c_double), pointer :: Alpha, Beta, Volume
  Contains
    procedure :: applyAction => kmc_applyAction
    procedure :: funTransitionRates => kmc_funTransitionRates
    procedure :: init => kmc_init
  End Type kmc_model

Contains
  Integer Function kmc_ode_fun(this, t, X, F) Result(info)
    class(kmc_ode), intent(in), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    ASSOCIATE(alpha => this%Alpha, beta => this%Beta, &
        cS => X(1), cI => X(2), cR => X(3))
      F(1) = -alpha*cS*cI
      F(3) = beta*cI
      F(2) = -F(1) - F(3)
      info = mlf_ODE_Continue
      If(cI == 0) info = mlf_ODE_StopTime
    END ASSOCIATE
  End Function kmc_ode_fun

  Integer Function kmc_init(this, Alpha, Beta, Volume, NIndiv, data_handler) Result(info)
    class(kmc_model), intent(inout), target :: this
    real(c_double), intent(in), optional :: Alpha, Beta, Volume
    integer(c_int64_t), intent(in), optional :: NIndiv(3)
    class(mlf_data_handler), intent(inout), optional :: data_handler
    type(mlf_step_numFields) :: numFields
    integer(c_int64_t) :: NCat
    CALL numFields%initFields(nRPar = 3, nRsc = 1)
    info = mlf_kmc_init(this, numFields, NActions = 2, data_handler = data_handler)
    If(info < 0) RETURN
    NCat = 3
    info = this%add_i64array(numFields, NCat, this%NIndiv, C_CHAR_"NIndiv", &
      data_handler = data_handler)
    If(info < 0) RETURN
    If(PRESENT(NIndiv)) this%NIndiv = NIndiv
    CALL this%addRPar(numFields, this%Alpha, "Alpha")
    CALL this%addRPar(numFields, this%Beta, "Beta")
    CALL this%addRPar(numFields, this%Volume, "Volume")
    If(PRESENT(Alpha)) this%Alpha = Alpha
    If(PRESENT(Beta)) this%Beta = Beta
    If(PRESENT(Volume)) this%Volume = Volume
  End Function kmc_init

  Integer Function kmc_funTransitionRates(this, Rates) Result(N)
    class(kmc_model), intent(inout), target :: this
    real(c_double), intent(out), target :: Rates(:)
    ASSOCIATE(alpha => this%Alpha, beta => this%Beta, Volume => this%Volume, &
        nS => REAL(this%NIndiv(1),8), nI => REAL(this%NIndiv(2),8))
      Rates(1) = alpha*nS*nI/Volume
      Rates(2) = beta*nI
    END ASSOCIATE
    N = 2
  End Function kmc_funTransitionRates

  Integer Function kmc_applyAction(this, id, Rate) Result(info)
    class(kmc_model), intent(inout), target :: this
    real(c_double), intent(in) :: Rate
    integer, intent(in) :: id
    SELECT CASE(id)
    Case(1)
      this%NIndiv(:) = this%NIndiv(:) + [-1,+1,0]
    Case(2)
      this%NIndiv(:) = this%NIndiv(:) + [0,-1,+1]
    END SELECT
    info = 0
  End Function kmc_applyAction

End Module test_kmc_model

