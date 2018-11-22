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

  Type, Public, Extends(mlf_hybrid_kmc_model) :: test_hybrid_kmc
    integer(c_int64_t), pointer :: NIndiv(:)
    real(c_double), pointer :: Alpha, Beta, Volume
  Contains
    procedure :: init => test_hybrid_init
    procedure :: applyAction => test_applyAction
    procedure :: funTransitionRates => test_funTransitionRates
    procedure :: evalOde => test_evalOde
  End Type test_hybrid_kmc
Contains
    Integer Function test_evalOde(this, t, X, F)
      class(test_hybrid_kmc), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Function test_evalOde
End Module test_hybrid_kmc_model

