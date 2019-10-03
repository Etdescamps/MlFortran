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

Module mlf_radius_sampler 
  ! This sampler is written for randomly choosing an element slighty above or below
  ! a reference value, in order to put some variability

  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_rand
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_1DRealSampler) :: mlf_radiusSampler
    real(c_double) :: dR, r0, D
    real(c_float) :: fdR, fr0, fD
  Contains
    procedure :: mlf_radiusSampler_initDouble, mlf_radiusSampler_initReal
    procedure :: random => mlf_radiusSampler_random
    procedure :: sample => mlf_radiusSampler_sample
    procedure :: sampleFloat => mlf_radiusSampler_sampleFloat
    generic :: init => mlf_radiusSampler_initDouble, mlf_radiusSampler_initReal
  End Type mlf_radiusSampler
Contains
  Integer Function mlf_radiusSampler_initDouble(this, r0, r1, D) Result(info)
    class(mlf_radiusSampler), intent(inout) :: this
    real(c_double), intent(in) :: r0, r1, D
    info = 0
    this%dR = r1 - r0; this%r0 = r0; this%D = D
    this%fdR = r1 - r0; this%fr0 = r0; this%fD = D
    Do While(REAL(this%fr0, 8) < r0) ! Select the first float supperior to r0
      this%fr0 = ieee_next_after(this%fr0, HUGE(1.0))
    End Do
    Do While(REAL(this%fr0 + this%fdR, 8) > r1) ! Select a float so that fr0 + fdR < r1
      this%fdR = ieee_next_after(this%fdR, -HUGE(1.0))
    End Do
  End Function mlf_radiusSampler_initDouble

  Integer Function mlf_radiusSampler_initReal(this, r0, r1, D) Result(info)
    class(mlf_radiusSampler), intent(inout) :: this
    real(c_float), intent(in) :: r0, r1, D
    info = 0
    this%dR = r1 - r0; this%r0 = r0; this%D = D
    this%fdR = r1 - r0; this%fr0 = r0; this%fD = D
  End Function mlf_radiusSampler_initReal

  Real(c_double) Function mlf_radiusSampler_random(this) Result(res)
    class(mlf_radiusSampler), intent(inout) :: this
    real(c_double) :: r
    CALL RANDOM_NUMBER(r)
    res = this%r0+this%dR*(r**this%D)
  End Function mlf_radiusSampler_random

  Subroutine mlf_radiusSampler_sample(this, X)
    class(mlf_radiusSampler), intent(inout) :: this
    real(c_double), intent(out) :: X(:)
    CALL RANDOM_NUMBER(X)
    X = this%r0+this%dR*(X**this%D)
  End Subroutine mlf_radiusSampler_sample

  Subroutine mlf_radiusSampler_sampleFloat(this, X)
    class(mlf_radiusSampler), intent(inout) :: this
    real(c_float), intent(out) :: X(:)
    CALL RANDOM_NUMBER(X)
    X = this%fr0+this%fdR*(X**this%fD)
  End Subroutine mlf_radiusSampler_sampleFloat
End Module mlf_radius_sampler 

