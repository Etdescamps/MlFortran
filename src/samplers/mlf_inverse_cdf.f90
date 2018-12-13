! Copyright (c) 2017-2018 Etienne Descamps
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

Module mlf_inverse_cdf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_rand
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_1DRealSampler) :: mlf_inverseCDFSampler
    real(c_double), allocatable :: yI(:), X(:), Y(:)
  Contains
    procedure :: init => mlf_inverseCDF_init
    procedure :: random => mlf_inverseCDF_random
    !procedure :: sample => mlf_inverseCDF_sample
  End Type mlf_inverseCDFSampler
Contains
  Integer Function mlf_inverseCDF_init(this, X, Y) Result(info)
    class(mlf_inverseCDFSampler), intent(inout) :: this
    real(c_double), intent(in) :: X(:), Y(:)
    integer :: i
    info = -1
    If(SIZE(X) < 1 .OR. SIZE(X) /= SIZE(Y)) RETURN
    this%X = X; this%Y = Y
    ALLOCATE(this%yI(SIZE(X)-1))
    ASSOCIATE(yI => this%yI)
      yI(1) = 0.5*(Y(1)+Y(2))*(X(2)-X(1))
      Do i = 2,SIZE(X)-1
        yI(i) = 0.5*(Y(i)+Y(i+1))*(X(i+1)-X(i))+yI(i-1)
      End Do
    END ASSOCIATE
    info = 0
  End Function mlf_inverseCDF_init

  Real(c_double) Function TriangleSample(revert) Result(res)
    logical, intent(in) :: revert
    CALL RANDOM_NUMBER(res)
    res = SQRT(res)
    If(revert) res = 1d0 - res
  End Function TriangleSample

  Real(c_double) Function mlf_inverseCDF_random(this) Result(res)
    class(mlf_inverseCDFSampler), intent(inout) :: this
    real(c_double) :: r
    integer(c_int) :: i
    ASSOCIATE(yI => this%yI, X => this%X, Y => this%Y)
      CALL RANDOM_NUMBER(r)
      r = r*yI(SIZE(yI))
      i = mlf_di_search(yI, r)
      If(i > 1) r = r - yI(i-1)
      If(r < MIN(Y(i), Y(i+1))*(X(i+1)-X(i))) Then
        CALL RANDOM_NUMBER(r)
        res = X(i)+r*(X(i+1)-X(i))
      Else
        res = X(i)+(X(i+1)-X(i))*TriangleSample(Y(i) > Y(i+1))
      Endif
    END ASSOCIATE
  End Function mlf_inverseCDF_random

  !Subroutine mlf_inverseCDFSampler(this, X)
  !  class(mlf_inverseCDFSampler), intent(inout) :: this
  !  real(c_double), intent(out) :: X(:)
  !End Subroutine mlf_inverseCDFSampler
End Module mlf_inverse_cdf

