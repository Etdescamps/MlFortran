! Copyright (c) 2019 Etienne Descamps
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

Module mlf_stat_fun
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Public :: Digamma, Trigamma

  real(c_double), parameter :: m_gamma = 0.5772156649015328606d0 ! Euler-Maschenori constant
  real(c_double), parameter :: zeta2 = (mlf_PI**2)/6d0 ! Zeta(2)
Contains
  ! Clean versions of Digamma and Trigamma that are elemental and use NaN as an error instead
  ! of using a supplementar argument
  Elemental Real(c_double) Function Digamma(x) Result(y)
    ! Based on Jose Bernado's Algorithms AS 103
    real(c_double), intent(in) :: x
    real(c_double), parameter :: S = 1d-6, C = 8.5d0
    real(c_double), parameter :: B(5) = [-1d0/12d0, 1d0/120d0, -1d0/252d0, &
      1d0/240d0, -1d0/132d0]
    real(c_double) :: z, r, r2
    If(x <= 0) Then
      ! Invalid input x
      ! TODO??: implement negative values =/= {0, -1, -2, ...}
      y = IEEE_VALUE(y, IEEE_QUIET_NAN)
      RETURN
    Endif
    If(x <= S) Then
      ! Approximation when x is small
      ! zeta2 term come from trigamma function at x -> 0:
      ! phi'(x) = sum_i>=0 1/(i+x)^2 ~= 1/x² + Zeta(2) + o(1)
      y = -m_gamma - 1d0/x + zeta2*x
      RETURN
    Endif
    z = x
    y = 0
    Do While(z <= C)
      ! Use Digamma(x+1) = Digamma(x) + 1/x
      y = y - 1d0/z
      z = z + 1d0
    End Do
    ! Approximation when x -> infty
    r = 1d0/z
    r2 = r*r
    y = y + LOG(z) - 0.5d0*r + r2*(B(1)+r2*(B(2)+r2*(B(3)+r2*(B(4)+r2*B(5)))))
  End Function Digamma

  Elemental Real(c_double) Function Trigamma(x) Result(y)
    ! Based on BE Schneider's Algorithms AS 121
    real(c_double), intent(in) :: x
    real(c_double), parameter :: S = 1d-4, C = 5d0
    real(c_double), parameter :: B(4) = [1d0/6d0, -1d0/30d0, 1d0/42d0, -1d0/30d0]
    real(c_double) :: z, r, r2
    If(x <= 0) Then
      ! Invalid input x
      y = IEEE_VALUE(y, IEEE_QUIET_NAN)
      RETURN
    Endif
    If(x <= S) Then
      ! Approximation when x is small
      ! phi'(x) = sum_i>=0 1/(i+x)^2 ~= 1/x² + Zeta(2) + o(1)
      y = x**(-2) + zeta2
    Endif
    z = x
    y = 0
    Do While(z <= C)
      ! Use psi'(x+1) = psi'(x) - 1/x²
      y = y + z**(-2)
      z = z + 1d0
    End Do
    r = 1d0/z
    r2 = r*r
    y = y + r * (1d0+0.5d0*r+r2*(B(1)+r2*(B(2)+r2*(B(3)+r2*B(4)))))
  End Function Trigamma
End Module mlf_stat_fun

