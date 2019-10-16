! Copyright (c) 2017-2019 Etienne Descamps
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

Module mlf_binomial
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_beta_dist
  IMPLICIT NONE
  PRIVATE

  Public :: BinomialCumulative, BinomialLikelyhood
Contains
  Real(c_double) Function BinomialCumulative(p, n, k) Result(val)
    real(c_double), intent(in) :: p
    integer(8) :: n, k
    val = IncompleteBeta(REAL(n-k,8), REAL(k+1,8), 1-p)
  End Function BinomialCumulative

  Real(c_double) Function BinomialLikelyhood(p, n, k) Result(val)
    real(c_double), intent(in) :: p
    integer(8) :: n, k
    real(c_double) :: vEq, vSup, vInf
    vEq = EXP(LOG_GAMMA(REAL(n+1,8)) - LOG_GAMMA(REAL(n-k+1,8)) &
             - LOG_GAMMA(REAL(k+1,8)))*p**k*(1-p)**(n-k)
    vInf = BinomialCumulative(p, n, k)
    vSup = 1d0 - vInf + vEq
    val = MIN(vInf, vSup)
  End Function BinomialLikelyhood
End Module mlf_binomial


