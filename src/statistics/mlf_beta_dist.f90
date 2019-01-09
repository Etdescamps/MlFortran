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

Module mlf_beta_dist
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_gamma_dist
  Use mlf_poly
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_beta_maxlikelihood
  Public :: LogBetaFunction, BetaFunction, RandomBeta
  Public :: IncompleteBeta, InverseIncompleteBeta, BetaDensity

Contains
  Real(c_double) Function RandomBeta(a, b) Result(y)
    real(c_double), intent(in) :: a, b
    real(c_double) :: R(2), X(2), S
    If(a <= 1d0 .AND. b <= 1d0) Then
      ! Johnk's algorithm
      S = HUGE(0d0)
      Do While(S > 1d0)
        CALL RANDOM_NUMBER(R)
        X = R**(1d0/[a, b])
        S = SUM(X)
      End Do
      If(S > 0) Then
        y = X(1)/S
      Else
        X = LOG(R)/[a,b]
        X = X - MAXVAL(X)
        S = SUM(EXP(X))
        y = EXP(X(1) - LOG(S))
      Endif
    Else
      ! Use gamma distribution sampling
      R(1) = RandomGamma(a)
      R(2) = RandomGamma(b)
      y = R(1)/SUM(R)
    Endif
  End Function RandomBeta

  Elemental Real(c_double) Function BetaDensity(a, b, x) Result(y)
    real(c_double), intent(in) :: a, b, x
    y = x**(a-1)*(1-x)**(b-1)/BetaFunction(a,b)
  End Function BetaDensity

  Elemental Real(c_double) Function LogBetaFunction(a, b) Result(y)
    real(c_double), intent(in) :: a, b
    y = LOG_GAMMA(a)+LOG_GAMMA(b)-LOG_GAMMA(a+b)
  End Function LogBetaFunction

  Elemental Real(c_double) Function BetaFunction(a, b) Result(y)
    real(c_double), intent(in) :: a, b
    y = EXP(LogBetaFunction(a, b))
  End Function BetaFunction

  Elemental Real(c_double) Function InverseIncompleteBeta(a, b, y) Result(x)
    real(c_double), intent(in) :: a, b, y
    real(c_double) :: F, dx, V
    integer, parameter :: nMax = 30
    integer :: i
    If(y <= 0 .OR. y >= 1) Then
      If(y == 0) Then
        x = 0
      Else If(y == 1) Then
        x = 1
      Else
        x = IEEE_VALUE(y, IEEE_QUIET_NAN)
      Endif
      RETURN
    Endif
    F = EXP(-LogBetaFunction(a,b))
    ! Evaluate first at y = x
    x = y
    Do i = 1, nMax
      V = IncompleteBeta(a, b, x) - y
      If(V == 0) RETURN
      dx = F*x**(a-1d0)*(1d0-x)**(b-1d0)
      If(V < 0) Then
        ! Root is located between x and 1
        x = x+(1d0-x)*(1d0-FindRoot2DInterp(1d0-y, V, -dx/(1d0-x)))
      Else
        ! Root is located between 0 and x
        x = x*FindRoot2DInterp(-y, V, dx/x)
      Endif
      If(ABS(V) < 1e-20) RETURN
    End Do
  End Function InverseIncompleteBeta

  Elemental Real(c_double) Function IncompleteBeta(a, b, x) Result(y)
    real(c_double), intent(in) :: a, b, x
    real(c_double) :: F
    If(x <= 0 .OR. x >= 1) Then
      If(x == 0) Then
        y = 0
      Else If(x == 1) Then
        y = 1
      Else
        y = IEEE_VALUE(y, IEEE_QUIET_NAN)
      Endif
      RETURN
    Endif
    F = EXP(a*LOG(x)+b*LOG(1d0-x)-LogBetaFunction(a,b))
    If(x <= (a+1d0)/(a+b+2d0)) Then
      y = F*BCF(a,b,x)
    Else
      y = 1d0-F*BCF(a,b,x)
    Endif
  Contains
    Pure Real(c_double) Function BCF(a, b, x) Result(y)
      real(c_double), intent(in) :: a, b, x
      integer, parameter :: nMax = 200
      integer :: i
      real(c_double) :: num, c, d
      real(c_double), parameter :: tinyX = 1d-30, eps = 1d-10
      y = 2d0; c = 2d0; d = 1d0
      Do i=1,nMax
        ! Proceed two iterations of the Lentz's algorithm for
        ! the evaluation of continued fractions
        ! Odd iteration
        num = -((a+i)*(a+b+i)*x)/((a+2d0*i)*(a+2d0*i+1))
        d = 1d0 + num*d
        If(ABS(d) < tinyX) d = tinyX
        d = 1d0/d
        c = 1d0+num/c
        If(ABS(c) < tinyX) c = tinyX
        y = y*c*d
        ! Even iteration
        num = (i*(b-i)*x)/((a+2d0*i-1d0)*(a+2d0*i))
        d = 1d0 + num*d
        If(ABS(d) < tinyX) d = tinyX
        d = 1d0/d
        c = 1d0+num/c
        If(ABS(c) < tinyX) c = tinyX
        y = y*c*d
        If(ABS(1d0-c*d) < eps) RETURN
      End Do
      ! Error: has not converge after nMax iterations
      y = IEEE_VALUE(y, IEEE_QUIET_NAN)
    End Function BCF
  End Function IncompleteBeta

  Integer Function mlf_beta_maxlikelihood(X, alpha, beta, a, b) Result(info)
    real(c_double), intent(in) :: X(:)
    real(c_double), intent(in), optional :: a, b
    real(c_double), intent(out) :: alpha, beta
    real(c_double) :: invAB, lGa, lGb, mu, var, invN, u
    real(c_double) :: dX(2), psi(3), psi2(3), haa, hab, hbb
    real(c_double) :: a0, b0, da, db
    integer, parameter :: nMax = 100
    real(c_double), parameter :: eps = 1d-9
    integer :: N, i
    info = -1
    N = SIZE(X)
    invN = 1d0/(REAL(N, KIND=8))
    mu = SUM(X)*invN
    If(PRESENT(a) .OR. PRESENT(b)) Then
      CALL InitOrDefault(a0, a, 0d0)
      CALL InitOrDefault(b0, b, 1d0)
      invAB = 1d0/(b0-a0)
      ! Compute the logarithm of the geometric mean of X-a and b-X
      lGa = SUM(LOG((X-a0)*invAB))*invN
      lGb = SUM(LOG((b0-X)*invAB))*invN
      ! Compute mean and variance of X
      var = SUM(((X-mu)*invAB)**2)/REAL(N-1,KIND=8)
      mu = (mu-a0)*invAB
    Else
      lGa = SUM(LOG(X))*invN
      lGb = SUM(LOG(1d0-X))*invN
      var = SUM((X-mu)**2)/REAL(N-1,KIND=8)
    Endif
    u = mu*(1-mu)
    If(var < u) Then
      ! Use the method of moments to estimate alpha and beta
      u = u/var-1 ! positive
      alpha = mu*u
      beta = (1-mu)*u
    Else
      ! N.L. Johnson and S. Kotz estimation from Continous Univariate Distribution Vol. 2
      u = 1d0/(1-EXP(lGa)-EXP(lGb))
      alpha = 0.5d0*(1d0+EXP(lGa)*u)
      beta  = 0.5d0*(1d0+EXP(lGb)*u)
    Endif
    ! Do a few iteration of the Newton's methods using the Hessian matrix
    Do i = 1, nMax
      ! Compute psi and psi'
      psi  = Digamma([alpha, beta, alpha+beta])
      psi2 = Trigamma([alpha, beta, alpha+beta])
      ! Compute the gradient dX and the hessian matrix H=[haa hab; hab hbb]
      dX = [lGa, lGb] - psi(1:2) + psi(3)
      haa = psi2(3)-psi2(1)
      hbb = psi2(3)-psi2(2)
      hab = psi2(3)
      u = haa*hbb-hab**2
      If(u == 0d0) RETURN
      u = 1d0/u
      ! Apply a Newton step on [alpha beta]
      da = -u*(hbb*dX(1)-hab*dX(2))
      db = -u*(haa*dX(2)-hab*dX(1))
      alpha = MAX(0d0, alpha + da)
      beta = MAX(0d0, beta + db)
      If(ABS(da) < eps*alpha .AND. ABS(db) < eps*beta) EXIT
    End Do
    info = 0
  End Function mlf_beta_maxlikelihood
End Module mlf_beta_dist

