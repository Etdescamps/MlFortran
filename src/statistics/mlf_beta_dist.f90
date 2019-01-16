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

  Public :: LogBetaFunction, BetaFunction, RandomBeta
  Public :: IncompleteBeta, InverseIncompleteBeta, BetaDensity
  Public :: InverseIncompleteBetaInterval, IntegrateIncompleteBeta
  Public :: MaxLikelihoodBeta, MaxLikelihoodBetaPriorBeta

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

  Real(c_double) Function InverseIncompleteBeta(a, b, y) Result(x)
    real(c_double), intent(in) :: a, b, y
    real(c_double) :: F, dx, V, xMin, xMax, vMin, vMax, dMin, dMax, d0, d1
    integer, parameter :: nMax = 20
    integer :: i, lastMin, lastMax
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
    xMin = 0d0; xMax = 1d0
    vMin = -y; vMax = 1d0-y
    dMin = 0d0; dMax = 0d0
    If(a <= 1) dMin = F*(0.5d0*EPSILON(0d0))**(a-1d0)*(1d0-0.5d0*EPSILON(0d0))**(b-1d0)
    If(.NOT. IEEE_IS_FINITE(dMin)) dMin = HUGE(1d0)
    If(b <= 1) dMax = F*(0.5d0*EPSILON(0d0))**(b-1d0)*(1d0-0.5d0*EPSILON(0d0))**(a-1d0)
    If(.NOT. IEEE_IS_FINITE(dMax)) dMax = HUGE(1d0)
    lastMin = 0; lastMax = 0
    Do i = 1, nMax
      If(lastMin-lastMax > 4) x = 0.5d0*(x+xMax)
      If(lastMax-lastMin > 4) x = 0.5d0*(x+xMin)
      V = IncompleteBeta(a, b, x) - y
      If(ABS(V) < 1d-30) RETURN
      dx = F*x**(a-1d0)*(1d0-x)**(b-1d0)
      If(ABS(V) < 1d-10*dx) Then
        ! Proceed a last Newton Raphson step before returning
        x = x -V/dx
        RETURN
      Endif
      If(V < 0) Then
        xMin = x; vMin = V; dMin = dx
        lastMin = i
      Else
        xMax = x; vMax = V; dMax = dx
        lastMax = i
      Endif
      d0 = MERGE(HUGE(1d0), dMin*(xMax-xMin)/(vMax-vMin), dMin == HUGE(1d0))
      d1 = MERGE(HUGE(1d0), dMax*(xMax-xMin)/(vMax-vMin), dMax == HUGE(1d0))
      !x = xMin + (xMax-xMin)*FindRoot2PDerivative(vMin, d0, vMax, d1)
      x = xMin + (xMax-xMin)*FindRoot3Bezier(d0, d1, -vMin/(vMax-vMin))
    End Do
  End Function InverseIncompleteBeta

  ! Compute the inverse CDF on a whole interval of values: used for integration
  Subroutine InverseIncompleteBetaInterval(a, b, X)
    real(c_double), intent(in)  :: a, b
    real(c_double), intent(out) :: X(:)
    real(c_double) :: Y(SIZE(X)), DY(SIZE(X)), V, F, dx
    real(c_double) :: M0(3), M1(3), d0, d1, t, Z
    real(c_double) :: xMin, xMax, vMin, vMax, dMin, dMax
    integer :: i, j, k, N, lastMin, lastMax
    integer, parameter :: nMax = 10
    N = SIZE(X)
    X = [(REAL(i-1, 8)/REAL(N-1, 8), i=1,N)]
    X(1) = EPSILON(1d0)
    X(N) = 1d0 - EPSILON(1d0)
    Y = IncompleteBeta(a, b, X)
    Y(1) = 0d0; Y(N) = 1d0
    F = EXP(-LogBetaFunction(a,b))
    DY = F*X**(a-1d0)*(1d0-X)**(b-1d0)
    WHERE(.NOT. IEEE_IS_FINITE(DY)) DY = HUGE(1d0)
    j = 1
    M0 = [0d0, 0d0, DY(1)]
    M1 = [X(1), Y(1), DY(1)]
    X(1) = 0d0; X(N) = 1d0
    Do i = 2, N-1
      Z = X(i)
      Do While(Z > Y(j))
        j = j+1
        M0 = M1
        M1 = [REAL(j-1, 8)/REAL(N-1, 8), Y(j), DY(j)]
      End Do
      lastMin = 0; lastMax = 0
      xMin = M0(1); xMax = M1(1)
      vMin = M0(2)-Z; vMax = M1(2)-Z
      dMin = M0(3); dMax = M1(3)
      Do k = 1, nMax
        d0 = MERGE(HUGE(1d0), dMin*(xMax-xMin)/(vMax-vMin), dMin == HUGE(1d0))
        d1 = MERGE(HUGE(1d0), dMax*(xMax-xMin)/(vMax-vMin), dMax == HUGE(1d0))
        t = xMin + (xMax-xMin)*FindRoot3Bezier(d0, d1, -vMin/(vMax-vMin))
        If(lastMin-lastMax > 4) t = 0.5d0*(t+xMax)
        If(lastMax-lastMin > 4) t = 0.5d0*(t+xMin)
        V = IncompleteBeta(a, b, t) - Z
        dx = F*t**(a-1d0)*(1d0-t)**(b-1d0)
        If(ABS(V) < MIN(1d-10*dx, 1d-20)) EXIT
        If(V < 0) Then
          xMin = t; vMin = V; dMin = dx
          lastMin = i
        Else
          xMax = t; vMax = V; dMax = dx
          lastMax = i
        Endif
      End Do
      M0 = [t, V+Z, dx]
      ! Proceed a last Newton Raphson step before returning
      If(k <= nMax) t = t - V/dx
      X(i) = t
    End Do
  End Subroutine InverseIncompleteBetaInterval

  ! Compute the probability that the distribution Beta(a,b) generates
  ! an element closer to z than the distribution whose CDF Interval is X
  ! (X generated by, for example, InverseIncompleteBetaInterval)
  Real(c_double) Function IntegrateIncompleteBeta(X, a, b, z) Result(y)
    real(c_double), intent(in) :: X(:), a, b, z
    real(c_double), parameter :: icoeff(3) =  [3d0/8d0, 7d0/6d0, 23d0/24d0]
    integer :: N
    N = SIZE(X)
    y = SUM(IC(X(4:N-3))) + DOT_PRODUCT(icoeff, IC(X(1:3))) &
      + DOT_PRODUCT(icoeff(3:1:-1), IC(X(N-2:N)))
  Contains
    Elemental Real(c_double) Function IC(t) Result(r)
      real(c_double), intent(in) :: t
      real(c_double) :: delta
      delta = ABS(t-z)
      If(t + delta < 1) Then 
        r = IncompleteBeta(a, b, t+delta)
      Else
        r = 1
      Endif
      If(t - delta > 0) r = r - IncompleteBeta(a, b, t-delta)
    End Function IC
  End Function IntegrateIncompleteBeta

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
    If(x < (a+1d0)/(a+b+2d0)) Then
      y = F*BCF(a,b,x)/a
    Else
      y = 1d0-F*BCF(b,a,1-x)/b
    Endif
  Contains
    Pure Real(c_double) Function BCF(a, b, x) Result(y)
      real(c_double), intent(in) :: a, b, x
      integer, parameter :: nMax = 200
      integer :: i
      real(c_double) :: num, c, d
      real(c_double), parameter :: tinyX = 1d-30, eps = 1d-15
      c = 1d0; d = MAX(1d0 - (a+b)*x/(a+1), tinyX)
      d = 1d0/d
      y = d
      Do i=1,nMax
        ! Proceed two iterations of the Lentz's algorithm for
        ! the evaluation of continued fractions
        ! Odd iteration
        num = i*(b-i)*x/((a+2*i-1)*(a+2*i))
        d = 1d0 + num*d
        If(ABS(d) < tinyX) d = tinyX
        d = 1d0/d
        c = 1d0+num/c
        If(ABS(c) < tinyX) c = tinyX
        y = y*c*d
        ! Even iteration
        num = -(a+i)*(a+b+i)*x/((a+2*i)*(a+2*i+1))
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

  Integer Function MaxLikelihoodBetaPriorBeta(X, alpha, beta, p_alpha, p_beta, a, b) Result(info)
    real(c_double), intent(in) :: X(:), p_alpha, p_beta
    real(c_double), intent(in), optional :: a, b
    real(c_double), intent(out) :: alpha, beta
    real(c_double) :: invAB, lGa, lGb, mu, var, invN, u
    real(c_double) :: dX(2), psi(3), psi2(3), haa, hab, hbb
    real(c_double) :: a0, b0, da, db, df, d2f, dY(2), y
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
      y = alpha/(alpha+beta)
      df = (p_alpha-1)/y - (p_beta-1)/(1-y)
      d2f = -(p_alpha-1)/y**2 - (p_beta-1)/(1-y)**2
      dY = 1/(alpha+beta)*[-y, 1-y]
      ! Compute the gradient dX and the hessian matrix H=[haa hab; hab hbb]
      dX = [lGa, lGb] - psi(1:2) + psi(3) + df*dY/N
      haa = psi2(3)-psi2(1) + d2f*dY(1)**2 + df*2*(1-y)/(alpha+beta)**2/N
      hbb = psi2(3)-psi2(2) + d2f*dY(2)**2 + df*2*y/(alpha+beta)**2/N
      hab = psi2(3) + d2f*dY(1)*dY(2) + df*(2-y)/(alpha+beta)**2/N
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
  End Function MaxLikelihoodBetaPriorBeta

  Integer Function MaxLikelihoodBeta(X, alpha, beta, W, a, b) Result(info)
    real(c_double), intent(in) :: X(:)
    real(c_double), intent(in), optional :: a, b, W(:)
    real(c_double), intent(out) :: alpha, beta
    real(c_double) :: invAB, lGa, lGb, mu, var, invN, u
    real(c_double) :: dX(2), psi(3), psi2(3), haa, hab, hbb
    real(c_double) :: a0, b0, da, db
    integer, parameter :: nMax = 100
    real(c_double), parameter :: eps = 1d-12
    integer :: N, i
    info = -1
    N = SIZE(X)
    If(PRESENT(W)) Then
      invN = 1d0/SUM(W)
      mu = SUM(W*X)*invN
    Else
      invN = 1d0/(REAL(N, KIND=8))
      mu = SUM(X)*invN
    Endif
    If(PRESENT(a) .OR. PRESENT(b)) Then
      CALL InitOrDefault(a0, a, 0d0)
      CALL InitOrDefault(b0, b, 1d0)
      invAB = 1d0/(b0-a0)
      If(PRESENT(W)) Then
        lGa = SUM(W*LOG((X-a0)*invAB))*invN
        lGb = SUM(W*LOG((b0-X)*invAB))*invN
      Else
        ! Compute the logarithm of the geometric mean of X-a and b-X
        lGa = SUM(LOG((X-a0)*invAB))*invN
        lGb = SUM(LOG((b0-X)*invAB))*invN
      Endif
      ! Compute mean and variance of X
      If(PRESENT(W)) Then
        var = SUM(W*((X-mu)*invAB)**2)*REAL(N,8)/REAL(N-1,8)*invN
      Else
        var = SUM(((X-mu)*invAB)**2)/REAL(N-1,8)
      Endif
      mu = (mu-a0)*invAB
    Else If(PRESENT(W)) Then
      lGa = SUM(W*LOG(X))*invN
      lGb = SUM(W*LOG(1d0-X))*invN
      ! Unbiased weighted estimate of variance
      var = SUM(W*(X-mu)**2)/(1/invN-SUM(W*W)*invN)
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
      If(MAX(ABS(da), ABS(dX(1))) < eps*alpha .AND. &
         MAX(ABS(db), ABS(dX(2))) < eps*beta) EXIT
    End Do
    info = 0
  End Function MaxLikelihoodBeta
End Module mlf_beta_dist

