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
  Use mlf_distribution
  Use mlf_2dplane_integration
  IMPLICIT NONE
  PRIVATE

  Public :: LogBetaFunction, BetaFunction, RandomBeta
  Public :: BetaDensity, IncompleteBeta, IntegrateIncompleteBeta
  Public :: QuantileTableBeta, QuantileBeta
  Public :: MaxLikelihoodBeta, MaxLikelihoodBetaPriorBeta
  Public :: BetaPriorConst, BetaPriorMaxLL

  real(c_double), parameter :: m_gamma = 0.5772156649015328606d0 ! Euler-Maschenori constant
  
  Type, Public, Extends(mlf_distributionWithQuantile_type) :: mlf_beta_distribution
    real(c_double) :: alpha, beta, a = 0d0, b = 1d0
  Contains
    procedure :: getParameters => Beta_getParameters
    procedure :: withParameters => Beta_withParameters
    procedure :: getStats => Beta_getStats
    procedure :: fitWithData => Beta_fitWithData
    procedure :: quantile => Beta_quantile
    procedure :: quantileTable => Beta_quantileTable
    procedure :: computeCDF => Beta_computeCDF
    procedure :: computePDF => Beta_computePDF
    procedure :: computeLogPDF => Beta_computeLogPDF
  End Type mlf_beta_distribution

  Type, Public, Extends(mlf_distribution_multivariate) :: mlf_beta_prior
    real(c_double) :: lambda, x0, y0
  Contains
    procedure :: getParameters => Beta_Prior_getParameters
    procedure :: withParameters => Beta_Prior_withParameters
    procedure :: fitWithData => Beta_Prior_fitWithData
    procedure :: computePDF => Beta_Prior_computePDF
    procedure :: computeLogPDF => Beta_Prior_computeLogPDF
  End Type mlf_beta_prior

  Type, Extends(mlf_2dplane_h_fun) :: Beta_prior_integrate
    real(c_double) :: lambda, x0, y0
  Contains
    procedure :: getValDer => Beta_prior_getValDer
  End Type Beta_prior_integrate
Contains
  Integer Function Beta_getParameters(this, X) Result(info)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    info = 0
    Select Case(SIZE(X))
    Case(2)
      X = [this%alpha, this%beta]
    Case(4)
      X = [this%alpha, this%beta, this%a, this%b]
    Case Default
      info = -1
    End Select
  End Function Beta_getParameters

  Integer Function Beta_withParameters(this, X) Result(info)
    class(mlf_beta_distribution), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    info = 0
    Select Case(SIZE(X))
    Case(2)
      this%alpha = X(1); this%beta = X(2)
      this%a = 0; this%b = 1
    Case(4)
      this%alpha = X(1); this%beta = X(2)
      this%a = X(3); this%b = X(4)
    Case Default
      info = -1
    End Select
  End Function Beta_withParameters

  Real(c_double) Function Beta_getStats(this, statType) Result(y)
    class(mlf_beta_distribution), intent(in) :: this
    integer, intent(in) :: statType
    ASSOCIATE(a => this%a, b => this%b, alpha => this%alpha, beta => this%beta)
      Select Case(statType)
      Case (mlf_dist_mode)
        y = a + (b-a)*(alpha - 1d0)/(alpha + beta - 2d0)
      Case (mlf_dist_median)
        y = this%quantile(0.5d0)
      Case (mlf_dist_mean)
        y = a + (b-a)*alpha/(alpha+beta)
      Case (mlf_dist_geom_mean)
        y = Digamma(alpha)-Digamma(alpha+beta)
      Case Default
        y = IEEE_VALUE(y, IEEE_QUIET_NAN)
      End Select
    END ASSOCIATE
  End Function Beta_getStats

  Integer Function Beta_fitWithData(this, Points, W, prior) Result(info)
    class(mlf_beta_distribution), intent(inout) :: this
    real(c_double), intent(in) :: Points(:)
    real(c_double), intent(in), optional :: W(:)
    class(mlf_distribution_abstract), optional, intent(in) :: prior
    If(PRESENT(prior)) Then
      Select Type(prior)
      Class is (mlf_beta_distribution)
        info = MaxLikelihoodBetaPriorBeta(Points, this%alpha, this%beta, &
          prior%alpha, prior%beta, this%a, this%b)
      Class Default
        info = -1
      End Select
    Else
      info = MaxLikelihoodBeta(Points, this%alpha, this%beta, W, this%a, this%b)
    Endif
  End Function Beta_fitWithData

  Subroutine Beta_quantileTable(this, X)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    CALL QuantileTableBeta(this%alpha, this%beta, X)
    X = this%a + (this%b-this%a)*X
  End Subroutine Beta_quantileTable

  Real(c_double) Function Beta_quantile(this, y) Result(x)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(in) :: y
    x = this%a + (this%b-this%a)*QuantileBeta(this%alpha, this%beta, y)
  End Function Beta_quantile

  Real(c_double) Function Beta_computeCDF(this, x) Result(y)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    If(x <= this%a) Then
      y = 0
    Else If(x >= this%b) Then
      y = 1
    Else
      y = IncompleteBeta(this%alpha, this%beta, (x-this%a)/(this%b-this%a))
    Endif
  End Function Beta_computeCDF

  Elemental Real(c_double) Function Beta_computePDF(this, x) Result(y)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(in) :: x 
    If(x >= this%a .AND. x <= this%b) Then
      y = BetaDensity(this%alpha, this%beta, (x-this%a)/(this%b-this%a))/(this%b-this%a)
    Else
      y = 0
    Endif
  End Function Beta_computePDF

  Elemental Real(c_double) Function Beta_computeLogPDF(this, x) Result(y)
    class(mlf_beta_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    If(x >= this%a .AND. x <= this%b) Then
      y = LogBetaDensity(this%alpha, this%beta, (x-this%a)/(this%b-this%a))-LOG(this%b-this%a)
    Else
      y = -HUGE(1d0)
    Endif
  End Function Beta_computeLogPDF

  Integer Function Beta_Prior_getParameters(this, X) Result(info)
    class(mlf_beta_prior), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    X = [this%lambda, this%x0, this%y0]
    info = 0
  End Function Beta_Prior_getParameters

  Integer Function Beta_Prior_withParameters(this, X) Result(info)
    class(mlf_beta_prior), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    this%lambda = X(1)
    this%x0     = X(2)
    this%y0     = X(3)
    info = 0
  End Function Beta_Prior_withParameters

  Real(c_double) Function Beta_Prior_computeLogPDF(this, X) Result(y)
    class(mlf_beta_prior), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    ASSOCIATE(x0 => this%x0, y0 => this%y0, lambda => this%lambda, &
        alpha => X(1), beta => X(2))
      y = lambda*(LOG_GAMMA(alpha+beta) - LOG_GAMMA(alpha) - LOG_GAMMA(beta)) &
        + alpha*LOG(x0) + beta*LOG(y0)
    END ASSOCIATE
  End Function Beta_Prior_computeLogPDF

  Real(c_double) Function Beta_Prior_computePDF(this, X) Result(y)
    class(mlf_beta_prior), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    y = EXP(this%computeLogPDF(X))
  End Function Beta_Prior_computePDF

  Integer Function Beta_Prior_fitWithData(this, Points, W, prior) Result(info)
    class(mlf_beta_prior), intent(inout) :: this
    real(c_double), intent(in) :: Points(:,:)
    real(c_double), intent(in), optional :: W(:)
    class(mlf_distribution_abstract), optional, intent(in) :: prior
    real(c_double) :: sLG, sAlpha, sBeta
    info = -1
    sLG = 0; sAlpha = 0; sBeta = 0
    If(PRESENT(W)) Then
      sLG = SUM(W*(LOG_GAMMA(Points(1,:) + Points(2,:)) &
                - LOG_GAMMA(Points(1,:)) - LOG_GAMMA(Points(2,:))))
      sAlpha = SUM(W*Points(1,:))
      sBeta = SUM(W*Points(2,:))
    Else
      sLG = SUM(LOG_GAMMA(Points(1,:) + Points(2,:)) &
                - LOG_GAMMA(Points(1,:)) - LOG_GAMMA(Points(2,:)))
      sAlpha = SUM(Points(1,:))
      sBeta = SUM(Points(2,:))
    Endif
  End Function Beta_Prior_fitWithData

  Pure Subroutine BetaPriorMaxLL(lambda, x0, y0, alpha, beta)
    real(c_double), intent(in) :: lambda, x0, y0
    real(c_double), intent(out) :: alpha, beta
    real(c_double) :: lX0, lY0, Fa
    real(c_double) :: X(2), F(2), H(3), Tr(3), invDetH
    integer :: i
    integer, parameter :: nMax = 50
    If(lambda <= 0) Then
      alpha = 0d0; beta = 0d0
      RETURN
    Endif
    alpha = IEEE_VALUE(alpha, IEEE_QUIET_NAN)
    beta  = IEEE_VALUE(beta,  IEEE_QUIET_NAN)
    ! We try to find the maximum likelihood of the distribution
    ! A(alpha, beta) ~ Beta(alpha, beta)^-lambda * x0^alpha * y0^beta

    lX0 = LOG(x0); lY0 = LOG(y0)

    ! Does not converge in this case: this function can go to infinity
    If(EXP(lX0/lambda) + EXP(lY0/lambda) >= 1) RETURN

    ! The maximum likelihood is found using a Newton's method
    ! We use the derivatives of this function A
    ! dA/dalpha(alpha,beta) = A(alpha,beta)*f1(alpha,beta)
    ! dA/dbeta(alpha,beta)  = A(alpha,beta)*f2(alpha,beta)
    ! with f1 and f2 defined as:
    ! f1(alpha, beta) = LOG(x0) - lambda*(Digamma(alpha) - Digamma(alpha+beta))
    ! f2(alpha, beta) = LOG(y0) - lambda*(Digamma(beta)  - Digamma(alpha+beta))

    ! Find a first approximation of alpha and beta:
    ! Using f1 = 0 = f2 we have Digamma(alpha)-Digamma(beta) = (log(x0)-log(y0))/lambda
    ! We use the approximation Digamma(alpha)-Digamma(beta) ~= log(alpha)-log(beta)
    ! => beta ~= alpha*U with U = (y0/x0)^1/lambda
    ! We try to find the root of the function
    ! u(alpha) = Digamma(alpha*(1+U))-1/2*(Digamma(alpha)+Digamma(alpha*U))+log(x0*y0)/2/lambda
    If(x0 < y0) Then
      Fa = EXP((lY0-lX0)/lambda)
      X(2) = EstimateRoot(Fa, lX0/lambda)
      X(1) = X(2)/Fa
    Else
      Fa = EXP((lX0-lY0)/lambda)
      X(1) = EstimateRoot(Fa, lY0/lambda)
      X(2) = X(1)/Fa
    Endif

    If(ANY(IEEE_IS_NAN(X))) X = [1d0, 1d0]
    ! Apply Newton's method on log likelihood for finding the roots
    Do i = 1, nMax
      F = lambda*(Digamma(X(1)+X(2)) - Digamma(X)) + [lX0, lY0]
      If(SUM(F*F) < MAX(1d-25, -1d-20*lambda*(lX0+lY0))) EXIT
      Tr = Trigamma([X(1), X(2), X(1)+X(2)])
      H = lambda*(Tr(3) - [Tr(1), Tr(2), 0d0])
      invDetH = 1d0/(H(1)*H(2)-H(3)**2)
      X = MAX(X - invDetH*[H(2)*F(1)-H(3)*F(2), -H(3)*F(1)+H(1)*F(2)], 0.25d0*X)
    End Do
    alpha = X(1)
    beta  = X(2)
  Contains
    Pure Real(c_double) Function UA(z, U)
      real(c_double), intent(in) :: z, U
      UA = Digamma(z*(1d0+U)) - Digamma(z)
    End Function UA

    Pure Real(c_double) Function DUA(z, U)
      real(c_double), intent(in) :: z, U
      DUA = (1d0+U)*Trigamma(z*(1d0+U)) - Trigamma(z)
    End Function DUA

    Pure Real(c_double) Function EstimateRoot(U, V) Result(z)
      real(c_double), intent(in) :: U, V
      real(c_double) :: z0, z1, a0, a1, a, da
      integer :: i
      integer, parameter :: nMax = 30
      z = IEEE_VALUE(z, IEEE_QUIET_NAN)
      a = UA(1d0, U) + V
      da = DUA(1d0, U)
      If(a*da > 0d0) Then
        z0 = 0.5d0; z1 = 1d0
        a0 = UA(z0, U) + V; a1 = a
        Do While(a0*a1 > 0d0)
          z1 = z0
          a1 = a0
          z0 = 0.5d0*z0
          a0 = UA(z0, U) + V
        End Do
      Else
        z0 = 1d0; z1 = 2d0
        a0 = a; a1 = UA(z1, U) + V
        Do While(a0*a1 > 0d0)
          z0 = z1
          a0 = a1
          z1 = 2d0*z1
          a1 = UA(z1, U) + V
          If(ABS(a1-a0)/z1 < 1d-20) RETURN
        End Do
      Endif
      z = z0+(z1-z0)*a0/(a0-a1)
      a = UA(z, U) + V
      Do i= 1, nMax
        If(ABS(a) < 1d-6) RETURN
        If(a*a0 > 0) Then
          a0 = a; z0 = z
        Else
          a1 = a; z1 = z
        Endif
        da = DUA(z, U)
        z = z - a/da
        If(da == 0 .OR. z < z0 .OR. z > z1) z = z0+(z1-z0)*a0/(a0-a1)
        a = UA(z, U) + V
      End Do
    End Function EstimateRoot
  End Subroutine BetaPriorMaxLL

  Subroutine Beta_prior_getValDer(this, x, y, D)
    class(Beta_prior_integrate), intent(inout) :: this
    real(c_double), intent(in) :: x, y
    type(mlf_2d_h_val), intent(out) :: D
    real(c_double) :: Dg(3), Tg(3), Ga, Gb, F
    ASSOCIATE(z => x+y, lambda => this%lambda, lX0 => LOG(this%x0), lY0 => LOG(this%y0))
      F = EXP(lambda*(LOG_GAMMA(z)-LOG_GAMMA(x)-LOG_GAMMA(y)) + lX0*x + lY0*y)
      Dg = Digamma([x, y, z])
      Tg = Trigamma([x, y, z])
      Ga = lambda*(Dg(3)-Dg(1)) + lX0
      Gb = lambda*(Dg(3)-Dg(2)) + lY0
      D%val = F
      D%der = F*[Ga, Gb]
      D%hes(1) = F*(Ga**2 + lambda*(Tg(3)-Tg(1)))
      D%hes(2) = F*(Gb**2 + lambda*(Tg(3)-Tg(2)))
      D%hes(3) = F*(Ga*Gb + lambda*Tg(3))
    END ASSOCIATE
  End Subroutine Beta_prior_getValDer

  Real(c_double) Function BetaPriorGetBorders(lambda, x0, y0, dx, b, N) Result(y)
    real(c_double), intent(in) :: lambda, x0, y0, dx, b
    integer, intent(in) :: N
    real(c_double) :: lY0, S, u, y0b, F, dxL
    integer :: i
    lY0 = LOG(y0)
    y0b = EXP(lY0*b)
    F = y0b/REAL(N+1,8)
    S = 0
    Do i = 1, N
      u = F*REAL(i)
      S = S + LOG_GAMMA(2*LOG(u)/lY0)*u
    End Do
    S = 2*S + LOG_GAMMA(b)*EXP(0.5d0*lY0*b)
    dxL = dx**(lambda+1d0)
    y = -y0b/lY0*dxL*(1d0/(lambda+1)+dx/(lambda+2)*(lambda*m_gamma+LOG(x0))) &
      + lambda/(lambda+2)*dxL*dx*(S-y0b*LOG_GAMMA(b))
  End Function BetaPriorGetBorders

  Real(c_double) Function BetaPriorConst2(lambda, x0, y0, eps, N, aMax, bMax) Result(y)
    real(c_double), intent(in) :: lambda, x0, y0, eps
    real(c_double), intent(in), optional :: aMax, bMax
    integer, intent(in) :: N
    real(c_double) :: alpha0, beta0, rN, rNij, invM, invM2, maxL, lX0, lY0
    real(c_double) :: dx, dy
    type(Beta_prior_integrate) :: fun
    If(PRESENT(aMax) .AND. PRESENT(bMax)) Then
      alpha0 = aMax; beta0 = bMax
    Else
      CALL BetaPriorMaxLL(lambda, x0, y0, alpha0, beta0)
    Endif
    fun%lambda = lambda; fun%x0 = x0; fun%y0 = y0
    dx = 1d-3*eps**(1d0/lambda); dy = dx
    y = BetaPriorGetBorders(lambda, x0, y0, dx, dy, N) &
      + BetaPriorGetBorders(lambda, y0, x0, dy, dx, N)
    y = y + fun%integrateOnPlane(eps, 20, x0 = x0, y0 = y0, xMin = dx, yMin = dy)
  End Function BetaPriorConst2

  Elemental Real(c_double) Function BetaPriorConst(lambda, x0, y0, N, aMax, bMax) Result(y)
    real(c_double), intent(in) :: lambda, x0, y0
    real(c_double), intent(in), optional :: aMax, bMax
    integer, intent(in) :: N
    real(c_double) :: alpha0, beta0, rN, rNij, invM, invM2, maxL, lX0, lY0
    real(c_double) :: S, Sa, alpha, beta, lGa, a, va, vb
    integer :: Nk, Nl, i, j, k, l, M
    If(PRESENT(aMax) .AND. PRESENT(bMax)) Then
      alpha0 = aMax; beta0 = bMax
    Else
      CALL BetaPriorMaxLL(lambda, x0, y0, alpha0, beta0)
    Endif
    Nk = CEILING(alpha0); Nl = CEILING(beta0)
    rN = REAL(N,8)
    M = MAX(CEILING(rN/(REAL(Nk,8)*REAL(Nl,8))), 8)
    invM = 1d0/REAL(M,8)
    invM2 = invM**2
    lX0 = LOG(x0); lY0 = LOG(y0)
    maxL = EXP(lambda*(LOG_GAMMA(alpha0+beta0)-LOG_GAMMA(alpha0)-LOG_GAMMA(beta0)) &
              + alpha0*lX0 + beta0*lY0)
    y = 0
    Do i = 1, M
      alpha = REAL(i,8)*invM
      lGa = LOG_GAMMA(alpha)
      S = 0
      Do j = 1, M
        beta = REAL(j,8)*invM
        a = EXP(lambda*(LOG_GAMMA(alpha+beta)-lGa-LOG_GAMMA(beta)) + alpha*lX0 + beta*lY0)
        Sa = 0
        va = a
        Do k = 1, MAX(Nk*4, Nk+20)
          vb = va
          Do l = 1, MAX(Nl*4, Nl+20)
            Sa = Sa + vb
            vb = vb*((alpha+beta+k+l)/(beta+l))**lambda*y0
            If(l > Nl .AND. vb < 1d-8*maxL) Then
              Sa = Sa + vb; EXIT
            Endif
          End Do
          va = va*((alpha+beta+k)/(alpha+k))**lambda*x0
          If(k > Nk .AND. va < 1d-8*maxL) EXIT
        End Do
        S = S + Sa
      End Do
      y = y + S*invM2
    End Do
  End Function BetaPriorConst

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

  Elemental Real(c_double) Function LogBetaDensity(a, b, x) Result(y)
    real(c_double), intent(in) :: a, b, x
    y = LOG(x)*(a-1)+LOG(1-x)*(b-1)-LogBetaFunction(a,b)
  End Function LogBetaDensity


  Elemental Real(c_double) Function LogBetaFunction(a, b) Result(y)
    real(c_double), intent(in) :: a, b
    y = LOG_GAMMA(a)+LOG_GAMMA(b)-LOG_GAMMA(a+b)
  End Function LogBetaFunction

  Elemental Real(c_double) Function BetaFunction(a, b) Result(y)
    real(c_double), intent(in) :: a, b
    y = EXP(LogBetaFunction(a, b))
  End Function BetaFunction

  Real(c_double) Function QuantileBeta(a, b, y) Result(x)
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
  End Function QuantileBeta

  ! Compute the inverse CDF on a whole interval of values: used for integration
  Subroutine QuantileTableBeta(a, b, X)
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
  End Subroutine QuantileTableBeta

  ! Compute the probability that the distribution Beta(a,b) generates
  ! an element closer to z than the distribution whose Quantile Table is X
  ! (X generated by, for example, QuantileTableBeta)
  Real(c_double) Function IntegrateIncompleteBeta(X, a, b, z) Result(y)
    real(c_double), intent(in) :: X(:), a, b, z
    real(c_double), parameter :: icoeff(3) =  [3d0/8d0, 7d0/6d0, 23d0/24d0]
    integer :: N
    N = SIZE(X)
    y = (SUM(IC(X(4:N-3))) + DOT_PRODUCT(icoeff, IC(X(1:3))) &
       + DOT_PRODUCT(icoeff(3:1:-1), IC(X(N-2:N))))/REAL(N, 8)
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

  ! Return the regularized incomplete beta function I_x(a,b)
  ! This algorithm is described in Numerical Recipes in C by Press
  ! (Chapter 6 Special Functions -> 6.4 Incomplete Beta Functions
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

