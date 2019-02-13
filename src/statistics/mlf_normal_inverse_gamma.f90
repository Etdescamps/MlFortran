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

Module mlf_normal_inverse_gamma
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_cfuns
  Use mlf_distribution
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_distribution_multivariate) :: mlf_distribution_normalInverseGamma
    real(c_double) :: mu, lambda, alpha, beta
  Contains
    procedure :: getParameters => normalInverseGamma_getParameters
    procedure :: withParameters => normalInverseGamma_withParameters
    procedure :: computePDF => normalInverseGamma_computePDF
    procedure :: computeLogPDF => normalInverseGamma_computeLogPDF
    procedure :: fitWithData => normalInverseGamma_fitWithData
  End Type mlf_distribution_normalInverseGamma

  Type, Public, Extends(mlf_distribution_univariate) :: mlf_posterior_normalInverseGamma
    real(c_double) :: mu, lambda, alpha, beta
  Contains
    procedure :: getParameters => posterior_normalInverseGamma_getParameters
    procedure :: withParameters => posterior_normalInverseGamma_withParameters
    procedure :: fitWithData => posterior_normalInverseGamma_fitWithData
    procedure :: computePDF => posterior_normalInverseGamma_computePDF
    procedure :: computeLogPDF => posterior_normalInverseGamma_computeLogPDF
    procedure :: getStats => posterior_normalInverseGamma_getStats
  End Type mlf_posterior_normalInverseGamma
Contains
  Integer Function normalInverseGamma_getParameters(this, X) Result(info)
    class(mlf_distribution_normalInverseGamma), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    info = 0
    X = [this%mu, this%lambda, this%alpha, this%beta]
  End Function normalInverseGamma_getParameters

  Integer Function normalInverseGamma_withParameters(this, X) Result(info)
    class(mlf_distribution_normalInverseGamma), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    info = 0
    this%mu = X(1); this%lambda = X(2); this%alpha = X(3); this%beta = X(4)
  End Function normalInverseGamma_withParameters

  Integer Function posterior_normalInverseGamma_getParameters(this, X) Result(info)
    class(mlf_posterior_normalInverseGamma), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    info = 0
    X = [this%mu, this%lambda, this%alpha, this%beta]
  End Function posterior_normalInverseGamma_getParameters

  Integer Function posterior_normalInverseGamma_withParameters(this, X) Result(info)
    class(mlf_posterior_normalInverseGamma), intent(inout) :: this
    real(c_double), intent(in) :: X(:)
    info = 0
    this%mu = X(1); this%lambda = X(2); this%alpha = X(3); this%beta = X(4)
  End Function posterior_normalInverseGamma_withParameters

  Elemental Real(c_double) Function posterior_normalInverseGamma_computePDF(this, x) Result(y)
    class(mlf_posterior_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: x
    real(c_double) :: a, b, l
    ASSOCIATE(alpha => this%alpha, beta => this%beta, lambda => this%lambda, mu => this%mu)
      a = alpha + 0.5d0
      l = lambda + 1d0
      b = beta + 0.5d0*lambda/(lambda+1)*(x-mu)**2
      y = SQRT(lambda)*beta**alpha*GAMMA(a)/(SQRT(2d0*mlf_PI*l)*b**a*GAMMA(alpha))
    END ASSOCIATE
  End Function posterior_normalInverseGamma_computePDF

  Elemental Real(c_double) Function posterior_normalInverseGamma_computeLogPDF(this, x) Result(y)
    class(mlf_posterior_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: x
    real(c_double) :: a, b, l
    ASSOCIATE(alpha => this%alpha, beta => this%beta, lambda => this%lambda, mu => this%mu)
      a = alpha + 0.5d0
      l = lambda + 1d0
      b = beta + 0.5d0*lambda/(lambda+1)*(x-mu)**2
      y = alpha*LOG(beta) - a*LOG(b) + LOG_GAMMA(a) - LOG_GAMMA(alpha) &
        + 0.5d0*(LOG(lambda) - LOG(l) - LOG(2d0*mlf_PI))
    END ASSOCIATE
  End Function posterior_normalInverseGamma_computeLogPDF

  Real(c_double) Function posterior_normalInverseGamma_getStats(this, statType) Result(y)
    class(mlf_posterior_normalInverseGamma), intent(in) :: this
    integer, intent(in) :: statType
    ASSOCIATE(alpha => this%alpha, beta => this%beta, lambda => this%lambda, mu => this%mu)
      Select Case(statType)
      Case (mlf_dist_mode, mlf_dist_median, mlf_dist_mean)
        y = mu
      Case Default
        y = IEEE_VALUE(y, IEEE_QUIET_NAN)
      End Select
    END ASSOCIATE
  End Function posterior_normalInverseGamma_getStats


  Integer Function posterior_normalInverseGamma_fitWithData(this, Points, W, prior) Result(info)
    class(mlf_posterior_normalInverseGamma), intent(inout) :: this
    real(c_double), intent(in) :: Points(:)
    real(c_double), intent(in), optional :: W(:)
    class(mlf_distribution_abstract), optional, intent(in) :: prior
    real(c_double) :: sN
    info = -1
    If(.NOT. PRESENT(prior)) RETURN
    Select Type(prior)
    Class is (mlf_distribution_normalInverseGamma)
      sN = REAL(SIZE(Points),8)
      this%alpha = prior%alpha + sN/2
      this%beta = prior%beta + 0.5d0*SUM(Points**2) + 0.5d0*prior%lambda*prior%mu**2 &
                - 0.5d0*(prior%lambda*prior%mu + SUM(Points))**2/(prior%lambda + sN)
      this%lambda = prior%lambda + sN
      this%mu = (prior%lambda*prior%mu + SUM(Points))/(prior%lambda + sN)
      info = 0
    End Select
  End Function posterior_normalInverseGamma_fitWithData

  Real(c_double) Function normalInverseGamma_computePDF(this, x) Result(y)
    class(mlf_distribution_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    ASSOCIATE(z => X(1), sigma => X(2), mu => this%mu, lambda => this%lambda, &
        alpha => this%alpha, beta => this%beta)
      y = SQRT(lambda)/(SQRT(2d0*mlf_PI)*sigma)*beta**alpha/GAMMA(alpha)*sigma**(-2*alpha+2) &
        * EXP(-(2d0*beta+lambda*(z-mu)**2)/(2d0*sigma**2))
    END ASSOCIATE
  End Function normalInverseGamma_computePDF

  Real(c_double) Function normalInverseGamma_computeLogPDF(this, X) Result(y)
    class(mlf_distribution_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    ASSOCIATE(z => X(1), sigma => X(2), mu => this%mu, lambda => this%lambda, &
        alpha => this%alpha, beta => this%beta)
      y = 0.5d0*(LOG(lambda) - LOG(2d0*mlf_PI)) - LOG(sigma) + alpha*LOG(beta) - LOG_GAMMA(alpha) &
        + (-2*alpha+2)*LOG(sigma) - (2d0*beta+lambda*(z-mu)**2)/(2d0*sigma**2)
    END ASSOCIATE
  End Function normalInverseGamma_computeLogPDF

  Integer Function normalInverseGamma_fitWithData(this, Points, W, prior) Result(info)
    class(mlf_distribution_normalInverseGamma), intent(inout) :: this
    real(c_double), intent(in) :: Points(:,:)
    real(c_double), intent(in), optional :: W(:)
    class(mlf_distribution_abstract), optional, intent(in) :: prior
    real(c_double) :: meanX, meanS2, varX, varS2, rN, invSumN
    ASSOCIATE(X => Points(1,:), sigma => Points(2,:), &
        mu => this%mu, lambda => this%lambda, alpha => this%alpha, beta => this%beta)
      rN = REAL(SIZE(X), 8)
      If(PRESENT(W)) Then
        invSumN = 1d0/SUM(W)
        meanX = SUM(X*W)*invSumN
        meanS2 = SUM(sigma**2*W)*invSumN
        varX = SUM(W*(X-meanX)**2)*rN/(rN-1)*invSumN
        varS2 = SUM(W*(sigma**2-meanS2)**2)*rN/(rN-1)*invSumN
      Else
        invSumN = 1d0/rN
        meanX = SUM(X)*invSumN
        meanS2 = SUM(sigma**2)*invSumN
        varX = SUM((X-meanX)**2)/(rN-1)
        varS2 = SUM((sigma**2-meanS2)**2)/(rN-1)
      Endif
      mu = meanX
      lambda = meanS2/varX
      alpha = 2d0 + meanS2**2/varS2
      beta = meanS2*(alpha-1)
    END ASSOCIATE
    info = 0
  End Function normalInverseGamma_fitWithData
End Module mlf_normal_inverse_gamma

