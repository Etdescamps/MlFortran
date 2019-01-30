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

Module mlf_student_t_dist
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_cfuns
  Use mlf_distribution
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_distribution_type) :: mlf_distribution_normalInverseGamma
    real(c_double) :: mu, lambda, alpha, beta
  Contains
    procedure :: computePDF => normalInverseGamma_computePDF
    procedure :: computeLogPDF => normalInverseGamma_computeLogPDF
    procedure :: fitWithData => normalInverseGamma_fitWithData
  End Type mlf_distribution_normalInverseGamma
Contains
  Real(c_double) Function normalInverseGamma_computePDF(this, x) Result(y)
    class(mlf_distribution_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    ASSOCIATE(z => X(1), sigma => X(2), mu => this%mu, lambda => this%lambda, &
        alpha => this%alpha, beta => this%beta)
      y = SQRT(lambda)/(SQRT(2d0*mlf_PI)*sigma)*beta**alpha/GAMMA(alpha)*sigma**(-2*alpha+2) &
        * EXP(-(2d0*beta+lambda*(z-mu)**2)/(2d0*sigma**2))
    END ASSOCIATE
  End Function normalInverseGamma_computePDF

  Real(c_double) Function normalInverseGamma_computeLogPDF(this, x) Result(y)
    class(mlf_distribution_normalInverseGamma), intent(in) :: this
    real(c_double), intent(in) :: X(:)
    ASSOCIATE(z => X(1), sigma => X(2), mu => this%mu, lambda => this%lambda, &
        alpha => this%alpha, beta => this%beta)
      y = 0.5d0*(LOG(lambda)-LOG(2d0*mlf_PI))-LOG(sigma)+alpha*LOG(beta)-LOG_GAMMA(alpha) &
        - (2*alpha+2)*LOG(sigma) - (2d0*beta+lambda*(z-mu)**2)/(2d0*sigma**2)
    END ASSOCIATE
  End Function normalInverseGamma_computeLogPDF

  Integer Function normalInverseGamma_fitWithData(this, Points, W) Result(info)
    class(mlf_distribution_normalInverseGamma), intent(inout) :: this
    real(c_double), intent(in) :: Points(:,:)
    real(c_double), intent(in), optional :: W(:)
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
End Module mlf_student_t_dist

